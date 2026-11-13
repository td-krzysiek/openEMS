// SPDX-License-Identifier: GPL-3.0-or-later
#include "FDTD/domain_performance.h"
#include <fstream>
#include <iterator>

struct Counts
{
	unsigned int allocated = 0, began = 0, ended = 0, read = 0, synchronized = 0;
	bool ready = true;
};
class FakeTimer : public GPU_DeviceTimer
{
public:
	explicit FakeTimer(Counts& counts) : m_Counts(counts)
	{
		++counts.allocated;
	}
	void Begin() override
	{
		++m_Counts.began;
	}
	void End() override
	{
		++m_Counts.ended;
	}
	bool ReadSeconds(double& value) override
	{
		++m_Counts.read;
		value = 0.125;
		return m_Counts.ready;
	}

private:
	Counts& m_Counts;
};
class FakeBackend : public GPU_Backend
{
public:
	Counts counts;
	bool supported = true;
	std::string GetName() const override
	{
		return "fake \"device\"\n";
	}
	bool Init(const Operator*) override
	{
		return true;
	}
	GPU_Backend* NewSubGridBackend() override
	{
		return nullptr;
	}
	void UpdateVoltages() override
	{
	}
	void UpdateCurrents() override
	{
	}
	void DownloadVoltages(ArrayLib::ArrayNIJK<float>&) override
	{
	}
	void DownloadCurrents(ArrayLib::ArrayNIJK<float>&) override
	{
	}
	void UploadVoltages(const ArrayLib::ArrayNIJK<float>&) override
	{
	}
	void UploadCurrents(const ArrayLib::ArrayNIJK<float>&) override
	{
	}
	void Synchronize() override
	{
		++counts.synchronized;
	}
	GPU_DeviceTimer* NewDeviceTimer() override
	{
		return supported ? new FakeTimer(counts) : nullptr;
	}
};
static void Require(bool value, const char* reason)
{
	if (!value)
	{
		throw std::runtime_error(reason);
	}
}
int main(int argc, char** argv)
{
	try
	{
		Require(argc == 2, "test needs a fresh output path");
		std::remove(argv[1]);
		FakeBackend a, b;
		std::vector<GPU_Backend*> devices = {&a, &b};
		openems::distributed::PartitionPlan plan(1, {9, 7, 5}, {0, 3, 9});
		{
			DomainPerformance profile(argv[1], plan, devices, 128, 256);
			for (uint64_t step = 0; step < 512; ++step)
			{
				profile.StartStep(step);
				for (size_t i = 0; i < devices.size(); ++i)
				{
					for (bool extension : {false, true})
					{
						auto* timing = profile.Begin(i, extension);
						profile.End(i, timing);
					}
					profile.Host(i, DomainPerformance::Pack, DomainPerformance::Now(), 280);
				}
				profile.FinishStep();
			}
			Require(a.counts.allocated == 2 && a.counts.began == 4 && a.counts.read == 4,
				"sampling or timer reuse failed");
			bool refused = false;
			try
			{
				DomainPerformance duplicate(argv[1], plan, devices, 128, 0);
			}
			catch (const std::runtime_error&)
			{
				refused = true;
			}
			Require(refused, "profile overwrote existing evidence");
			// A device which is not ready is reported as unavailable, never waited for.
			a.counts.ready = false;
			profile.StartStep(512);
			profile.End(0, profile.Begin(0, false));
			profile.FinishStep();
		}
		Require(a.counts.synchronized == 0 && b.counts.synchronized == 0,
			"profiler added device synchronization");
		std::ifstream stream(argv[1]);
		const std::string text((std::istreambuf_iterator<char>(stream)), {});
		Require(text.find("\"step\":256") != std::string::npos, "warmup sample absent");
		Require(text.find("\"step\":384") != std::string::npos, "second sample absent");
		Require(text.find("\"main_device_seconds\":0.125") != std::string::npos,
			"wrong device interval total");
		Require(text.find("\"main_device_seconds\":null") != std::string::npos,
			"unready device must not look like zero compute time");
		Require(text.find("\"type\":\"closed\",\"samples\":3") != std::string::npos,
			"buffered output was not finalized");
		// Detail labels must reuse the aggregate profiler's exact timer workload.
		Counts reference_counts;
		for (bool details : {false, true})
		{
			FakeBackend first, second;
			const std::string path = std::string(argv[1]) + (details ? ".detail" : ".aggregate");
			std::remove(path.c_str());
			{
				DomainPerformance profile(path, plan, {&first, &second}, 128, 0, details);
				const auto upml = profile.RegisterExtension("UPML \"region\"\n");
				const auto rlc = profile.RegisterExtension("Lumped RLC");
				for (uint64_t step : {uint64_t(0), uint64_t(128)})
				{
					profile.StartStep(step);
					for (size_t domain = 0; domain < 2; ++domain)
					{
						profile.End(domain, profile.Begin(domain, false));
						for (auto id : {upml, rlc})
						{
							for (auto hook : {DomainPerformance::PreVoltage, DomainPerformance::PostCurrent})
							{
								profile.End(domain, profile.BeginExtension(domain, id, hook));
							}
						}
					}
					profile.FinishStep();
				}
				if (details)
				{
					first.counts.ready = false;
					profile.StartStep(256);
					profile.End(0, profile.BeginExtension(0, upml, DomainPerformance::PostCurrent));
					profile.FinishStep();
				}
			}
			Require(first.counts.synchronized == 0 && second.counts.synchronized == 0,
				"extension detail added a device wait");
			if (!details)
			{
				reference_counts = first.counts;
			}
			else
			{
				Require(first.counts.allocated == reference_counts.allocated &&
					first.counts.began == reference_counts.began + 1 &&
					first.counts.read == reference_counts.read + 1,
					"extension labels added timers/events beyond the explicit unready control");
				std::ifstream detail_stream(path);
				const std::string detail_text((std::istreambuf_iterator<char>(detail_stream)), {});
				Require(detail_text.find("\"name\":\"Lumped RLC\",\"calls\":2,\"missing_device_intervals\":0,\"device_seconds\":0.25") != std::string::npos,
					"extension totals were not attributed to their own instance");
				Require(detail_text.find("\"name\":\"post_current\",\"calls\":1,\"missing_device_intervals\":1,\"device_seconds\":null") != std::string::npos,
					"unready extension interval must remain missing, not zero");
			}
		}
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}
