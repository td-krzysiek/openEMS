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
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}
