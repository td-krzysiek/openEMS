// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_DOMAIN_PERFORMANCE_H
#define OPENEMS_DOMAIN_PERFORMANCE_H

#include "distributed/domain_plan.h"
#include "gpu_backend.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <array>

// Optional observations only: never changes ownership or waits for a timer.
// Header-only so the standalone fake-device test exercises the actual writer.
class DomainPerformance
{
public:
	using Clock = std::chrono::steady_clock;
	using Time = Clock::time_point;
	enum ExtensionHook
	{
		PreVoltage, PostVoltage, ApplyVoltage,
		PreCurrent, PostCurrent, ApplyCurrent, SynchronizeHook, ExtensionHookCount
	};
	static constexpr size_t NoExtension = size_t(-1);
	enum HostMetric
	{
		UpdateWait,
		Pack,
		Unpack,
		UnpackWait,
		Boundary,
		HostMetricCount
	};
	struct Interval
	{
		std::unique_ptr<GPU_DeviceTimer> timer;
		bool extension = false;
		bool ended = false;
		size_t extension_id = NoExtension;
		ExtensionHook hook = PreVoltage;
		double submit_seconds = 0;
		Time started;
	};

	static std::unique_ptr<DomainPerformance> FromEnvironment(
		const openems::distributed::PartitionPlan& plan, const std::vector<GPU_Backend*>& devices)
	{
		const char* path = std::getenv("OPENEMS_DOMAIN_PROFILE");
		if (!path || !*path)
		{
			return nullptr;
		}
		try
		{
			return std::unique_ptr<DomainPerformance>(new DomainPerformance(path, plan, devices,
				Setting("OPENEMS_DOMAIN_PROFILE_INTERVAL", 4096, 128),
				Setting("OPENEMS_DOMAIN_PROFILE_WARMUP", 4096, 0),
				Setting("OPENEMS_DOMAIN_PROFILE_EXTENSIONS", 0, 0) != 0));
		}
		catch (const std::exception& error)
		{
			std::cerr << "Domain profiling disabled: " << error.what() << std::endl;
			return nullptr;
		}
	}

	DomainPerformance(const std::string& path, const openems::distributed::PartitionPlan& plan,
		const std::vector<GPU_Backend*>& devices, uint64_t interval, uint64_t warmup,
		bool extension_detail = false)
		: m_Devices(devices), m_Data(devices.size()), m_Interval(interval), m_Warmup(warmup),
		  m_ExtensionDetail(extension_detail)
	{
		if (!interval || devices.size() != plan.Domains().size())
		{
			throw std::invalid_argument("invalid domain profile configuration");
		}
		// Refuse to truncate a previous run's evidence. No file is opened when disabled.
		m_File = std::fopen(path.c_str(), "wx");
		if (!m_File)
		{
			throw std::runtime_error("cannot exclusively create domain profile: " + path);
		}
		std::setvbuf(m_File, nullptr, _IOFBF, 65536);
		std::ostringstream out;
		out.imbue(std::locale::classic());
		out << "{\"type\":\"header\",\"schema\":\"openems.domain-performance.v1\","
			<< "\"solver_version\":" << Quote(GIT_VERSION) << ",\"interval_steps\":" << interval
			<< ",\"warmup_steps\":" << warmup << ",\"repartition_during_run\":false,\"grid\":["
			<< plan.Shape().nx << ',' << plan.Shape().ny << ',' << plan.Shape().nz
			<< "],\"halo_mode\":" << Quote(std::getenv("OPENEMS_DOMAIN_HALOS") ?
			    std::getenv("OPENEMS_DOMAIN_HALOS") : "serial")
			<< ",\"host_buffer_mode\":" << Quote(std::getenv("OPENEMS_DOMAIN_HOST_BUFFERS") ?
			    std::getenv("OPENEMS_DOMAIN_HOST_BUFFERS") : "copy") << ",\"domains\":[";
		for (size_t i = 0; i < devices.size(); ++i)
		{
			const auto& domain = plan.Domains()[i];
			out << (i ? "," : "") << "{\"id\":" << i
				<< ",\"backend\":" << Quote(devices[i]->GetName()) << ",\"begin\":" << domain.begin
				<< ",\"end\":" << domain.end
				<< ",\"owned_nodes\":" << uint64_t(domain.end - domain.begin) * plan.Shape().Plane()
				<< '}';
		}
		out << "],\"extension_detail\":" << (m_ExtensionDetail ? "true" : "false")
			<< ",\"device_time_scope\":\"main_updates_and_native_extension_hooks_only\","
			<< "\"host_wait_is_gpu_idle\":false,\"added_device_waits\":false}";
		Write(out.str());
		if (m_File)
		{
			Flush();
		}
	}
	~DomainPerformance()
	{
		if (m_File)
		{
			Write("{\"type\":\"closed\",\"samples\":" + std::to_string(m_Samples) + "}");
			if (m_File)
			{
				std::fclose(m_File);
			}
		}
	}
	bool Active() const
	{
		return m_Active;
	}
	//! Labels are registered once at setup; existing event pairs are reused unchanged.
	size_t RegisterExtension(const std::string& name)
	{
		if (!m_ExtensionDetail)
		{
			return NoExtension;
		}
		m_ExtensionNames.push_back(name);
		for (auto& data : m_Data)
		{
			data.extensions.emplace_back();
		}
		return m_ExtensionNames.size() - 1;
	}
	Interval* BeginExtension(size_t domain, size_t extension, ExtensionHook hook)
	{
		auto* sample = Begin(domain, true);
		if (sample && extension != NoExtension)
		{
			sample->extension_id = extension;
			sample->hook = hook;
		}
		return sample;
	}
	static Time Now()
	{
		return Clock::now();
	}
	static double Seconds(Time since)
	{
		return std::chrono::duration<double>(Now() - since).count();
	}
	void StartStep(uint64_t step)
	{
		m_Step = step;
		m_Active = m_File && step >= m_Warmup && (step - m_Warmup) % m_Interval == 0;
		if (!m_Active)
		{
			return;
		}
		m_Exchange = 0;
		for (auto& data : m_Data)
		{
			data.used = 0;
			data.main_submit = data.extension_submit = 0;
			data.sent = data.received = 0;
			for (auto& extension : data.extensions)
			{
				extension = {};
			}
			for (double& value : data.host)
			{
				value = 0;
			}
		}
	}
	Interval* Begin(size_t domain, bool extension)
	{
		if (!m_Active)
		{
			return nullptr;
		}
		auto& data = m_Data.at(domain);
		if (data.used == data.intervals.size())
		{
			data.intervals.emplace_back();
			try
			{
				data.intervals.back().timer.reset(m_Devices.at(domain)->NewDeviceTimer());
			}
			catch (const std::exception& error)
			{
				Warn(error.what());
			}
		}
		auto& sample = data.intervals[data.used++];
		sample.extension = extension;
		sample.ended = false;
		sample.extension_id = NoExtension;
		sample.submit_seconds = 0;
		if (sample.timer)
		{
			try
			{
				sample.timer->Begin();
			}
			catch (const std::exception& error)
			{
				sample.timer.reset();
				Warn(error.what());
			}
		}
		sample.started = Now();
		return &sample;
	}
	void End(size_t domain, Interval* sample)
	{
		if (!sample)
		{
			return;
		}
		auto& data = m_Data.at(domain);
		sample->submit_seconds = Seconds(sample->started);
		(sample->extension ? data.extension_submit : data.main_submit) += sample->submit_seconds;
		if (sample->timer)
		{
			try
			{
				sample->timer->End();
				sample->ended = true;
			}
			catch (const std::exception& error)
			{
				Warn(error.what());
			}
		}
	}
	void Host(size_t domain, HostMetric metric, Time started, uint64_t bytes = 0)
	{
		if (m_Active)
		{
			auto& data = m_Data.at(domain);
			data.host[metric] += Seconds(started);
			if (metric == Pack)
			{
				data.sent += bytes;
			}
			if (metric == Unpack)
			{
				data.received += bytes;
			}
		}
	}
	void Exchange(Time started)
	{
		if (m_Active)
		{
			m_Exchange += Seconds(started);
		}
	}
	void ExchangeSeconds(double elapsed)
	{
		if (m_Active)
		{
			m_Exchange += elapsed;
		}
	}
	void FinishStep()
	{
		if (!m_Active)
		{
			return;
		}
		std::ostringstream out;
		out.imbue(std::locale::classic());
		out << std::setprecision(17) << "{\"type\":\"sample\",\"step\":" << m_Step
			<< ",\"epoch_milliseconds\":"
			<< std::chrono::duration_cast<std::chrono::milliseconds>(
				   std::chrono::system_clock::now().time_since_epoch())
				   .count()
			<< ",\"exchange_wall_seconds\":" << m_Exchange << ",\"domains\":[";
		for (size_t i = 0; i < m_Data.size(); ++i)
		{
			auto& data = m_Data[i];
			double main = 0, extension = 0;
			size_t missing = 0;
			for (size_t j = 0; j < data.used; ++j)
			{
				auto& sample = data.intervals[j];
				double seconds = 0;
				bool ready = false;
				try
				{
					ready = sample.ended && sample.timer && sample.timer->ReadSeconds(seconds) &&
						std::isfinite(seconds) && seconds >= 0;
				}
				catch (const std::exception& error)
				{
					Warn(error.what());
				}
				if (!ready)
				{
					++missing;
				}
				else
				{
					(sample.extension ? extension : main) += seconds;
				}
				if (sample.extension_id != NoExtension)
				{
					auto& detail = data.extensions.at(sample.extension_id).at(sample.hook);
					++detail.calls;
					detail.submit += sample.submit_seconds;
					if (ready)
					{
						detail.device += seconds;
					}
					else
					{
						++detail.missing;
					}
				}
			}
			out << (i ? "," : "") << "{\"id\":" << i
				<< ",\"device_timing_complete\":" << (!missing && data.used ? "true" : "false")
				<< ",\"missing_device_intervals\":" << missing
				<< ",\"interval_count\":" << data.used
				<< ",\"main_device_seconds\":" << (missing || !data.used ? "null" : Number(main))
				<< ",\"extension_device_seconds\":"
				<< (missing || !data.used ? "null" : Number(extension))
				<< ",\"main_submit_host_seconds\":" << data.main_submit
				<< ",\"extension_submit_host_seconds\":" << data.extension_submit
				<< ",\"update_wait_host_seconds\":" << data.host[UpdateWait]
				<< ",\"pack_host_seconds\":" << data.host[Pack]
				<< ",\"unpack_host_seconds\":" << data.host[Unpack]
				<< ",\"unpack_wait_host_seconds\":" << data.host[UnpackWait]
				<< ",\"boundary_host_seconds\":" << data.host[Boundary]
				<< ",\"halo_sent_bytes\":" << data.sent
				<< ",\"halo_received_bytes\":" << data.received;
			if (m_ExtensionDetail)
			{
				WriteExtensions(out, data);
			}
			out << '}';
		}
		out << "]}";
		Write(out.str());
		if (++m_Samples % 8 == 0 && m_File)
		{
			Flush();
		}
		m_Active = false;
	}

private:
	struct HookData
	{
		double device = 0, submit = 0;
		size_t calls = 0, missing = 0;
	};
	struct Data
	{
		std::vector<Interval> intervals;
		size_t used = 0;
		double main_submit = 0, extension_submit = 0;
		double host[HostMetricCount] = {};
		uint64_t sent = 0, received = 0;
		std::vector<std::array<HookData, ExtensionHookCount>> extensions;
	};
	void WriteExtensions(std::ostream& out, const Data& data) const
	{
		static const char* names[] = {"pre_voltage", "post_voltage", "apply_voltage",
			"pre_current", "post_current", "apply_current", "synchronize"};
		out << ",\"extensions\":[";
		for (size_t id = 0; id < data.extensions.size(); ++id)
		{
			HookData total;
			for (const auto& hook : data.extensions[id])
			{
				total.device += hook.device;
				total.submit += hook.submit;
				total.calls += hook.calls;
				total.missing += hook.missing;
			}
			out << (id ? "," : "") << "{\"id\":" << id << ",\"name\":" << Quote(m_ExtensionNames[id]);
			auto values = [&](const HookData& hook)
			{
				out << ",\"calls\":" << hook.calls << ",\"missing_device_intervals\":" << hook.missing
					<< ",\"device_seconds\":" << (!hook.calls || hook.missing ? "null" : Number(hook.device))
					<< ",\"submit_host_seconds\":" << hook.submit;
			};
			values(total);
			out << ",\"hooks\":[";
			for (size_t hook = 0; hook < ExtensionHookCount; ++hook)
			{
				out << (hook ? "," : "") << "{\"name\":" << Quote(names[hook]);
				values(data.extensions[id][hook]);
				out << '}';
			}
			out << "]}";
		}
		out << ']';
	}
	static uint64_t Setting(const char* name, uint64_t fallback, uint64_t minimum)
	{
		const char* raw = std::getenv(name);
		if (!raw)
		{
			return fallback;
		}
		const std::string value(raw);
		if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
		{
			throw std::invalid_argument(std::string("invalid ") + name);
		}
		const auto number = std::stoull(value);
		if (number < minimum)
		{
			throw std::invalid_argument(std::string("invalid ") + name);
		}
		return number;
	}
	static std::string Number(double value)
	{
		std::ostringstream out;
		out.imbue(std::locale::classic());
		out << std::setprecision(17) << value;
		return out.str();
	}
	static std::string Quote(const std::string& value)
	{
		std::ostringstream out;
		out.imbue(std::locale::classic());
		out << '"';
		for (unsigned char c : value)
		{
			if (c == '"' || c == '\\')
			{
				out << '\\' << char(c);
			}
			else if (c < 32)
			{
				out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << unsigned(c);
			}
			else
			{
				out << char(c);
			}
		}
		out << '"';
		return out.str();
	}
	void Warn(const std::string& reason)
	{
		if (!m_Warned)
		{
			std::cerr << "Domain profiling warning: " << reason << std::endl;
			m_Warned = true;
		}
	}
	void Write(const std::string& value)
	{
		if (m_File &&
			(std::fwrite(value.data(), 1, value.size(), m_File) != value.size() ||
				std::fputc('\n', m_File) == EOF))
		{
			Warn("output failed; simulation continues without further profile writes");
			std::fclose(m_File);
			m_File = nullptr;
		}
	}
	void Flush()
	{
		if (m_File && std::fflush(m_File) != 0)
		{
			Warn("profile flush failed; simulation continues without further profile writes");
			std::fclose(m_File);
			m_File = nullptr;
		}
	}
	std::vector<GPU_Backend*> m_Devices;
	std::vector<Data> m_Data;
	std::vector<std::string> m_ExtensionNames;
	FILE* m_File = nullptr;
	uint64_t m_Interval, m_Warmup, m_Step = 0, m_Samples = 0;
	bool m_Active = false, m_Warned = false;
	bool m_ExtensionDetail = false;
	double m_Exchange = 0;
};
#endif
