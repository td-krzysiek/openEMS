/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/cartesiansubgrid_advance.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using Grid = CartesianSubgrid;

static void Require(bool valid, const char *message)
{
	if (!valid)
	{
		throw std::runtime_error(message);
	}
}
template <class Action> static void Reject(Action action)
{
	try
	{
		action();
	}
	catch (const std::logic_error &)
	{
		return;
	}
	throw std::runtime_error("invalid compact advance was accepted");
}

enum class Mode { Normal, Missing, UnsupportedDevice, FailMagnetic, FailObservation };
struct Counters
{
	unsigned int builds = 0, factories = 0, cpuElectric = 0, cpuMagnetic = 0;
	unsigned int deviceElectric = 0, deviceMagnetic = 0, observations = 0;
	std::vector<unsigned int> hooks;
};

//! A transparent test recurrence, not a new physics implementation.
class CountedDevice : public GPU_CartesianSubgrid
{
  public:
	CountedDevice(std::unique_ptr<GPU_CartesianSubgrid> real, Counters &counts, Mode mode)
		: m_Real(std::move(real)), m_Counts(counts), m_Mode(mode)
	{
	}
	void Step(double dt, double drive) override
	{
		UpdateElectric(dt, drive);
		UpdateMagnetic(dt);
	}
	void UpdateElectric(double dt, double drive) override
	{
		++m_Counts.deviceElectric;
		m_Real->UpdateElectric(dt, drive);
	}
	void UpdateMagnetic(double dt) override
	{
		++m_Counts.deviceMagnetic;
		if (m_Mode == Mode::FailMagnetic)
		{
			throw std::runtime_error("injected device magnetic failure");
		}
		m_Real->UpdateMagnetic(dt);
	}
	void ReadFields(Grid &grid) override { m_Real->ReadFields(grid); }
	void RequireFixedTimestep(double dt) override { m_Real->RequireFixedTimestep(dt); }
	GPU_CartesianSubgridFields Fields() const override { return m_Real->Fields(); }
	std::array<double, 2> ReadLedger() override { return m_Real->ReadLedger(); }
	std::size_t StorageBytes() const override { return m_Real->StorageBytes(); }
  private:
	std::unique_ptr<GPU_CartesianSubgrid> m_Real;
	Counters &m_Counts;
	Mode m_Mode;
};

class CountedAdvance : public Engine_Ext_CartesianSubgrid_Advance
{
  public:
	CountedAdvance(Operator_Extension *op, Grid &grid, Counters &counts, Mode mode)
		: Engine_Ext_CartesianSubgrid_Advance(op), m_Grid(grid), m_Counts(counts), m_Mode(mode)
	{
	}
	void DoPreVoltageUpdates() override { m_Counts.hooks.push_back(0); }
	void DoPostVoltageUpdates() override { m_Counts.hooks.push_back(2); }
	void Apply2Voltages() override { m_Counts.hooks.push_back(3); }
	void DoPreCurrentUpdates() override { m_Counts.hooks.push_back(4); }
	void DoPostCurrentUpdates() override { m_Counts.hooks.push_back(6); }
	void Apply2Current() override { m_Counts.hooks.push_back(7); }
	void UpdateElectric(double dt, double drive) override
	{
		++m_Counts.cpuElectric;
		m_Counts.hooks.push_back(1);
		m_Grid.UpdateElectric(dt, drive);
	}
	void UpdateMagnetic(double dt) override
	{
		++m_Counts.cpuMagnetic;
		m_Counts.hooks.push_back(5);
		if (m_Mode == Mode::FailMagnetic)
		{
			throw std::runtime_error("injected CPU magnetic failure");
		}
		m_Grid.UpdateMagnetic(dt);
	}
	std::unique_ptr<GPU_CartesianSubgrid> CreateDevice(GPU_Backend &backend) override
	{
		++m_Counts.factories;
		if (m_Mode == Mode::UnsupportedDevice)
		{
			return {};
		}
		return std::make_unique<CountedDevice>(backend.CreateCartesianSubgrid(m_Grid), m_Counts, m_Mode);
	}
	void ReadDeviceDiagnostics(GPU_CartesianSubgrid &) override
	{
		++m_Counts.observations;
		if (m_Mode == Mode::FailObservation)
		{
			throw std::runtime_error("injected device observation failure");
		}
	}
  private:
	Grid &m_Grid;
	Counters &m_Counts;
	Mode m_Mode;
};

class CountedOperator : public Operator_Ext_CartesianSubgrid_Advance
{
  public:
	CountedOperator(Operator_CartesianSubgrid *op, Counters &counts, Mode mode)
		: Operator_Ext_CartesianSubgrid_Advance(op), m_Compact(op), m_Counts(counts), m_Mode(mode)
	{
	}
	bool BuildExtension() override
	{
		++m_Counts.builds;
		return true;
	}
	Engine_Ext_CartesianSubgrid_Advance *CreateEngineExtention() override
	{
		if (m_Mode == Mode::Missing)
		{
			return nullptr;
		}
		m_Eng_Ext = new CountedAdvance(this, m_Compact->Grid(), m_Counts, m_Mode);
		return static_cast<Engine_Ext_CartesianSubgrid_Advance *>(m_Eng_Ext);
	}
	std::string GetExtensionName() const override { return "test-only compact advance"; }
  private:
	Operator_CartesianSubgrid *m_Compact;
	Counters &m_Counts;
	Mode m_Mode;
};

struct FixtureOptions
{
	Mode mode = Mode::Normal;
	bool duplicate = false, absorb = false, provider = true;
};

struct Fixture
{
	ContinuousStructure geometry;
	Excitation signal;
	Counters counts;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	Fixture(bool device, const SubgridDeviceTarget &target, FixtureOptions options = {})
	{
		geometry.GetGrid()->SetDeltaUnit(.001);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= (options.absorb && n == 0 ? 14u : 8u); ++i)
			{
				geometry.GetGrid()->AddDiscLine(n, i);
			}
		}
		signal.SetupDiracPulse(30e9);
		op.reset(Operator_CartesianSubgrid::New({2, 2, 2, 6, 6, 6}));
		op->SetEngineType(device ? Engine::GPU : Engine::BASIC);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&geometry), "advance fixture geometry rejected");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		if (options.absorb)
		{
			bc[1] = 3;
			op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(op.get(), {0, 4, 0, 0, 0, 0}));
		}
		op->SetBoundaryCondition(bc);
		if (options.provider)
		{
			op->AddExtension(new CountedOperator(op.get(), counts, options.mode));
		}
		if (options.duplicate)
		{
			op->AddExtension(new CountedOperator(op.get(), counts, options.mode));
		}
		Require(op->CalcECOperator() == 0, "advance fixture setup failed");
		signal.SetupDiracPulse(30e9);
		Require(signal.buildExcitationSignal(1024) && signal.GetLength() > 0 &&
			std::any_of(signal.GetVoltageSignal(), signal.GetVoltageSignal() + signal.GetLength(),
				[](double value) { return value != 0; }),
			"advance control requires a real nonzero built excitation");
		op->Grid().SetElectric([](unsigned int n, const Grid::Point &p)
		{
			return std::sin(1234 * p[0] + 441 * p[1] - 731 * p[2] + n);
		});
		std::size_t source = 0;
		double distance = 1;
		for (std::size_t i = 0; i < op->Grid().EdgeCount(); ++i)
		{
			const auto point = op->Grid().EdgePosition(i);
			const double d = std::abs(point[0] - .004) + std::abs(point[1] - .004) + std::abs(point[2] - .004);
			if (d < distance)
			{
				distance = d;
				source = i;
			}
		}
		op->Grid().SetSources({{source, .01}});
		if (device)
		{
			std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
			Require(bool(backend), "requested advance backend unavailable");
			target.Verify(*backend);
			op->SetGPUBackend(backend.release());
		}
	}
};

static void Control(bool device, const SubgridDeviceTarget &target, bool sourceOnly = false)
{
	Fixture fixture(device, target);
	if (sourceOnly)
	{
		fixture.op->Grid().SetElectric([](unsigned int, const Grid::Point &) { return 0; });
		fixture.op->Grid().SetMagnetic([](unsigned int, const Grid::Point &) { return 0; });
	}
	auto expected = fixture.op->Grid();
	std::unique_ptr<Engine> engine(fixture.op->CreateEngine());
	constexpr unsigned int steps = 512;
	const double dt = fixture.op->GetTimestep();
	for (unsigned int i = 0; i < steps; ++i)
	{
		const double drive = i < fixture.signal.GetLength() ? fixture.signal.GetVoltageSignal()[i] : 0;
		expected.Step(dt, drive);
	}
	engine->IterateTS(steps / 2);
	engine->IterateTS(steps / 2);
	const auto a = expected.DeviceData(), b = fixture.op->Grid().DeviceData();
	double delta = 0, scale = 0;
	for (const auto pair : {std::make_pair(&a.electric, &b.electric), std::make_pair(&a.magnetic, &b.magnetic)})
	{
		for (std::size_t i = 0; i < pair.first->size(); ++i)
		{
			delta = std::max(delta, std::abs(double(pair.first->at(i)) - pair.second->at(i)));
			scale = std::max(scale, std::abs(double(pair.first->at(i))));
		}
	}
	Require(scale > 0 && delta / scale < 2e-4, "advance boundary changed fields or double-advanced a timestep");
	const auto &c = fixture.counts;
	Require(c.builds == 1, "operator half was not built exactly once");
	if (device)
	{
		Require(c.factories == 1 && c.deviceElectric == steps && c.deviceMagnetic == steps &&
			c.observations == 2 && c.cpuElectric == 0 && c.cpuMagnetic == 0 && c.hooks.empty(),
			"device recurrence used CPU updates or missed observation");
	}
	else
	{
		Require(c.cpuElectric == steps && c.cpuMagnetic == steps && c.factories == 0 &&
			c.deviceElectric == 0 && c.hooks.size() == steps * 8, "CPU advance call count differs");
		for (std::size_t i = 0; i < c.hooks.size(); ++i)
		{
			Require(c.hooks[i] == i % 8, "CPU pre/advance/post/apply order differs");
		}
	}
	std::cout << "ADVANCE source_only=" << sourceOnly << " device=" << device << " backend=" << target.backend
		<< " field_error=" << delta / scale << std::endl;
}

static void AbsorberControl(bool device, const SubgridDeviceTarget &target)
{
	FixtureOptions options;
	options.absorb = true;
	options.provider = false;
	Fixture baseline(device, target, options);
	options.provider = true;
	Fixture advanced(device, target, options);
	Require(baseline.op->GetTimestep() == advanced.op->GetTimestep(), "PML advance timestep differs");
	std::unique_ptr<Engine> first(baseline.op->CreateEngine()), second(advanced.op->CreateEngine());
	first->IterateTS(256);
	second->IterateTS(256);
	const auto a = baseline.op->Grid().DeviceData(), b = advanced.op->Grid().DeviceData();
	Require(a.electric == b.electric && a.magnetic == b.magnetic,
		"transparent advance changed PML half-step ordering");
	Require(device ? advanced.counts.deviceElectric == 256 && advanced.counts.cpuElectric == 0 :
		advanced.counts.cpuElectric == 256, "PML wrapper bypassed advance provider");
	std::cout << "PASS: transparent advance preserves existing PML fields device=" << device << std::endl;
}

static void Guards(const SubgridDeviceTarget &target)
{
	for (bool device : {false, true})
	{
		for (Mode mode : {Mode::Normal, Mode::Missing})
		{
			Fixture fixture(device, target, {mode, mode == Mode::Normal});
			Reject([&] { std::unique_ptr<Engine> engine(fixture.op->CreateEngine()); });
		}
		Fixture fixture(device, target, {Mode::FailMagnetic});
		std::unique_ptr<Engine> engine(fixture.op->CreateEngine());
		bool failed = false;
		try
		{
			engine->IterateTS(1);
		}
		catch (const std::runtime_error &)
		{
			failed = true;
		}
		Require(failed, "injected advance failure was ignored");
		Reject([&] { engine->IterateTS(1); });
		Reject([&] { engine->GetVolt(0, 2, 2, 2); });
		Reject([&] { engine->GetCurr(0, 2, 2, 2); });
	}
	Fixture unavailable(true, target, {Mode::UnsupportedDevice});
	Reject([&] { std::unique_ptr<Engine> engine(unavailable.op->CreateEngine()); });
	Require(unavailable.counts.cpuElectric == 0, "unsupported device silently fell back to CPU");
	Fixture observation(true, target, {Mode::FailObservation});
	std::unique_ptr<Engine> engine(observation.op->CreateEngine());
	try
	{
		engine->IterateTS(1);
		throw std::logic_error("observation failure was ignored");
	}
	catch (const std::runtime_error &)
	{
	}
	Reject([&] { engine->IterateTS(1); });
	std::cout << "PASS: duplicate/missing/device/failure/observation guards" << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "usage: test_cartesian_subgrid_device_advance [--metal|--cuda]");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromFlag(argv[1]) : SubgridDeviceTarget{};
		Control(false, target);
		Control(false, target, true);
		Control(true, target);
		Control(true, target, true);
		AbsorberControl(false, target);
		AbsorberControl(true, target);
		Guards(target);
		std::cout << "PASS: paired compact advance boundary, not conductivity implementation" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}
