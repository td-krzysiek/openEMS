/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid_runtime.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_conductivity.h"
#include "subgrid_device_target.h"
#include "openems.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <thread>

using Grid = CartesianSubgrid;
static void Require(bool value, const char *message)
{
	if (!value)
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
	throw std::runtime_error("invalid native domain option was accepted");
}

struct Fixture
{
	ContinuousStructure geometry;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	std::unique_ptr<Engine> engine;
	Engine_Ext_CartesianSubgrid_PML *pml = nullptr;
	Fixture(const SubgridDeviceTarget &target, unsigned int axis, bool loss,
		double unit, double offset, bool domains)
	{
		geometry.GetGrid()->SetDeltaUnit(unit);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= 14; ++i)
			{
				geometry.GetGrid()->AddDiscLine(n, offset + i);
			}
		}
		signal.SetupDiracPulse(30e9);
		op.reset(Operator_CartesianSubgrid::New({offset+6, offset+6, offset+6,
			offset+8, offset+8, offset+8}));
		op->SetEngineType(Engine::GPU);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&geometry), "native domain geometry rejected");
		int boundaries[6] = {3, 3, 3, 3, 3, 3};
		op->SetBoundaryCondition(boundaries);
		op->SetTimestepFactor(.8);
		op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(op.get(), {4, 4, 4, 4, 4, 4}));
		if (loss)
		{
			op->AddExtension(new Operator_Ext_CartesianSubgrid_Conductivity(op.get(),
				[unit](const Grid::Point &lo, const Grid::Point &hi)
				{
					bool inside = true;
					for (unsigned int n = 0; n < 3; ++n)
					{
						const double p = .5 * (lo[n] + hi[n]) / unit;
						inside = inside && p >= 6 && p < 8;
					}
					return inside ? Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{.02, 50} :
						Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{};
				}));
		}
		Require(op->CalcECOperator() == 0, "native domain operator failed");
		signal.SetupDiracPulse(30e9);
		Require(signal.buildExcitationSignal(1024), "native domain waveform failed");
		op->Grid().SetElectric([unit](unsigned int n, const Grid::Point &p)
		{
			return (n == 1 ? 1. : -.5) * std::exp(-std::pow((p[0]/unit-7)/2, 2) -
				std::pow((p[1]/unit-7)/2, 2) - std::pow((p[2]/unit-7)/2, 2));
		});
		if (domains)
		{
			std::vector<std::unique_ptr<GPU_Backend>> backends;
			for (std::size_t k = 0; k < 2; ++k)
			{
				target.SelectDomain(k);
				backends.emplace_back(GPU_Backend::New(target.backend));
				Require(bool(backends.back()), "native worker unavailable");
				target.Verify(*backends.back());
			}
			op->SetDomainExecution(axis, {offset+6.1}, std::move(backends));
		}
		else
		{
			target.SelectDomain(0);
			std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
			Require(bool(backend), "native single backend unavailable");
			target.Verify(*backend);
			op->SetGPUBackend(backend.release());
		}
		engine.reset(op->CreateEngine());
		for (std::size_t k = 0; k < engine->GetExtensionCount(); ++k)
		{
			if (auto *extension = dynamic_cast<Engine_Ext_CartesianSubgrid_PML *>(engine->GetExtension(k)))
			{
				pml = extension;
			}
		}
		Require(pml, "native PML mirror missing");
	}
};

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool loss, bool domains)
{
	Fixture legacy(target, axis, loss, .001, 0, domains), sparse(target, axis, loss, .001, 0, domains);
	auto &grid = sparse.op->Grid();
	const auto &a = grid.BaseLines();
	std::vector<Grid::Integral> functions;
	for (unsigned int component = 0; component < 3; ++component)
	{
		for (double fraction : {.25, .5, .75})
		{
			const Grid::Point p{a[0][6] + fraction * (a[0][8] - a[0][6]), a[1][7], a[2][7]};
			functions.push_back(grid.ElectricSample(component, p));
			functions.push_back(grid.MagneticSample(component, p));
		}
	}
	functions.push_back(grid.VoltagePath({a[0][6], a[1][7], a[2][7]}, {a[0][8], a[1][7], a[2][7]}));
	functions.push_back(grid.CurrentContour(2, {a[0][6], a[1][6], a[2][7]}, {a[0][8], a[1][8], a[2][7]}));
	functions.push_back(functions.front());
	Grid::Integral bad;
	bad.weights = {{grid.EdgeCount(), 1}};
	Reject([&] { sparse.op->RegisterObservation(bad); });
	for (std::size_t i = 0; i < functions.size(); ++i)
	{
		Require(sparse.op->RegisterObservation(functions[i]) == i, "native registration changed order");
	}
	sparse.op->SealObservations();
	Reject([&] { sparse.op->SealObservations(); });
	Reject([&] { sparse.op->RegisterObservation(functions.front()); });
	Reject([&] { sparse.op->EvaluateObservation(functions.size()); });
	const auto initial = sparse.op->NativeObservationInfo();
	Require(initial.sealed && initial.device && initial.fullFieldReads == 0 && initial.selectedReadbacks == 1,
		"native sealing downloaded complete fields");
	double sampleError = 0, energyError = 0;
	const unsigned int stride[3] = {1, 1, 1};
	for (unsigned int batch = 0; batch < 4; ++batch)
	{
		legacy.engine->IterateTS(16);
		sparse.engine->IterateTS(16);
		for (std::size_t i = 0; i < functions.size(); ++i)
		{
			const auto expected = legacy.op->ObservationGrid().Evaluate(functions[i]);
			const auto actual = sparse.op->EvaluateObservation(i);
			sampleError = std::max(sampleError, std::abs(expected - actual) / std::max(std::abs(expected), 1e-25));
		}
		double e = 0, h = 0, es = 0, hs = 0;
		legacy.engine->CalcFastEnergy(stride, e, h);
		sparse.engine->CalcFastEnergy(stride, es, hs);
		energyError = std::max({energyError, std::abs(e - es) / std::max(std::abs(e), 1e-300),
			std::abs(h - hs) / std::max(std::abs(h), 1e-300)});
		const auto once = sparse.op->NativeObservationInfo();
		double repeatedE = 0, repeatedH = 0;
		sparse.engine->CalcFastEnergy(stride, repeatedE, repeatedH);
		const auto twice = sparse.op->NativeObservationInfo();
		Require(repeatedE == es && repeatedH == hs && once.energyReadbacks == twice.energyReadbacks &&
			once.energyDownloadBytes == twice.energyDownloadBytes, "same-epoch energy read was not cached");
		Require(twice.fullFieldReads == 0 && twice.selectedReadbacks == batch + 2 && twice.energyReadbacks == batch + 1,
			"sparse U/I/energy observation refreshed complete native fields");
	}
	Require(sampleError < 2e-12 && energyError < 2e-6, "native sparse fields/energy differ from eager mirrors");
	const auto selected = sparse.op->NativeObservationInfo();
	Require(selected.energyDownloadBytes > 0 && selected.energyDownloadBytes <= (domains ? 2 : 1) * 4 * 24,
		"native energy scalar transfer accounting failed");
	const auto x = sparse.op->ObservationGrid().DeviceStateData(), y = legacy.op->ObservationGrid().DeviceStateData();
	Require(x.electric == y.electric && x.magnetic == y.magnetic, "explicit native full mirror has stale fields");
	sparse.op->ObservationGrid();
	Require(sparse.op->NativeObservationInfo().fullFieldReads == 1, "same-epoch full mirror refreshed twice");
	for (unsigned int pass = 0; pass < 2; ++pass)
	{
		legacy.engine->IterateTS(1);
		sparse.engine->IterateTS(1);
		std::vector<std::thread> readers;
		std::array<double, 4> values{};
		for (unsigned int k = 0; k < 4; ++k)
		{
			readers.emplace_back([&, k]
			{
				const unsigned int point[3] = {6, 7, 7};
				values[k] = sparse.engine->GetVolt(1, point) + sparse.engine->GetCurr(2, point);
			});
		}
		for (auto &reader : readers)
		{
			reader.join();
		}
		const unsigned int point[3] = {6, 7, 7};
		const auto expected = legacy.engine->GetVolt(1, point) + legacy.engine->GetCurr(2, point);
		Require(std::all_of(values.begin(), values.end(), [&](double value) { return value == expected; }),
			"concurrent native fallback reads changed fields");
		Require(sparse.op->NativeObservationInfo().fullFieldReads == pass + 2,
			"native fallback failed to refresh exactly once per epoch");
	}
	sparse.op->InvalidateObservations();
	Reject([&] { sparse.op->ObservationGrid(); });
	Reject([&] { sparse.op->EvaluateObservation(0); });
	std::cout << "NATIVE_SPARSE_OBSERVATIONS axis=" << axis << " loss=" << loss << " domains=" << domains
		<< " steps=66 functionals=" << functions.size() << " sample_error=" << sampleError
		<< " energy_error=" << energyError << " sparse_batches=4 full_reads_before_fallback=0 fallback_reads=3"
		<< " energy_bytes=" << selected.energyDownloadBytes << " host_bytes=" << selected.hostObservationBytes << '\n';
}
int main(int argc, char **argv)
{
	try
	{
		Require(argc <= 2, "unexpected native sparse observation argument");
		const auto target = argc == 1 ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromDomainFlag(argv[1]);
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			for (bool loss : {false, true})
			{
				for (bool domains : {false, true})
				{
					Control(target, axis, loss, domains);
				}
			}
		}
		std::cout << "PASS: sealed native sparse probes, device energy and fresh concurrent full fallback\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}
