/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef TEST_SUBGRID_NATIVE_GATHER_FIXTURE_H
#define TEST_SUBGRID_NATIVE_GATHER_FIXTURE_H
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid_runtime.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/engine.h"
#include "FDTD/engine_interface_fdtd.h"
#include "FDTD/subgrid_native_field_gather.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_conductivity.h"
#include "subgrid_device_target.h"
#include <limits>
#include "openems.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <iostream>
#include <memory>
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
	throw std::runtime_error("invalid native field gather operation accepted");
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
				geometry.GetGrid()->AddDiscLine(n, offset + i + .01 * i * i);
			}
		}
		signal.SetupDiracPulse(30e9);
		const auto coordinate = [offset](double k) { return offset + k + .01 * k * k; };
		op.reset(Operator_CartesianSubgrid::New({coordinate(6), coordinate(6), coordinate(6),
			coordinate(8), coordinate(8), coordinate(8)}));
		Grid::Axes children;
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i < 14; ++i)
			{
				children[n].push_back(coordinate(i));
				children[n].push_back(coordinate(i) + .3 * (coordinate(i + 1) - coordinate(i)));
			}
			children[n].push_back(coordinate(14));
		}
		op->SetChildLines(children);
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
				[unit, coordinate, offset](const Grid::Point &lo, const Grid::Point &hi)
				{
					bool inside = true;
					for (unsigned int n = 0; n < 3; ++n)
					{
						const double p = .5 * (lo[n] + hi[n]) / unit;
						inside = inside && p >= coordinate(6) - offset && p < coordinate(8) - offset;
					}
					return inside ? Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{.02, 50} :
						Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{};
				}));
		}
		Require(op->CalcECOperator() == 0, "native domain operator failed");
		signal.SetupDiracPulse(30e9);
		Require(signal.buildExcitationSignal(1024), "native domain waveform failed");
		op->Grid().SetElectric([unit, coordinate, offset](unsigned int n, const Grid::Point &p)
		{
			return (n == 1 ? 1. : -.5) * std::exp(-std::pow((p[0]/unit-coordinate(7)+offset)/2, 2) -
				std::pow((p[1]/unit-coordinate(7)+offset)/2, 2) - std::pow((p[2]/unit-coordinate(7)+offset)/2, 2));
		});
		op->Grid().SetSources({{op->Grid().EdgeCount() / 2, .01}});
		if (domains)
		{
			std::vector<std::unique_ptr<GPU_Backend>> backends;
			for (std::size_t k = 0; k < 2; ++k)
			{
				target.SelectDomain(k);
				backends.emplace_back(GPU_Backend::New(target.backend));
				Require(bool(backends.back()), "native worker unavailable");
				target.Verify(*backends.back());
				std::cout << "GATHER_WORKER " << k << " " << backends.back()->GetName() << std::endl;
			}
			op->SetDomainExecution(axis, {coordinate(6.1)}, std::move(backends));
		}
		else
		{
			target.SelectDomain(0);
			std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
			Require(bool(backend), "native single backend unavailable");
			target.Verify(*backend);
			std::cout << "GATHER_SINGLE " << backend->GetName() << std::endl;
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

#endif
