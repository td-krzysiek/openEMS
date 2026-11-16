/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/engine.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "tools/useful.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

using Loss = Operator_Ext_CartesianSubgrid_Conductivity;
static void Require(bool ok, const char *message)
{
	if (!ok)
	{
		throw std::runtime_error(message);
	}
}
template <class F> static void Reject(F call)
{
	try
	{
		call();
	}
	catch (const std::logic_error &)
	{
		return;
	}
	throw std::runtime_error("invalid conductivity input/phase was admitted");
}
struct Fixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	Fixture()
	{
		csx.GetGrid()->SetDeltaUnit(.001);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= 4; ++i)
			{
				csx.GetGrid()->AddDiscLine(n, i);
			}
		}
		signal.SetupDiracPulse(30e9);
		op.reset(Operator_CartesianSubgrid::New({1, 1, 1, 3, 3, 3}));
		op->SetExcitationSignal(&signal);
		op->SetEngineType(Engine::BASIC);
		Require(op->SetGeometryCSX(&csx), "guard fixture geometry rejected");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		op->SetBoundaryCondition(bc);
		Require(op->CalcECOperator() == 0, "guard fixture setup rejected");
	}
};
static Loss::Sampler Constant(double electric, double magnetic)
{
	return [=](const CartesianSubgrid::Point &, const CartesianSubgrid::Point &)
	{
		return Loss::Coefficients{electric, magnetic};
	};
}
int main()
{
	try
	{
		SetThreadLimit(1);
		Fixture f;
		const auto surface = []
		{
			Loss::SurfaceAdmittance value;
			value.start = {.001, .001, .001};
			value.stop = {.0015, .001, .001};
			value.conductance = .01;
			value.poles = {{25, 1e-10}};
			return value;
		};
		unsigned int count = 0;
		const auto rejectSurface = [&](Loss::SurfaceAdmittance value)
		{
			Reject([&]
			{
				Loss bad(f.op.get(), Constant(0, 0));
				bad.SetSurfaceAdmittances({value});
				bad.BuildExtension();
			});
			++count;
		};
		for (double invalid : {-1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
		{
			auto value = surface();
			value.conductance = invalid;
			rejectSurface(value);
			value = surface();
			value.poles[0].resistance = invalid;
			rejectSurface(value);
			value = surface();
			value.poles[0].inductance = invalid;
			rejectSurface(value);
		}
		auto value = surface();
		value.poles[0].inductance = 0;
		rejectSurface(value);
		value = surface();
		value.stop = value.start;
		rejectSurface(value);
		value = surface();
		value.start = {0, 0, 0};
		value.stop = {.001, 0, 0};
		rejectSurface(value);
		value = surface();
		value.start = {.001, .001, .001};
		value.stop = {.003, .001, .001};
		rejectSurface(value); // A nonlocal path escapes the existing mass graph.
		value = surface();
		value.poles[0] = {1e308, 1e-308};
		rejectSurface(value);
		std::cout << "PASS: surface scalar/path/mass-graph/rate guards count=" << count << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}
