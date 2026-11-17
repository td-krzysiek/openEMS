/* SPDX-License-Identifier: GPL-3.0-or-later */
// Read actual ordinary coefficients and compact upload metadata; no private API.
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPropLumpedElement.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/engine.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

using Grid = CartesianSubgrid;
static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}
struct Fixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator> op;
	Grid::Axes lines;
	Grid::Point start{}, stop{};
	Fixture(bool compact, bool graded, unsigned int direction, bool reverse)
	{
		csx.GetGrid()->SetDeltaUnit(.001);
		for (unsigned int n = 0; n < 3; ++n)
		{
			const double base[5] = {0, graded ? .5 : 1, graded ? 1.5 : 2, graded ? 3 : 3, 4};
			for (unsigned int k = 0; k < 4; ++k)
			{
				lines[n].push_back(base[k]);
				lines[n].push_back(.5 * (base[k] + base[k + 1]));
			}
			lines[n].push_back(4);
			for (double value : lines[n])
			{
				csx.GetGrid()->AddDiscLine(n, value);
			}
			start[n] = lines[n][n == direction ? 1 : 4];
			stop[n] = lines[n][n == direction ? 7 : 4];
		}
		auto *resistor = new CSPropLumpedElement(csx.GetParameterSet());
		resistor->SetResistance(50);
		resistor->SetDirection(direction);
		resistor->SetCaps(false);
		resistor->SetLEtype(CSPropLumpedElement::PARALLEL);
		csx.AddProperty(resistor);
		auto *box = new CSPrimBox(csx.GetParameterSet(), resistor);
		for (unsigned int n = 0; n < 3; ++n)
		{
			box->SetCoord(2 * n, reverse ? stop[n] : start[n]);
			box->SetCoord(2 * n + 1, reverse ? start[n] : stop[n]);
		}
		signal.SetupDiracPulse(30e9);
		if (compact)
		{
			op.reset(Operator_CartesianSubgrid::New({.5, .5, .5, 3, 3, 3}, Grid::Index{1, 1, 1}));
			op->SetTimestepFactor(.5);
		}
		else
		{
			op.reset(Operator::New(1));
		}
		op->SetEngineType(Engine::BASIC);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&csx), "line diagnostic geometry failed");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		op->SetBoundaryCondition(bc);
	}
};
static void Control(unsigned int direction, bool graded, bool reverse)
{
	Fixture compact(true, graded, direction, reverse), ordinary(false, graded, direction, reverse);
	Require(compact.op->CalcECOperator() == 0, "native line diagnostic preparation failed");
	const double dt = compact.op->GetTimestep();
	ordinary.op->SetTimestep(dt);
	Require(ordinary.op->CalcECOperator() == 0 && ordinary.op->GetTimestepValid(), "ordinary line diagnostic preparation failed");
	auto &grid = dynamic_cast<Operator_CartesianSubgrid *>(compact.op.get())->Grid();
	Grid::Point a = compact.start, b = compact.stop;
	for (unsigned int n = 0; n < 3; ++n)
	{
		a[n] *= .001;
		b[n] *= .001;
	}
	const auto segments = grid.LocalEdgePathSegments(reverse ? b : a, reverse ? a : b);
	const auto upload = grid.DeviceSetupMetadata();
	Require(segments.size() == 6 && upload.rootConductance.size() == 6, "native diagnostic line footprint differs");
	double error = 0;
	for (std::size_t k = 0; k < segments.size(); ++k)
	{
		unsigned int pos[3];
		for (unsigned int n = 0; n < 3; ++n)
		{
			const double x = std::min(segments[k].first[n], segments[k].second[n]) / .001;
			auto found = std::min_element(ordinary.lines[n].begin(), ordinary.lines[n].end(), [&](double u, double v)
			{
				return std::abs(u - x) < std::abs(v - x);
			});
			pos[n] = unsigned(found - ordinary.lines[n].begin());
		}
		const double expected = (1 - ordinary.op->GetVV(direction, pos)) / ordinary.op->GetVI(direction, pos);
		const double actual = double(upload.rootConductance[k]) * upload.rootConductance[k];
		error = std::max(error, std::abs(actual - expected) / expected);
		std::cout << std::setprecision(12) << "LINE_EDGE axis=" << direction << " graded=" << graded
			<< " reverse=" << reverse << " segment=" << k << " ordinary_g=" << expected
			<< " native_g=" << actual << std::endl;
	}
	Require(error < 2e-5, "native line conductance differs from actual ordinary coefficients");
	std::cout << "LINE_NORMALIZATION axis=" << direction << " graded=" << graded
		<< " reverse=" << reverse << " segments=" << segments.size()
		<< " relative_coefficient_error=" << error << " fdtd_steps=0" << std::endl;
}
int main()
{
	try
	{
		SetThreadLimit(1);
		for (unsigned int direction = 0; direction < 3; ++direction)
		{
			for (bool graded : {false, true})
			{
				for (bool reverse : {false, true})
				{
					Control(direction, graded, reverse);
				}
			}
		}
		std::cout << "PASS: native line normalization matches actual ordinary coefficients on uniform and graded axes" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}
