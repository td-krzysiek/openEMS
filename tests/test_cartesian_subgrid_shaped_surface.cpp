/* SPDX-License-Identifier: GPL-3.0-or-later */
#if __has_include("CSPropLossyMetal.h")
#include "FDTD/operator_cartesiansubgrid.h"
#include "subgrid_native_surface_metadata.h"
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPrimLinPoly.h"
#include "CSPrimCylinder.h"
#include "CSPrimCylindricalShell.h"
#include "CSPropLossyMetal.h"
#include "CSPropMaterial.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/engine.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using Grid = CartesianSubgrid;
using Point = std::array<double, 2>;
static constexpr double pi = 3.14159265358979323846;
static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}
enum class Shape
{
	Concave, Curve, Cylinder, Barrel, Bowtie, Backtrack, Flat, Sheet, Cut, RepeatedClosing
};
static bool Inside(Shape shape, double x, double y)
{
	const double r = std::hypot(x - 2, y - 2);
	if (shape == Shape::Concave)
	{
		return x > 1.25 && x < 2.75 && y > 1.25 && y < 2.75 && (x < 1.75 || y < 1.75);
	}
	if (shape == Shape::Curve)
	{
		return y > 2 && r > .4 && r < .7;
	}
	return r < .7 && (shape != Shape::Barrel || r > .4);
}
static std::vector<Point> Contour(Shape shape)
{
	if (shape == Shape::RepeatedClosing)
	{
		auto points = Contour(Shape::Concave);
		points.push_back(points.front());
		points.push_back(points.front());
		return points;
	}
	if (shape == Shape::Concave)
	{
		return {{1.25, 1.25}, {2.75, 1.25}, {2.75, 1.75}, {1.75, 1.75}, {1.75, 2.75}, {1.25, 2.75}};
	}
	if (shape == Shape::Bowtie)
	{
		return {{1.3, 1.3}, {2.7, 2.7}, {1.3, 2.7}, {2.7, 1.3}};
	}
	if (shape == Shape::Backtrack)
	{
		return {{1.3, 1.3}, {2., 1.3}, {1.6, 1.3}, {2.7, 2.7}, {1.3, 2.7}};
	}
	if (shape == Shape::Flat)
	{
		return {{1.3, 1.3}, {2., 2.}, {2.7, 2.7}};
	}
	std::vector<Point> points;
	for (unsigned int i = 0; i <= 32; ++i)
	{
		const double angle = pi * i / 32;
		points.push_back({2 + .7 * std::cos(angle), 2 + .7 * std::sin(angle)});
	}
	for (unsigned int i = 0; i <= 32; ++i)
	{
		const double angle = pi * (32 - i) / 32;
		points.push_back({2 + .4 * std::cos(angle), 2 + .4 * std::sin(angle)});
	}
	return points;
}
struct Settings
{
	Shape shape = Shape::Concave;
	unsigned int normal = 2;
	double thickness = 18e-6, unit = .001, offset = 0;
	bool reverse = false, bulk = false, closed = false;
};
struct Result
{
	std::vector<double> fields;
	std::vector<bool> occupied;
	double loss = 0, balance = 0, volume = 0;
	std::size_t faces = 0;
};
static Result Run(const Settings &settings, bool boxes, const SubgridDeviceTarget *target)
{
	const auto axis = settings.normal, u = (axis + 1) % 3, v = (axis + 2) % 3;
	const auto draw = [&](double mm)
	{
		return (mm + settings.offset) * .001 / settings.unit;
	};
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(settings.unit);
	for (unsigned int n : {u, v})
	{
		for (unsigned int i = 0; i <= 4; ++i)
		{
			csx.GetGrid()->AddDiscLine(n, draw(i));
		}
	}
	const double bottom = .77178;
	const double top = bottom + settings.thickness / .001;
	for (double z : {0., .5, bottom, top, 1.5, 2.})
	{
		csx.GetGrid()->AddDiscLine(axis, draw(z));
	}
	auto *metal = new CSPropLossyMetal(csx.GetParameterSet());
	metal->SetConductivity(5.8e7);
	csx.AddProperty(metal);
	const auto addBox = [&](CSProperties *property, double x0, double y0, double x1, double y1,
		double z0, double z1, int priority)
	{
		auto *box = new CSPrimBox(csx.GetParameterSet(), property);
		box->SetPriority(priority);
		box->SetCoord(2 * u, draw(x0));
		box->SetCoord(2 * u + 1, draw(x1));
		box->SetCoord(2 * v, draw(y0));
		box->SetCoord(2 * v + 1, draw(y1));
		box->SetCoord(2 * axis, draw(z0));
		box->SetCoord(2 * axis + 1, draw(z1));
	};
	std::vector<Point> points;
	CSPrimLinPoly *polygon = nullptr;
	if (boxes)
	{
		for (unsigned int i = 0; i < 8; ++i)
		{
			for (unsigned int j = 0; j < 8; ++j)
			{
				const double x = 1 + .25 * i, y = 1 + .25 * j;
				if (Inside(settings.shape, x + .125, y + .125))
				{
					addBox(metal, x, y, x + .25, y + .25, bottom, top, 10);
				}
			}
		}
	}
	else if (settings.shape == Shape::Cylinder || settings.shape == Shape::Barrel)
	{
		CSPrimCylinder *cylinder = nullptr;
		if (settings.shape == Shape::Barrel)
		{
			auto *shell = new CSPrimCylindricalShell(csx.GetParameterSet(), metal);
			shell->SetShellWidth(.3 * .001 / settings.unit);
			cylinder = shell;
		}
		else
		{
			cylinder = new CSPrimCylinder(csx.GetParameterSet(), metal);
		}
		cylinder->SetPriority(10);
		cylinder->SetRadius((settings.shape == Shape::Barrel ? .55 : .7) * .001 / settings.unit);
		for (unsigned int n = 0; n < 3; ++n)
		{
			cylinder->SetCoord(2 * n, draw(n == axis ? bottom : 2));
			cylinder->SetCoord(2 * n + 1, draw(n == axis ? top : 2));
		}
	}
	else
	{
		polygon = new CSPrimLinPoly(csx.GetParameterSet(), metal);
		polygon->SetPriority(10);
		polygon->SetNormDir(int(axis));
		polygon->SetElevation(draw(settings.shape == Shape::Cut ? bottom + .003 : settings.reverse ? top : bottom));
		polygon->SetLength(settings.shape == Shape::Sheet ? 0 :
			(settings.reverse ? -1 : 1) * settings.thickness / settings.unit);
		points = Contour(settings.shape);
		if (settings.closed)
		{
			points.push_back(points.front());
		}
		if (settings.reverse)
		{
			std::reverse(points.begin(), points.end());
		}
		for (const auto &point : points)
		{
			polygon->AddCoord(draw(point[0]));
			polygon->AddCoord(draw(point[1]));
		}
	}
	if (settings.bulk)
	{
		auto *substrate = new CSPropMaterial(csx.GetParameterSet());
		substrate->SetEpsilon(3.556);
		substrate->SetKappa(.02);
		csx.AddProperty(substrate);
		addBox(substrate, 0, 0, 4, 4, .5, bottom, 0);
	}
	Require(csx.Update().empty(), "CSXCAD rejected shaped fixture before native preparation");
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::array<double, 6> region;
	std::array<unsigned int, 3> ratio{4, 4, 4};
	for (unsigned int n = 0; n < 3; ++n)
	{
		region[n] = draw(n == axis ? .5 : 1);
		region[n + 3] = draw(n == axis ? 1.5 : 3);
	}
	ratio[axis] = 1;
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New(std::vector<double>(region.begin(), region.end()),
		ratio));
	op->SetEngineType(target ? Engine::GPU : Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "shaped-conductor geometry rejected");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(bc);
	Require(op->CalcECOperator() == 0, "shaped-conductor preparation failed");
	if (polygon)
	{
		Require(polygon->GetQtyCoords() == points.size(), "native contour vertex count changed");
		for (std::size_t i = 0; i < points.size(); ++i)
		{
			Require(polygon->GetCoord(int(2 * i)) == draw(points[i][0]) &&
				polygon->GetCoord(int(2 * i + 1)) == draw(points[i][1]), "native contour was rounded or replaced");
		}
	}
	auto &grid = op->Grid();
	Result result;
	result.faces = grid.OpaqueBoundaryFaces().size();
	std::size_t count = 0, bore = 0;
	for (std::size_t leaf = 0; leaf < grid.CellCount(); ++leaf)
	{
		const auto lo = grid.LeafStart(leaf), length = grid.LeafLengths(leaf);
		const double x = (lo[u] + .5 * length[u]) / .001;
		const double y = (lo[v] + .5 * length[v]) / .001;
		const double z = (lo[axis] + .5 * length[axis]) / .001;
		const bool expected = z > bottom && z < top && Inside(settings.shape, x, y);
		Require(grid.LeafIsOpaque(leaf) == expected, "native shaped occupancy differs from analytic centre classification");
		result.occupied.push_back(expected);
		if (expected)
		{
			++count;
			result.volume += length[0] * length[1] * length[2];
		}
		if (settings.shape == Shape::Barrel && z > bottom && z < top && std::hypot(x - 2, y - 2) < .4)
		{
			Require(!grid.LeafIsOpaque(leaf), "barrel bore was filled with conductor");
			++bore;
		}
	}
	Require(count && result.faces, "shaped conductor lost its occupied cells or exposed faces");
	Require(settings.shape != Shape::Barrel || bore > 0, "barrel fixture did not resolve its bore");
	grid.SetMagnetic([&](unsigned int n, const Grid::Point &p)
	{
		return n == u ? .002 * std::exp(-std::pow((p[axis] - .0011) / .0004, 2)) : 0;
	});
	if (target)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
		Require(bool(backend), "shaped-conductor backend unavailable");
		target->Verify(*backend);
		op->SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> engine(op->CreateEngine());
	auto *surface = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(0));
	Require(surface && engine->GetExtensionCount() == 1, "shaped conductor lost its single paired provider");
	const double dt = op->GetTimestep(), initial = grid.Energy(dt);
	Require(initial > 0, "shaped fixture has no incident energy");
	for (unsigned int step = 0; step < 256; ++step)
	{
		engine->IterateTS(1);
		result.balance = std::max(result.balance, std::abs(grid.Energy(dt) + surface->AdditionalStoredEnergy() +
			surface->ElectricLoss() + surface->MagneticLoss() + surface->SurfaceLoss() - initial) / initial);
		result.loss = surface->SurfaceLoss() / initial;
		if (step % 32 == 0)
		{
			for (std::size_t leaf = 0; leaf < grid.CellCount(); ++leaf)
			{
				const auto lo = grid.LeafStart(leaf), length = grid.LeafLengths(leaf);
				auto point = lo;
				for (unsigned int n = 0; n < 3; ++n)
				{
					point[n] += .5 * length[n];
				}
				for (unsigned int n = 0; n < 3; ++n)
				{
					result.fields.push_back(grid.Electric(n, point));
					result.fields.push_back(grid.Magnetic(n, point) * 376.730313668);
				}
			}
		}
	}
	Require(result.balance < (target ? 2e-4 : 2e-9), "shaped conductor energy ledger failed");
	Require(result.loss > 1e-8, "shaped conductor has no measurable Ohmic loss");
	std::cout << "SHAPED_SURFACE shape=" << int(settings.shape) << " boxes=" << boxes << " device=" << bool(target)
		<< " axis=" << axis << " closed=" << settings.closed << " reverse=" << settings.reverse
		<< " cells=" << count << " faces=" << result.faces << " volume=" << result.volume
		<< " loss=" << result.loss << " balance=" << result.balance << std::endl;
	return result;
}
static double Error(const Result &a, const Result &b)
{
	Require(a.occupied == b.occupied && a.faces == b.faces && a.fields.size() == b.fields.size(),
		"shaped/decomposed topology or observation count differs");
	double peak = 0, delta = 0;
	for (std::size_t i = 0; i < a.fields.size(); ++i)
	{
		Require(std::isfinite(b.fields[i]), "nonfinite shaped-conductor field");
		peak = std::max(peak, std::abs(a.fields[i]));
		delta = std::max(delta, std::abs(a.fields[i] - b.fields[i]));
	}
	return delta / std::max(peak, 1e-300);
}
int main(int argc, char **argv)
{
	if (!SubgridHasNativeSurfaceMetadata<CSPropLossyMetal>::value)
	{
		std::cout << "SKIP: CSXCAD lacks volume LossyMetal metadata" << std::endl;
		return 77;
	}
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "expected optional --metal or --cuda");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromFlag(argv[1]) : SubgridDeviceTarget{};
		for (const Settings settings : {Settings{}, Settings{Shape::Curve}, Settings{Shape::Cylinder},
			Settings{Shape::Barrel, 2, .00027178}, Settings{Shape::Curve, 2, 35e-6, .001, 0, true, true},
			Settings{Shape::Concave, 0}, Settings{Shape::Curve, 1, 18e-6, 1e-6, 10},
			Settings{Shape::Concave, 2, 18e-6, .001, 0, true, true, true}})
		{
			const auto shaped = Run(settings, false, nullptr), decomposed = Run(settings, true, nullptr);
			const auto device = Run(settings, false, &target);
			const double exact = Error(decomposed, shaped), error = Error(shaped, device);
			Require(exact == 0 && decomposed.loss == shaped.loss && decomposed.volume == shaped.volume,
				"shaped conductor differs from independent occupied-box decomposition");
			Require(error < 2e-4, "shaped conductor CPU/device fields disagree");
			std::cout << "SHAPED_COMPARE shape=" << int(settings.shape) << " exact_error=" << exact
				<< " device_error=" << error << std::endl;
		}
		for (Shape shape : {Shape::Bowtie, Shape::Backtrack, Shape::Flat, Shape::Sheet, Shape::Cut, Shape::RepeatedClosing})
		{
			bool rejected = false;
			try
			{
				Run(Settings{shape}, false, nullptr);
			}
			catch (const std::invalid_argument &error)
			{
				rejected = true;
				std::cout << "SHAPED_REJECT shape=" << int(shape) << " message=" << error.what() << std::endl;
			}
			Require(rejected, "invalid shaped conductor was admitted");
		}
		std::cout << "PASS: native shaped conductor occupancy, exact decomposition, fields and energy" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}
#else
#include <iostream>
int main()
{
	std::cout << "SKIP: CSXCAD has no volume LossyMetal property" << std::endl;
	return 77;
}
#endif
