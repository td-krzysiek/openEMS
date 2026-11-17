/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_material_quadrature.h"
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/engine.h"
#include "CSPrimLinPoly.h"
#include "CSPrimCylinder.h"
#include "CSPropMaterial.h"
#include "ContinuousStructure.h"
#include "tools/constants.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "subgrid_device_target.h"
#include "subgrid_fp32_oracle.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace Q = CartesianSubgridMaterialQuadrature;

static void Require(bool ok, const char *message)
{
	if (!ok)
	{
		throw std::runtime_error(message);
	}
}

static void Quadrature()
{
	const Q::Values background{1, 1, 0, 0}, slab{3.556, 1.7, .01, .02};
	for (unsigned int count : {0, 2, 4, 8})
	{
		unsigned int calls = 0;
		const auto value = Q::Average({0, 0, 0}, {1, 1, 1}, count, [&](const Q::Point &)
		{
			++calls;
			return slab;
		});
		Require(value == slab && calls == (count ? count * count * count : 1),
			"constant leaf identity or explicit quadrature cardinality changed");
	}
	double lastTriangle = 1, lastBore = 1;
	for (unsigned int count : {2, 4, 8})
	{
		// Independent analytic diagonal half-volume; its interface goes through
		// quadrature nodes, exposing a finite sampling error rather than hiding it.
		const auto triangle = Q::Average({0, 0, 0}, {1, 1, 1}, count, [&](const Q::Point &p)
		{
			return p[0] + p[1] <= 1 ? slab : background;
		});
		const double fraction = (triangle[0] - 1) / (slab[0] - 1);
		const double error = std::abs(fraction - .5);
		Require(error < lastTriangle && error <= .5 / count + 1e-14,
			"diagonal material sampling failed its analytic convergence bound");
		lastTriangle = error;
		// A higher-priority cylindrical air bore cuts the dielectric. Resolve
		// the bore at every point, never average primitive contributions first.
		const auto bore = Q::Average({0, 0, 0}, {1, 1, 1}, count, [&](const Q::Point &p)
		{
			return std::hypot(p[0] - .5, p[1] - .5) < .35 ? background : slab;
		});
		const double air = (slab[0] - bore[0]) / (slab[0] - 1);
		const double boreError = std::abs(air - std::acos(-1.) * .35 * .35);
		// The 2x2 nodes miss this particular circle; subsequent levels improve.
		Require(boreError <= lastBore, "priority bore volume does not converge");
		lastBore = boreError;
		Require(std::abs(bore[2] - .01 * (1 - air)) < 1e-15 &&
			std::abs(bore[3] - .02 * (1 - air)) < 1e-15,
			"conductivity and permittivity used different native priority fractions");
		std::cout << "MATERIAL_QUADRATURE n=" << count << " diagonal_error=" << error
			<< " priority_bore_error=" << boreError << '\n';
	}
	Require(lastBore < .011 && lastTriangle <= .062500000001,
		"finest analytic volume approximation exceeded its stated bound");
	unsigned int guards = 0;
	for (unsigned int bad : {1, 3, 16, 1024})
	{
		try
		{
			Q::Average({0, 0, 0}, {1, 1, 1}, bad, [&](const Q::Point &) { return slab; });
		}
		catch (const std::invalid_argument &)
		{
			++guards;
		}
	}
	Require(guards == 4, "invalid quadrature counts were accepted");
	for (unsigned int component = 0; component < 4; ++component)
	{
		for (double bad : {-1., std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
		{
			bool rejected = false;
			try
			{
				Q::Average({0, 0, 0}, {1, 1, 1}, 2, [&](const Q::Point &)
				{
					auto value = slab;
					value[component] = bad;
					return value;
				});
			}
			catch (const std::invalid_argument &)
			{
				rejected = true;
			}
			Require(rejected, "invalid material sample entered the positive mass/loss coefficients");
		}
	}
}

static void PassiveCore()
{
	using Grid = CartesianSubgrid;
	Grid::Axes axes;
	for (auto &axis : axes)
	{
		for (unsigned int n = 0; n <= 6; ++n)
		{
			axis.push_back(.001 * n);
		}
	}
	auto sample = [](const Grid::Point &lo, const Grid::Point &hi)
	{
		const auto average = Q::Average(lo, hi, 4, [](const Q::Point &p)
		{
			const bool inside = std::hypot(p[0] - .003, p[1] - .003) < .00135;
			return Q::Values{EPS0 * (inside ? 3.556 : 1), MUE0, 0, 0};
		});
		return Grid::Medium{average[0], average[1]};
	};
	Grid grid(axes, EPS0, MUE0,
		std::vector<Grid::Region>{{{{2, 2, 2}, {4, 4, 4}}, {2, 2, 2}}}, false, sample);
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return std::sin(517 * p[0] + 943 * p[1] - 617 * p[2] + n);
	});
	const double dt = .7 * grid.StableTimestep(), initial = grid.Energy(dt);
	double drift = 0;
	for (unsigned int k = 0; k < 2000; ++k)
	{
		grid.Step(dt);
		if (k % 50 == 0)
		{
			drift = std::max(drift, std::abs(grid.Energy(dt) / initial - 1));
		}
	}
	Require(drift < 1e-9, "sampled positive material broke the lossless energy invariant");
	std::cout << "SAMPLED_MATERIAL_ENERGY steps=2000 relative_drift=" << drift << '\n';
}

static void DeviceFields(CartesianSubgrid &native, const SubgridDeviceTarget &target, bool domains)
{
	using Grid = CartesianSubgrid;
	auto cpu = native;
	cpu.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return std::sin(517 * p[0] + 943 * p[1] - 617 * p[2] + n);
	});
	cpu.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return .001 * std::cos(733 * p[0] - 113 * p[1] + 631 * p[2] + n);
	});
	auto device = cpu, oracle = cpu;
	std::vector<std::unique_ptr<GPU_Backend>> owners;
	std::vector<GPU_Backend *> backends;
	for (std::size_t k = 0; k < (domains ? 3 : 1); ++k)
	{
		target.SelectDomain(k);
		owners.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(owners.back()), "shaped-material device backend unavailable");
		target.Verify(*owners.back());
		std::cout << "MATERIAL_DEVICE owner=" << k << " backend=" << owners.back()->GetName() << std::endl;
		backends.push_back(owners.back().get());
	}
	std::unique_ptr<GPU_CartesianSubgrid> execution;
	if (backends.size() == 1)
	{
		execution = backends[0]->CreateCartesianSubgrid(device);
	}
	else
	{
		const auto &lines = device.BaseLines()[0];
		auto plan = std::make_shared<CartesianSubgridDomainPlan>(device, 0,
			std::vector<double>{lines[2], lines[4]});
		execution = CreateCartesianSubgridDomains(plan, backends);
	}
	auto reference = CreateCartesianSubgridReference(oracle);
	Require(bool(execution) && bool(reference), "shaped-material device execution missing");
	const double dt = .7 * cpu.StableTimestep();
	for (unsigned int k = 0; k < 256; ++k)
	{
		cpu.Step(dt);
		execution->Step(dt, 0);
		reference->Step(dt, 0);
	}
	execution->ReadFields(device);
	reference->ReadFields(oracle);
	const auto a = cpu.DeviceStateData(), b = device.DeviceStateData(), c = oracle.DeviceStateData();
	const auto error = [](const std::vector<float> &x, const std::vector<float> &y)
	{
		Require(x.size() == y.size(), "shaped-material field cardinality changed");
		double delta = 0, scale = 0;
		for (std::size_t k = 0; k < x.size(); ++k)
		{
			Require(std::isfinite(x[k]) && std::isfinite(y[k]), "nonfinite shaped-material field");
			delta = std::max(delta, std::abs(double(x[k]) - y[k]));
			scale = std::max(scale, std::abs(double(x[k])));
		}
		return delta / std::max(scale, 1e-300);
	};
	const double doubleError = std::max(error(a.electric, b.electric), error(a.magnetic, b.magnetic));
	const double floatError = std::max(error(c.electric, b.electric), error(c.magnetic, b.magnetic));
	Require(doubleError < 2e-4 && floatError < 2e-6,
		"shaped-material device fields exceeded retained double/FP32 bounds");
	std::cout << "NATIVE_SHAPED_DEVICE steps=256 workers=" << backends.size()
		<< " double_error=" << doubleError << " fp32_error=" << floatError << std::endl;
}

static void NativeShape(unsigned int axis, double unit, bool reverse, const SubgridDeviceTarget &target, bool domains)
{
	using Grid = CartesianSubgrid;
	ContinuousStructure csx;
	const Grid::Point origin{31, -23, 17};
	const unsigned int a = (axis + 1) % 3, b = (axis + 2) % 3;
	const auto draw = [&](unsigned int n, double value) { return (origin[n] + value) * .001 / unit; };
	csx.GetGrid()->SetDeltaUnit(unit);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int k = 0; k <= 6; ++k)
		{
			csx.GetGrid()->AddDiscLine(n, draw(n, k));
		}
	}
	const auto add = [&](bool air)
	{
		auto *material = new CSPropMaterial(csx.GetParameterSet());
		material->SetEpsilon(air ? 1 : 3.556);
		material->SetMue(air ? 1 : 1.7);
		csx.AddProperty(material);
		if (air)
		{
			auto *cylinder = new CSPrimCylinder(csx.GetParameterSet(), material);
			cylinder->SetPriority(9);
			cylinder->SetRadius(.6 * .001 / unit);
			for (unsigned int n = 0; n < 3; ++n)
			{
				cylinder->SetCoord(2 * n, draw(n, n == axis ? 1 : 3));
				cylinder->SetCoord(2 * n + 1, draw(n, n == axis ? 5 : 3));
			}
		}
		else
		{
			auto *polygon = new CSPrimLinPoly(csx.GetParameterSet(), material);
			polygon->SetPriority(1);
			polygon->SetNormDir(int(axis));
			polygon->SetElevation(draw(axis, reverse ? 5 : 1));
			polygon->SetLength((reverse ? -4 : 4) * .001 / unit);
			const double points[8][2] = {{1.6, 1}, {4.4, 1}, {5, 1.6}, {5, 4.4},
				{4.4, 5}, {1.6, 5}, {1, 4.4}, {1, 1.6}};
			for (unsigned int k = 0; k < 8; ++k)
			{
				const auto &point = points[reverse ? 7 - k : k];
				polygon->AddCoord(draw(a, point[0]));
				polygon->AddCoord(draw(b, point[1]));
			}
		}
	};
	// Change property insertion order while retaining explicit CAD priority.
	add(reverse);
	add(!reverse);
	Require(csx.Update().empty(), "native shaped-material CSX input is invalid");
	std::vector<double> region;
	for (double endpoint : {2., 4.})
	{
		for (unsigned int n = 0; n < 3; ++n)
		{
			region.push_back(draw(n, endpoint));
		}
	}
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> strict(Operator_CartesianSubgrid::New(region));
	strict->SetEngineType(Engine::BASIC);
	strict->SetExcitationSignal(&signal);
	bool rejected = false;
	try
	{
		strict->SetGeometryCSX(&csx);
	}
	catch (const std::invalid_argument &error)
	{
		rejected = std::string(error.what()).find("axis-aligned volume boxes") != std::string::npos;
	}
	Require(rejected, "strict default silently admitted shaped material");
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New(region));
	op->SetMaxLeaves(4000000);
	op->SetMaterialSamples(4);
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx) && op->CalcECOperator() == 0,
		"explicit native shaped-material preparation failed");
	unsigned int partial = 0, air = 0;
	double volume = 0;
	for (std::size_t leaf = 0; leaf < op->Grid().CellCount(); ++leaf)
	{
		const auto lo = op->Grid().LeafStart(leaf), length = op->Grid().LeafLengths(leaf);
		Grid::Point hi;
		for (unsigned int n = 0; n < 3; ++n)
		{
			hi[n] = lo[n] + length[n];
		}
		// Independent analytic octagon/bore membership. This does not call
		// CSXCAD or the native material sampler to construct the expected data.
		const auto expected = Q::Average(lo, hi, 4, [&](const Q::Point &p)
		{
			const double u = p[a] * 1000, v = p[b] * 1000, z = p[axis] * 1000;
			const bool slab = z > 1 && z < 5 && u > 1 && u < 5 && v > 1 && v < 5 &&
				u + v > 2.6 && u + v < 9.4 && u - v < 3.4 && v - u < 3.4;
			const bool bore = z > 1 && z < 5 && std::hypot(u - 3, v - 3) < .6;
			return Q::Values{EPS0 * (slab && !bore ? 3.556 : 1),
				MUE0 * (slab && !bore ? 1.7 : 1), 0, 0};
		});
		const auto actual = op->Grid().LeafMedium(leaf);
		Require(std::abs(actual.epsilon / expected[0] - 1) < 2e-13 &&
			std::abs(actual.permeability / expected[1] - 1) < 2e-13,
			"native material priority, axis, units or finite-leaf fraction differs from analytic CAD");
		const double fraction = (actual.epsilon / EPS0 - 1) / 2.556;
		partial += fraction > 1e-8 && fraction < 1 - 1e-8;
		air += fraction < 1e-8;
		volume += fraction * length[0] * length[1] * length[2] * 1e9;
	}
	const double exact = (16 - 4 * .6 * .6 / 2 - std::acos(-1.) * .6 * .6) * 4;
	Require(partial > 0 && air > 0 && std::abs(volume / exact - 1) < .02,
		"native shape did not retain a resolved partial boundary and priority bore");
	DeviceFields(op->Grid(), target, domains);
	std::cout << "NATIVE_SHAPED_MATERIAL axis=" << axis << " unit=" << unit << " reversed=" << reverse
		<< " partial=" << partial << " volume_relative_error=" << std::abs(volume / exact - 1) << '\n';
}

int main(int argc, char **argv)
{
	try
	{
		Require(argc <= 2, "unexpected shaped-material device arguments");
		const bool referenceDomains = argc == 2 && std::string(argv[1]) == "--reference-domains";
		const auto target = argc == 1 || referenceDomains ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromDomainFlag(argv[1]);
		const bool domains = referenceDomains || target.physicalDevices > 1;
		Quadrature();
		PassiveCore();
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			NativeShape(axis, .001, false, target, domains);
			NativeShape(axis, 1e-6, true, target, domains);
		}
		std::cout << "PASS: explicit material quadrature, priority fractions and positive energy\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}
