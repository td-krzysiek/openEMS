/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_candidate_index.h"
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPrimCylinder.h"
#include "CSPrimLinPoly.h"
#include "CSPropMaterial.h"
#include "CSPropMetal.h"
#include <iostream>
#include <random>
#include <stdexcept>

namespace
{
using Index = CartesianSubgridCandidateIndex<CSPrimitives *>;
using Point = Index::Point;
using Bounds = Index::Bounds;

void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

void AddBox(ContinuousStructure &csx, CSProperties *property, int priority, const Bounds &bounds)
{
	auto *box = new CSPrimBox(csx.GetParameterSet(), property);
	box->SetPriority(priority);
	for (unsigned int n = 0; n < bounds.size(); ++n)
	{
		box->SetCoord(n, bounds[n]);
	}
}
}

int main()
{
	try
	{
		ContinuousStructure csx;
		csx.SetCoordInputType(CARTESIAN);
		csx.SetDrawingTolerance(.0001);
		auto *dielectric = new CSPropMaterial(csx.GetParameterSet());
		auto *metal = new CSPropMetal(csx.GetParameterSet());
		csx.AddProperty(dielectric);
		csx.AddProperty(metal);
		std::vector<CSPrimLinPoly *> polygons;
		for (unsigned int i = 0; i < 64; ++i)
		{
			const double x = 2 * (i % 8), y = 2 * (i / 8);
			AddBox(csx, dielectric, int(i % 7), {x, x + .8, y, y + .8, -.2, .4});
			auto *polygon = new CSPrimLinPoly(csx.GetParameterSet(), metal);
			polygon->SetNormDir(2);
			polygon->SetElevation(.25);
			polygon->SetLength(i % 2 ? .1 : -.1);
			polygon->SetPriority(5);
			for (const Point &point : {Point{0, 0, 0}, {1, 0, 0}, {1, .4, 0},
				{.4, .4, 0}, {.4, 1, 0}, {0, 1, 0}})
			{
				polygon->AddCoord(x + point[0]);
				polygon->AddCoord(y + point[1]);
			}
			polygons.push_back(polygon);
			auto *cylinder = new CSPrimCylinder(csx.GetParameterSet(), metal);
			const Bounds ends{x + 1.4, x + 1.4, y + 1.4, y + 1.4, -.2, .6};
			for (unsigned int n = 0; n < ends.size(); ++n)
			{
				cylinder->SetCoord(n, ends[n]);
			}
			cylinder->SetRadius(.15);
			cylinder->SetPriority(int(i % 11));
		}
		AddBox(csx, dielectric, -100, {-1, 17, -1, 17, -1, 2});
		AddBox(csx, dielectric, 100, {0, .2, 0, .2, .2, .4});
		AddBox(csx, metal, 5, {0, .3, 0, .3, .2, .4});
		AddBox(csx, metal, 10, {2, 4, 2, 4, 0, 0});
		AddBox(csx, metal, 10, {50, 51, 50, 51, 50, 51});
		Require(csx.Update().empty(), "actual CSXCAD fixture did not update");
		const Bounds domain{-1, 17, -1, 17, -1, 2};
		std::size_t queries = 0, inspected = 0;
		for (const auto type : {CSProperties::MATERIAL,
			static_cast<CSProperties::PropertyType>(CSProperties::MATERIAL | CSProperties::METAL)})
		{
			const auto all = csx.GetAllPrimitives(true, type);
			std::vector<Bounds> bounds;
			for (auto *primitive : all)
			{
				Bounds box{};
				const bool exact = primitive->GetBoundBox(box.data());
				Require(!primitive->HasTransform(), "fixture unexpectedly transformed");
				Require(exact || dynamic_cast<CSPrimPolygon *>(primitive), "fixture bounds not conservative");
				bounds.push_back(box);
			}
			Index index(all, bounds, domain, csx.GetDrawingTolerance());
			Require(index.Indexed(), "actual CAD candidates did not enable indexing");
			auto check = [&](const Point &point)
			{
				CSPrimitives *expected = nullptr, *actual = nullptr;
				auto *property = csx.GetPropertyByCoordPriority(point.data(), all, true, &expected);
				const auto &candidates = index.At(point.data());
				auto *selected = csx.GetPropertyByCoordPriority(point.data(), candidates, true, &actual);
				Require(property == selected && expected == actual,
					"actual CSXCAD property or tied-priority primitive winner changed");
				inspected += candidates.size();
				++queries;
			};
			std::mt19937 random(700170);
			std::uniform_real_distribution<double> x(-1, 17), z(-.3, .7);
			for (unsigned int i = 0; i < 20000; ++i)
			{
				check({x(random), x(random), z(random)});
			}
			for (const Bounds &box : bounds)
			{
				for (unsigned int corner = 0; corner < 8; ++corner)
				{
					Point point;
					for (unsigned int n = 0; n < 3; ++n)
					{
						point[n] = box[2 * n + ((corner >> n) & 1)];
					}
					check(point);
					for (unsigned int n = 0; n < 3; ++n)
					{
						auto near = point;
						near[n] += .5 * csx.GetDrawingTolerance();
						check(near);
					}
				}
			}
			for (auto *polygon : polygons)
			{
				// Concave missing corner, signed extrusion and actual native
				// tolerance are tested through CSXCAD rather than an AABB oracle.
				const double px = polygon->GetCoord(0), py = polygon->GetCoord(1);
				for (const double elevation : {.14995, .15, .24995, .25, .25005, .35, .35005})
				{
					check({px + .8, py + .8, elevation});
					check({px + .2, py + .2, elevation});
				}
			}
		}
		std::cout << "NATIVE_CANDIDATES exact_cad_queries=" << queries
			<< " priority_primitive_identity=1 concave_signed_extrusion=1 average_candidates="
			<< double(inspected) / queries << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}
