/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "openems.h"
#include "ContinuousStructure.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/engine.h"
#include "tools/constants.h"
#include "tools/useful.h"
#include "tinyxml.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

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
	bool rejected = false;
	try
	{
		action();
	}
	catch (const std::invalid_argument &)
	{
		rejected = true;
	}
	Require(rejected, "invalid child-coordinate API request was accepted");
}

static std::string Serialized(openEMS &solver)
{
	TiXmlDocument document;
	Require(solver.Write2XML(&document), "child-coordinate serialization failed");
	TiXmlPrinter printer;
	document.Accept(&printer);
	return printer.CStr();
}

static std::string ChildDefinition(const std::string &xml)
{
	TiXmlDocument document;
	document.Parse(xml.c_str());
	Require(!document.Error(), "serialized API fixture is malformed");
	auto *element = document.RootElement()->FirstChildElement("FDTD")->FirstChildElement("CartesianSubgridChildLines");
	Require(element, "serialized child-coordinate element is missing");
	TiXmlPrinter printer;
	element->Accept(&printer);
	return printer.CStr();
}

static void Serialization()
{
	openEMS solver;
	solver.SetCSX(new ContinuousStructure());
	solver.SetGaussExcite(13e9, 9e9);
	solver.SetupCartesianSubgrid({1, 1, 1, 3, 3, 3});
	const std::vector<std::vector<double>> axes{{0, .13, 1, 1.7, 2}, {}, {-1, -.1, 0, .7, 1}};
	solver.SetupCartesianSubgridChildLines(axes);
	const std::string saved = Serialized(solver);
	TiXmlDocument parsed;
	parsed.Parse(saved.c_str());
	Require(!parsed.Error(), "generated child-coordinate XML is malformed");
	openEMS restored;
	restored.SetCSX(new ContinuousStructure());
	Require(restored.Parse_XML_FDTDSetup(parsed.RootElement()->FirstChildElement("FDTD")),
		"child-coordinate XML parsing failed");
	Require(ChildDefinition(Serialized(restored)) == ChildDefinition(saved),
		"child-coordinate XML does not roundtrip exactly");
	Reject([&] { solver.SetupCartesianSubgridChildLines({{0, 1}, {}}); });
	Reject([&] { solver.SetupCartesianSubgridChildLines({{0}, {}, {}}); });
	Reject([&] { solver.SetupCartesianSubgridChildLines({{0, 1, 1}, {}, {}}); });
	Reject([&] { solver.SetupCartesianSubgridChildLines({{}, {0, std::numeric_limits<double>::infinity()}, {}}); });
	Require(Serialized(solver) == saved, "invalid child-coordinate request changed stored configuration");
	for (const std::string &contents : {
		"<CartesianSubgridChildLines/>",
		"<CartesianSubgridChildLines Q=\"0,1\"/>",
		"<CartesianSubgridChildLines X=\"0,,1\"/>",
		"<CartesianSubgridChildLines X=\"0,1,\"/>",
		"<CartesianSubgridChildLines X=\"0,1suffix\"/>",
		"<CartesianSubgridChildLines X=\"0,1\"><Unexpected/></CartesianSubgridChildLines>",
		"<CartesianSubgridChildLines X=\"0,1\"/><CartesianSubgridChildLines X=\"0,1\"/>"})
	{
		TiXmlDocument bad;
		bad.Parse(("<FDTD>" + contents + "</FDTD>").c_str());
		Require(!bad.Error(), "negative API fixture is not valid XML");
		Reject([&] { openEMS rejected; rejected.Parse_XML_FDTDSetup(bad.RootElement()); });
	}
	solver.SetupCartesianSubgridChildLines({{}, {}, {}});
	Require(Serialized(solver).find("CartesianSubgridChildLines") == std::string::npos,
		"clearing child coordinates leaves the XML element");
	std::cout << "CHILD_API roundtrip=1 invalid_requests=11 clear=1" << std::endl;
}

static void Native(double unit)
{
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(unit);
	Grid::Axes base, child, physicalBase, physicalChild;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int k = 0; k <= 6; ++k)
		{
			const double value = 1.7 * n + (.001 / unit) * (k + .04 * k * k);
			base[n].push_back(value);
			csx.GetGrid()->AddDiscLine(n, value);
		}
		for (std::size_t k = 0; k + 1 < base[n].size(); ++k)
		{
			child[n].push_back(base[n][k]);
			child[n].push_back(base[n][k] + (.27 + .07 * n) * (base[n][k + 1] - base[n][k]));
		}
		child[n].push_back(base[n].back());
		for (double value : base[n])
		{
			physicalBase[n].push_back((value - base[n][0]) * unit);
		}
		for (double value : child[n])
		{
			physicalChild[n].push_back((value - base[n][0]) * unit);
		}
	}
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New(
		{base[0][2], base[1][2], base[2][2], base[0][4], base[1][4], base[2][4]}));
	op->SetChildLines(child);
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "native child-coordinate geometry rejected");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(bc);
	Require(op->CalcECOperator() == 0, "native child-coordinate preparation failed");
	Grid manual(physicalBase, EPS0, MUE0, {{{{2, 2, 2}, {4, 4, 4}}, {2, 2, 2}}},
		false, {}, {}, {}, {}, Grid::CurlLayout::Patterns, physicalChild);
	Require(op->Grid().StableTimestep() == manual.StableTimestep(), "native coordinate conversion changes timestep");
	Require(op->Grid().CellCount() == manual.CellCount(), "native coordinate conversion changes leaf count");
	for (std::size_t leaf = 0; leaf < manual.CellCount(); ++leaf)
	{
		Require(op->Grid().LeafStart(leaf) == manual.LeafStart(leaf) &&
			op->Grid().LeafLengths(leaf) == manual.LeafLengths(leaf), "native conversion changes exact geometry");
	}
	const auto initial = [](unsigned int n, const Grid::Point &p)
	{
		return std::sin(731 * p[(n + 1) % 3]) + .2 * std::cos(451 * p[(n + 2) % 3]);
	};
	op->Grid().SetElectric(initial);
	manual.SetElectric(initial);
	std::unique_ptr<Engine> engine(op->CreateEngine());
	Require(engine->IterateTS(256), "native child-coordinate stepping failed");
	for (unsigned int k = 0; k < 256; ++k)
	{
		manual.Step(op->GetTimestep());
	}
	const auto native = op->Grid().DeviceStateData(), expected = manual.DeviceStateData();
	Require(native.electric == expected.electric && native.magnetic == expected.magnetic,
		"native adapter changes exact manual-core fields");
	bool late = false;
	try
	{
		op->SetChildLines({});
	}
	catch (const std::logic_error &)
	{
		late = true;
	}
	Require(late, "prepared operator accepts coordinate mutation");
	std::cout << "CHILD_API unit=" << unit << " exact_geometry=1 exact_fields=1 late_change_rejected=1" << std::endl;
}

int main()
{
	try
	{
		SetThreadLimit(1);
		Serialization();
		Native(.001);
		Native(1.0);
		std::cout << "PASS: child-coordinate XML/API, SI drawing conversion and native/manual fields" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}
