// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Krzysztof Rosinski

#include "FDTD/operator.h"
#include "FDTD/engine.h"
#include "FDTD/excitation.h"
#include "FDTD/extensions/operator_extension.h"
#include "FDTD/extensions/operator_ext_upml.h"
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPrimLinPoly.h"
#include "CSPropExcitation.h"
#include "CSPropMaterial.h"
#include "CSPropMetal.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
void Require(bool condition, const char* message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

void Environment(const char* value)
{
#ifdef _WIN32
	_putenv_s("OPENEMS_GEOMETRY_INDEX", value ? value : "");
#else
	if (value)
	{
		setenv("OPENEMS_GEOMETRY_INDEX", value, 1);
	}
	else
	{
		unsetenv("OPENEMS_GEOMETRY_INDEX");
	}
#endif
}

struct RestoreEnvironment
{
	bool existed;
	std::string value;
	RestoreEnvironment()
		: existed(std::getenv("OPENEMS_GEOMETRY_INDEX") != nullptr),
		  value(existed ? std::getenv("OPENEMS_GEOMETRY_INDEX") : "")
	{
	}
	~RestoreEnvironment()
	{
		Environment(existed ? value.c_str() : nullptr);
	}
};

class TestOperator : public Operator
{
public:
	TestOperator()
	{
		Init();
		SetNumThreads(2);
	}
	bool HasIndex() const
	{
		return static_cast<bool>(m_GeometryLookup);
	}
	bool fail = false;
protected:
	bool Calc_EC() override
	{
		if (fail)
		{
			throw std::runtime_error("intentional setup failure");
		}
		return Operator::Calc_EC();
	}
};

void AddBox(ContinuousStructure& csx, CSProperties* property, int priority,
	const std::array<double, 6>& bounds)
{
	auto* box = new CSPrimBox(csx.GetParameterSet(), property);
	box->SetPriority(priority);
	for (int n=0; n<6; ++n)
	{
		box->SetCoord(n, bounds[n]);
	}
}

struct Result
{
	double timestep;
	std::vector<float> coefficients;
	std::vector<float> fields;
};

class LookupExtension : public Operator_Extension
{
public:
	LookupExtension(TestOperator& op, ContinuousStructure& csx, bool& seen,
		bool indexed, bool fail)
		: Operator_Extension(&op), m_Operator(op), m_CSX(csx), m_Seen(seen),
		  m_Indexed(indexed), m_Fail(fail)
	{
	}
	bool BuildExtension() override
	{
		m_Seen = true;
		Require(m_Operator.HasIndex() == m_Indexed, "index not retained through extension construction");
		const double bounds[6] = {3.5, 4.5, 3.5, 4.5, 0, 8};
		const double point[3] = {4, 4, 2.25};
		for (auto type : {CSProperties::MATERIAL, CSProperties::EXCITATION,
			static_cast<CSProperties::PropertyType>(CSProperties::MATERIAL | CSProperties::METAL)})
		{
			const auto actual = m_Operator.GetPrimitivesBoundBox(8, 8, -1, type);
			const auto legacy = m_CSX.GetPrimitivesByBoundBox(bounds, true, type);
			Require(m_CSX.GetPropertyByCoordPriority(point, actual, false) ==
				m_CSX.GetPropertyByCoordPriority(point, legacy, false), "extension material winner changed");
		}
		if (m_Fail)
		{
			throw std::runtime_error("intentional extension failure");
		}
		return false; // Test hook only; no engine extension.
	}
private:
	TestOperator& m_Operator;
	ContinuousStructure& m_CSX;
	bool& m_Seen;
	bool m_Indexed;
	bool m_Fail;
};

Result Run(const char* option, int fail=0)
{
	Environment(option);
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(.001);
	for (int axis=0; axis<3; ++axis)
	{
		for (int i=0; i<=16; ++i)
		{
			csx.GetGrid()->AddDiscLine(axis, i*.5);
		}
	}
	auto* dielectric = new CSPropMaterial(csx.GetParameterSet());
	dielectric->SetEpsilon(3.5);
	dielectric->SetKappa(.002);
	auto* metal = new CSPropMetal(csx.GetParameterSet());
	auto* voidMaterial = new CSPropMaterial(csx.GetParameterSet());
	csx.AddProperty(dielectric);
	csx.AddProperty(metal);
	csx.AddProperty(voidMaterial);
	AddBox(csx, dielectric, 0, {{.75, 7.25, .75, 7.25, .75, 4.75}});
	AddBox(csx, metal, 5, {{1, 7, 1, 7, 1, 1}});
	AddBox(csx, voidMaterial, 10, {{3, 4, 3, 4, .75, 1.25}});
	AddBox(csx, dielectric, 5, {{2, 5, 2, 5, 3, 6}});
	auto* polygon = new CSPrimLinPoly(csx.GetParameterSet(), metal);
	polygon->SetNormDir(2);
	polygon->SetElevation(2.1);
	polygon->SetLength(.25);
	polygon->SetPriority(6);
	for (const auto& point : {std::array<double, 2>{{2, 2}}, {{6, 2}}, {{6, 3}},
		{{3, 3}}, {{3, 6}}, {{2, 6}}})
	{
		polygon->AddCoord(point[0]);
		polygon->AddCoord(point[1]);
	}
	auto* source = new CSPropExcitation(csx.GetParameterSet());
	csx.AddProperty(source);
	source->SetExcitType(0);
	AddBox(csx, source, 0, {{4, 4, 4, 4, 2, 2.5}});
	Require(csx.Update().empty(), "fixture must update");
	Excitation excitation;
	excitation.SetupGaussianPulse(1e9, 1e9);
	TestOperator op;
	op.SetEngineType(Engine::BASIC);
	op.SetExcitationSignal(&excitation);
	Require(op.SetGeometryCSX(&csx), "fixture grid must load");
	op.fail = fail == 1;
	bool extensionSeen = false;
	op.AddExtension(new LookupExtension(op, csx, extensionSeen,
		option && std::string(option) != "0", fail == 2));
	const int boundaries[6] = {3, 3, 3, 3, 3, 3};
	const unsigned int widths[6] = {2, 2, 2, 2, 2, 2};
	Require(Operator_Ext_UPML::Create_UPML(&op, boundaries, widths, ""), "UPML fixture creation failed");
	try
	{
		Require(op.CalcECOperator() == 0, "setup must succeed");
		Require(!fail, "intentional setup failure was not reached");
	}
	catch (const std::runtime_error& error)
	{
		Require(!op.HasIndex(), "setup exception retained borrowed index");
		if ((fail == 1 && std::string(error.what()) == "intentional setup failure") ||
			(fail == 2 && std::string(error.what()) == "intentional extension failure"))
		{
			return {};
		}
		throw;
	}
	Require(!op.HasIndex(), "index must be released before FDTD");
	Require(extensionSeen, "extension lookup check never ran");
	// Subsequent lookups must see geometry edits rather than stale tile lists.
	AddBox(csx, metal, 100, {{0, 8, 0, 8, 0, 8}});
	Require(csx.Update().empty(), "post-setup geometry must update");
	double bounds[6] = {3.5, 4.5, 3.5, 4.5, 0, 8};
	Require(op.GetPrimitivesBoundBox(8, 8, -1, CSProperties::METAL) ==
		csx.GetPrimitivesByBoundBox(bounds, true, CSProperties::METAL),
		"post-setup query is stale");
	Result result;
	result.timestep = op.GetTimestep();
	for (unsigned int n=0; n<3; ++n)
	{
		for (unsigned int x=0; x<17; ++x)
		{
			for (unsigned int y=0; y<17; ++y)
			{
				for (unsigned int z=0; z<17; ++z)
				{
					result.coefficients.insert(result.coefficients.end(), {
						op.GetVV(n,x,y,z), op.GetVI(n,x,y,z),
						op.GetII(n,x,y,z), op.GetIV(n,x,y,z)});
				}
			}
		}
	}
	std::unique_ptr<Engine> engine(op.CreateEngine());
	engine->SetVolt(0, 8, 8, 8, 1);
	Require(engine->IterateTS(48), "short field control must run");
	for (unsigned int n=0; n<3; ++n)
	{
		for (unsigned int x=0; x<17; ++x)
		{
			for (unsigned int y=0; y<17; ++y)
			{
				for (unsigned int z=0; z<17; ++z)
				{
					result.fields.push_back(engine->GetVolt(n,x,y,z));
					result.fields.push_back(engine->GetCurr(n,x,y,z));
				}
			}
		}
	}
	return result;
}

void Compare(const Result& a, const Result& b)
{
	Require(std::abs(a.timestep-b.timestep) <= 1e-15*a.timestep, "timestep differs");
	Require(a.coefficients == b.coefficients, "coefficient selection changed");
	Require(a.fields.size() == b.fields.size(), "field size differs");
	float error = 0, magnitude = 0;
	for (std::size_t i=0; i<a.fields.size(); ++i)
	{
		Require(std::isfinite(a.fields[i]) && std::isfinite(b.fields[i]), "nonfinite field");
		error = std::max(error, std::abs(a.fields[i]-b.fields[i]));
		magnitude = std::max(magnitude, std::abs(a.fields[i]));
	}
	Require(magnitude > 0, "field control must be nontrivial");
	Require(error <= 1e-7*magnitude, "field control exceeds relative 1e-7 tolerance");
	std::cout << "Coefficient arrays identical; field max relative error " << error/magnitude << '\n';
}
}

int main()
{
	RestoreEnvironment restore;
	try
	{
		const auto legacy = Run(nullptr);
		Compare(legacy, Run("0"));
#if __has_include("CSPrimitiveIndex.h")
		Compare(legacy, Run("8,8,1"));
		Compare(legacy, Run("4,4,4"));
		Run("4,4,1", true);
		Run("4,4,1", 2);
#endif
		bool rejected = false;
		try
		{
			Run("not-a-grid");
		}
		catch (const std::exception&)
		{
			rejected = true;
		}
		Require(rejected, "invalid/unsupported index configuration must fail");
	}
	catch (const std::exception& error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
	std::cout << "Geometry index integration tests passed\n";
}
