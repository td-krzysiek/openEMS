/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "openems.h"
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPropMetal.h"
#include <iostream>
#include <stdexcept>
#include <string>

#if __has_include("CSPropLossyMetal.h")
#include "CSPropLossyMetal.h"

class GuardSolver : public openEMS
{
  public:
	GuardSolver(unsigned int choice, bool local, bool surface, bool primitive)
	{
		const EngineType choices[] = {EngineType_Basic, EngineType_SSE, EngineType_CPU, EngineType_GPU_Reference};
		if (choice >= 4)
		{
			throw std::invalid_argument("unexpected engine choice");
		}
		m_engine = choices[choice];
		SetNumberOfThreads(1);
		auto *csx = new ContinuousStructure();
		SetCSX(csx);
		csx->GetGrid()->SetDeltaUnit(.001);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= 4; ++i)
			{
				csx->GetGrid()->AddDiscLine(n, i);
			}
		}
		CSProperties *property = nullptr;
		if (surface)
		{
			auto *metal = new CSPropLossyMetal(csx->GetParameterSet());
			metal->SetConductivity(5.8e7);
			property = metal;
		}
		else
		{
			property = new CSPropMetal(csx->GetParameterSet());
		}
		csx->AddProperty(property);
		if (primitive)
		{
			auto *box = new CSPrimBox(csx->GetParameterSet(), property);
			for (unsigned int n = 0; n < 3; ++n)
			{
				box->SetCoord(2 * n, 1.5);
				box->SetCoord(2 * n + 1, 2.5);
			}
		}
		if (local)
		{
			SetupCartesianSubgrid({1, 1, 1, 3, 3, 3});
		}
	}

	bool PrepareOperator()
	{
		return SetupOperator();
	}
};

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

int main()
{
	try
	{
		unsigned int rejected = 0;
		for (unsigned int engine = 0; engine < 4; ++engine)
		{
			for (bool primitive : {false, true})
			{
				GuardSolver ordinary(engine, false, true, primitive);
				bool failed = false;
				try
				{
					ordinary.PrepareOperator();
				}
				catch (const std::invalid_argument &error)
				{
					failed = std::string(error.what()).find("LossyMetal requires the native Cartesian subgrid path") != std::string::npos;
				}
				Require(failed, "ordinary LossyMetal was accepted without its conductor extension");
				++rejected;
			}
		}
		unsigned int controls = 0;
		for (unsigned int engine : {0u, 3u})
		{
			GuardSolver native(engine, true, true, true);
			Require(native.PrepareOperator(), "supported native LossyMetal operator was rejected");
			GuardSolver pec(engine, false, false, true);
			Require(pec.PrepareOperator(), "ordinary legacy PEC operator was rejected");
			controls += 2;
		}
		std::cout << "SURFACE_CAPABILITY rejected=" << rejected << " supported_operator_controls=" << controls
			<< " no_field_or_time_advance=1\n";
		std::cout << "PASS: unsupported ordinary LossyMetal rejects before backend/operator creation\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}
#else
int main()
{
	std::cout << "SKIP: optional CSXCAD LossyMetal property is not available\n";
	return 0;
}
#endif
