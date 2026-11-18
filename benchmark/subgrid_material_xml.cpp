/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "CSPropMaterial.h"
#include "FDTD/engine.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "tests/subgrid_device_target.h"
#include "tools/constants.h"
#include "tools/useful.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

using Grid = CartesianSubgrid;

static void Require(bool condition, const std::string &message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

int main(int argc, char **argv)
{
	try
	{
		Require(argc == 6, "expected input.xml, case index, magnetic seed, backend flag and expected rejection");
		SetThreadLimit(1);
		const unsigned int index = unsigned(std::stoul(argv[2]));
		Require(index < 5, "invalid material control index");
		const bool magnetic = std::string(argv[3]) == "1";
		const std::string flag = argv[4], rejection = argv[5];
		Require(flag == "--cpu" || flag == "--reference" || flag == "--metal" || flag == "--cuda", "invalid backend");
		const bool device = flag != "--cpu";
		SubgridDeviceTarget target;
		if (device && flag != "--reference")
		{
			target = SubgridDeviceTarget::FromFlag(flag);
		}
		const std::array<double, 3> eps = index == 1 ? std::array<double, 3>{2, 2, 2} :
			std::array<double, 3>{2, 3, 5};
		const std::array<double, 3> mu = index == 1 ? std::array<double, 3>{1.2, 1.2, 1.2} :
			std::array<double, 3>{1.2, 1.8, 2.5};
		const std::array<double, 3> electric = index == 0 ? std::array<double, 3>{0, .3, .9} :
			index == 2 || index == 3 ? std::array<double, 3>{.4, 0, 0} : std::array<double, 3>{0, 0, 0};
		const std::array<double, 3> magneticLoss = index == 1 || index == 3 ?
			std::array<double, 3>{0, 2500, 7000} : std::array<double, 3>{0, 0, 0};
		ContinuousStructure csx;
		const auto readError = csx.ReadFromXML(argv[1]);
		Require(!readError || !*readError, std::string("native XML read error: ") + (readError ? readError : ""));
		const double unit = index % 2 ? 1e-6 : .001;
		Require(csx.GetGrid()->GetDeltaUnit() == unit, "Python XML changed drawing units");
		const std::array<double, 3> origin{31, -23, 17};
		std::vector<double> region;
		for (double endpoint : {1., 3.})
		{
			for (unsigned int n = 0; n < 3; ++n)
			{
				region.push_back((origin[n] + endpoint) * .001 / unit);
			}
		}
		Excitation signal;
		signal.SetupDiracPulse(30e9);
		std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New(region));
		op->SetEngineType(device ? Engine::GPU : Engine::BASIC);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&csx), "native Python XML geometry rejected");
		int boundary[6];
		std::fill(boundary, boundary + 6, magnetic ? 0 : 1);
		op->SetBoundaryCondition(boundary);
		op->SetTimestepFactor(.4);
		try
		{
			Require(op->CalcECOperator() == 0, "native Python XML operator failed without diagnostic");
		}
		catch (const std::invalid_argument &error)
		{
			Require(rejection != "none" && std::string(error.what()).find(rejection) != std::string::npos,
				std::string("unexpected native rejection: ") + error.what());
			std::cout << "PASS: PYTHON_XML_REJECT case=" << index << " reason=" << rejection << '\n';
			return 0;
		}
		Require(rejection == "none", "invalid Python XML material accepted");
		for (std::size_t leaf = 0; leaf < op->Grid().CellCount(); ++leaf)
		{
			const auto material = op->Grid().LeafMedium(leaf);
			for (unsigned int n = 0; n < 3; ++n)
			{
				Require(material.Epsilon(n) == eps[n] * EPS0 && material.Permeability(n) == mu[n] * MUE0,
					"Python XML material differs from independent SI expectation");
			}
		}
		const auto seed = [](unsigned int n, const Grid::Point &) { return .2 + .3 * n; };
		if (magnetic)
		{
			op->Grid().SetMagnetic(seed);
		}
		else
		{
			op->Grid().SetElectric(seed);
		}
		if (device)
		{
			std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
			Require(bool(backend), "explicit Python XML backend missing");
			target.Verify(*backend);
			op->SetGPUBackend(backend.release());
		}
		std::unique_ptr<Engine> engine(op->CreateEngine());
		const bool hasLoss = index != 4;
		auto *loss = hasLoss ? dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(0)) : nullptr;
		Require(hasLoss ? bool(loss) : op->GetNumberOfExtentions() == 0, "Python XML loss provider omitted/added");
		const double dt = op->GetTimestep(), initial = op->Grid().Energy(dt);
		Require(engine->IterateTS(100), "Python XML field evolution failed");
		double error = 0;
		for (unsigned int n = 0; n < 3; ++n)
		{
			const double rate = magnetic ? magneticLoss[n] / (mu[n] * MUE0) : electric[n] / (eps[n] * EPS0);
			const double a = .5 * dt * rate, expected = std::pow((1 - a) / (1 + a), 100);
			const Grid::Point point{.0015, .0015, .0015};
			const double actual = magnetic ? op->Grid().Magnetic(n, point) : op->Grid().Electric(n, point);
			error = std::max(error, std::abs(actual / (.2 + .3 * n) - expected));
		}
		const double balance = std::abs(op->Grid().Energy(dt) +
			(loss ? loss->ElectricLoss() + loss->MagneticLoss() : 0) - initial) / initial;
		Require(error < (device ? 2e-4 : 2e-10) && balance < (device ? 2e-4 : 2e-10),
			"Python XML analytic constitutive decay or loss ledger differs");
		std::cout << std::setprecision(17) << "PASS: PYTHON_XML_NATIVE case=" << index << " magnetic=" << magnetic
			<< " backend=" << flag << " leaves=" << op->Grid().CellCount() << " timestep=" << dt
			<< " error=" << error << " balance=" << balance << '\n';
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}
