/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/engine.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_conductivity.h"
#include "tools/useful.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

int main(int argc, char **argv)
{
	try
	{
		if (argc != 3)
		{
			throw std::invalid_argument("usage: preparation_benchmark bulk|surface count");
		}
		const std::string mode = argv[1];
		const auto count = std::stoull(argv[2]);
		if ((mode != "bulk" && mode != "surface") || count == 0 ||
			(mode == "bulk" && (count < 4 || count > 48)) || (mode == "surface" && count > 250000))
		{
			throw std::invalid_argument("bounded diagnostic count exceeded");
		}
		using Grid = CartesianSubgrid;
		using Loss = Operator_Ext_CartesianSubgrid_Conductivity;
		ContinuousStructure csx;
		csx.GetGrid()->SetDeltaUnit(.001);
		const auto cells = mode == "bulk" ? count : 4;
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (std::size_t i = 0; i <= cells; ++i)
			{
				csx.GetGrid()->AddDiscLine(n, double(i));
			}
		}
		Excitation signal;
		signal.SetupDiracPulse(30e9);
		std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New(
			{1, 1, 1, 2, 2, 2}, Grid::Index{2, 2, 2}));
		op->SetEngineType(Engine::BASIC);
		op->SetExcitationSignal(&signal);
		if (!op->SetGeometryCSX(&csx))
		{
			throw std::runtime_error("diagnostic geometry failed");
		}
		int boundary[6] = {0, 0, 0, 0, 0, 0};
		op->SetBoundaryCondition(boundary);
		auto *loss = new Loss(op.get(), [=](const Grid::Point &, const Grid::Point &)
		{
			return Loss::Coefficients{mode == "bulk" ? .02 : 0., 0.};
		});
		if (mode == "surface")
		{
			Loss::SurfaceAdmittance surface;
			surface.start = {.001, .0015, .001};
			surface.stop = {.0015, .0015, .001};
			surface.conductance = .01;
			for (unsigned int pole = 0; pole < 7; ++pole)
			{
				surface.poles.push_back({10. * (pole + 1), 1e-10 * (pole + 1)});
			}
			// Repeated passive paths stress storage, not a physical metal fit.
			loss->SetSurfaceAdmittances(std::vector<Loss::SurfaceAdmittance>(count, surface));
		}
		op->AddExtension(loss);
		if (op->CalcECOperator() != 0)
		{
			throw std::runtime_error("diagnostic operator preparation failed");
		}
		auto &grid = op->Grid();
		signal.SetupDiracPulse(30e9);
		if (!signal.buildExcitationSignal(32))
		{
			throw std::runtime_error("diagnostic excitation preparation failed");
		}
		grid.SetElectric([](unsigned int axis, const Grid::Point &point)
		{
			return 1e-5 * std::sin(731 * point[0] + 423 * point[1] - 511 * point[2] + axis);
		});
		grid.SetMagnetic([](unsigned int axis, const Grid::Point &point)
		{
			return 1e-8 * std::cos(933 * point[0] - 331 * point[1] + 211 * point[2] + axis);
		});
		std::unique_ptr<Engine> engine(op->CreateEngine());
		engine->IterateTS(16);
		std::cout << "BENCHMARK " << mode << " count=" << count << " cells=" << grid.CellCount()
			<< " E=" << grid.EdgeCount() << " H=" << grid.FaceCount() << std::setprecision(17)
			<< " energy=" << grid.Energy(op->GetTimestep());
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			std::cout << " e" << axis << '=' << grid.Electric(axis, {.00125, .00125, .00125})
				<< " h" << axis << '=' << grid.Magnetic(axis, {.00125, .00125, .00125});
		}
		std::cout << '\n';
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}
