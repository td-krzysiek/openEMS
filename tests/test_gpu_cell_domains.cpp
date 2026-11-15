// SPDX-License-Identifier: GPL-3.0-or-later
#include "openems.h"
#include "FDTD/operator.h"
#include "FDTD/cell_model.h"
#include "FDTD/gpu_cell_model.h"
#include "FDTD/gpu_backend_reference.h"
#include "FDTD/gpu_ext.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>

class Simulation : public openEMS
{
public:
	const Operator* GetOperator() const { return FDTD_Op; }
};

static void Require(bool value, const char* message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

int main(int argc, char** argv)
{
	try
	{
		Require(argc == 2, "expected a public XML fixture");
		Simulation sim;
		sim.SetLibraryArguments({"engine=basic", "numThreads=1"});
		Require(sim.ParseFDTDSetup(argv[1]) && sim.SetupFDTD() == 0, "fixture setup");
		const Operator* op = sim.GetOperator();
		unsigned int lines[3];
		for (int n = 0; n < 3; ++n)
		{
			lines[n] = op->GetNumberOfLines(n, true);
		}
		CellModel cm(op);
		cm.Analyse(true);
		cm.Build(2, 0, 1);
		const GPU_CellModel model(cm, lines);
		const size_t plane = size_t(lines[1]) * lines[2];
		for (unsigned int begin : {0u, 1u, lines[0] / 2, lines[0] - 1})
		{
			const unsigned int end = std::min(begin + 3, lines[0]);
			const GPU_CellModel slice(model, begin, end);
			Require(slice.numLines[0] == end - begin + 2, "storage includes two halos");
			for (int f = 0; f < 2; ++f)
			{
				const auto& src = model.field[f];
				const auto& dst = slice.field[f];
				uint64_t slots = 0;
				Require(dst.index.size() == size_t(end - begin + 2) * plane, "local index size");
				Require(dst.poles.size() == src.poles.size(), "preserve every ADE pole");
				for (unsigned int x = 0; x < slice.numLines[0]; ++x)
				{
					const unsigned int gx = unsigned(std::max(0, std::min(int(lines[0]) - 1, int(begin) + int(x) - 1)));
					for (size_t i = 0; i < plane; ++i)
					{
						const size_t a = size_t(x) * plane + i, b = size_t(gx) * plane + i;
						Require(dst.index[a] == src.index[b], "unchanged node combination");
						if (!dst.slot.empty())
						{
							const bool state = x > 0 && x + 1 < slice.numLines[0] && src.slot[b] != UINT32_MAX;
							Require(dst.slot[a] == (state ? uint32_t(slots++) : UINT32_MAX), "owned state only, dense local slots");
						}
					}
					for (int n = 0; n < 3; ++n)
					{
						Require(std::equal(dst.gxy[n].begin() + size_t(x) * lines[1],
						    dst.gxy[n].begin() + size_t(x + 1) * lines[1],
						    src.gxy[n].begin() + size_t(gx) * lines[1]), "local mesh factors");
					}
				}
				Require(dst.numSlots == slots, "ADE allocation excludes halo histories");
			}
			GPU_Backend_Reference backend;
			Require(backend.InitDomain(op, model, begin, end), "domain init");
			for (unsigned int c = 0; c < 3; ++c)
			{
				for (unsigned int x = begin; x < end; ++x)
				{
					const uint64_t global = (uint64_t(c) * lines[0] + x) * plane + plane - 1;
					const uint64_t local = (uint64_t(c) * slice.numLines[0] + x - begin + 1) * plane + plane - 1;
					Require(backend.LocalFlat(global) == local, "component and halo stride remapping");
				}
			}
			GPU_PointList points;
			const std::vector<unsigned int> x = {lines[0] - 1, begin, begin, 0, end};
			points.Init(x, &backend);
			const auto values = points.Sorted(std::vector<int>{0, 1, 2, 3, 4});
			Require(values.size() == points.Size(), "point data size follows owned count");
			for (auto n : points.Order())
			{
				Require(backend.OwnsX(x[n]), "extensions never own a halo point");
			}
		}
		bool rejected = false;
		try
		{
			GPU_CellModel bad(model, 1, 1);
		}
		catch (const std::invalid_argument&)
		{
			rejected = true;
		}
		Require(rejected, "empty domain rejected");
		std::cout << "CellModel domain packing: PASS" << std::endl;
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}
