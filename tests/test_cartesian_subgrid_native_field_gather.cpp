/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "subgrid_native_gather_fixture.h"

struct Dump
{
	std::unique_ptr<Engine_Interface_FDTD> direct, selected;
	std::unique_ptr<Engine_Field_Gather> gather;
	ArrayLib::ArrayNIJK<float> values;
	std::vector<std::complex<float>> sums, expectedSums;
	bool magnetic = false;
	unsigned int mode = 0;
};

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool loss, bool domains)
{
	const double unit = axis == 1 ? .0007 : .001, offset = axis == 2 ? -3 : 0;
	Fixture legacy(target, axis, loss, unit, offset, domains), sparse(target, axis, loss, unit, offset, domains);
	std::vector<unsigned int> positions[3]{{0, 5, 6, 7, 8, 13, 14}, {0, 5, 6, 7, 8, 13, 14}, {0, 5, 6, 7, 8, 13, 14}};
	unsigned int *p[3] = {positions[0].data(), positions[1].data(), positions[2].data()};
	const unsigned int counts[3] = {7, 7, 7};
	const std::size_t lines = 49, nodes = 343;
	std::vector<std::unique_ptr<Dump>> dumps;
	// A whole over-budget candidate must leave no partly registered samples.
	Engine_Interface_FDTD guard(sparse.op.get());
	const auto before = sparse.op->RegisterObservation({});
	sparse.op->SetFieldGatherBudget(1);
	std::unique_ptr<Engine_Field_Gather> refused(guard.CreateFieldGather(false, counts, p));
	Require(!refused && sparse.op->RegisterObservation({}) == before + 1,
		"over-budget gather left partial registration");
	sparse.op->SetFieldGatherBudget(std::size_t(256) << 20);
	positions[0][0] = 15;
	Reject([&] { std::unique_ptr<Engine_Field_Gather> bad(guard.CreateFieldGather(false, counts, p)); });
	positions[0][0] = 0;
	Reject([&] { std::unique_ptr<Engine_Field_Gather> bad(sparse.op->CreateFieldGather(nullptr, false, counts, p)); });
	const unsigned int empty[3] = {0, 7, 7};
	Require(guard.CreateFieldGather(false, empty, p) == nullptr, "empty gather accepted");
	for (unsigned int mode = 0; mode < 3; ++mode)
	{
		for (bool magnetic : {false, true})
		{
			auto dump = std::make_unique<Dump>();
			dump->direct = std::make_unique<Engine_Interface_FDTD>(legacy.op.get());
			dump->selected = std::make_unique<Engine_Interface_FDTD>(sparse.op.get());
			dump->mode = mode;
			dump->magnetic = magnetic;
			for (auto *interface : {dump->direct.get(), dump->selected.get()})
			{
				interface->SetInterpolationType(static_cast<Engine_Interface_Base::InterpolationType>(mode));
			}
			dump->gather.reset(dump->selected->CreateFieldGather(magnetic, counts, p));
			Require(dynamic_cast<Field_Gather_Subgrid *>(dump->gather.get()), "prepared native gather missing");
			dump->values.Init("Selected", std::array<unsigned int, 3>{7, 7, 7});
			dump->sums.resize(3 * nodes);
			dump->expectedSums.resize(3 * nodes);
			dumps.push_back(std::move(dump));
		}
	}
	std::vector<std::unique_ptr<Dump>> zeroDumps;
	for (bool magnetic : {false, true})
	{
		auto dump = std::make_unique<Dump>();
		dump->direct = std::make_unique<Engine_Interface_FDTD>(legacy.op.get());
		dump->selected = std::make_unique<Engine_Interface_FDTD>(sparse.op.get());
		dump->magnetic = magnetic;
		dump->mode = magnetic ? Engine_Interface_Base::NODE_INTERPOLATE : Engine_Interface_Base::CELL_INTERPOLATE;
		for (auto *interface : {dump->direct.get(), dump->selected.get()})
		{
			interface->SetInterpolationType(static_cast<Engine_Interface_Base::InterpolationType>(dump->mode));
		}
		unsigned int one[3] = {1, 1, 1};
		unsigned int corner[3] = {magnetic ? 0u : 14u, magnetic ? 0u : 14u, magnetic ? 0u : 14u};
		unsigned int *coordinates[3] = {&corner[0], &corner[1], &corner[2]};
		dump->gather.reset(dump->selected->CreateFieldGather(magnetic, one, coordinates));
		auto *gather = dynamic_cast<Field_Gather_Subgrid *>(dump->gather.get());
		Require(gather && gather->region.Size() == 0, "zero-stencil fixture did not create an empty region");
		dump->values.Init("Zero", std::array<unsigned int, 3>{1, 1, 1});
		zeroDumps.push_back(std::move(dump));
	}
	Reject([&] { sparse.op->SetFieldGatherBudget(0); });
	sparse.op->SealObservations();
	Require(guard.CreateFieldGather(false, counts, p) == nullptr, "late gather registered after sealing");
	double fieldError = 0, spectrumError = 0;
	for (unsigned int batch = 0; batch < 4; ++batch)
	{
		legacy.engine->IterateTS(16);
		sparse.engine->IterateTS(16);
		for (auto &dump : dumps)
		{
			std::array<std::exception_ptr, 4> errors;
			std::vector<std::thread> readers;
			for (std::size_t worker = 0; worker < 4; ++worker)
			{
				readers.emplace_back([&, worker]
				{
					try
					{
						dump->gather->Evaluate(worker * lines / 4, (worker + 1) * lines / 4, dump->values);
					}
					catch (...)
					{
						errors[worker] = std::current_exception();
					}
				});
			}
			for (auto &reader : readers)
			{
				reader.join();
			}
			for (const auto &error : errors)
			{
				if (error)
				{
					std::rethrow_exception(error);
				}
			}
			auto *gather = dynamic_cast<Field_Gather_Subgrid *>(dump->gather.get());
			Require(gather->Readbacks() == batch + 1, "concurrent gather refreshed more than once per epoch");
			dump->gather->Evaluate(0, lines, dump->values);
			Require(gather->Readbacks() == batch + 1, "same-epoch gather was not cached");
			double delta = 0, scale = 0;
			const auto weight = std::polar(.1f, static_cast<float>(-2 * std::acos(-1.) * 24.125e9 *
				dump->direct->GetTime(dump->magnetic)));
			for (unsigned int i = 0; i < counts[0]; ++i)
			{
				for (unsigned int j = 0; j < counts[1]; ++j)
				{
					for (unsigned int k = 0; k < counts[2]; ++k)
					{
						const unsigned int point[3] = {positions[0][i], positions[1][j], positions[2][k]};
						double expected[3];
						if (dump->magnetic)
						{
							dump->direct->GetHField(point, expected);
						}
						else
						{
							dump->direct->GetEField(point, expected);
						}
						for (unsigned int n = 0; n < 3; ++n)
						{
							const float reference = static_cast<float>(expected[n]), actual = dump->values(n, i, j, k);
							Require(std::isfinite(reference) && std::isfinite(actual), "nonfinite gathered field");
							delta = std::max(delta, std::abs(double(reference) - actual));
							scale = std::max(scale, std::abs(double(reference)));
							const auto index = n * nodes + (i * counts[1] + j) * counts[2] + k;
							dump->sums[index] += actual * weight;
							dump->expectedSums[index] += reference * weight;
						}
					}
				}
			}
			Require(scale > 0, "native gather fixture never excited fields");
			fieldError = std::max(fieldError, delta / scale);
		}
		for (auto &dump : zeroDumps)
		{
			dump->gather->Evaluate(0, 1, dump->values);
			unsigned int corner[3] = {dump->magnetic ? 0u : 14u, dump->magnetic ? 0u : 14u, dump->magnetic ? 0u : 14u};
			double expected[3];
			if (dump->magnetic)
			{
				dump->direct->GetHField(corner, expected);
			}
			else
			{
				dump->direct->GetEField(corner, expected);
			}
			for (unsigned int n = 0; n < 3; ++n)
			{
				Require(expected[n] == 0 && dump->values(n, 0, 0, 0) == 0,
					"empty boundary gather did not preserve the direct zero field");
			}
			Require(dynamic_cast<Field_Gather_Subgrid *>(dump->gather.get())->Readbacks() == batch + 1,
				"empty region cache epoch changed");
		}
		Require(sparse.op->NativeObservationInfo().fullFieldReads == 0,
			"prepared field gather downloaded complete fields");
	}
	for (auto &dump : dumps)
	{
		double delta = 0, scale = 0;
		for (std::size_t i = 0; i < dump->sums.size(); ++i)
		{
			delta = std::max(delta, double(std::abs(dump->sums[i] - dump->expectedSums[i])));
			scale = std::max(scale, double(std::abs(dump->expectedSums[i])));
		}
		spectrumError = std::max(spectrumError, delta / std::max(scale, 1e-300));
		Reject([&] { dump->gather->Evaluate(0, lines + 1, dump->values); });
		const float invalidSource = 0;
		Reject([&] { dump->gather->Evaluate(0, lines, dump->values, &invalidSource); });
	}
	Require(fieldError < 2e-6 && spectrumError < 2e-6, "native gathered field/DFT differs from direct full fallback");
	sparse.op->ObservationGrid();
	Require(sparse.op->NativeObservationInfo().fullFieldReads == 1, "explicit complete fallback missing");
	dumps.front()->gather->Evaluate(0, lines, dumps.front()->values);
	Require(sparse.op->NativeObservationInfo().fullFieldReads == 1, "gather repeated a complete fallback");
	std::cout << "NATIVE_FIELD_GATHER axis=" << axis << " loss=" << loss << " domains=" << domains
		<< " modes=3 fields=2 batches=4 threads=4 nodes=" << nodes << " field_error=" << fieldError
		<< " dft_error=" << spectrumError << " full_reads_before_fallback=0 fallback_reads=1 zero_regions=2\n";
}

int main(int argc, char **argv)
{
	try
	{
		Require(argc <= 2, "unexpected native gather argument");
		const auto target = argc == 1 ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromDomainFlag(argv[1]);
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			for (bool loss : {false, true})
			{
				for (bool domains : {false, true})
				{
					Control(target, axis, loss, domains);
				}
			}
		}
		std::cout << "PASS: prepared native E/H fields, bounded tables, concurrent gathers and complex DFT\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}
