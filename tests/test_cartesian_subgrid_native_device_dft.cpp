/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid_runtime.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/engine.h"
#include "FDTD/engine_interface_fdtd.h"
#include "FDTD/subgrid_native_field_gather.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_conductivity.h"
#include "subgrid_device_target.h"
#include <limits>
#include "openems.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <thread>

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
	try
	{
		action();
	}
	catch (const std::logic_error &)
	{
		return;
	}
	throw std::runtime_error("invalid native field gather operation accepted");
}

struct Fixture
{
	ContinuousStructure geometry;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	std::unique_ptr<Engine> engine;
	Engine_Ext_CartesianSubgrid_PML *pml = nullptr;
	Fixture(const SubgridDeviceTarget &target, unsigned int axis, bool loss,
		double unit, double offset, bool domains)
	{
		geometry.GetGrid()->SetDeltaUnit(unit);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= 14; ++i)
			{
				geometry.GetGrid()->AddDiscLine(n, offset + i + .01 * i * i);
			}
		}
		signal.SetupDiracPulse(30e9);
		const auto coordinate = [offset](double k) { return offset + k + .01 * k * k; };
		op.reset(Operator_CartesianSubgrid::New({coordinate(6), coordinate(6), coordinate(6),
			coordinate(8), coordinate(8), coordinate(8)}));
		Grid::Axes children;
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i < 14; ++i)
			{
				children[n].push_back(coordinate(i));
				children[n].push_back(coordinate(i) + .3 * (coordinate(i + 1) - coordinate(i)));
			}
			children[n].push_back(coordinate(14));
		}
		op->SetChildLines(children);
		op->SetEngineType(Engine::GPU);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&geometry), "native domain geometry rejected");
		int boundaries[6] = {3, 3, 3, 3, 3, 3};
		op->SetBoundaryCondition(boundaries);
		op->SetTimestepFactor(.8);
		op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(op.get(), {4, 4, 4, 4, 4, 4}));
		if (loss)
		{
			op->AddExtension(new Operator_Ext_CartesianSubgrid_Conductivity(op.get(),
				[unit, coordinate, offset](const Grid::Point &lo, const Grid::Point &hi)
				{
					bool inside = true;
					for (unsigned int n = 0; n < 3; ++n)
					{
						const double p = .5 * (lo[n] + hi[n]) / unit;
						inside = inside && p >= coordinate(6) - offset && p < coordinate(8) - offset;
					}
					return inside ? Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{.02, 50} :
						Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{};
				}));
		}
		Require(op->CalcECOperator() == 0, "native domain operator failed");
		signal.SetupDiracPulse(30e9);
		Require(signal.buildExcitationSignal(1024), "native domain waveform failed");
		op->Grid().SetElectric([unit, coordinate, offset](unsigned int n, const Grid::Point &p)
		{
			return (n == 1 ? 1. : -.5) * std::exp(-std::pow((p[0]/unit-coordinate(7)+offset)/2, 2) -
				std::pow((p[1]/unit-coordinate(7)+offset)/2, 2) - std::pow((p[2]/unit-coordinate(7)+offset)/2, 2));
		});
		op->Grid().SetSources({{op->Grid().EdgeCount() / 2, .01}});
		if (domains)
		{
			std::vector<std::unique_ptr<GPU_Backend>> backends;
			for (std::size_t k = 0; k < 2; ++k)
			{
				target.SelectDomain(k);
				backends.emplace_back(GPU_Backend::New(target.backend));
				Require(bool(backends.back()), "native worker unavailable");
				target.Verify(*backends.back());
				std::cout << "GATHER_WORKER " << k << " " << backends.back()->GetName() << std::endl;
			}
			op->SetDomainExecution(axis, {coordinate(6.1)}, std::move(backends));
		}
		else
		{
			target.SelectDomain(0);
			std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
			Require(bool(backend), "native single backend unavailable");
			target.Verify(*backend);
			std::cout << "GATHER_SINGLE " << backend->GetName() << std::endl;
			op->SetGPUBackend(backend.release());
		}
		engine.reset(op->CreateEngine());
		for (std::size_t k = 0; k < engine->GetExtensionCount(); ++k)
		{
			if (auto *extension = dynamic_cast<Engine_Ext_CartesianSubgrid_PML *>(engine->GetExtension(k)))
			{
				pml = extension;
			}
		}
		Require(pml, "native PML mirror missing");
	}
};


struct Dump
{
	std::unique_ptr<Engine_Interface_FDTD> direct, selected;
	std::unique_ptr<Engine_Field_Gather> gather;
	std::vector<std::unique_ptr<ArrayLib::ArrayNIJK<std::complex<float>>>> fields;
	std::vector<std::vector<std::complex<float>>> expected;
	bool magnetic = false;
	unsigned int mode = 0;
	int id = -1;
};

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool loss, bool domains)
{
	const double unit = axis == 1 ? .0007 : .001, offset = axis == 2 ? -3 : 0;
	Fixture legacy(target, axis, loss, unit, offset, domains), native(target, axis, loss, unit, offset, domains);
	std::vector<unsigned int> positions[3]{{0, 5, 6, 7, 8, 13, 14}, {0, 5, 6, 7, 8, 13, 14}, {0, 5, 6, 7, 8, 13, 14}};
	unsigned int *p[3] = {positions[0].data(), positions[1].data(), positions[2].data()};
	const unsigned int counts[3] = {7, 7, 7};
	constexpr std::size_t nodes = 343;
	std::vector<std::unique_ptr<Dump>> dumps;
	Engine_Interface_FDTD guard(native.op.get());
	const auto empty = native.op->RegisterObservation({});
	Require(guard.CreateFieldDFT(nullptr, 3) == -1, "null/native DFT fallback changed");
	std::vector<ArrayLib::ArrayNIJK<std::complex<float>> *> outputs;
	Require(!guard.ReadFieldDFT(99, outputs), "unknown DFT ID read succeeded");
	Reject([&] { guard.AccumulateFieldDFT(99, {{1, 0}}); });
	// Whole budget fallback must retain the gather as a genuine TD/host reader.
	auto fallback = std::unique_ptr<Engine_Field_Gather>(guard.CreateFieldGather(false, counts, p));
	native.op->SetFieldDFTBudget(1);
	Require(guard.CreateFieldDFT(fallback.get(), 3) == -1 && native.op->FieldDFTStorageBytes() == 0,
		"budget fallback uploaded or partly registered DFT");
	native.op->SetFieldDFTBudget(std::size_t(256) << 20);
	const bool hostFallback = axis == 0 && !loss && !domains;
	if (!hostFallback)
	{
		Require(guard.CreateFieldDFT(fallback.get(), 1) >= 0, "retry after whole budget fallback failed");
	}
	for (unsigned int mode = 0; mode < 3; ++mode)
	{
		for (bool magnetic : {false, true})
		{
			auto dump = std::make_unique<Dump>();
			dump->direct = std::make_unique<Engine_Interface_FDTD>(legacy.op.get());
			dump->selected = std::make_unique<Engine_Interface_FDTD>(native.op.get());
			dump->mode = mode;
			dump->magnetic = magnetic;
			for (auto *interface : {dump->direct.get(), dump->selected.get()})
			{
				interface->SetInterpolationType(static_cast<int>(mode));
			}
			dump->gather.reset(dump->selected->CreateFieldGather(magnetic, counts, p));
			Require(dynamic_cast<Field_Gather_Subgrid *>(dump->gather.get()), "native DFT gather missing");
			Require(dump->selected->CreateFieldDFT(dump->gather.get(), 0) == -1, "zero-frequency DFT accepted");
			Reject([&] { guard.CreateFieldDFT(dump->gather.get(), 3); });
			dump->id = dump->selected->CreateFieldDFT(dump->gather.get(), 3);
			Require(dump->id >= 0, "native DFT unavailable");
			Reject([&] { dump->selected->CreateFieldDFT(dump->gather.get(), 3); });
			Reject([&] { dump->selected->AccumulateFieldDFT(dump->id, {{1, 0}, {1, 0}, {1, 0}}); });
			for (unsigned int frequency = 0; frequency < 3; ++frequency)
			{
				auto field = std::make_unique<ArrayLib::ArrayNIJK<std::complex<float>>>();
				field->Init("NativeDFT", std::array<unsigned int, 3>{7, 7, 7});
				dump->fields.push_back(std::move(field));
			}
			dump->expected.assign(3, std::vector<std::complex<float>>(3 * nodes));
			dumps.push_back(std::move(dump));
		}
	}
	std::vector<std::unique_ptr<Dump>> zeros;
	for (bool magnetic : {false, true})
	{
		auto dump = std::make_unique<Dump>();
		dump->selected = std::make_unique<Engine_Interface_FDTD>(native.op.get());
		dump->selected->SetInterpolationType(magnetic ? Engine_Interface_Base::NODE_INTERPOLATE : Engine_Interface_Base::CELL_INTERPOLATE);
		unsigned int corner[3] = {magnetic ? 0u : 14u, magnetic ? 0u : 14u, magnetic ? 0u : 14u};
		unsigned int *coords[3] = {&corner[0], &corner[1], &corner[2]};
		unsigned int one[3] = {1, 1, 1};
		dump->gather.reset(dump->selected->CreateFieldGather(magnetic, one, coords));
		Require(dynamic_cast<Field_Gather_Subgrid *>(dump->gather.get())->region.Size() == 0, "empty DFT support fixture differs");
		dump->id = dump->selected->CreateFieldDFT(dump->gather.get(), 1);
		Require(dump->id >= 0, "empty native DFT rejected");
		auto field = std::make_unique<ArrayLib::ArrayNIJK<std::complex<float>>>();
		field->Init("EmptyDFT", std::array<unsigned int, 3>{1, 1, 1});
		dump->fields.push_back(std::move(field));
		zeros.push_back(std::move(dump));
	}
	Reject([&] { native.op->SetFieldDFTBudget(0); });
	native.op->SealObservations();
	Require(native.op->EvaluateObservation(empty) == 0, "empty host placeholder changed");
	const auto configured = native.op->NativeObservationInfo();
	Require(configured.fullFieldReads == 0 && configured.selectedMagneticStates == 0,
		"device-only H recipes remain in host registry");
	Require(hostFallback ? configured.selectedElectricStates > 0 : configured.selectedElectricStates == 0,
		"DFT-only E support remains active, or budget fallback lost its live support");
	double maximum = 0;
	for (unsigned int batch = 0; batch < 8; ++batch)
	{
		legacy.engine->IterateTS(8);
		native.engine->IterateTS(8);
		for (auto &dump : dumps)
		{
			std::vector<std::complex<float>> weights;
			for (unsigned int frequency = 0; frequency < 3; ++frequency)
			{
				const double phase = -2 * std::acos(-1.) * (23.5e9 + .625e9 * frequency) * dump->direct->GetTime(dump->magnetic);
				// Include the actual very small timestep-scaled FD weight, plus a
				// different schedule/phase for H through GetTime(dualTime).
				const float scale = static_cast<float>(2 * legacy.op->GetTimestep() * 8);
				weights.push_back(std::polar(scale, static_cast<float>(phase)));
			}
			const auto before = native.op->NativeObservationInfo();
			dump->selected->AccumulateFieldDFT(dump->id, weights);
			Require(native.op->NativeObservationInfo().fieldDFTDownloadBytes == before.fieldDFTDownloadBytes,
				"native DFT sample downloaded field/sums");
			Reject([&] { dump->selected->AccumulateFieldDFT(dump->id, weights); });
			Reject([&] { dump->selected->AccumulateFieldDFT(dump->id, {}); });
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
						for (unsigned int frequency = 0; frequency < 3; ++frequency)
						{
							for (unsigned int n = 0; n < 3; ++n)
							{
								dump->expected[frequency][n * nodes + (i * 7 + j) * 7 + k] +=
									static_cast<float>(expected[n]) * weights[frequency];
							}
						}
					}
				}
			}
			outputs.clear();
			for (auto &field : dump->fields)
			{
				outputs.push_back(field.get());
			}
			Require(dump->selected->ReadFieldDFT(dump->id, outputs), "native DFT snapshot unavailable");
			for (unsigned int frequency = 0; frequency < 3; ++frequency)
			{
				double delta = 0, scale = 0;
				for (std::size_t k = 0; k < 3 * nodes; ++k)
				{
					delta = std::max(delta, double(std::abs(dump->fields[frequency]->data(k) - dump->expected[frequency][k])));
					scale = std::max(scale, double(std::abs(dump->expected[frequency][k])));
				}
				Require(scale > 0, "native DFT control has no signal");
				maximum = std::max(maximum, delta / scale);
			}
		}
		for (auto &dump : zeros)
		{
			const auto before = native.op->NativeObservationInfo().fieldDFTDownloadBytes;
			dump->selected->AccumulateFieldDFT(dump->id, {{1, 0}});
			outputs = {dump->fields.front().get()};
			Require(dump->selected->ReadFieldDFT(dump->id, outputs), "empty DFT snapshot unavailable");
			Require(native.op->NativeObservationInfo().fieldDFTDownloadBytes == before, "empty native DFT transferred states");
			for (std::size_t k = 0; k < dump->fields.front()->size(); ++k)
			{
				Require(dump->fields.front()->data(k) == std::complex<float>{}, "empty native DFT is nonzero");
			}
		}
		Require(native.op->NativeObservationInfo().fullFieldReads == 0, "DFT-only processing read complete fields");
	}
	Require(maximum < 2e-5, "native device DFT differs from direct host complex arrays");
	// Genuine host fallback remains current; it must not read DFT placeholders.
	ArrayLib::ArrayNIJK<float> live;
	live.Init("Fallback", std::array<unsigned int, 3>{7, 7, 7});
	if (hostFallback)
	{
		fallback->Evaluate(0, 49, live);
		Require(native.op->NativeObservationInfo().fullFieldReads == 0, "budget fallback lost prepared selected support");
	}
	// Ad-hoc evaluation of an already device-only gather is explicit full fallback.
	dumps.front()->gather->Evaluate(0, 49, live);
	Require(native.op->NativeObservationInfo().fullFieldReads == 1, "device-only TD request returned stale placeholders");
	auto &sample = dumps.front();
	unsigned int point[3] = {positions[0][3], positions[1][3], positions[2][3]};
	double reference[3];
	sample->direct->GetEField(point, reference);
	for (unsigned int n = 0; n < 3; ++n)
	{
		Require(std::abs(double(live(n, 3, 3, 3)) - reference[n]) < 2e-6 * std::max(1., std::abs(reference[n])),
			"device-only explicit TD fallback is stale");
	}
	outputs = {sample->fields.front().get()};
	Reject([&] { sample->selected->ReadFieldDFT(sample->id, outputs); });
	outputs.assign(3, nullptr);
	Reject([&] { sample->selected->ReadFieldDFT(sample->id, outputs); });
	outputs = {sample->fields[0].get(), sample->fields[1].get(), sample->fields[2].get()};
	sample->selected->SetInterpolationType(Engine_Interface_Base::NODE_INTERPOLATE);
	Reject([&] { sample->selected->ReadFieldDFT(sample->id, outputs); });
	sample->selected->SetInterpolationType(Engine_Interface_Base::NO_INTERPOLATION);
	std::cout << "NATIVE_DEVICE_DFT axis=" << axis << " loss=" << loss << " domains=" << domains
		<< " modes=3 fields=2 frequencies=3 batches=8 nodes=343 complex_error=" << maximum
		<< " no_full_reads_before_fallback=1 host_E_support=" << configured.selectedElectricStates
		<< " host_H_support=0 empty_groups=2\n";
}

int main(int argc, char **argv)
{
	try
	{
		Require(argc <= 2, "unexpected native DFT argument");
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
		std::cout << "PASS: native device DFT arrays, filtered host selection, empty support and live fallback\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}
