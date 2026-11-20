/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/engine.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "tools/useful.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

using Grid = CartesianSubgrid;
using Loss = Operator_Ext_CartesianSubgrid_Conductivity;
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;
static void Require(bool valid, const char *message)
{
	if (!valid)
	{
		throw std::runtime_error(message);
	}
}
static Matrix Dense(std::size_t rows, std::size_t columns)
{
	return Matrix(rows, Vector(columns));
}
static double Dot(const Vector &a, const Vector &b)
{
	return std::inner_product(a.begin(), a.end(), b.begin(), 0.);
}
static Vector Multiply(const Matrix &matrix, const Vector &value)
{
	Vector result(matrix.size());
	for (std::size_t i = 0; i < matrix.size(); ++i)
	{
		result[i] = Dot(matrix[i], value);
	}
	return result;
}
static Matrix Factor(Matrix matrix)
{
	for (std::size_t i = 0; i < matrix.size(); ++i)
	{
		for (std::size_t j = 0; j <= i; ++j)
		{
			for (std::size_t k = 0; k < j; ++k)
			{
				matrix[i][j] -= matrix[i][k] * matrix[j][k];
			}
			Require(std::isfinite(matrix[i][j]) && (i != j || matrix[i][j] > 0), "dense oracle matrix is not SPD");
			matrix[i][j] = i == j ? std::sqrt(matrix[i][j]) : matrix[i][j] / matrix[j][j];
		}
	}
	return matrix;
}
static Vector Solve(const Matrix &lower, Vector rhs)
{
	for (std::size_t i = 0; i < rhs.size(); ++i)
	{
		for (std::size_t j = 0; j < i; ++j)
		{
			rhs[i] -= lower[i][j] * rhs[j];
		}
		rhs[i] /= lower[i][i];
	}
	for (std::size_t i = rhs.size(); i-- > 0;)
	{
		for (std::size_t j = i + 1; j < rhs.size(); ++j)
		{
			rhs[i] -= lower[j][i] * rhs[j];
		}
		rhs[i] /= lower[i][i];
	}
	return rhs;
}

//! Independent physical path/volume assembly; no extension's matrices/factors.
struct Oracle
{
	Matrix mass, conductance, curl;
	Vector magneticRate, electric, magnetic, sourceScale;
	std::vector<Grid::Integral> rootProbe;
	Oracle(const Grid &g, const Loss::Sampler &sample)
	{
		const auto ne = g.EdgeCount(), nh = g.FaceCount();
		Require(ne < 1024, "keep the independent dense oracle bounded");
		mass = Dense(ne, ne);
		conductance = Dense(ne, ne);
		curl = Dense(nh, ne);
		magneticRate.resize(nh);
		sourceScale.resize(ne);
		rootProbe.resize(ne);
		for (std::size_t c = 0; c < g.CellCount(); ++c)
		{
			const auto lo = g.LeafStart(c), length = g.LeafLengths(c);
			auto hi = lo;
			for (unsigned int n = 0; n < 3; ++n)
			{
				hi[n] += length[n];
			}
			const auto medium = g.LeafMedium(c);
			const auto sigma = sample(lo, hi);
			const double volume = length[0] * length[1] * length[2];
			for (unsigned int n = 0; n < 3; ++n)
			{
				const auto a = (n + 1) % 3, b = (n + 2) % 3;
				for (unsigned int i = 0; i < 2; ++i)
				{
					for (unsigned int j = 0; j < 2; ++j)
					{
						auto start = lo;
						start[a] += i * length[a];
						start[b] += j * length[b];
						auto stop = start;
						stop[n] += length[n];
						const auto line = g.VoltagePath(start, stop);
						if (line.weights.size() == 1)
						{
							const auto term = line.weights.front();
							rootProbe[term.first] = line;
							sourceScale[term.first] = length[n] / term.second;
						}
						for (const auto &x : line.weights)
						{
							for (const auto &y : line.weights)
							{
								const double weight = volume / (4 * length[n] * length[n]) * x.second * y.second;
								mass[x.first][y.first] += medium.epsilon * weight;
								conductance[x.first][y.first] += sigma.electric * weight;
							}
						}
					}
				}
				const double area = length[a] * length[b];
				const double reluctance = std::sqrt(volume / (2 * medium.permeability * area * area));
				for (unsigned int side = 0; side < 2; ++side)
				{
					const auto face = 6 * c + 2 * n + side;
					magneticRate[face] = sigma.magnetic / medium.permeability;
					auto edge = [&](unsigned int direction, unsigned int cross, unsigned int high, double sign)
					{
						auto start = lo;
						start[n] += side * length[n];
						start[cross] += high * length[cross];
						auto stop = start;
						stop[direction] += length[direction];
						for (const auto &term : g.VoltagePath(start, stop).weights)
						{
							curl[face][term.first] += sign * reluctance * term.second;
						}
					};
					edge(b, a, 1, 1);
					edge(b, a, 0, -1);
					edge(a, b, 1, -1);
					edge(a, b, 0, 1);
				}
			}
		}
		for (std::size_t i = 0; i < ne; ++i)
		{
			Require(rootProbe[i].weights.size() == 1 && sourceScale[i] > 0, "oracle cannot observe an independent root");
			Require(std::abs(mass[i][i] - 1) < 2e-13, "public line basis does not reproduce normalized mass");
		}
		electric = ReadElectric(g);
		magnetic = ReadMagnetic(g);
	}
	Vector ReadElectric(const Grid &g) const
	{
		Vector result(g.EdgeCount());
		for (std::size_t i = 0; i < result.size(); ++i)
		{
			result[i] = g.Evaluate(rootProbe[i]) / rootProbe[i].weights.front().second;
		}
		return result;
	}
	Vector ReadMagnetic(const Grid &g) const
	{
		Vector result(g.FaceCount());
		for (std::size_t c = 0; c < g.CellCount(); ++c)
		{
			const auto start = g.LeafStart(c), length = g.LeafLengths(c);
			const double normalization = std::sqrt(g.LeafMedium(c).permeability * length[0] * length[1] * length[2] / 2);
			for (unsigned int n = 0; n < 3; ++n)
			{
				auto a = start, b = start;
				for (unsigned int axis = 0; axis < 3; ++axis)
				{
					a[axis] += length[axis] * (axis == n ? 1. / 3 : .5);
					b[axis] += length[axis] * (axis == n ? 2. / 3 : .5);
				}
				const double ha = g.Magnetic(n, a), hb = g.Magnetic(n, b);
				result[6 * c + 2 * n] = (2 * ha - hb) * normalization;
				result[6 * c + 2 * n + 1] = (2 * hb - ha) * normalization;
			}
		}
		return result;
	}
};

static void Control(double scale, bool ports, bool soft, unsigned int steps, bool dc, const SubgridDeviceTarget *target, bool electric = true, bool magnetic = true, bool xyz = false)
{
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(.001);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= 4; ++i)
		{
			csx.GetGrid()->AddDiscLine(n, i);
		}
	}
	const Loss::Sampler sample = [=](const Grid::Point &, const Grid::Point &stop)
	{
		return Loss::Coefficients{electric ? scale * (stop[0] <= .002 ? .02 : .06) : 0, magnetic ? scale * (stop[1] <= .002 ? 100 : 250) : 0};
	};
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New({1, 1, 1, 3, 3, 3}, xyz ? Grid::Index{2, 2, 2} : Grid::Index{2, 2, 1}));
	op->SetEngineType(target ? Engine::GPU : Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "conductivity fixture geometry failed");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(bc);
	auto *loss = new Loss(op.get(), sample);
	op->AddExtension(loss);
	Require(op->CalcECOperator() == 0, "conductivity fixture setup failed");
	if (dc)
	{
		signal.SetupCustomExcite("1", 30e9, 30e9);
	}
	else
	{
		signal.SetupDiracPulse(30e9);
	}
	Require(signal.buildExcitationSignal(steps) && signal.GetLength() > 0, "conductivity test drive was not built");
	Require(std::any_of(signal.GetVoltageSignal(), signal.GetVoltageSignal() + signal.GetLength(),
		[](double value) { return value != 0; }), "conductivity test drive is identically zero");
	if (dc)
	{
		Require(signal.GetLength() >= steps && std::all_of(signal.GetVoltageSignal(),
			signal.GetVoltageSignal() + steps, [](double value) { return value == 1; }),
			"DC control requires a sustained unit drive for the entire window");
	}
	auto &grid = op->Grid();
	const double dt = op->GetTimestep();
	std::vector<Grid::Resistor> resistors;
	if (ports)
	{
		resistors = {{{.002, .002, .001}, {.002, .002, .003}, 50, .002},
			{{.001, .001, .002}, {.003, .001, .002}, 75, -.001}};
		grid.SetResistors(resistors, dt);
	}
	grid.SetElectric([](unsigned int n, const Grid::Point &p) { return std::sin(1731 * p[0] + 933 * p[1] - 671 * p[2] + n); });
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p) { return std::cos(633 * p[0] - 731 * p[1] + 411 * p[2] + n) / 376.730313668; });
	Oracle oracle(grid, sample);
	std::size_t first = 0, second = 0;
	double commutator = 0, conductivityScale = 0;
	for (std::size_t i = 0; i < oracle.mass.size(); ++i)
	{
		for (std::size_t j = 0; j < oracle.mass.size(); ++j)
		{
			if (i != j && oracle.mass[i][j] > 0)
			{
				first = i;
				second = j;
			}
			double value = 0;
			for (std::size_t k = 0; k < oracle.mass.size(); ++k)
			{
				value += oracle.mass[i][k] * oracle.conductance[k][j] - oracle.conductance[i][k] * oracle.mass[k][j];
			}
			commutator = std::max(commutator, std::abs(value));
			conductivityScale = std::max(conductivityScale, std::abs(oracle.conductance[i][j]));
		}
	}
	Require(first != second && (!scale || !electric || commutator > conductivityScale * 1e-4), "conductivity fixture lacks noncommuting interface coupling");
	std::vector<std::pair<std::size_t, double>> sources;
	if (soft)
	{
		sources = {{first, 1.2}, {second, -.7}, {first, -.2}};
		grid.SetSources(sources);
	}
	Require(loss->BuildExtension(), "conductivity port refactor failed");
	if (target)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
		Require(bool(backend), "requested conductivity backend unavailable");
		target->Verify(*backend);
		op->SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> engine(op->CreateEngine());
	auto *actual = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(0));
	Require(actual, "conductivity engine pair is missing");
	Matrix portRows;
	Vector current;
	for (const auto &port : resistors)
	{
		Vector row(grid.EdgeCount());
		for (const auto &term : grid.VoltagePath(port.start, port.stop).weights)
		{
			row[term.first] = term.second / std::sqrt(port.resistance);
		}
		portRows.push_back(std::move(row));
		current.push_back(port.current * std::sqrt(port.resistance));
	}
	Matrix totalG = oracle.conductance, implicit = oracle.mass;
	for (const auto &row : portRows)
	{
		for (std::size_t i = 0; i < row.size(); ++i)
		{
			for (std::size_t j = 0; j < row.size(); ++j)
			{
				totalG[i][j] += row[i] * row[j];
			}
		}
	}
	for (std::size_t i = 0; i < implicit.size(); ++i)
	{
		for (std::size_t j = 0; j < implicit.size(); ++j)
		{
			implicit[i][j] += .5 * dt * totalG[i][j];
		}
	}
	const auto factor = Factor(implicit);
	const double initial = grid.Energy(dt);
	double electricLoss = 0, magneticLoss = 0, portLoss = 0, sourceWork = 0;
	double largestError = 0, peak = 0, balance = 0, lossError = 0;
	for (unsigned int step = 0; step < steps; ++step)
	{
		const double drive = step < signal.GetLength() ? signal.GetVoltageSignal()[step] : 0;
		Vector transpose(grid.EdgeCount(), 0);
		for (std::size_t f = 0; f < grid.FaceCount(); ++f)
		{
			for (std::size_t i = 0; i < grid.EdgeCount(); ++i)
			{
				transpose[i] += oracle.curl[f][i] * oracle.magnetic[f];
			}
		}
		const auto ge = Multiply(totalG, oracle.electric);
		auto rhs = transpose;
		for (std::size_t i = 0; i < rhs.size(); ++i)
		{
			rhs[i] -= ge[i];
			for (std::size_t p = 0; p < portRows.size(); ++p)
			{
				rhs[i] += portRows[p][i] * current[p] * drive;
			}
		}
		const auto increment = Solve(factor, rhs);
		auto next = oracle.electric, midpoint = oracle.electric;
		for (std::size_t i = 0; i < next.size(); ++i)
		{
			next[i] += dt * increment[i];
			midpoint[i] = .5 * (oracle.electric[i] + next[i]);
		}
		electricLoss += dt * Dot(midpoint, Multiply(oracle.conductance, midpoint));
		for (std::size_t p = 0; p < portRows.size(); ++p)
		{
			const double voltage = Dot(portRows[p], midpoint);
			portLoss += dt * voltage * voltage;
			sourceWork += dt * voltage * current[p] * drive;
		}
		for (const auto &source : sources)
		{
			const double kick = drive * source.second * oracle.sourceScale[source.first];
			const auto massField = Multiply(oracle.mass, next);
			sourceWork += kick * (massField[source.first] + .5 * kick - .5 * dt * transpose[source.first]);
			next[source.first] += kick;
		}
		oracle.electric = std::move(next);
		const auto curl = Multiply(oracle.curl, oracle.electric);
		for (std::size_t f = 0; f < oracle.magnetic.size(); ++f)
		{
			const double previous = oracle.magnetic[f], rate = oracle.magneticRate[f], a = .5 * dt * rate;
			oracle.magnetic[f] = ((1 - a) * previous - dt * curl[f]) / (1 + a);
			const double mid = .5 * (previous + oracle.magnetic[f]);
			magneticLoss += dt * rate * mid * mid;
		}
		engine->IterateTS(1);
		if (target && step % 16 == 0)
		{
			const auto before = std::make_pair(actual->ElectricLoss(), actual->MagneticLoss());
			const double resistor = grid.ResistorLoss(), work = grid.SourceWork();
			engine->GetVolt(0, 1, 1, 1);
			engine->GetCurr(0, 1, 1, 1);
			Require(before == std::make_pair(actual->ElectricLoss(), actual->MagneticLoss()) &&
				resistor == grid.ResistorLoss() && work == grid.SourceWork(),
				"repeated observation changes the accumulated conductivity/source ledger");
		}
		for (const auto &pair : {std::make_pair(oracle.ReadElectric(grid), oracle.electric),
			std::make_pair(oracle.ReadMagnetic(grid), oracle.magnetic)})
		{
			for (std::size_t i = 0; i < pair.first.size(); ++i)
			{
				peak = std::max(peak, std::abs(pair.second[i]));
				largestError = std::max(largestError, std::abs(pair.first[i] - pair.second[i]));
			}
		}
		const double scaleEnergy = std::max(initial, std::abs(sourceWork));
		balance = std::max(balance, std::abs(grid.Energy(dt) + grid.ResistorLoss() + actual->ElectricLoss() +
			actual->MagneticLoss() - grid.SourceWork() - initial) / scaleEnergy);
		lossError = std::max({lossError, std::abs(actual->ElectricLoss() - electricLoss) / scaleEnergy,
			std::abs(actual->MagneticLoss() - magneticLoss) / scaleEnergy,
			std::abs(grid.ResistorLoss() - portLoss) / scaleEnergy,
			std::abs(grid.SourceWork() - sourceWork) / scaleEnergy});
	}
	const double tolerance = target ? 2e-4 : 2e-10;
	Require(largestError / peak < tolerance && balance < tolerance && lossError < tolerance,
		"conductivity fields/loss/source balance differ from independent dense physical oracle");
	std::cout << "CONDUCTIVITY xyz=" << xyz << " electric=" << electric << " magnetic=" << magnetic << " device=" << bool(target) << " scale=" << scale << " ports=" << ports << " soft=" << soft << " steps=" << steps
		<< " dc=" << dc << " source_work=" << grid.SourceWork() << " field_error=" << largestError / peak << " balance=" << balance << " loss_error=" << lossError
		<< " commutator=" << commutator / std::max(conductivityScale, 1e-300) << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		SubgridDeviceTarget selected;
		const SubgridDeviceTarget *target = nullptr;
		if (argc > 1)
		{
			const std::string flag = argv[1];
			selected = flag == "--reference" ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromFlag(flag);
			target = &selected;
		}
		Control(0, false, false, 256, false, target);
		Control(0, true, true, 256, false, target);
		Control(1, false, false, 256, false, target);
		Control(1, true, false, 256, false, target);
		Control(1, true, true, 256, false, target);
		Control(1, false, true, 512, false, target);
		Control(10000, false, true, 512, false, target);
		Control(1, true, true, 2048, true, target);
		Control(10000, true, true, 256, false, target);
		Control(1, true, true, 256, false, target, true, false);
		Control(1, true, true, 256, false, target, false, true);
		Control(1, true, true, 256, false, target, true, true, true);
		std::cout << "PASS: bounded conductivity oracle controls; not native CAD/material or production qualification" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}
