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

// Independent monolithic midpoint solve in (E, sqrt(L)*I). The production
// extension eliminates branch states and uses sparse mass-graph factors.
struct FullSystem
{
	Matrix lu;
	std::vector<std::size_t> pivot;
	explicit FullSystem(Matrix a) : lu(std::move(a)), pivot(lu.size())
	{
		for (std::size_t k = 0; k < lu.size(); ++k)
		{
			std::size_t p = k;
			for (std::size_t i = k + 1; i < lu.size(); ++i)
			{
				if (std::abs(lu[i][k]) > std::abs(lu[p][k]))
				{
					p = i;
				}
			}
			Require(std::isfinite(lu[p][k]) && lu[p][k] != 0, "monolithic surface matrix is singular");
			pivot[k] = p;
			std::swap(lu[k], lu[p]);
			for (std::size_t i = k + 1; i < lu.size(); ++i)
			{
				lu[i][k] /= lu[k][k];
				for (std::size_t j = k + 1; j < lu.size(); ++j)
				{
					lu[i][j] -= lu[i][k] * lu[k][j];
				}
			}
		}
	}
	Vector Solve(Vector rhs) const
	{
		for (std::size_t k = 0; k < rhs.size(); ++k)
		{
			std::swap(rhs[k], rhs[pivot[k]]);
		}
		for (std::size_t i = 0; i < rhs.size(); ++i)
		{
			for (std::size_t j = 0; j < i; ++j)
			{
				rhs[i] -= lu[i][j] * rhs[j];
			}
		}
		for (std::size_t i = rhs.size(); i-- > 0;)
		{
			for (std::size_t j = i + 1; j < rhs.size(); ++j)
			{
				rhs[i] -= lu[i][j] * rhs[j];
			}
			rhs[i] /= lu[i][i];
		}
		return rhs;
	}
};

static void SurfaceControl(bool bulk, bool ports, bool soft, bool pureInductor, unsigned int poles, double rateScale = 1, bool measuredFit = false, const SubgridDeviceTarget *target = nullptr)
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
	const Loss::Sampler sample = [=](const Grid::Point &, const Grid::Point &)
	{
		return Loss::Coefficients{bulk ? .02 : 0, bulk ? 100. : 0};
	};
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New({1, 1, 1, 2, 2, 2}, Grid::Index{2, 2, 2}));
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "surface fixture geometry failed");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(bc);
	Loss::SurfaceAdmittance surface;
	surface.start = {.001, .0015, .001};
	surface.stop = {.0015, .0015, .001};
	surface.conductance = pureInductor ? 0 : .01;
	for (unsigned int i = 0; i < poles; ++i)
	{
		surface.poles.push_back({pureInductor ? 0. : rateScale * 10. * (i + 1), (i + 1) * 1e-10});
	}
	if (measuredFit)
	{
		// Unchanged PR242 FitFace(5.8e7 S/m,18um,Rq=0,30GHz,8 poles).
		// A constitutive path control, not native metal-face qualification.
		surface.conductance = 4.895575978021397;
		surface.poles = {
			{0.0030629169056125664, 1.3275887376533706e-11},
			{0.005772562048392223, 6.814047666831367e-12},
			{0.05564158049785245, 1.7887276818984933e-11},
			{0.027854729560116007, 2.4386637893263613e-12},
			{0.09780046600837815, 2.3318581086566667e-12},
			{0.10409923028703612, 6.759536657507219e-13},
			{0.5860379233540944, 1.0363425815679262e-12},
		};
		poles = static_cast<unsigned int>(surface.poles.size());
	}
	auto *loss = new Loss(op.get(), sample);
	loss->SetSurfaceAdmittances({surface});
	op->AddExtension(loss);
	Require(op->CalcECOperator() == 0, "surface fixture setup failed");
	signal.SetupDiracPulse(30e9);
	Require(signal.buildExcitationSignal(512) && signal.GetLength() == 2 && signal.GetVoltageSignal()[0] == 0 && signal.GetVoltageSignal()[1] == 1,
		"surface fixture requires a real built unit Dirac drive");
	auto &g = op->Grid();
	const double dt = op->GetTimestep();
	const auto line = g.VoltagePath(surface.start, surface.stop);
	Require(line.weights.size() > 1, "surface fixture misses a constrained interface trace");
	std::vector<Grid::Resistor> resistors;
	if (ports)
	{
		resistors = {{{.001, .001, .001}, {.001, .001, .002}, 50, 1e-5}};
		g.SetResistors(resistors, dt);
	}
	g.SetElectric([](unsigned int n, const Grid::Point &p) { return std::sin(731 * p[0] + 423 * p[1] - 511 * p[2] + n); });
	g.SetMagnetic([](unsigned int n, const Grid::Point &p) { return std::cos(933 * p[0] - 331 * p[1] + 211 * p[2] + n) / 376.730313668; });
	Oracle oracle(g, sample);
	const auto ne = g.EdgeCount();
	std::vector<std::pair<std::size_t, double>> sources;
	if (soft)
	{
		const auto edge = line.weights.front().first;
		sources = {{edge, .5}, {edge, -.1}, {line.weights.back().first, -.2}};
		g.SetSources(sources);
	}
	Require(loss->BuildExtension(), "surface/port refactor failed");
	if (target)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
		Require(bool(backend), "surface backend unavailable");
		target->Verify(*backend);
		op->SetEngineType(Engine::GPU);
		op->SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> engine(op->CreateEngine());
	auto *actual = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(0));
	Require(actual, "surface constitutive engine missing");
	bool rejected = false;
	try
	{
		loss->SetSurfaceAdmittances({});
	}
	catch (const std::invalid_argument &)
	{
		rejected = true;
	}
	Require(rejected, "prepared surface operator allowed coefficient mutation");
	Matrix physicalG = oracle.conductance, totalG;
	for (const auto &a : line.weights)
	{
		for (const auto &b : line.weights)
		{
			physicalG[a.first][b.first] += surface.conductance * a.second * b.second;
		}
	}
	totalG = physicalG;
	Vector portRow(ne, 0);
	if (ports)
	{
		for (const auto &term : g.VoltagePath(resistors[0].start, resistors[0].stop).weights)
		{
			portRow[term.first] = term.second / std::sqrt(50.);
		}
		for (std::size_t i = 0; i < ne; ++i)
		{
			for (std::size_t j = 0; j < ne; ++j)
			{
				totalG[i][j] += portRow[i] * portRow[j];
			}
		}
	}
	Matrix branch = Dense(poles, ne), lhs = Dense(ne + poles, ne + poles);
	Vector rate(poles), state(poles, 0);
	for (std::size_t i = 0; i < ne; ++i)
	{
		for (std::size_t j = 0; j < ne; ++j)
		{
			lhs[i][j] = oracle.mass[i][j] + .5 * dt * totalG[i][j];
		}
	}
	for (unsigned int p = 0; p < poles; ++p)
	{
		rate[p] = surface.poles[p].resistance / surface.poles[p].inductance;
		lhs[ne + p][ne + p] = 1 + .5 * dt * rate[p];
		for (const auto &term : line.weights)
		{
			const double k = term.second / std::sqrt(surface.poles[p].inductance);
			branch[p][term.first] = k;
			lhs[term.first][ne + p] = .5 * dt * k;
			lhs[ne + p][term.first] = -.5 * dt * k;
		}
	}
	FullSystem solve(lhs);
	const double initial = g.Energy(dt);
	double peak = 0, fieldError = 0, balance = 0, ledgerError = 0;
	double electricLoss = 0, surfaceLoss = 0, magneticLoss = 0, portLoss = 0, work = 0;
	for (unsigned int step = 0; step < 512; ++step)
	{
		const double drive = step < signal.GetLength() ? signal.GetVoltageSignal()[step] : signal.GetVoltageSignal()[0];
		Vector transpose(ne, 0);
		for (std::size_t f = 0; f < g.FaceCount(); ++f)
		{
			for (std::size_t i = 0; i < ne; ++i)
			{
				transpose[i] += oracle.curl[f][i] * oracle.magnetic[f];
			}
		}
		auto rhs = Multiply(totalG, oracle.electric);
		for (std::size_t i = 0; i < ne; ++i)
		{
			rhs[i] = dt * (transpose[i] - rhs[i] + portRow[i] * std::sqrt(50.) * 1e-5 * drive);
			for (unsigned int p = 0; p < poles; ++p)
			{
				rhs[i] -= dt * branch[p][i] * state[p];
			}
		}
		for (unsigned int p = 0; p < poles; ++p)
		{
			rhs.push_back(dt * (Dot(branch[p], oracle.electric) - rate[p] * state[p]));
		}
		const auto delta = solve.Solve(rhs);
		auto next = oracle.electric, midpoint = next;
		for (std::size_t i = 0; i < ne; ++i)
		{
			next[i] += delta[i];
			midpoint[i] = .5 * (oracle.electric[i] + next[i]);
		}
		for (unsigned int p = 0; p < poles; ++p)
		{
			const double mid = state[p] + .5 * delta[ne + p];
			surfaceLoss += dt * rate[p] * mid * mid;
			state[p] += delta[ne + p];
		}
		electricLoss += dt * Dot(midpoint, Multiply(physicalG, midpoint));
		const double voltage = Dot(portRow, midpoint);
		portLoss += dt * voltage * voltage;
		work += dt * voltage * std::sqrt(50.) * 1e-5 * drive;
		for (const auto &source : sources)
		{
			const double kick = drive * source.second * oracle.sourceScale[source.first];
			work += kick * (Multiply(oracle.mass, next)[source.first] + .5 * kick - .5 * dt * transpose[source.first]);
			next[source.first] += kick;
		}
		oracle.electric = std::move(next);
		const auto curl = Multiply(oracle.curl, oracle.electric);
		for (std::size_t f = 0; f < g.FaceCount(); ++f)
		{
			const double previous = oracle.magnetic[f], a = .5 * dt * oracle.magneticRate[f];
			oracle.magnetic[f] = ((1 - a) * previous - dt * curl[f]) / (1 + a);
			const double mid = .5 * (previous + oracle.magnetic[f]);
			magneticLoss += dt * oracle.magneticRate[f] * mid * mid;
		}
		engine->IterateTS(1);
		if (target && step % 16 == 0)
		{
			const auto before = std::make_pair(actual->SurfaceLoss(), actual->AdditionalStoredEnergy());
			engine->GetVolt(0, 1, 1, 1);
			engine->GetCurr(0, 1, 1, 1);
			Require(before == std::make_pair(actual->SurfaceLoss(), actual->AdditionalStoredEnergy()),
				"repeated observation changes surface histories or ledgers");
		}
		for (const auto &pair : {std::make_pair(oracle.ReadElectric(g), oracle.electric),
			std::make_pair(oracle.ReadMagnetic(g), oracle.magnetic)})
		{
			for (std::size_t i = 0; i < pair.first.size(); ++i)
			{
				peak = std::max(peak, std::abs(pair.second[i]));
				fieldError = std::max(fieldError, std::abs(pair.first[i] - pair.second[i]));
			}
		}
		const double stored = .5 * Dot(state, state), scale = std::max(initial, std::abs(work));
		balance = std::max(balance, std::abs(g.Energy(dt) + actual->AdditionalStoredEnergy() + actual->SurfaceLoss() +
			actual->ElectricLoss() + actual->MagneticLoss() + g.ResistorLoss() - g.SourceWork() - initial) / scale);
		ledgerError = std::max({ledgerError, std::abs(actual->AdditionalStoredEnergy() - stored) / scale,
			std::abs(actual->SurfaceLoss() - surfaceLoss) / scale,
			std::abs(actual->ElectricLoss() - electricLoss) / scale,
			std::abs(actual->MagneticLoss() - magneticLoss) / scale,
			std::abs(g.ResistorLoss() - portLoss) / scale, std::abs(g.SourceWork() - work) / scale});
		unsigned int nodes[3] = {3, 3, 3};
		double e = 0, h = 0;
		Require(engine->CalcFastEnergy(nodes, e, h), "surface stored-energy observation missing");
		const double expected = 2 * (g.ElectricEnergy() + g.MagneticEnergy() + actual->AdditionalStoredEnergy());
		Require(std::abs(EPS0 * e + MUE0 * h - expected) < 1e-13 * std::max(expected, initial),
			"termination energy omits the surface histories");
	}

	std::cout << "SURFACE device=" << bool(target) << " bulk=" << bulk << " ports=" << ports << " soft=" << soft << " pureL=" << pureInductor
		<< " fit=" << measuredFit << " rate_scale=" << rateScale << " poles=" << poles << " trace=" << line.weights.size() << " field_error=" << fieldError / peak
		<< " balance=" << balance << " ledger_error=" << ledgerError << std::endl;
	const double tolerance = target ? 2e-4 : 2e-10;
	Require(fieldError / peak < tolerance && balance < tolerance && ledgerError < tolerance,
		"surface recurrence differs from monolithic physical oracle");
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "expected optional --reference, --metal or --cuda");
		SubgridDeviceTarget selected;
		const SubgridDeviceTarget *target = nullptr;
		if (argc == 2)
		{
			selected = std::string(argv[1]) == "--reference" ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromFlag(argv[1]);
			target = &selected;
		}
		SurfaceControl(false, false, false, false, 1, 1, false, target);
		SurfaceControl(false, false, false, true, 1, 1, false, target);
		SurfaceControl(true, false, false, false, 8, 1, false, target);
		SurfaceControl(true, true, true, false, 8, 1, false, target);
		SurfaceControl(false, true, true, false, 0, 1, false, target);
		SurfaceControl(true, true, true, false, 8, 1e6, false, target);
		SurfaceControl(true, true, true, false, 7, 1, true, target);
		std::cout << "PASS: surface states, full monolithic oracle and stored-energy/Ohmic ledgers" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}
