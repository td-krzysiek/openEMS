/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_device_precision.h"
#include "FDTD/kernels/subgrid_update.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

static void Require(bool condition, const char *message)
{
	if (!condition)
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
	catch (const std::invalid_argument &)
	{
		return;
	}
	throw std::runtime_error("invalid/unsafe factor conversion was accepted");
}


static void CurlCancellation()
{
	const std::uint32_t start[] = {0, 3}, column[] = {0, 1, 2}, base[] = {0}, pattern[] = {0};
	const double exact[] = {100000001, 3.141592653589793, -100000000};
	const float field[] = {1, .75f, 1};
	float high[3], low[3], output = 0, naive = 0;
	double expected = 0;
	for (unsigned int i = 0; i < 3; ++i)
	{
		high[i] = static_cast<float>(exact[i]);
		low[i] = static_cast<float>(exact[i] - double(high[i]));
		expected += exact[i] * field[i];
		naive += high[i] * field[i];
	}
	Op_SubgridCurl operation{start, column, base, pattern, high, low, field, &output, 1, 0, 0};
	operation(0);
	Require(std::abs(double(output) - expected) < 5e-7, "split curl lost cancelled coefficient residuals");
	Require(std::abs(double(naive) - expected) > 1, "cancellation negative control no longer discriminates");
	operation.additive = 1;
	operation.scale = -.25f;
	const float before = output;
	operation(0);
	Require(std::abs(double(output) - (double(before) - .25 * expected)) < 5e-7,
		"compensated additive magnetic curl differs");
}


static void TimestepPrecision()
{
	std::size_t improved = 0;
	for (unsigned int i = 1; i < 10000; ++i)
	{
		const double timestep = 1e-12 * (1 + double(i % 101) / 107);
		const float high = static_cast<float>(timestep);
		const float low = static_cast<float>(timestep - double(high));
		const float previous = float(std::sin(double(i))), increment = float(1e12 * std::cos(double(i)));
		const float expected = float(double(previous) + timestep * double(increment));
		const float actual = subgrid_field_increment(previous, high, low, increment);
		const float naive = previous + high * increment;
		const double bound = 2 * std::numeric_limits<float>::epsilon() *
			std::max(std::abs(double(previous)), std::abs(timestep * increment));
		Require(std::abs(double(actual) - expected) <= bound, "split timestep field update lost precision");
		if (std::abs(double(actual) - expected) < std::abs(double(naive) - expected))
		{
			++improved;
		}
	}
	Require(improved > 1000, "timestep negative control no longer discriminates");
}

int main()
{
	try
	{
		CurlCancellation();
		TimestepPrecision();
		using namespace CartesianSubgridFloat;
		TailInfo stateInfo;
		const auto state = State({1, -.25, 1e-100, -1e-41, 0}, stateInfo);
		Require(state == std::vector<float>({1, -.25f, 0, 0, 0}) && stateInfo.terms == 2 &&
			stateInfo.relativeBound < 2e-14, "initial state FTZ bound/normals mishandled");
		Reject([&] { State({1e-40}, stateInfo); });
		Reject([&] { State({1e-25, 1e-39, 1e-39}, stateInfo); });
		Reject([&] { State({1e50}, stateInfo); });
		Reject([&] { State({std::numeric_limits<double>::quiet_NaN()}, stateInfo); });
		TailInfo info;
		const auto normal = Factor({1, 0, .8, .1}, 2, 1, info);
		Require(normal == std::vector<float>({1, 0, .8f, .1f}) && info.terms == 0,
				"normal FP32 coefficients changed");
		const auto tail = Factor({1, 0, .8, 1e-40}, 2, 1, info);
		Require(tail[0] == 1 && tail[2] == .8f && tail[3] == 0 && info.terms == 1,
				"tiny off-diagonal fill or positive diagonal mishandled");
		Require(std::abs(info.l1 / 1e-40 - 1) < 1e-12 &&
				std::abs(info.relativeBound / 1.25e-40 - 1) < 1e-12, "wrong dropped perturbation bound");
		TailInfo minimumInfo;
		const double smallest = std::numeric_limits<float>::min();
		const auto minimum = Factor({1, 0, 1, smallest}, 2, 1, minimumInfo);
		Require(minimum[3] == smallest && minimumInfo.terms == 0, "normal minimum was dropped");
		Reject([&] { Factor({1e-40}, 1, 0, info); });
		Reject([&] { Factor({1, 0, 2e-38, 1e-40}, 2, 1, info); });
		Reject([&] { Factor({1, 0, -1, .1}, 2, 1, info); });
		Reject([&] { Factor({1, 0, 1, std::numeric_limits<double>::infinity()}, 2, 1, info); });
		Reject([&] { Factor({1, 0, 1, 1e50}, 2, 1, info); });
		Reject([&] { Factor({}, 0, 0, info); });
		Reject([&] { Factor({1}, 1, 1, info); });
		Reject([&] { Factor({1}, std::numeric_limits<std::size_t>::max(), 1, info); });
		std::cout << "PASS: explicit factor FTZ tails, unchanged normals/diagonals and perturbation guards" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}
