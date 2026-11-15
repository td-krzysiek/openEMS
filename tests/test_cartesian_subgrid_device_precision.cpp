/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_device_precision.h"
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

int main()
{
	try
	{
		using namespace CartesianSubgridFloat;
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
