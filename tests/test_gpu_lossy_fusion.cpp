// SPDX-License-Identifier: GPL-3.0-or-later
#include "FDTD/gpu_lossy_fusion.h"
#include <iostream>

int main()
{
	int failures = 0;
	auto check = [&](bool pass, const char* message) {
		if (!pass) { std::cerr << message << '\n'; ++failures; }
	};
	unsigned int x[] = {0, 0, 1}, y[] = {0, 1, 0}, z[] = {2, 0, 0};
	unsigned int* pos[] = {x, y, z};
	std::vector<GPU_LossyBranch> branches(8, {3, pos, true, false, false});
	check(GPU_CanFuseLossy(true, branches), "ordered eight-branch surface accepted");
	check(!GPU_CanFuseLossy(false, branches), "generic Lorentz extension rejected");
	for (int mode = 0; mode < 4; ++mode)
	{
		auto invalid = branches;
		if (mode == 0) invalid[3].current = true;
		if (mode == 1) invalid[3].lorentz = true;
		if (mode == 2) invalid[3].voltage = false;
		if (mode == 3) invalid[3].count = 2;
		check(!GPU_CanFuseLossy(true, invalid), "incompatible history rejected");
	}
	unsigned int other_z[] = {3, 0, 0};
	unsigned int* other[] = {x, y, other_z};
	branches[4].pos = other;
	check(!GPU_CanFuseLossy(true, branches), "different branch positions rejected");
	branches[4].pos = pos;
	y[1] = 0; z[1] = 2;
	check(!GPU_CanFuseLossy(true, branches), "duplicate positions rejected");
	z[1] = 1;
	check(!GPU_CanFuseLossy(true, branches), "unsorted positions rejected");
	y[1] = 1; z[1] = 0;
	branches.pop_back();
	check(!GPU_CanFuseLossy(true, branches), "unsupported branch count rejected");
	check(GPU_ParseLossyFusion(nullptr) == GPU_LossyFusion::Off, "default off");
	check(GPU_ParseLossyFusion("apply") == GPU_LossyFusion::Apply, "apply mode");
	check(GPU_ParseLossyFusion("all") == GPU_LossyFusion::All, "all mode");
	bool rejected = false;
	try { GPU_ParseLossyFusion("typo"); } catch (const std::runtime_error&) { rejected = true; }
	check(rejected, "unknown policy rejected");
	return failures ? 1 : 0;
}
