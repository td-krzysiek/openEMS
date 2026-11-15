// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_GPU_LOSSY_FUSION_H
#define OPENEMS_GPU_LOSSY_FUSION_H

#include <array>
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <vector>

enum class GPU_LossyFusion { Off, Apply, All };

inline GPU_LossyFusion GPU_ParseLossyFusion(const char* value)
{
	if (!value || std::strcmp(value, "off") == 0)
		return GPU_LossyFusion::Off;
	if (std::strcmp(value, "apply") == 0)
		return GPU_LossyFusion::Apply;
	if (std::strcmp(value, "all") == 0)
		return GPU_LossyFusion::All;
	throw std::runtime_error("OPENEMS_GPU_LOSSY_FUSION must be off, apply or all");
}

struct GPU_LossyBranch
{
	size_t count;
	unsigned int** pos;
	bool voltage, current, lorentz;
};

// LossyMetal emits one strictly ordered position list, shared by all eight
// voltage-only branches. Prove that contract before sharing field accesses.
inline bool GPU_CanFuseLossy(bool lossy, const std::vector<GPU_LossyBranch>& branches)
{
	if (!lossy || branches.size() != 8 || branches[0].count == 0)
		return false;
	const size_t count = branches[0].count;
	for (const auto& branch : branches)
	{
		if (branch.count != count || !branch.voltage || branch.current || branch.lorentz || !branch.pos)
			return false;
		for (unsigned int n = 0; n < 3; ++n)
			if (!branch.pos[n] || std::memcmp(branch.pos[n], branches[0].pos[n], count*sizeof(unsigned int)))
				return false;
	}
	for (size_t i = 1; i < count; ++i)
	{
		std::array<unsigned int, 3> previous, current;
		for (unsigned int n = 0; n < 3; ++n)
		{
			previous[n] = branches[0].pos[n][i-1];
			current[n] = branches[0].pos[n][i];
		}
		if (!(previous < current))
			return false;
	}
	return true;
}

#endif
