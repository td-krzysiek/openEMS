// SPDX-License-Identifier: GPL-3.0-or-later
#include "FDTD/gpu_coeff_sets.h"
#include <cstring>
#include <iostream>
#include <stdexcept>

int main()
{
	int failures = 0;
	auto check = [&](bool pass, const char* message) {
		if (!pass) { std::cerr << message << '\n'; ++failures; }
	};
	check(!GPU_UseRawMainCoefficients(nullptr), "default auto");
	check(!GPU_UseRawMainCoefficients("auto"), "explicit auto");
	check(GPU_UseRawMainCoefficients("raw"), "explicit raw");
	bool rejected = false;
	try { GPU_UseRawMainCoefficients("typo"); } catch (const std::runtime_error&) { rejected = true; }
	check(rejected, "invalid policy rejected");
	for (size_t distinct : {size_t(4), size_t(65537)})
	{
		GPU_CoeffSets sets;
		auto values = [&](size_t i, float* out) {
			for (unsigned int c = 0; c < 12; ++c) out[c] = float(i%distinct) + float(c)/16.f;
		};
		const size_t count = distinct*4;
		check(GPU_FindSets(count, 12, values, sets), "compressible test table");
		check(sets.mode == (distinct < 65537 ? 1u : 2u), "16/32 bit index selection retained");
		for (size_t i = 0; i < count; ++i)
		{
			float raw[12]; values(i, raw);
			check(std::memcmp(raw, sets.table.data()+12*sets.index[i], sizeof(raw)) == 0,
			      "compressed coefficients retain raw FP32 bits");
		}
	}
	return failures ? 1 : 0;
}
