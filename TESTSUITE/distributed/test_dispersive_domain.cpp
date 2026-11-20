// SPDX-License-Identifier: GPL-3.0-or-later
#include "FDTD/operator.h"
#include "FDTD/gpu_dispersive_domain.h"
#include <stdexcept>

static void Require(bool condition)
{
	if (!condition)
	{
		throw std::runtime_error("dispersive domain packing regression");
	}
}

int main()
{
	unsigned int x[] = { 0, 3, 4, 7, 8 };
	unsigned int y[] = { 1, 2, 1, 2, 1 };
	unsigned int z[] = { 2, 1, 2, 1, 2 };
	unsigned int* coordinates[] = { x, y, z };
	float c0[] = { 1, 2, 3, 4, 5 };
	float c1[] = { 6, 7, 8, 9, 10 };
	float c2[] = { 11, 12, 13, 14, 15 };
	float* coefficients[] = { c0, c1, c2 };
	GPU_DomainLayout whole;
	GPU_DispersiveDomain original(whole, 5, coordinates, 3, 4);
	Require(original.source == std::vector<size_t>({ 0, 1, 2, 3, 4 }));
	Require(original.position[2] == (4 * 3 + 1) * 4 + 2);
	Require(original.Coefficients(coefficients) == std::vector<float>({ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }));
	GPU_DomainLayout slice;
	slice.enabled = true;
	slice.begin = 4;
	slice.end = 8;
	GPU_DispersiveDomain local(slice, 5, coordinates, 3, 4);
	Require(local.source == std::vector<size_t>({ 2, 3 }));
	Require(local.position == std::vector<size_t>({ (1 * 3 + 1) * 4 + 2, (4 * 3 + 2) * 4 + 1 }));
	Require(local.Coefficients(coefficients) == std::vector<float>({ 3, 4, 8, 9, 13, 14 }));
	Require(local.Coefficients<float>(nullptr) == std::vector<float>(6, 0));
	slice.begin = 1;
	slice.end = 3;
	GPU_DispersiveDomain empty(slice, 5, coordinates, 3, 4);
	Require(empty.source.empty() && empty.position.empty());
	Require(empty.Coefficients(coefficients).empty());
	return 0;
}
