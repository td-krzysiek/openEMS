// SPDX-License-Identifier: GPL-3.0-or-later
#include "tools/lossymetal_extent_cache.h"
#include <iostream>

int main()
{
	LossyMetalExtentCache cache;
	LossyMetalExtentCache::Run result;
	int failures = 0;
	auto check = [&](bool ok, const char* message) {
		if (!ok) {
			std::cerr << message << '\n';
			++failures;
		}
	};
	check(!cache.Find(3, 1, result), "empty cache");
	cache.Remember({3, 7, 7, 0.125}, 1);
	check(cache.Find(3, 1, result) && result.boundary == 7 && result.tail == 0.125,
		"inclusive first node and original boundary tail");
	check(cache.Find(7, 1, result), "inclusive last node");
	check(!cache.Find(2, 1, result) && !cache.Find(8, 1, result), "outside known run");
	check(!cache.Find(5, -1, result), "directions are independent");
	cache.Remember({1, 3, 7, 0.125}, 1);
	check(cache.Find(1, 1, result) && result.first == 1 && result.last == 7,
		"joining a newly scanned prefix to a known tail");
	cache.Remember({10, 13, 13, 0.25}, 1);
	check(!cache.Find(9, 1, result) && cache.Find(10, 1, result) && result.boundary == 13,
		"disjoint bodies retain the intervening gap");
	cache.Remember({1, 7, 1, 0.5}, -1);
	check(cache.Find(7, -1, result) && result.boundary == 1 && result.tail == 0.5,
		"negative direction stores its own boundary");
	for (int i = 0; i < 80; ++i) {
		cache.Remember({100 + 3*i, 101 + 3*i, 101 + 3*i, 0}, 1);
	}
	check(cache.Find(337, 1, result) && result.boundary == 338,
		"bounded cache retains the newest run after eviction");
	check(!cache.Find(3, 1, result) && cache.Find(3, -1, result),
		"eviction only drops reusable metadata in the selected direction");
	return failures ? 1 : 0;
}
