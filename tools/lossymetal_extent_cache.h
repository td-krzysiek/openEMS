// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LOSSYMETAL_EXTENT_CACHE_H
#define LOSSYMETAL_EXTENT_CACHE_H

#include <algorithm>
#include <vector>

// Sparse known runs of metal grid nodes. The cached tail is the unchanged
// InsideLength() result from the last inside node to the first outside node.
class LossyMetalExtentCache
{
public:
	struct Run
	{
		int first, last, boundary;
		double tail;
	};

	bool Find(int index, int sign, Run& result) const
	{
		const auto& runs = m_Runs[sign > 0];
		auto it = std::lower_bound(runs.begin(), runs.end(), index,
			[](const Run& run, int value) {
				return run.last < value;
			});
		if (it == runs.end() || index < it->first) {
			return false;
		}
		result = *it;
		return true;
	}

	void Remember(Run value, int sign)
	{
		auto& runs = m_Runs[sign > 0];
		// Bounding entry count limits retained metadata, never changes a result.
		if (runs.size() >= 64) {
			runs.clear();
		}
		auto it = std::lower_bound(runs.begin(), runs.end(), value.first,
			[](const Run& run, int first) {
				return run.last < first - 1;
			});
		while (it != runs.end() && it->first <= value.last + 1) {
			if (it->boundary != value.boundary || it->tail != value.tail) {
				++it;
				continue;
			}
			value.first = std::min(value.first, it->first);
			value.last = std::max(value.last, it->last);
			it = runs.erase(it);
		}
		it = std::lower_bound(runs.begin(), runs.end(), value.first,
			[](const Run& run, int first) {
				return run.first < first;
			});
		runs.insert(it, value);
	}

private:
	std::vector<Run> m_Runs[2];
};
#endif
