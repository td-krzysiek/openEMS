/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_device_curl.h"
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/kernels/subgrid_update.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <string>

static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

static void Check(const std::vector<std::uint32_t> &start, const std::vector<std::uint32_t> &column,
	const std::vector<float> &high, const std::vector<float> &low, std::size_t input, bool arbitrary = false)
{
	const CartesianSubgridDeviceCurl patterns(start, column, high, low, input, arbitrary);
	for (std::size_t row = 0; row + 1 < start.size(); ++row)
	{
		const std::size_t p = patterns.pattern[row], first = patterns.start[p];
		Require(patterns.start[p + 1] - first == start[row + 1] - start[row], "row size changed");
		for (std::size_t k = start[row]; k < start[row + 1]; ++k)
		{
			const std::size_t entry = first + k - start[row];
			Require(patterns.base[row] + patterns.offset[entry] == column[k], "column changed");
			Require(std::memcmp(&patterns.weight[entry], &high[k], sizeof(float)) == 0, "high changed");
			Require(std::memcmp(&patterns.residual[entry], &low[k], sizeof(float)) == 0, "low changed");
		}
	}
	std::vector<float> field(input), expanded(start.size() - 1), compressed(expanded.size());
	std::vector<std::uint32_t> zero(expanded.size(), 0), identity(expanded.size());
	for (std::size_t i = 0; i < field.size(); ++i)
	{
		field[i] = float(std::sin(1.731 * i) + .31 * std::cos(.87 * i));
	}
	for (std::size_t row = 0; row < identity.size(); ++row)
	{
		identity[row] = static_cast<std::uint32_t>(row);
		expanded[row] = compressed[row] = float(std::cos(.451 * row));
	}
	Op_SubgridCurl original{start.data(), column.data(), zero.data(), identity.data(),
		high.data(), low.data(), field.data(), expanded.data(), -.375f, 1e-9f, 0};
	Op_SubgridCurl candidate{patterns.start.data(), patterns.offset.data(), patterns.base.data(),
		patterns.pattern.data(), patterns.weight.data(), patterns.residual.data(), field.data(),
		compressed.data(), -.375f, 1e-9f, 0};
	for (std::uint32_t additive : {0, 1})
	{
		original.additive = candidate.additive = additive;
		for (std::size_t row = 0; row < expanded.size(); ++row)
		{
			original(row);
			candidate(row);
			Require(std::memcmp(&expanded[row], &compressed[row], sizeof(float)) == 0,
				"patterned point operation differs from expanded gather");
		}
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
	throw std::runtime_error("malformed pattern input was accepted");
}

int main()
{
	try
	{
		std::cout << std::setprecision(10);
		for (const std::string name : {"uniform", "local", "graded"})
		{
			CartesianSubgrid::Axes axes;
			for (unsigned int n = 0; n < 3; ++n)
			{
				for (unsigned int i = 0; i <= 12; ++i)
				{
					axes[n].push_back(.001 * i * (1 + .1 * n) + (name == "graded" ? 1e-6 * i * i : 0));
				}
			}
			CartesianSubgrid grid(axes, 8.85418781762e-12, 1.2566370614359173e-6,
				{{{{3, 3, 3}, {9, 9, 9}}, name == "uniform" ? CartesianSubgrid::Index{1, 1, 1} : CartesianSubgrid::Index{2, 2, 2}}});
			const auto data = grid.DeviceData();
			Check(data.curlStart, data.curlColumn, data.curlWeight, data.curlLow, data.electric.size());
			Check(data.transposeStart, data.transposeColumn, data.transposeWeight, data.transposeLow, data.magnetic.size());
			const CartesianSubgridDeviceCurl curl(data.curlStart, data.curlColumn, data.curlWeight,
				data.curlLow, data.electric.size());
			const CartesianSubgridDeviceCurl transpose(data.transposeStart, data.transposeColumn,
				data.transposeWeight, data.transposeLow, data.magnetic.size());
			const std::size_t expanded = (data.curlStart.size() + data.transposeStart.size() +
				data.curlColumn.size() + data.transposeColumn.size()) * sizeof(std::uint32_t) +
				(data.curlWeight.size() + data.curlLow.size() + data.transposeWeight.size() + data.transposeLow.size()) * sizeof(float);
			std::cout << name << " expanded_bytes=" << expanded << " patterned_bytes="
				<< curl.DeviceBytes() + transpose.DeviceBytes() << " curl_patterns=" << curl.start.size() - 1
				<< " transpose_patterns=" << transpose.start.size() - 1 << std::endl;
		}
		// Repeated translated rows, empty rows and signed-zero residual identity.
		const std::vector<std::uint32_t> start{0, 2, 4, 4, 6, 8};
		const std::vector<std::uint32_t> column{0, 2, 3, 5, 0, 2, 8, 10};
		const std::vector<float> high{1, -3, 1, -3, 1, -3, 1, -3};
		const std::vector<float> low{0, .125f, 0, .125f, -0.0f, .125f, 0, .125f};
		Check(start, column, high, low, 11);
		const CartesianSubgridDeviceCurl p(start, column, high, low, 11);
		Require(p.pattern[0] == p.pattern[1] && p.pattern[1] == p.pattern[4], "translated reuse lost");
		Require(p.pattern[0] != p.pattern[3], "signed-zero residual was quantized");
		// Owned-first domain mapping can make columns non-monotone. Preserve
		// exact term order, including repeats, without unsigned offset wrap.
		Check({0, 3, 6}, {5, 1, 3, 15, 11, 13}, {1, 1e8f, -1e8f, 1, 1e8f, -1e8f},
			{0, 0, 0, 0, 0, 0}, 16, true);
		Check({0, 3}, {5, 1, 5}, {1, -3, .5f}, {0, .125f, 0}, 6, true);
		Reject([] { CartesianSubgridDeviceCurl({0, 2}, {2, 0}, {1, 1}, {0, 0}, 2, true); });
		Check({0}, {}, {}, {}, 0);
		Check({0, 0, 0}, {}, {}, {}, 0);
		Reject([] { CartesianSubgridDeviceCurl({}, {}, {}, {}, 0); });
		Reject([] { CartesianSubgridDeviceCurl({1}, {}, {}, {}, 0); });
		Reject([] { CartesianSubgridDeviceCurl({0, 1}, {}, {}, {}, 0); });
		Reject([] { CartesianSubgridDeviceCurl({0, 1, 0}, {}, {}, {}, 0); });
		Reject([] { CartesianSubgridDeviceCurl({0, 1}, {0}, {}, {0}, 1); });
		Reject([] { CartesianSubgridDeviceCurl({0, 1}, {1}, {1}, {0}, 1); });
		Reject([] { CartesianSubgridDeviceCurl({0, 2}, {1, 0}, {1, 1}, {0, 0}, 2); });
		Reject([] { CartesianSubgridDeviceCurl({0, 2}, {0, 0}, {1, 1}, {0, 0}, 1); });
		Reject([] { CartesianSubgridDeviceCurl({0, 1}, {0}, {1e-40f}, {0}, 1); });
		Reject([] { CartesianSubgridDeviceCurl({0, 1}, {0}, {1}, {1e-40f}, 1); });
		Reject([] { CartesianSubgridDeviceCurl({0, 1}, {0}, {std::numeric_limits<float>::infinity()}, {0}, 1); });
		Reject([] { CartesianSubgridDeviceCurl({0, 1}, {0}, {1}, {std::numeric_limits<float>::quiet_NaN()}, 1); });
		Reject([] { CartesianSubgridDeviceCurl({0}, {}, {}, {}, std::size_t(1) << 32); });
		std::mt19937 random(42);
		std::vector<std::uint32_t> starts{0}, columns;
		std::vector<float> highs, lows;
		for (unsigned int row = 0; row < 2000; ++row)
		{
			const unsigned int origin = random() % 90, count = random() % 10;
			for (unsigned int k = 0; k < count; ++k)
			{
				columns.push_back(origin + k);
				highs.push_back(float(int(random() % 20) - 10));
				lows.push_back(float(int(random() % 20) - 10) * 1e-8f);
			}
			starts.push_back(static_cast<std::uint32_t>(columns.size()));
		}
		Check(starts, columns, highs, lows, 100);
		std::cout << "PASS: exact device row patterns and input guards" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}
