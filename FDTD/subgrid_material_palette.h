/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_MATERIAL_PALETTE_H
#define CARTESIAN_SUBGRID_MATERIAL_PALETTE_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

//! Exact immutable leaf media. Repeated values share their original double bits.
//! A bounded setup dictionary falls back to dense storage for arbitrary media.
template <class Medium> class CartesianSubgridMaterialPalette
{
  public:
	static constexpr std::size_t PaletteLimit = 4096;

	void Reserve(std::size_t count)
	{
		if (m_Size || m_Sealed || m_Expected)
		{
			throw std::logic_error("CartesianSubgrid: material storage reservation is setup-only");
		}
		m_Expected = count;
	}

	void Append(const Medium &value)
	{
		if (m_Sealed || m_Size >= m_Expected)
		{
			throw std::logic_error("CartesianSubgrid: material storage is sealed or full");
		}
		if (m_DenseMode)
		{
			m_Dense.push_back(value);
			++m_Size;
			return;
		}
		const auto key = Bits(value);
		auto found = m_Dictionary.find(key);
		if (found == m_Dictionary.end())
		{
			if (m_Palette.size() == PaletteLimit)
			{
				// A fully varying sampled medium must not retain an unbounded map.
				// Release the temporary IDs immediately after exact expansion.
				m_Dense.reserve(m_Expected);
				for (std::size_t i = 0; i < m_Size; ++i)
				{
					m_Dense.push_back((*this)[i]);
				}
				std::vector<std::uint32_t>().swap(m_IDs);
				std::vector<Medium>().swap(m_Palette);
				std::map<Key, std::uint32_t>().swap(m_Dictionary);
				m_DenseMode = true;
				m_Dense.push_back(value);
				++m_Size;
				return;
			}
			const auto id = static_cast<std::uint32_t>(m_Palette.size());
			m_Palette.push_back(value);
			found = m_Dictionary.emplace(key, id).first;
			if (id == 1)
			{
				// Homogeneous grids need neither one ID per leaf nor dense media.
				m_IDs.reserve(m_Expected);
				m_IDs.assign(m_Size, 0);
			}
		}
		if (!m_IDs.empty() || m_Palette.size() > 1)
		{
			m_IDs.push_back(found->second);
		}
		++m_Size;
	}

	void Finish()
	{
		if (m_Sealed || m_Size != m_Expected)
		{
			throw std::logic_error("CartesianSubgrid: incomplete or already sealed material storage");
		}
		std::map<Key, std::uint32_t>().swap(m_Dictionary);
		m_Sealed = true;
	}

	const Medium &operator[](std::size_t leaf) const
	{
		return m_DenseMode ? m_Dense[leaf] : m_Palette[m_IDs.empty() ? 0 : m_IDs[leaf]];
	}
	const Medium &at(std::size_t leaf) const
	{
		if (leaf >= m_Size)
		{
			throw std::out_of_range("CartesianSubgrid: material leaf is out of range");
		}
		return (*this)[leaf];
	}
	std::size_t Size() const { return m_Size; }
	std::size_t UniqueCount() const { return m_DenseMode ? 0 : m_Palette.size(); }
	bool Dense() const { return m_DenseMode; }
	std::size_t StorageBytes() const
	{
		return m_IDs.capacity() * sizeof(std::uint32_t) +
			(m_Palette.capacity() + m_Dense.capacity()) * sizeof(Medium);
	}

  private:
	using Key = std::array<std::uint64_t, 8>;
	static Key Bits(const Medium &value)
	{
		Key key;
		static_assert(sizeof(key[0]) == sizeof(value.epsilon), "64-bit material epsilon required");
		static_assert(sizeof(key[1]) == sizeof(value.permeability), "64-bit material mu required");
		std::memcpy(&key[0], &value.epsilon, sizeof(key[0]));
		std::memcpy(&key[1], &value.permeability, sizeof(key[1]));
		for (unsigned int n = 0; n < 3; ++n)
		{
			std::memcpy(&key[2 + n], &value.epsilonDiagonal[n], sizeof(key[0]));
			std::memcpy(&key[5 + n], &value.permeabilityDiagonal[n], sizeof(key[0]));
		}
		return key;
	}
	std::vector<Medium> m_Palette, m_Dense;
	std::vector<std::uint32_t> m_IDs;
	std::map<Key, std::uint32_t> m_Dictionary;
	std::size_t m_Size = 0, m_Expected = 0;
	bool m_Sealed = false, m_DenseMode = false;
};
#endif
