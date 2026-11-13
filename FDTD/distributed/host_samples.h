// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_DISTRIBUTED_HOST_SAMPLES_H
#define OPENEMS_DISTRIBUTED_HOST_SAMPLES_H

#include <memory>
#include <stdexcept>
#include <vector>

namespace openems
{
namespace distributed
{

// Contiguous FP32 wire storage. Borrowing is explicit at both endpoints and
// keeps the native staging owner alive; owned vectors are never implicitly lent.
class HostSamples
{
  public:
	explicit HostSamples(size_t count, float value = 0)
		: m_Owned(count, value), m_Data(m_Owned.data()), m_Size(count)
	{
	}
	HostSamples(size_t count, float *data, std::shared_ptr<void> owner)
		: m_Data(data), m_Size(count), m_Owner(std::move(owner))
	{
		if (!data || !m_Owner || !count)
		{
			throw std::invalid_argument("invalid borrowed host halo storage");
		}
	}
	HostSamples(const HostSamples &) = delete;
	HostSamples &operator=(const HostSamples &) = delete;
	size_t size() const
	{
		return m_Size;
	}
	float *data()
	{
		return m_Data;
	}
	const float *data() const
	{
		return m_Data;
	}
	float *begin()
	{
		return m_Data;
	}
	float *end()
	{
		return m_Data + m_Size;
	}
	float &at(size_t index)
	{
		if (index >= m_Size)
		{
			throw std::out_of_range("host halo sample index");
		}
		return m_Data[index];
	}
	const std::shared_ptr<void> &Owner() const
	{
		return m_Owner;
	}
	void AllowBorrowing()
	{
		m_AllowBorrowing = true;
	}
	void RestoreOwned()
	{
		if (m_Owned.size() != m_Size)
		{
			throw std::logic_error("external-only halo cannot receive an owned copy");
		}
		m_Data = m_Owned.data();
		m_Owner.reset();
	}
	bool Adopt(const HostSamples &source)
	{
		if (!m_AllowBorrowing || !source.m_Owner)
		{
			return false;
		}
		if (source.size() != size())
		{
			throw std::invalid_argument("borrowed host halo size mismatch");
		}
		m_Data = source.m_Data;
		m_Owner = source.m_Owner;
		return true;
	}

  private:
	std::vector<float> m_Owned;
	float *m_Data;
	size_t m_Size;
	std::shared_ptr<void> m_Owner;
	bool m_AllowBorrowing = false;
};

} // namespace distributed
} // namespace openems
#endif
