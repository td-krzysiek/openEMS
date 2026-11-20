/* SPDX-License-Identifier: GPL-3.0-or-later */
#pragma once

#include <type_traits>
#include <utility>

template <class T, class = void>
struct SubgridHasNativeSurfaceMetadata : std::false_type
{
};

template <class T>
struct SubgridHasNativeSurfaceMetadata<T, std::void_t<
	decltype(std::declval<T &>().GetConductivity()),
	decltype(std::declval<T &>().GetRq()),
	decltype(std::declval<T &>().GetRoughnessModel())>> : std::true_type
{
};
