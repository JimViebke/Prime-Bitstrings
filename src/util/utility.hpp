#pragma once

#include <limits>

#include "stdint.h"

namespace mbp::util
{
	long long current_time_in_ms();

	void print_as_bits(const uint64_t n);
	void print_as_bits(const uint32_t n);

	__forceinline bool upper_32_bits_match(const uint64_t a, const uint64_t b)
	{
		return (a >> 32) == (b >> 32);
	}

	namespace util_detail
	{
		template<size_t n_bits>
		consteval auto narrowest_uint_for_n_bits()
		{
			if constexpr (n_bits <= 8) return uint8_t(0);
			else if constexpr (n_bits <= 16) return uint16_t(0);
			else if constexpr (n_bits <= 32) return uint32_t(0);
			else if constexpr (n_bits <= 64) return uint64_t(0);
		}

		template<uint64_t val>
		consteval auto narrowest_uint_for_val()
		{
			if constexpr (val <= std::numeric_limits<uint8_t>::max()) return uint8_t(0);
			else if constexpr (val <= std::numeric_limits<uint16_t>::max()) return uint16_t(0);
			else if constexpr (val <= std::numeric_limits<uint32_t>::max()) return uint32_t(0);
			else return uint64_t(0);
		}
	}

	template<size_t n_bits>
	using narrowest_uint_for_n_bits = decltype(util_detail::narrowest_uint_for_n_bits<n_bits>());

	template<uint64_t val>
	using narrowest_uint_for_val = decltype(util_detail::narrowest_uint_for_val<val>());

	uint64_t hash(uint64_t x);

	template<typename T>
	__forceinline T min(T a, T b)
	{
		return b < a ? b : a;
	}

}
