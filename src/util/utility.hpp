#pragma once

#include <chrono>

#include "direct.h"

namespace mbp::util
{
	inline auto current_time_in_ms()
	{
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	}

	inline auto current_time_in_us()
	{
		return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	}

	void create_folder(const std::string& path)
	{
		std::ignore = _mkdir(path.c_str());
	}

	std::string to_bits(const uint64_t n)
	{
		std::stringstream ss;
		ss << std::bitset<8>(uint8_t(n >> (8 * 7))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 6))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 5))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 4))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 3))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 2))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 1))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 0)));
		return ss.str();
	}
	std::string to_bits(const uint32_t n)
	{
		std::stringstream ss;
		ss << std::bitset<8>(uint8_t(n >> (8 * 3))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 2))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 1))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 0)));
		return ss.str();
	}

	__forceinline bool upper_32_bits_match(const size_t a, const size_t b)
	{
		constexpr size_t mask = size_t(-1) << 32;
		return (a & mask) == (b & mask);
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

	// via https://stackoverflow.com/a/12996028/2924233
	uint64_t hash(uint64_t x)
	{
		x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
		x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
		x = x ^ (x >> 31);
		return x;
	}

	template<typename T>
	__forceinline T min(T a, T b)
	{
		return b < a ? b : a;
	}

	template<typename T>
	__forceinline T min(T a, T b, T c)
	{
		T x = min(a, b);
		return c < x ? c : x;
	}

	template<typename T>
	__forceinline T min(T a, T b, T c, T d)
	{
		T x = min(a, b);
		T y = min(c, d);
		return y < x ? y : x;
	}
}
