#pragma once

#include <chrono>

#include "direct.h"

namespace mbp
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

	// Out and in shall be aligned on a multiple of 256 bits.
	// Size can be anything at least as large as 4 * sizeof(__m256i)
	__forceinline void vectorized_copy(__m256i* out,
									   const __m256i* in,
									   size_t size_in_bytes)
	{
		const size_t leftover_bytes = size_in_bytes % (4 * sizeof(__m256i));
		const __m256i* const aligned_end = in + ((size_in_bytes - leftover_bytes) / sizeof(__m256i));

		for (; in != aligned_end; in += 4, out += 4)
		{
			const __m256i ymm0 = _mm256_load_si256(in + 0);
			const __m256i ymm1 = _mm256_load_si256(in + 1);
			const __m256i ymm2 = _mm256_load_si256(in + 2);
			const __m256i ymm3 = _mm256_load_si256(in + 3);
			_mm256_store_si256(out + 0, ymm0);
			_mm256_store_si256(out + 1, ymm1);
			_mm256_store_si256(out + 2, ymm2);
			_mm256_store_si256(out + 3, ymm3);
		}

		// Move both pointers to the (unaligned) end of the data
		in = (__m256i*)(((uint8_t*)in) + leftover_bytes);
		out = (__m256i*)(((uint8_t*)out) + leftover_bytes);
		// Move one "stride" back
		in -= 4;
		out -= 4;
		// Ignore overlap and write the final unaligned stride
		const __m256i ymm0 = _mm256_loadu_si256(in + 0);
		const __m256i ymm1 = _mm256_loadu_si256(in + 1);
		const __m256i ymm2 = _mm256_loadu_si256(in + 2);
		const __m256i ymm3 = _mm256_loadu_si256(in + 3);
		_mm256_storeu_si256(out + 0, ymm0);
		_mm256_storeu_si256(out + 1, ymm1);
		_mm256_storeu_si256(out + 2, ymm2);
		_mm256_storeu_si256(out + 3, ymm3);
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
}
