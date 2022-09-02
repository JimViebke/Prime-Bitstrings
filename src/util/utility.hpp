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

	// Horizontally add two vectors together
	__forceinline auto vcl_hadd2_x(const __m256i a, const __m256i b)
	{
		const __m256i sum_a = _mm256_sad_epu8(a, _mm256_setzero_si256());
		const __m256i sum_b = _mm256_sad_epu8(b, _mm256_setzero_si256());
		const __m256i sum1 = _mm256_add_epi16(sum_a, sum_b);
		const __m256i sum2 = _mm256_shuffle_epi32(sum1, 2);
		const __m256i sum3 = _mm256_add_epi16(sum1, sum2);
		const __m128i sum4 = _mm256_extracti128_si256(sum3, 1);
		const __m128i sum5 = _mm_add_epi16(_mm256_castsi256_si128(sum3), sum4);
		return _mm_cvtsi128_si32(sum5);
	}

	// Convert a 32-bit bitmask to a 32-byte/256-bit bitmask
	__forceinline __m256i expand_bits_to_bytes(uint32_t x)
	{
		__m256i xbcast = _mm256_set1_epi32(x);    // we only use the low 32bits of each lane, but this is fine with AVX2

		// Each byte gets the source byte containing the corresponding bit
		__m256i shufmask = _mm256_set_epi64x(0x0303030303030303, 0x0202020202020202, 0x0101010101010101, 0x0000000000000000);
		__m256i shuf = _mm256_shuffle_epi8(xbcast, shufmask);

		__m256i andmask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);  // every 8 bits -> 8 bytes, pattern repeats.
		__m256i isolated_inverted = _mm256_andnot_si256(shuf, andmask);

		// this is the extra step: compare each byte == 0 to produce 0 or -1
		return _mm256_cmpeq_epi8(isolated_inverted, _mm256_setzero_si256());
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
