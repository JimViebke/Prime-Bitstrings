#pragma once

#include <immintrin.h>

namespace mbp
{
	using uint128_t = __m128i;
	using uint256_t = __m256i;
}

namespace mbp::util
{
	// Horizontally add a vector. Adapted from Agner Fog's Vector Class Library.
	__forceinline uint64_t vcl_hadd_x(const uint256_t a)
	{
		const uint256_t sum1 = _mm256_sad_epu8(a, _mm256_setzero_si256());
		const uint256_t sum2 = _mm256_shuffle_epi32(sum1, 0b01'01'01'10);
		const uint256_t sum3 = _mm256_add_epi16(sum1, sum2);
		const uint128_t sum4 = _mm256_extracti128_si256(sum3, 1);
		const uint128_t sum5 = _mm_add_epi16(_mm256_castsi256_si128(sum3), sum4);
		return _mm_cvtsi128_si64(sum5);
	}

	// Horizontally add two vectors together. Adapted from the above.
	__forceinline uint64_t vcl_hadd2_x(const uint256_t a, const uint256_t b)
	{
		const uint256_t sum_a = _mm256_sad_epu8(a, _mm256_setzero_si256());
		const uint256_t sum_b = _mm256_sad_epu8(b, _mm256_setzero_si256());
		const uint256_t sum1 = _mm256_add_epi16(sum_a, sum_b);
		const uint256_t sum2 = _mm256_shuffle_epi32(sum1, 0b01'01'01'10);
		const uint256_t sum3 = _mm256_add_epi16(sum1, sum2);
		const uint128_t sum4 = _mm256_extracti128_si256(sum3, 1);
		const uint128_t sum5 = _mm_add_epi16(_mm256_castsi256_si128(sum3), sum4);
		return _mm_cvtsi128_si64(sum5);
	}

	// Convert a 32-bit bitmask to a 32-byte/256-bit bitmask. Via Stack Overflow.
	__forceinline uint256_t expand_bits_to_bytes(uint32_t val)
	{
		const uint256_t shuffle_mask = _mm256_set_epi64x(0x0303030303030303, 0x0202020202020202, 0x0101010101010101, 0x0000000000000000);
		const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

		uint256_t ymm0 = _mm256_set1_epi32(val);
		ymm0 = _mm256_shuffle_epi8(ymm0, shuffle_mask);
		ymm0 = _mm256_andnot_si256(ymm0, and_mask);
		return _mm256_cmpeq_epi8(ymm0, _mm256_setzero_si256());
	}

	// Convert a 16-bit bitmask to a 32-byte/256-bit bitmask. Adapted from the above.
	__forceinline uint256_t expand_16_bits_to_bytes(uint16_t val)
	{
		const uint256_t shuffle_mask = _mm256_set_epi64x(0x0101010101010101, 0x0000000000000000, 0x0101010101010101, 0x0000000000000000);
		const uint256_t and_mask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);

		uint256_t ymm0 = _mm256_set1_epi16(val);
		ymm0 = _mm256_shuffle_epi8(ymm0, shuffle_mask);
		ymm0 = _mm256_andnot_si256(ymm0, and_mask);
		return _mm256_cmpeq_epi8(ymm0, _mm256_setzero_si256());
	}

	namespace detail
	{
		void setw_wrapper(const size_t n);

		void print(const char c);
		void print(const char* c);
		void print(const size_t s);
	}

	template<typename scalar_t, typename vector_t>
	void print_vector_as(const vector_t& data)
	{
		static_assert(std::is_unsigned<scalar_t>());
		
		using namespace detail;

		const scalar_t* ptr = (const scalar_t*)&data;

		print('[');

		for (size_t i = 0; i < sizeof(vector_t) / sizeof(scalar_t); ++i)
		{
			if (i % 8 == 0 && i > 0) print("  "); // padding

			print(' ');
			setw_wrapper(3);
			if (ptr[i] != 0)
				print(size_t(ptr[i]));
			else // write a dot in place of a 0
				print('.');
		}

		print("]\n");
	}

	size_t vector_count_ones(const uint8_t* data, const size_t size);
}
