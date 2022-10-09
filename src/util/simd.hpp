#pragma once

// Hacky way to hide the contents of zmmintrin.h
// For some reason MSVC always includes it from immintrin.h
#ifndef _ZMMINTRIN_H_INCLUDED
#define _ZMMINTRIN_H_INCLUDED
#endif

#include "immintrin.h"
#include "types.hpp"

#include <bitset>
#include <iostream>

namespace mbp::util
{

	// Out and in must be aligned on a multiple of 32 bytes.
	// Size must be at least 4 * 32 bytes.
	__forceinline void vectorized_copy(uint256_t* out,
									   const uint256_t* in,
									   size_t size_in_bytes)
	{
		const size_t leftover_bytes = size_in_bytes % (32ull * 4);
		const uint256_t* const aligned_end = in + ((size_in_bytes - leftover_bytes) / 32);

		for (; in != aligned_end; in += 4, out += 4)
		{
			const uint256_t ymm0 = _mm256_load_si256(in + 0);
			const uint256_t ymm1 = _mm256_load_si256(in + 1);
			const uint256_t ymm2 = _mm256_load_si256(in + 2);
			const uint256_t ymm3 = _mm256_load_si256(in + 3);
			_mm256_store_si256(out + 0, ymm0);
			_mm256_store_si256(out + 1, ymm1);
			_mm256_store_si256(out + 2, ymm2);
			_mm256_store_si256(out + 3, ymm3);
		}

		// Move both pointers to the (unaligned) end of the data
		in = (uint256_t*)(((uint8_t*)in) + leftover_bytes);
		out = (uint256_t*)(((uint8_t*)out) + leftover_bytes);
		// Move one "stride" back
		in -= 4;
		out -= 4;
		// Ignore overlap and write the final unaligned stride
		const uint256_t ymm0 = _mm256_loadu_si256(in + 0);
		const uint256_t ymm1 = _mm256_loadu_si256(in + 1);
		const uint256_t ymm2 = _mm256_loadu_si256(in + 2);
		const uint256_t ymm3 = _mm256_loadu_si256(in + 3);
		_mm256_storeu_si256(out + 0, ymm0);
		_mm256_storeu_si256(out + 1, ymm1);
		_mm256_storeu_si256(out + 2, ymm2);
		_mm256_storeu_si256(out + 3, ymm3);
	}

	// Out and in must be aligned on a multiple of 32 bytes.
	// Size must be at least 4 * 32 bytes.
	template<size_t n_bytes>
	__forceinline void vectorized_copy_n(uint256_t* out,
										 const uint256_t* in)
	{
		constexpr size_t leftover_bytes = n_bytes % (32 * 4);
		const uint256_t* const aligned_end = in + ((n_bytes - leftover_bytes) / 32);

		for (; in != aligned_end; in += 4, out += 4)
		{
			const uint256_t ymm0 = _mm256_load_si256(in + 0);
			const uint256_t ymm1 = _mm256_load_si256(in + 1);
			const uint256_t ymm2 = _mm256_load_si256(in + 2);
			const uint256_t ymm3 = _mm256_load_si256(in + 3);
			_mm256_store_si256(out + 0, ymm0);
			_mm256_store_si256(out + 1, ymm1);
			_mm256_store_si256(out + 2, ymm2);
			_mm256_store_si256(out + 3, ymm3);
		}

		// <32*4 bytes left. Generate instructions for 0-3 aligned copies.

		if constexpr (leftover_bytes >= 32 * 1)
		{
			_mm256_store_si256(out++, _mm256_load_si256(in++));
		}

		if constexpr (leftover_bytes >= 32 * 2)
		{
			_mm256_store_si256(out++, _mm256_load_si256(in++));
		}

		if constexpr (leftover_bytes >= 32 * 3)
		{
			_mm256_store_si256(out++, _mm256_load_si256(in++));
		}

		// <32 bytes left. Generate instructions for 0-1 unaligned copies.

		constexpr size_t adjust = 32 - (leftover_bytes % 32);
		if constexpr (adjust > 0)
		{
			in = (uint256_t*)(((uint8_t*)in) - adjust);
			out = (uint256_t*)(((uint8_t*)out) - adjust);

			_mm256_storeu_si256(out, _mm256_loadu_si256(in));
		}
	}

	// Horizontally add a vector
	__forceinline uint64_t vcl_hadd_x(const uint256_t a)
	{
		const uint256_t sum1 = _mm256_sad_epu8(a, _mm256_setzero_si256());
		const uint256_t sum2 = _mm256_shuffle_epi32(sum1, 0b01'01'01'10);
		const uint256_t sum3 = _mm256_add_epi16(sum1, sum2);
		const uint128_t sum4 = _mm256_extracti128_si256(sum3, 1);
		const uint128_t sum5 = _mm_add_epi16(_mm256_castsi256_si128(sum3), sum4);
		return _mm_cvtsi128_si64(sum5);
	}

	// Horizontally add two vectors together
	__forceinline uint64_t vcl_hadd2_x(const uint256_t a, const uint256_t b)
	{
		const uint256_t sum_a = _mm256_sad_epu8(a, _mm256_setzero_si256()); // add adjacent bytes
		const uint256_t sum_b = _mm256_sad_epu8(b, _mm256_setzero_si256());
		const uint256_t sum1 = _mm256_add_epi16(sum_a, sum_b);
		const uint256_t sum2 = _mm256_shuffle_epi32(sum1, 0b01'01'01'10);
		const uint256_t sum3 = _mm256_add_epi16(sum1, sum2);
		const uint128_t sum4 = _mm256_extracti128_si256(sum3, 1);
		const uint128_t sum5 = _mm_add_epi16(_mm256_castsi256_si128(sum3), sum4);
		return _mm_cvtsi128_si64(sum5);
	}

	// Convert a 32-bit bitmask to a 32-byte/256-bit bitmask
	__forceinline uint256_t expand_bits_to_bytes(uint32_t x)
	{
		uint256_t xbcast = _mm256_set1_epi32(x);    // we only use the low 32bits of each lane, but this is fine with AVX2

		// Each byte gets the source byte containing the corresponding bit
		uint256_t shufmask = _mm256_set_epi64x(0x0303030303030303, 0x0202020202020202, 0x0101010101010101, 0x0000000000000000);
		uint256_t shuf = _mm256_shuffle_epi8(xbcast, shufmask);

		uint256_t andmask = _mm256_set1_epi64x(0x80'40'20'10'08'04'02'01);  // every 8 bits -> 8 bytes, pattern repeats.
		uint256_t isolated_inverted = _mm256_andnot_si256(shuf, andmask);

		// this is the extra step: compare each byte == 0 to produce 0 or -1
		return _mm256_cmpeq_epi8(isolated_inverted, _mm256_setzero_si256());
	}

	void print_vector256_u8(const uint256_t& data)
	{
		for (size_t i = 0; i < 32; ++i)
		{
			if (i % 8 == 0) // padding
				std::cout << "  ";

			// write a dot in place of a 0
			if (data.m256i_u8[i] != 0)
				std::cout << std::setw(3) << size_t(data.m256i_u8[i]);
			else
				std::cout << std::setw(3) << '.';
		}

		std::cout << '\n';
	}
	void print_vector256_u32(const uint256_t& data)
	{
		for (size_t i = 0; i < 8; ++i)
		{
			if (i % 4 == 0) // padding
				std::cout << "  ";

			// write a dot in place of a 0
			if (data.m256i_u32[i] != 0)
				std::cout << ' ' << std::setw(4) << size_t(data.m256i_u32[i]);
			else
				std::cout << ' ' << std::setw(4) << '.';
		}

		std::cout << '\n';
	}
	template<typename scalar_t, typename vector_t>
	void print_vector_as(const vector_t& data)
	{
		static_assert(std::is_unsigned<scalar_t>());

		const scalar_t* ptr = (const scalar_t*)&data;

		std::cout << '[';

		for (size_t i = 0; i < sizeof(vector_t) / sizeof(scalar_t); ++i)
		{
			if (i % 8 == 0 && i > 0) std::cout << "  "; // padding

			std::cout << ' ' << std::setw(3);
			if (ptr[i] != 0)
				std::cout << size_t(ptr[i]);
			else // write a dot in place of a 0
				std::cout << '.';
		}

		std::cout << "]\n";
	}



	size_t vector_count_ones(const uint8_t* data, const size_t size)
	{
		const uint8_t* const end = data + size;
		const uint256_t* const aligned_end = (uint256_t*)(end - (size % 32));
		const uint256_t* in = (uint256_t*)data;

		uint256_t vector_sum = _mm256_setzero_si256();

		while (in != aligned_end)
		{
			// vertically sum in 32-element blocks
			uint256_t inner_sum = _mm256_setzero_si256();
			for (size_t i = 0; i < 255 && in != aligned_end; ++i, ++in)
				inner_sum = _mm256_adds_epu8(inner_sum, *in);

			// horizontally sum to 4x uint64s
			inner_sum = _mm256_sad_epu8(inner_sum, _mm256_setzero_si256());
			// add to running total
			vector_sum = _mm256_add_epi64(vector_sum, inner_sum);
		}

		size_t sum = vector_sum.m256i_u64[0] +
			vector_sum.m256i_u64[1] +
			vector_sum.m256i_u64[2] +
			vector_sum.m256i_u64[3];

		// sum last few elements
		for (uint8_t* ptr = (uint8_t*)in; ptr < end; ++ptr)
			sum += *ptr;

		return sum;
	}



	const char* find_avx2(const char* b, const char* e, char c)
	{
		const char* i = b;
		uint256_t q = _mm256_set1_epi8(c);
		for (; i + 32 < e; i += 32)
		{
			uint256_t x = _mm256_lddqu_si256(
				reinterpret_cast<const uint256_t*>(i));
			uint256_t r = _mm256_cmpeq_epi8(x, q);
			unsigned int z = (unsigned int)_mm256_movemask_epi8(r);
			if (z)
				return i + std::countr_zero(z); // or, i + __builtin_ffs(z) - 1;
		}

		for (; i < e; ++i)
			if (*i == c)
				return i;

		return e;
	}



	const char* find_sse(const char* b, const char* e, char c)
	{
		const char* i = b;
		uint128_t q = _mm_set1_epi8(c);
		for (; i + 16 < e; i += 16)
		{
			uint128_t x = _mm_lddqu_si128(
				reinterpret_cast<const uint128_t*>(i));
			uint128_t r = _mm_cmpeq_epi8(x, q);
			unsigned int z = _mm_movemask_epi8(r);
			if (z)
				return i + std::countr_zero(z); // or, i + __builtin_ffs(z) - 1;
		}

		for (; i < e; ++i)
			if (*i == c)
				return i;

		return e;
	}

}
