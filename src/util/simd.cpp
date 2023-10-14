
#include <iomanip>
#include <iostream>

#include "simd.hpp"

namespace mbp::util
{
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

		size_t sum = _mm256_extract_epi64(vector_sum, 0) +
			_mm256_extract_epi64(vector_sum, 1) +
			_mm256_extract_epi64(vector_sum, 2) +
			_mm256_extract_epi64(vector_sum, 3);

		// sum last few elements
		for (uint8_t* ptr = (uint8_t*)in; ptr < end; ++ptr)
			sum += *ptr;

		return sum;
	}

	namespace detail
	{
		void setw_wrapper(const size_t n) { std::cout << std::setw(n); }

		void print(const char c) { std::cout << c; }
		void print(const char* c) { std::cout << c; }
		void print(const size_t s) { std::cout << s; }
	}

}
