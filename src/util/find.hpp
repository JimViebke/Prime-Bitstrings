
#include <bitset>

#include "immintrin.h"

namespace mbp::util
{

	const char* find_avx2(const char* b, const char* e, char c)
	{
		const char* i = b;
		__m256i q = _mm256_set1_epi8(c);
		for (; i + 32 < e; i += 32)
		{
			__m256i x = _mm256_lddqu_si256(
				reinterpret_cast<const __m256i*>(i));
			__m256i r = _mm256_cmpeq_epi8(x, q);
			unsigned int z = (unsigned int)_mm256_movemask_epi8(r);
			if (z)
				return i + std::countr_zero(z); // or, i + __builtin_ffs(z) - 1;
		}

		for (; i < e; ++i)
			if (*i == c)
				return i;

		return e;
	}

	// _tzcnt_u32 instead of countr_zero may save an instruction

	const char* find_sse(const char* b, const char* e, char c)
	{
		const char* i = b;
		__m128i q = _mm_set1_epi8(c);
		for (; i + 16 < e; i += 16)
		{
			__m128i x = _mm_lddqu_si128(
				reinterpret_cast<const __m128i*>(i));
			__m128i r = _mm_cmpeq_epi8(x, q);
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
