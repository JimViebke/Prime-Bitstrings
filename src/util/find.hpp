
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

}
