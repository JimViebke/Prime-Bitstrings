#pragma once

#include <stdint.h>

namespace pk
{
	// Prime library
	// Copyright 2014 Peter Knight
	// This code is released under GPLv2 license.

	constexpr uint64_t mulMod(uint64_t a, uint64_t b, uint64_t m)
	{
		/* Calculate ab (mod m)
		**
		** Decompose into sum of powers of 2 (mod m)
		*/
		uint64_t r = 0;
		while (b) {
			if (b & 1) {
				uint64_t r2 = r + a;
				if (r2 < r) r2 -= m; // Correct for an overflow
				r = r2 % m;
			}
			b >>= 1;
			if (b) {
				uint64_t a2 = a + a;
				if (a2 < a) a2 -= m; // Correct for an overflow
				a = a2 % m;
			}
		}
		return r;
	}

	constexpr uint64_t powMod(uint64_t a, uint64_t b, uint64_t m)
	{
		/* Calculate a^b (mod m)
		**
		** Decomposes into product of squares (mod m)
		*/
		uint64_t r = 1;
		while (b) {
			if (b & 1) r = mulMod(r, a, m);
			b >>= 1;
			if (b) a = mulMod(a, a, m);
		}
		return r;
	}
}
