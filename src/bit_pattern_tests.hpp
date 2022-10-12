#pragma once

#include "config.hpp"
#include "math/math.hpp"

namespace mbp
{
	tests_are_inlined size_t* gcd_test(size_t* input,
									   const size_t* const candidates_end)
	{
		/*
		A candidate's alternating bit sum must share a greatest common divisor
		with a product of primes of 1.

		That is, GCD(diff, product_of_primes) must equal 1, given:
				diff = absolute value of ((# even bits) - (# odd bits))
				product_of_primes = 3 * 5 * 7 * 11 * 13

		The absolute value of the alternating bitsum of a 64-bit integer is in
		range 0 - 32 inclusive. See build_gcd_lookup() for how we implement this.
		*/

		constexpr size_t gcd_lookup = build_gcd_lookup();

		const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) & 0b11);

		size_t* output = input;

		for (; input < candidates_end_rounded; input += 4)
		{
			// To avoid a branch, always write each candidate to the output, then
			// use a bittest + cmove to conditionally increment the output pointer
			const size_t c0 = *input;
			*output = c0;
			const size_t c0_pca = pop_count(c0 & 0xAAAAAAAAAAAAAAAA);
			const size_t c0_pcb = pop_count(c0 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (c0_pca + 32 - c0_pcb))) ++output;

			const size_t c1 = *(input + 1);
			*output = c1;
			const size_t c1_pca = pop_count(c1 & 0xAAAAAAAAAAAAAAAA);
			const size_t c1_pcb = pop_count(c1 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (c1_pca + 32 - c1_pcb))) ++output;

			const size_t c2 = *(input + 2);
			*output = c2;
			const size_t c2_pca = pop_count(c2 & 0xAAAAAAAAAAAAAAAA);
			const size_t c2_pcb = pop_count(c2 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (c2_pca + 32 - c2_pcb))) ++output;

			const size_t c3 = *(input + 3);
			*output = c3;
			const size_t c3_pca = pop_count(c3 & 0xAAAAAAAAAAAAAAAA);
			const size_t c3_pcb = pop_count(c3 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (c3_pca + 32 - c3_pcb))) ++output;
		}

		// handle last few elements
		for (; input < candidates_end; ++input)
		{
			const size_t c = *input;
			*output = c;
			const size_t c_pca = pop_count(c & 0xAAAAAAAAAAAAAAAA);
			const size_t c_pcb = pop_count(c & 0x5555555555555555);
			if (gcd_lookup & (1ull << (c_pca + 32 - c_pcb))) ++output;
		}

		return output;
	}
}
