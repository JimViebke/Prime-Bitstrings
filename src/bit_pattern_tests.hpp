#pragma once

#include "config.hpp"
#include "math/math.hpp"

namespace mbp
{
	tests_are_inlined size_t* prime_popcount_test(size_t* input,
												  const size_t* const candidates_end)
	{
		constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();

		const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) & 0b11);

		// Alternate between two output ptrs to split the workload into two dependency chains
		size_t* output_a = input;
		size_t* output_b = input + 1;

		for (; input < candidates_end_rounded; input += 4)
		{
			const size_t c0 = *input;
			const size_t pc_c0 = pop_count(c0);
			*output_a = c0;
			if (tiny_primes_lookup & (1ull << pc_c0)) output_a += 2;

			const size_t c1 = *(input + 1);
			const size_t pc_c1 = pop_count(c1);
			*output_b = c1;
			if (tiny_primes_lookup & (1ull << pc_c1)) output_b += 2;

			const size_t c2 = *(input + 2);
			const size_t pc_c2 = pop_count(c2);
			*output_a = c2;
			if (tiny_primes_lookup & (1ull << pc_c2)) output_a += 2;

			const size_t c3 = *(input + 3);
			const size_t pc_c3 = pop_count(c3);
			*output_b = c3;
			if (tiny_primes_lookup & (1ull << pc_c3)) output_b += 2;
		}

		// We've effectively created two interleaved arrays in one.
		// When the output ptrs end up nonadjacent, move data from a->b or b->a
		// until our data is contiguous again.

		// Scenario 1
		while (output_a > output_b + 1)
		{
			output_a -= 2;
			*output_b = *output_a;
			output_b += 2;
		}

		// Scenario 2
		while (output_a < output_b - 1)
		{
			output_b -= 2;
			*output_a = *output_b;
			output_a += 2;
		}

		// keep the smaller of the two
		size_t* output = (output_a < output_b) ? output_a : output_b;

		// handle last few elements
		for (; input < candidates_end; ++input)
		{
			const size_t c0 = *input;
			const size_t pc_c0 = pop_count(c0);
			*output = c0;
			if (tiny_primes_lookup & (1ull << pc_c0)) ++output;
		}

		return output;
	}

	tests_are_inlined size_t* gcd_test(size_t* input,
									   const size_t* const candidates_end)
	{
		constexpr size_t gcd_lookup = build_gcd_lookup();

		const size_t* const candidates_end_rounded = candidates_end - (size_t(candidates_end - input) & 0b11);

		size_t* output = input;

		for (; input < candidates_end_rounded; input += 4)
		{
			const size_t n0 = *input;
			*output = n0;
			const size_t n0_pca = pop_count(n0 & 0xAAAAAAAAAAAAAAAA);
			const size_t n0_pcb = pop_count(n0 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n0_pca + 32 - n0_pcb))) ++output;

			const size_t n1 = *(input + 1);
			*output = n1;
			const size_t n1_pca = pop_count(n1 & 0xAAAAAAAAAAAAAAAA);
			const size_t n1_pcb = pop_count(n1 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n1_pca + 32 - n1_pcb))) ++output;

			const size_t n2 = *(input + 2);
			*output = n2;
			const size_t n2_pca = pop_count(n2 & 0xAAAAAAAAAAAAAAAA);
			const size_t n2_pcb = pop_count(n2 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n2_pca + 32 - n2_pcb))) ++output;

			const size_t n3 = *(input + 3);
			*output = n3;
			const size_t n3_pca = pop_count(n3 & 0xAAAAAAAAAAAAAAAA);
			const size_t n3_pcb = pop_count(n3 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n3_pca + 32 - n3_pcb))) ++output;
		}

		// handle last few elements
		for (; input < candidates_end; ++input)
		{
			const size_t n0 = *input;
			*output = n0;
			const size_t n0_pca = pop_count(n0 & 0xAAAAAAAAAAAAAAAA);
			const size_t n0_pcb = pop_count(n0 & 0x5555555555555555);
			if (gcd_lookup & (1ull << (n0_pca + 32 - n0_pcb))) ++output;
		}

		return output;
	}
}
