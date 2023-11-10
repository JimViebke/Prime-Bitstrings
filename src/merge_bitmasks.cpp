
#include <memory>
#include <iostream>

#include "math/math.hpp"
#include "merge_bitmasks.hpp"
#include "sieve.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/types.hpp"

namespace mbp
{
	std::vector<bit_array<pow_2_16>> build_popcounts_lookup()
	{
		constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();

		std::vector<bit_array<pow_2_16>> lookup;

		// subtract 2 to normalize outer popcount of 2-48 to idx 0-46 (47 elements)
		lookup.reserve(47);

		// for each outer popcount
		for (size_t outer_pc = 2; outer_pc <= 48; ++outer_pc)
		{
			const size_t shifted_primes_lookup = tiny_primes_lookup >> outer_pc;
			uint64_t* bit_array_ptr = (uint64_t*)lookup.emplace_back().data();

			// inner bits == bits 1 through 16 (not 0 through 15)

			for (size_t inner_bits = 0; inner_bits < pow_2_16; /* increment below */)
			{
				uint64_t chunk = 0;
				for (size_t j = 0; j < 64; ++j, ++inner_bits)
				{
					const size_t bit = (shifted_primes_lookup >> pop_count(inner_bits)) & 1;
					chunk |= (bit << j);
				}
				*bit_array_ptr = chunk;
				++bit_array_ptr;
			}
		}

		return lookup;
	}

	std::vector<bit_array<pow_2_16>> build_gcd_lookup()
	{
		/*
		We have 48 "outer" bits (47 upper and 1 lower) and 16 "inner" bits (bits 1 through 16; not 0 through 15)

		With 48 outer bits, we have 24 even and 24 odd bits:
			even bits can be 1 through 24	(bottom bit is always set)
			odd  bits can be 0 through 24
			Added together, (even - odd) yields range: -23 to +24

		Normalize outer alternating bitsum from -23,+24 to 0,47 by adding 23
		*/

		// contains a 1 bit at indexes that share a GCD of 1 with a product of primes
		constexpr size_t tiny_gcd_lookup = []() consteval {
			size_t val = 0;
			for (size_t i = 0; i < 32; ++i)
			{
				val |= size_t(gcd(i, size_t(3 * 5 * 7 * 11 * 13)) == 1ull) << i;
			}
			return val;
		}();

		std::vector<bit_array<pow_2_16>> lookup;
		lookup.reserve(48);

		// for each outer alternating bitsum
		for (int outer_abs = -23; outer_abs <= 24; ++outer_abs)
		{
			uint64_t* bit_array_ptr = (uint64_t*)lookup.emplace_back().data();

			for (size_t inner_bits = 0; inner_bits < pow_2_16; /* increment below */)
			{
				// combine the next 64 writes
				uint64_t chunk = 0;
				for (size_t j = 0; j < 64; ++j, ++inner_bits)
				{
					// calculate the alternating bitsum of the inner bits, add outer_abs
					const auto even_pc = pop_count(inner_bits & 0xAAAAAAAAAAAAAAAA); // We're generating a lookup for bits 1-16,
					const auto odd_pc = pop_count(inner_bits & 0x5555555555555555); //  so the even/odd masks are swapped
					const auto alternating_bitsum = (even_pc - odd_pc) + outer_abs;

					const size_t bit = (tiny_gcd_lookup >> abs(alternating_bitsum)) & 1;
					chunk |= (bit << j);
				}

				*bit_array_ptr = chunk;
				++bit_array_ptr;
			}
		}

		return lookup;
	}

}
