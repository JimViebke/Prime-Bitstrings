#pragma once

#include "config.hpp"
#include "math/math.hpp"

namespace mbp
{
	constexpr size_t pow_2_16 = 1ull << 16;

	std::vector<std::vector<bit_array<pow_2_16>>> build_popcounts_lookup()
	{
		constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();

		std::vector<std::vector<bit_array<pow_2_16>>> pc_lookup;
		pc_lookup.reserve(8);

		// for each offset
		for (size_t bit_offset = 0; bit_offset < 8; ++bit_offset)
		{
			auto& offset_lookup = pc_lookup.emplace_back();
			// subtract 2 to normalize outer popcount of 2-48 to idx 0-46 (47 elements)
			offset_lookup.reserve(47);

			// for each outer popcount
			for (size_t outer_pc = 2; outer_pc <= 48; ++outer_pc)
			{
				const size_t shifted_primes_lookup = tiny_primes_lookup >> outer_pc;
				uint64_t* bit_array_ptr = (uint64_t*)offset_lookup.emplace_back().data();

				// inner bits == bits 1 through 16 (not 0 through 15)

				// Start inner_bits at 0-7 instead of 0. This drops off the first 0-7 (unused) values,
				// and has the same effect as shifting this part of the lookup left by 0-7 bits.
				for (size_t inner_bits = bit_offset; inner_bits < pow_2_16; /* increment below */)
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
		}

		return pc_lookup;
	}
	// [offset 0-7][outer pc 0-46][bitstring inner bits]
	const std::vector<std::vector<bit_array<pow_2_16>>> pc_lookup = build_popcounts_lookup();

	std::vector<std::vector<bit_array<pow_2_16>>> build_gcd_lookup()
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

		std::vector<std::vector<bit_array<pow_2_16>>> gcd_lookup;
		gcd_lookup.reserve(8);

		// for each offset
		for (size_t bit_offset = 0; bit_offset < 8; ++bit_offset)
		{
			auto& offset_lookup = gcd_lookup.emplace_back();
			offset_lookup.reserve(48);

			// for each outer alternating bitsum
			for (int outer_abs = -23; outer_abs <= 24; ++outer_abs)
			{
				uint64_t* bit_array_ptr = (uint64_t*)offset_lookup.emplace_back().data();

				// Start inner_bits at 0-7 instead of 0. This drops off the first 0-7 (unused) values,
				// and has the same effect as shifting this part of the lookup left by 0-7 bits.
				for (size_t inner_bits = bit_offset; inner_bits < pow_2_16; /* increment below */)
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
		}

		return gcd_lookup;
	}
	// [offset 0-7][outer abs 0-47][bitstring inner bits]
	const std::vector<std::vector<bit_array<pow_2_16>>> gcd_lookup = build_gcd_lookup();

}
