#pragma once

#include <memory>

#include "config.hpp"
#include "math/math.hpp"
#include "trial_division/multibase_div_tests.hpp"

namespace mbp
{
	constexpr size_t pow_2_16 = 1ull << 16;

	std::vector<bit_array<pow_2_16>> build_popcounts_lookup()
	{
		constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();

		std::vector<bit_array<pow_2_16>> pc_lookup;

		// subtract 2 to normalize outer popcount of 2-48 to idx 0-46 (47 elements)
		pc_lookup.reserve(47);

		// for each outer popcount
		for (size_t outer_pc = 2; outer_pc <= 48; ++outer_pc)
		{
			const size_t shifted_primes_lookup = tiny_primes_lookup >> outer_pc;
			uint64_t* bit_array_ptr = (uint64_t*)pc_lookup.emplace_back().data();

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

		return pc_lookup;
	}
	// [outer popcount 0-46][bitstring inner bits]
	const std::vector<bit_array<pow_2_16>> pc_lookup = build_popcounts_lookup();

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

		std::vector<bit_array<pow_2_16>> gcd_lookup;
		gcd_lookup.reserve(48);

		// for each outer alternating bitsum
		for (int outer_abs = -23; outer_abs <= 24; ++outer_abs)
		{
			uint64_t* bit_array_ptr = (uint64_t*)gcd_lookup.emplace_back().data();

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

		return gcd_lookup;
	}
	// [outer alternating bitsum 0-47][bitstring inner bits]
	const std::vector<bit_array<pow_2_16>> gcd_lookup = build_gcd_lookup();

	template<size_t base, size_t prime>
	std::unique_ptr<std::array<bit_array<pow_2_16>, prime>> build_bit_pattern_filter_for()
	{
		using namespace div_test::detail;

		using lookup_t = decltype(build_bit_pattern_filter_for<base, prime>())::element_type;

		std::unique_ptr<lookup_t> lookup = std::make_unique<lookup_t>();

		// Generate a small lookup mapping a bit index to its remainder
		// This is size 16+1 so we can lookup with rems[bit_index] instead of rems[bit_index % period]
		std::array<size_t, 16 + 1> rems{};
		for (size_t i = 0; i < 16ull + 1; ++i)
		{
			rems[i] = pk::powMod(base, i, prime);
		}

		for (size_t outer_rem = 0; outer_rem < prime; ++outer_rem)
		{
			bit_array<pow_2_16>& bit_array = (*lookup)[outer_rem];

			for (size_t i = 0; i < pow_2_16; ++i)
			{
				// our 16 bits represent bits 1-16, not 0-15
				const size_t bitstring = i << 1;
				size_t sum = outer_rem; // start with the remainder of the outer 48 bits (0 through prime-1)

				// for each bit (1-16) in the bitstring
				for (size_t bit_idx = 1; bit_idx <= 16; ++bit_idx)
				{
					// if the bit is set
					if ((bitstring >> bit_idx) & 1)
					{
						// add that bit's remainder to the sum
						sum += rems[bit_idx];
					}
				}

				// keep the candidate if the remainder is non-zero
				if (sum % prime != 0)
				{
					bit_array.set_bit(i);
				}
			}
		}

		return std::move(lookup);
	}
	// [outer rem 0-(prime - 1)][bitstring inner bits]
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 5>> b3m5_lookup = build_bit_pattern_filter_for<3, 5>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 7>> b3m7_lookup = build_bit_pattern_filter_for<3, 7>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b3m13_lookup = build_bit_pattern_filter_for<3, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 7>> b4m7_lookup = build_bit_pattern_filter_for<4, 7>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b4m13_lookup = build_bit_pattern_filter_for<4, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 17>> b4m17_lookup = build_bit_pattern_filter_for<4, 17>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 7>> b5m7_lookup = build_bit_pattern_filter_for<5, 7>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b5m13_lookup = build_bit_pattern_filter_for<5, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b8m13_lookup = build_bit_pattern_filter_for<8, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b9m13_lookup = build_bit_pattern_filter_for<9, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b10m13_lookup = build_bit_pattern_filter_for<10, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 17>> b13m17_lookup = build_bit_pattern_filter_for<13, 17>();
}
