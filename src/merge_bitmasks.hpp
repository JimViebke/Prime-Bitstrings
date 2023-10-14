#pragma once

#include <memory>

#include "config.hpp"
#include "math/math.hpp"
#include "sieve.hpp"
#include "trial_division/multibase_div_tests.hpp"

namespace mbp
{
	constexpr size_t pow_2_16 = 1ull << 16;

	static const sieve_container static_sieve = prime_sieve::generate_static_sieve();

	std::vector<bit_array<pow_2_16>> build_popcounts_lookup();
	// [outer popcount 0-46][bitstring inner bits]
	static const std::vector<bit_array<pow_2_16>> pc_lookup = build_popcounts_lookup();

	std::vector<bit_array<pow_2_16>> build_gcd_lookup();
	// [outer alternating bitsum 0-47][bitstring inner bits]
	static const std::vector<bit_array<pow_2_16>> gcd_lookup = build_gcd_lookup();

	template<size_t base, size_t prime>
	std::unique_ptr<std::array<bit_array<pow_2_16>, prime>> build_bitmask_for()
	{
		using namespace div_test::detail;

		using lookup_t = decltype(build_bitmask_for<base, prime>())::element_type;

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
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 5>> b3m5_lookup = build_bitmask_for<3, 5>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 7>> b3m7_lookup = build_bitmask_for<3, 7>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b3m13_lookup = build_bitmask_for<3, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 7>> b4m7_lookup = build_bitmask_for<4, 7>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b4m13_lookup = build_bitmask_for<4, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 17>> b4m17_lookup = build_bitmask_for<4, 17>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 7>> b5m7_lookup = build_bitmask_for<5, 7>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b5m13_lookup = build_bitmask_for<5, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b8m13_lookup = build_bitmask_for<8, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b9m13_lookup = build_bitmask_for<9, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 13>> b10m13_lookup = build_bitmask_for<10, 13>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 17>> b13m17_lookup = build_bitmask_for<13, 17>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 11>> b3m11_lookup = build_bitmask_for<3, 11>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 11>> b4m11_lookup = build_bitmask_for<4, 11>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 11>> b5m11_lookup = build_bitmask_for<5, 11>();
	static const std::unique_ptr<std::array<bit_array<pow_2_16>, 11>> b9m11_lookup = build_bitmask_for<9, 11>();



	// 64 bits == (47 high bits) + (16 "inner" bits) + (1 low bit (always set))
	constexpr uint64_t bits_1_16_mask = 0xFFFFull << 1;
	constexpr uint64_t outer_48_bits_mask = ~bits_1_16_mask;

	__forceinline size_t pc_lookup_idx(const uint64_t number)
	{
		return pop_count(number & outer_48_bits_mask) - 2; // -2 to normalize popcount 2-48 to idx 0-46
	}

	__forceinline size_t gcd_lookup_idx(const uint64_t number)
	{
		constexpr uint64_t even_mask = outer_48_bits_mask & 0x5555555555555555;
		constexpr uint64_t odd_mask = outer_48_bits_mask & 0xAAAAAAAAAAAAAAAA;

		const auto outer_even_pc = pop_count(number & even_mask);
		const auto outer_odd_pc = pop_count(number & odd_mask);
		return outer_even_pc - outer_odd_pc + 23; // +23 to normalize -23,24 to 0,47
	}

	// out_ptr is used to set in_ptr
	template<size_t pass>
	__forceinline void set_lookup_ptrs(const uint8_t* out_ptr,
									   const size_t sieve_offset,
									   const uint64_t next_number,
									   const uint8_t*& in_ptr,
									   const uint8_t*& lookup_1_ptr,
									   const uint8_t*& lookup_2_ptr,
									   const uint8_t*& lookup_3_ptr,
									   const uint8_t*& lookup_4_ptr,
									   const uint8_t*& lookup_5_ptr,
									   const uint8_t*& lookup_6_ptr,
									   const uint8_t*& lookup_7_ptr)
	{
		using namespace div_test;
		using namespace div_test::detail;

		const uint64_t bits_1_16 = (next_number & bits_1_16_mask) >> 1;
		const size_t byte_offset = bits_1_16 / 8;

		if constexpr (pass == 1)
		{
			const size_t pc_idx = pc_lookup_idx(next_number);
			const size_t gcd_idx = gcd_lookup_idx(next_number);

			in_ptr = static_sieve.data() + sieve_offset; // the offset into the sieve is also the offset into the static sieve
			lookup_1_ptr = pc_lookup[pc_idx].data() + byte_offset;
			lookup_2_ptr = gcd_lookup[gcd_idx].data() + byte_offset;

			constexpr uint64_t bitmask = bitmask_for<3, 5>::val;
			static_assert(bitmask == bitmask_for<4, 17>::val);
			static_assert(bitmask == bitmask_for<5, 13>::val);
			static_assert(bitmask == bitmask_for<8, 13>::val);
			static_assert(bitmask == bitmask_for<13, 17>::val);
			static_assert(period_of<bitmask>() == 4);

			const size_t pc_0 = pop_count(next_number & ((bitmask << 0) & outer_48_bits_mask));
			size_t b3_sum = pc_0;
			size_t b4_sum = pc_0;
			size_t b5_sum = pc_0;
			size_t b8_sum = pc_0;
			size_t b13_sum = pc_0;
			const size_t pc_1 = pop_count(next_number & ((bitmask << 1) & outer_48_bits_mask));
			b3_sum += pc_1 * pow_mod<3, 1, 5>::rem;
			b4_sum += pc_1 * pow_mod<4, 1, 17>::rem;
			b5_sum += pc_1 * pow_mod<5, 1, 13>::rem;
			b8_sum += pc_1 * pow_mod<8, 1, 13>::rem;
			b13_sum += pc_1 * pow_mod<13, 1, 17>::rem;
			const size_t pc_2 = pop_count(next_number & ((bitmask << 2) & outer_48_bits_mask));
			b3_sum += pc_2 * pow_mod<3, 2, 5>::rem;
			b4_sum += pc_2 * pow_mod<4, 2, 17>::rem;
			b5_sum += pc_2 * pow_mod<5, 2, 13>::rem;
			b8_sum += pc_2 * pow_mod<8, 2, 13>::rem;
			b13_sum += pc_2 * pow_mod<13, 2, 17>::rem;
			const size_t pc_3 = pop_count(next_number & ((bitmask << 3) & outer_48_bits_mask));
			b3_sum += pc_3 * pow_mod<3, 3, 5>::rem;
			b4_sum += pc_3 * pow_mod<4, 3, 17>::rem;
			b5_sum += pc_3 * pow_mod<5, 3, 13>::rem;
			b8_sum += pc_3 * pow_mod<8, 3, 13>::rem;
			b13_sum += pc_3 * pow_mod<13, 3, 17>::rem;

			//__assume(b3_sum <= get_sum_of_rems<5, in_base<3>>(outer_48_bits_mask));
			//__assume(b4_sum <= get_sum_of_rems<17, in_base<4>>(outer_48_bits_mask));
			//__assume(b5_sum <= get_sum_of_rems<13, in_base<5>>(outer_48_bits_mask));
			//__assume(b8_sum <= get_sum_of_rems<13, in_base<8>>(outer_48_bits_mask));
			//__assume(b13_sum <= get_sum_of_rems<17, in_base<13>>(outer_48_bits_mask));

			lookup_3_ptr = (*b3m5_lookup)[b3_sum % 5].data() + byte_offset;
			lookup_4_ptr = (*b4m17_lookup)[b4_sum % 17].data() + byte_offset;
			lookup_5_ptr = (*b5m13_lookup)[b5_sum % 13].data() + byte_offset;
			lookup_6_ptr = (*b8m13_lookup)[b8_sum % 13].data() + byte_offset;
			lookup_7_ptr = (*b13m17_lookup)[b13_sum % 17].data() + byte_offset;
		}
		else // all passes except for pass 1
		{
			in_ptr = out_ptr; // wherever we're writing to in the sieve is also the location we have to read from
		}

		if constexpr (pass == 2)
		{
			{
				constexpr uint64_t bitmask = bitmask_for<3, 13>::val;
				static_assert(bitmask == bitmask_for<4, 7>::val);
				static_assert(bitmask == bitmask_for<9, 13>::val);
				static_assert(period_of<bitmask>() == 3);

				const size_t pc_0 = pop_count(next_number & ((bitmask << 0) & outer_48_bits_mask));
				size_t b3_sum = pc_0;
				size_t b4_sum = pc_0;
				size_t b9_sum = pc_0;
				const size_t pc_1 = pop_count(next_number & ((bitmask << 1) & outer_48_bits_mask));
				b3_sum += pc_1 * pow_mod<3, 1, 13 >::rem;
				b4_sum += pc_1 * pow_mod<4, 1, 7>::rem;
				b9_sum += pc_1 * pow_mod<9, 1, 13>::rem;
				const size_t pc_2 = pop_count(next_number & ((bitmask << 2) & outer_48_bits_mask));
				b3_sum += pc_2 * pow_mod<3, 2, 13>::rem;
				b4_sum += pc_2 * pow_mod<4, 2, 7>::rem;
				b9_sum += pc_2 * pow_mod<9, 2, 13>::rem;

				// __assume(b3_sum <= get_sum_of_rems<13, in_base<3>>(outer_48_bits_mask));
				// __assume(b4_sum <= get_sum_of_rems<7, in_base<4>>(outer_48_bits_mask));
				// __assume(b9_sum <= get_sum_of_rems<13, in_base<9>>(outer_48_bits_mask));

				lookup_1_ptr = (*b3m13_lookup)[b3_sum % 13].data() + byte_offset;
				lookup_2_ptr = (*b4m7_lookup)[b4_sum % 7].data() + byte_offset;
				lookup_3_ptr = (*b9m13_lookup)[b9_sum % 13].data() + byte_offset;
			}

			{
				constexpr uint64_t bitmask = bitmask_for<3, 7>::val;
				static_assert(bitmask == bitmask_for<4, 13>::val);
				static_assert(bitmask == bitmask_for<5, 7>::val);
				static_assert(bitmask == bitmask_for<10, 13>::val);
				static_assert(period_of<bitmask>() == 6);

				const size_t pc_0 = pop_count(next_number & ((bitmask << 0) & outer_48_bits_mask));
				size_t b3_sum = pc_0;
				size_t b4_sum = pc_0;
				size_t b5_sum = pc_0;
				size_t b10_sum = pc_0;
				const size_t pc_1 = pop_count(next_number & ((bitmask << 1) & outer_48_bits_mask));
				b3_sum += pc_1 * pow_mod<3, 1, 7>::rem;
				b4_sum += pc_1 * pow_mod<4, 1, 13>::rem;
				b5_sum += pc_1 * pow_mod<5, 1, 7>::rem;
				b10_sum += pc_1 * pow_mod<10, 1, 13>::rem;
				const size_t pc_2 = pop_count(next_number & ((bitmask << 2) & outer_48_bits_mask));
				b3_sum += pc_2 * pow_mod<3, 2, 7>::rem;
				b4_sum += pc_2 * pow_mod<4, 2, 13>::rem;
				b5_sum += pc_2 * pow_mod<5, 2, 7>::rem;
				b10_sum += pc_2 * pow_mod<10, 2, 13>::rem;
				const size_t pc_3 = pop_count(next_number & ((bitmask << 3) & outer_48_bits_mask));
				b3_sum += pc_3 * pow_mod<3, 3, 7>::rem;
				b4_sum += pc_3 * pow_mod<4, 3, 13>::rem;
				b5_sum += pc_3 * pow_mod<5, 3, 7>::rem;
				b10_sum += pc_3 * pow_mod<10, 3, 13>::rem;
				const size_t pc_4 = pop_count(next_number & ((bitmask << 4) & outer_48_bits_mask));
				b3_sum += pc_4 * pow_mod<3, 4, 7>::rem;
				b4_sum += pc_4 * pow_mod<4, 4, 13>::rem;
				b5_sum += pc_4 * pow_mod<5, 4, 7>::rem;
				b10_sum += pc_4 * pow_mod<10, 4, 13>::rem;
				const size_t pc_5 = pop_count(next_number & ((bitmask << 5) & outer_48_bits_mask));
				b3_sum += pc_5 * pow_mod<3, 5, 7>::rem;
				b4_sum += pc_5 * pow_mod<4, 5, 13>::rem;
				b5_sum += pc_5 * pow_mod<5, 5, 7>::rem;
				b10_sum += pc_5 * pow_mod<10, 5, 13>::rem;

				// __assume(b3_sum <= get_sum_of_rems<7, in_base<3>>(outer_48_bits_mask));
				// __assume(b4_sum <= get_sum_of_rems<13, in_base<4>>(outer_48_bits_mask));
				// __assume(b5_sum <= get_sum_of_rems<7, in_base<5>>(outer_48_bits_mask));
				// __assume(b10_sum <= get_sum_of_rems<13, in_base<10>>(outer_48_bits_mask));

				lookup_4_ptr = (*b3m7_lookup)[b3_sum % 7].data() + byte_offset;
				lookup_5_ptr = (*b4m13_lookup)[b4_sum % 13].data() + byte_offset;
				lookup_6_ptr = (*b5m7_lookup)[b5_sum % 7].data() + byte_offset;
				lookup_7_ptr = (*b10m13_lookup)[b10_sum % 13].data() + byte_offset;
			}
		}

		if constexpr (pass == 3)
		{
			constexpr uint64_t bitmask = bitmask_for<3, 11>::val;
			static_assert(bitmask == bitmask_for<4, 11>::val);
			static_assert(bitmask == bitmask_for<5, 11>::val);
			static_assert(bitmask == bitmask_for<9, 11>::val);
			static_assert(period_of<bitmask>() == 5);

			const size_t pc_0 = pop_count(next_number & ((bitmask << 0) & outer_48_bits_mask));
			size_t b3_sum = pc_0;
			size_t b4_sum = pc_0;
			size_t b5_sum = pc_0;
			size_t b9_sum = pc_0;
			const size_t pc_1 = pop_count(next_number & ((bitmask << 1) & outer_48_bits_mask));
			b3_sum += pc_1 * pow_mod<3, 1, 11>::rem;
			b4_sum += pc_1 * pow_mod<4, 1, 11>::rem;
			b5_sum += pc_1 * pow_mod<5, 1, 11>::rem;
			b9_sum += pc_1 * pow_mod<9, 1, 11>::rem;
			const size_t pc_2 = pop_count(next_number & ((bitmask << 2) & outer_48_bits_mask));
			b3_sum += pc_2 * pow_mod<3, 2, 11>::rem;
			b4_sum += pc_2 * pow_mod<4, 2, 11>::rem;
			b5_sum += pc_2 * pow_mod<5, 2, 11>::rem;
			b9_sum += pc_2 * pow_mod<9, 2, 11>::rem;
			const size_t pc_3 = pop_count(next_number & ((bitmask << 3) & outer_48_bits_mask));
			b3_sum += pc_3 * pow_mod<3, 3, 11>::rem;
			b4_sum += pc_3 * pow_mod<4, 3, 11>::rem;
			b5_sum += pc_3 * pow_mod<5, 3, 11>::rem;
			b9_sum += pc_3 * pow_mod<9, 3, 11>::rem;
			const size_t pc_4 = pop_count(next_number & ((bitmask << 4) & outer_48_bits_mask));
			b3_sum += pc_4 * pow_mod<3, 4, 11>::rem;
			b4_sum += pc_4 * pow_mod<4, 4, 11>::rem;
			b5_sum += pc_4 * pow_mod<5, 4, 11>::rem;
			b9_sum += pc_4 * pow_mod<9, 4, 11>::rem;

			// __assume(b3_sum <= get_sum_of_rems<11, in_base<3>>(outer_48_bits_mask));
			// __assume(b4_sum <= get_sum_of_rems<11, in_base<4>>(outer_48_bits_mask));
			// __assume(b5_sum <= get_sum_of_rems<11, in_base<5>>(outer_48_bits_mask));
			// __assume(b9_sum <= get_sum_of_rems<11, in_base<9>>(outer_48_bits_mask));

			lookup_1_ptr = (*b3m11_lookup)[b3_sum % 11].data() + byte_offset;
			lookup_2_ptr = (*b4m11_lookup)[b4_sum % 11].data() + byte_offset;
			lookup_3_ptr = (*b5m11_lookup)[b5_sum % 11].data() + byte_offset;
			lookup_4_ptr = (*b9m11_lookup)[b9_sum % 11].data() + byte_offset;
		}
	}

	constexpr size_t last_pass = 3;

	template<size_t pass>
	__forceinline void merge_one_block(uint8_t*& out,
									   const size_t sieve_offset,
									   uint64_t& number,
									   const uint8_t*& in,
									   const uint8_t*& lookup_1_ptr,
									   const uint8_t*& lookup_2_ptr,
									   const uint8_t*& lookup_3_ptr,
									   const uint8_t*& lookup_4_ptr,
									   const uint8_t*& lookup_5_ptr,
									   const uint8_t*& lookup_6_ptr,
									   const uint8_t*& lookup_7_ptr,
									   size_t& sieve_popcount)
	{
		const size_t elements_to_rollover = (pow_2_16 - ((number & bits_1_16_mask) >> 1)) % pow_2_16; // map 65,536 -> 0

		uint64_t mask = *(uint64_t*)in;
		uint64_t bit_patterns_mask{};
		if constexpr (pass < 3)
			bit_patterns_mask = (*(uint64_t*)lookup_1_ptr &
								 *(uint64_t*)lookup_2_ptr &
								 *(uint64_t*)lookup_3_ptr &
								 *(uint64_t*)lookup_4_ptr &
								 *(uint64_t*)lookup_5_ptr &
								 *(uint64_t*)lookup_6_ptr &
								 *(uint64_t*)lookup_7_ptr);
		else
			bit_patterns_mask = (*(uint64_t*)lookup_1_ptr &
								 *(uint64_t*)lookup_2_ptr &
								 *(uint64_t*)lookup_3_ptr &
								 *(uint64_t*)lookup_4_ptr);

		in += sizeof(uint64_t);
		lookup_1_ptr += sizeof(uint64_t);
		lookup_2_ptr += sizeof(uint64_t);
		lookup_3_ptr += sizeof(uint64_t);
		lookup_4_ptr += sizeof(uint64_t);
		if constexpr (pass < 3)
		{
			lookup_5_ptr += sizeof(uint64_t);
			lookup_6_ptr += sizeof(uint64_t);
			lookup_7_ptr += sizeof(uint64_t);
		}

		if (elements_to_rollover >= 64) // copy another block using lookup data
		{
			mask &= bit_patterns_mask;
		}
		else // elements_to_rollover is 0 to 63
		{
			const uint8_t v{}; // we don't want set_lookup_ptrs() to modify the input ptr
			auto dummy_ptr = &v;

			set_lookup_ptrs<pass>(out, 0, number + (elements_to_rollover * 2), dummy_ptr,
								  lookup_1_ptr, lookup_2_ptr, lookup_3_ptr, lookup_4_ptr, lookup_5_ptr, lookup_6_ptr, lookup_7_ptr);

			uint64_t new_mask{};
			if constexpr (pass < 3)
				new_mask = (*(uint64_t*)lookup_1_ptr &
							*(uint64_t*)lookup_2_ptr &
							*(uint64_t*)lookup_3_ptr &
							*(uint64_t*)lookup_4_ptr &
							*(uint64_t*)lookup_5_ptr &
							*(uint64_t*)lookup_6_ptr &
							*(uint64_t*)lookup_7_ptr) << elements_to_rollover;
			else
				new_mask = (*(uint64_t*)lookup_1_ptr &
							*(uint64_t*)lookup_2_ptr &
							*(uint64_t*)lookup_3_ptr &
							*(uint64_t*)lookup_4_ptr) << elements_to_rollover;

			const uint64_t select_from_old = (1ull << elements_to_rollover) - 1; // up to and including rollover
			const uint64_t select_from_new = ~select_from_old;

			mask &= ((bit_patterns_mask & select_from_old) |
					 (new_mask & select_from_new));

			// (re)set
			set_lookup_ptrs<pass>(out, sieve_offset + sizeof(uint64_t), number + (64ull * 2), dummy_ptr,
								  lookup_1_ptr, lookup_2_ptr, lookup_3_ptr, lookup_4_ptr, lookup_5_ptr, lookup_6_ptr, lookup_7_ptr);
		}

		*(uint64_t*)out = mask;
		out += sizeof(uint64_t);

		if constexpr (pass == last_pass)
		{
			sieve_popcount += pop_count(mask);
		}

		number += (64ull * 2);
	}



	template<size_t pass>
	inline_toggle size_t merge_bitmasks(uint64_t number, sieve_container& sieve)
	{
		constexpr static uint8_t static_pc_shuf_lookup[32] = {
			0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
			0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
		constexpr static uint64_t static_nybble_mask[4] = {
			0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F };
		// while out != aligned_end:
		// - iterate by 32 bytes until we hit aligned_end, or a rollover of bits 1-16, whichever happens first
		// - if we broke before aligned_end, we hit a rollover
		//   - handle the next 32 bytes separately, then loop
		// handle cleanup

		size_t sieve_popcount = 0;

		uint8_t* out = sieve.data();

		const uint8_t* in{}; // either the static sieve (pass 1) or the sieve (every other pass)
		const uint8_t* lookup_1_ptr{}; // pointers into our large bitmask lookups
		const uint8_t* lookup_2_ptr{};
		const uint8_t* lookup_3_ptr{};
		const uint8_t* lookup_4_ptr{};
		const uint8_t* lookup_5_ptr{};
		const uint8_t* lookup_6_ptr{};
		const uint8_t* lookup_7_ptr{};

		set_lookup_ptrs<pass>(out, 0, number, in,
							  lookup_1_ptr, lookup_2_ptr, lookup_3_ptr, lookup_4_ptr, lookup_5_ptr, lookup_6_ptr, lookup_7_ptr);

		constexpr size_t bytes_per_step = sizeof(uint256_t);

		const uint8_t* const aligned_end = out + ((sieve_container::size() / 256) * bytes_per_step);

		for (; out != aligned_end; )
		{
			// The hot loop handles 256 elements (256 bits) per iteration.
			// How many times can we do this before reaching the lookups' ends?
			// Calculate this up front so the hot loop can run using a simpler condition.

			constexpr size_t elements_per_step = bytes_per_step * 8;

			const size_t elements_to_rollover = pow_2_16 - ((number & bits_1_16_mask) >> 1);

			const size_t n_steps = std::min(elements_to_rollover / elements_per_step,
											(aligned_end - out) / bytes_per_step);

			uint256_t nybble_mask{};
			uint256_t pc_shuf_lookup{};
			uint256_t pc{};
			if constexpr (pass == last_pass)
			{
				nybble_mask = _mm256_loadu_si256((uint256_t*)static_nybble_mask);
				pc_shuf_lookup = _mm256_loadu_si256((uint256_t*)static_pc_shuf_lookup);
			}

			for (size_t offset = 0; offset < n_steps * bytes_per_step; offset += bytes_per_step)
			{
				uint256_t data_0;
				if constexpr (pass == 1)
					data_0 = _mm256_loadu_si256((uint256_t*)(in + offset));
				else
					data_0 = _mm256_loadu_si256((uint256_t*)(out + offset));

				const uint256_t data_1 = _mm256_loadu_si256((uint256_t*)(lookup_1_ptr + offset));
				const uint256_t data_2 = _mm256_loadu_si256((uint256_t*)(lookup_2_ptr + offset));
				const uint256_t data_3 = _mm256_loadu_si256((uint256_t*)(lookup_3_ptr + offset));
				const uint256_t data_4 = _mm256_loadu_si256((uint256_t*)(lookup_4_ptr + offset));
				uint256_t data_5{};
				uint256_t data_6{};
				uint256_t data_7{};
				if constexpr (pass < 3)
				{
					data_5 = _mm256_loadu_si256((uint256_t*)(lookup_5_ptr + offset));
					data_6 = _mm256_loadu_si256((uint256_t*)(lookup_6_ptr + offset));
					data_7 = _mm256_loadu_si256((uint256_t*)(lookup_7_ptr + offset));
				}

				if (offset % (2 * bytes_per_step) == 0) // every other iteration
				{
					constexpr size_t prefetch_offset = 2ull * 64;

					_mm_prefetch((char*)in + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_1_ptr + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_2_ptr + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_3_ptr + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_4_ptr + offset + prefetch_offset, _MM_HINT_T0);
					if constexpr (pass < 3)
					{
						_mm_prefetch((char*)lookup_5_ptr + offset + prefetch_offset, _MM_HINT_T0);
						_mm_prefetch((char*)lookup_6_ptr + offset + prefetch_offset, _MM_HINT_T0);
						_mm_prefetch((char*)lookup_7_ptr + offset + prefetch_offset, _MM_HINT_T0);
					}

					// Demote input and output cache lines. Both will be reused, but some mask data will be reused sooner.
					// In passes >1, input and output pointers are the same, so [in] only needs to be demoted in pass 1.

					//if constexpr (pass == 1)
					//{
					//	_mm_cldemote(in + offset - 64);
					//}

					//_mm_cldemote(out + offset - 64);
				}

				uint256_t merged_data{};
				const uint256_t m1 = _mm256_and_si256(data_0, data_1);
				const uint256_t m2 = _mm256_and_si256(data_2, data_3);
				if constexpr (pass < 3)
				{
					const uint256_t m3 = _mm256_and_si256(data_4, data_5);
					const uint256_t m4 = _mm256_and_si256(data_6, data_7);
					merged_data = _mm256_and_si256(_mm256_and_si256(m1, m2),
												   _mm256_and_si256(m3, m4));
				}
				else
				{
					merged_data = _mm256_and_si256(_mm256_and_si256(m1, m2),
												   data_4);
				}

				_mm256_storeu_si256((uint256_t*)(out + offset), merged_data);

				if constexpr (pass == last_pass)
				{
					// data -> nybbles
					const uint256_t nybbles_lo = _mm256_and_si256(merged_data, nybble_mask);
					const uint256_t nybbles_hi = _mm256_and_si256(_mm256_srli_epi64(merged_data, 4), nybble_mask);
					// nybbles -> 8-bit pcs
					uint256_t pc_256 = _mm256_add_epi8(_mm256_shuffle_epi8(pc_shuf_lookup, nybbles_lo),
													   _mm256_shuffle_epi8(pc_shuf_lookup, nybbles_hi));
					// 8-bit pcs -> 64-bit pcs
					pc_256 = _mm256_sad_epu8(pc_256, _mm256_setzero_si256());
					// 64-bit pcs -> running total
					pc = _mm256_add_epi64(pc, pc_256);
				}
			} // end hot copy/merge loop

			out += n_steps * bytes_per_step;
			in += n_steps * bytes_per_step;
			lookup_1_ptr += n_steps * bytes_per_step;
			lookup_2_ptr += n_steps * bytes_per_step;
			lookup_3_ptr += n_steps * bytes_per_step;
			lookup_4_ptr += n_steps * bytes_per_step;
			if constexpr (pass < 3)
			{
				lookup_5_ptr += n_steps * bytes_per_step;
				lookup_6_ptr += n_steps * bytes_per_step;
				lookup_7_ptr += n_steps * bytes_per_step;
			}

			number += n_steps * elements_per_step * 2; // * 2 because the sieve only contains odd numbers

			if constexpr (pass == last_pass)
			{
				alignas(sizeof(uint256_t)) uint64_t buf[4]{};
				_mm256_storeu_si256((uint256_t*)buf, pc);

				sieve_popcount += buf[0];
				sieve_popcount += buf[1];
				sieve_popcount += buf[2];
				sieve_popcount += buf[3];
			}

			// if we are not at aligned_end, we are at a rollover
			if (out != aligned_end)
			{
				// handle 4 64-bit blocks
				for (size_t block = 0; block < 4; ++block)
				{
					merge_one_block<pass>(out, out - sieve.data(), number, in,
										  lookup_1_ptr, lookup_2_ptr, lookup_3_ptr, lookup_4_ptr, lookup_5_ptr, lookup_6_ptr, lookup_7_ptr,
										  sieve_popcount);
				}
			}
		} // end main copy/merge loop

		// Less than 256 elements left. Generate instructions for 0-4 final copies.
		constexpr size_t leftover_elements = sieve_container::size() % 256;

		if constexpr (leftover_elements > 64ull * 0)
		{
			merge_one_block<pass>(out, out - sieve.data(), number, in,
								  lookup_1_ptr, lookup_2_ptr, lookup_3_ptr, lookup_4_ptr, lookup_5_ptr, lookup_6_ptr, lookup_7_ptr,
								  sieve_popcount);
		}
		if constexpr (leftover_elements > 64ull * 1)
		{
			merge_one_block<pass>(out, out - sieve.data(), number, in,
								  lookup_1_ptr, lookup_2_ptr, lookup_3_ptr, lookup_4_ptr, lookup_5_ptr, lookup_6_ptr, lookup_7_ptr,
								  sieve_popcount);
		}
		if constexpr (leftover_elements > 64ull * 2)
		{
			merge_one_block<pass>(out, out - sieve.data(), number, in,
								  lookup_1_ptr, lookup_2_ptr, lookup_3_ptr, lookup_4_ptr, lookup_5_ptr, lookup_6_ptr, lookup_7_ptr,
								  sieve_popcount);
		}
		if constexpr (leftover_elements > 64ull * 3)
		{
			merge_one_block<pass>(out, out - sieve.data(), number, in,
								  lookup_1_ptr, lookup_2_ptr, lookup_3_ptr, lookup_4_ptr, lookup_5_ptr, lookup_6_ptr, lookup_7_ptr,
								  sieve_popcount);
		}

		return sieve_popcount;
	}


	void merge_bitmasks_one_by_one(uint64_t number,
								   std::array<sieve_container, prime_sieve::steps>& sieves,
								   std::array<size_t, prime_sieve::steps>& sieve_popcounts);
}
