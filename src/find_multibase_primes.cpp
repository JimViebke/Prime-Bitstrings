
#include <iomanip>
#include <iostream>
#include <sstream>

#include "bit_pattern_tests.hpp"
#include "find_multibase_primes.hpp"
#include "hardcoded_div_tests.hpp"
#include "io/io.hpp"
#include "sieve.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/simd.hpp"
#include "util/types.hpp"

namespace mbp
{
	static const sieve_container static_sieve = prime_sieve::generate_static_sieve();
	static std::unique_ptr<std::array<sieve_container, prime_sieve::steps>> sieves =
		std::make_unique<decltype(sieves)::element_type>();
	static std::array<size_t, prime_sieve::steps> sieve_popcounts{};

	namespace detail
	{
		// 64 bits == (47 high bits) + (16 "inner" bits) + (1 low bit (always set))
		constexpr uint64_t bits_1_16_mask = 0xFFFFull << 1;
		constexpr uint64_t outer_48_bits_mask = ~bits_1_16_mask;
	}

	__forceinline size_t pc_lookup_idx(const uint64_t number)
	{
		return pop_count(number & detail::outer_48_bits_mask) - 2; // -2 to normalize popcount 2-48 to idx 0-46
	}

	__forceinline size_t gcd_lookup_idx(const uint64_t number)
	{
		constexpr uint64_t even_mask = detail::outer_48_bits_mask & 0x5555555555555555;
		constexpr uint64_t odd_mask = detail::outer_48_bits_mask & 0xAAAAAAAAAAAAAAAA;

		const auto outer_even_pc = pop_count(number & even_mask);
		const auto outer_odd_pc = pop_count(number & odd_mask);
		return outer_even_pc - outer_odd_pc + 23; // +23 to normalize -23,24 to 0,47
	}

	template<size_t base, size_t prime>
	__forceinline size_t get_lookup_idx_for(const uint64_t number)
	{
		using namespace detail;
		using namespace div_test;
		using namespace div_test::detail;

		constexpr size_t max_sum_of_rems = get_sum_of_rems<prime, in_base<base>>(outer_48_bits_mask);

		const size_t sum = get_sum_of_rems<prime, in_base<base>>(number & outer_48_bits_mask);

		__assume(sum <= max_sum_of_rems);
		return sum % prime;
	}

	// out_ptr is used to set ss_ptr in pass 1, and sieve_data in pass 2
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
		using namespace detail;
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
			static_assert(period_of<bitmask>::val == 4);

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

			__assume(b3_sum <= get_sum_of_rems<5, in_base<3>>(outer_48_bits_mask));
			__assume(b4_sum <= get_sum_of_rems<17, in_base<4>>(outer_48_bits_mask));
			__assume(b5_sum <= get_sum_of_rems<13, in_base<5>>(outer_48_bits_mask));
			__assume(b8_sum <= get_sum_of_rems<13, in_base<8>>(outer_48_bits_mask));
			__assume(b13_sum <= get_sum_of_rems<17, in_base<13>>(outer_48_bits_mask));

			lookup_3_ptr = (*b3m5_lookup)[b3_sum % 5].data() + byte_offset;
			lookup_4_ptr = (*b4m17_lookup)[b4_sum % 17].data() + byte_offset;
			lookup_5_ptr = (*b5m13_lookup)[b5_sum % 13].data() + byte_offset;
			lookup_6_ptr = (*b8m13_lookup)[b8_sum % 13].data() + byte_offset;
			lookup_7_ptr = (*b13m17_lookup)[b13_sum % 17].data() + byte_offset;
		}

		if constexpr (pass == 2)
		{
			in_ptr = out_ptr; // wherever we're writing to in the sieve is also the location we have to read from

			{
				constexpr uint64_t bitmask = bitmask_for<3, 13>::val;
				static_assert(bitmask == bitmask_for<4, 7>::val);
				static_assert(bitmask == bitmask_for<9, 13>::val);
				static_assert(period_of<bitmask>::val == 3);

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

				__assume(b3_sum <= get_sum_of_rems<13, in_base<3>>(outer_48_bits_mask));
				__assume(b4_sum <= get_sum_of_rems<7, in_base<4>>(outer_48_bits_mask));
				__assume(b9_sum <= get_sum_of_rems<13, in_base<9>>(outer_48_bits_mask));

				lookup_1_ptr = (*b3m13_lookup)[b3_sum % 13].data() + byte_offset;
				lookup_2_ptr = (*b4m7_lookup)[b4_sum % 7].data() + byte_offset;
				lookup_3_ptr = (*b9m13_lookup)[b9_sum % 13].data() + byte_offset;
			}

			{
				constexpr uint64_t bitmask = bitmask_for<3, 7>::val;
				static_assert(bitmask == bitmask_for<4, 13>::val);
				static_assert(bitmask == bitmask_for<5, 7>::val);
				static_assert(bitmask == bitmask_for<10, 13>::val);
				static_assert(period_of<bitmask>::val == 6);

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

				__assume(b3_sum <= get_sum_of_rems<7, in_base<3>>(outer_48_bits_mask));
				__assume(b4_sum <= get_sum_of_rems<13, in_base<4>>(outer_48_bits_mask));
				__assume(b5_sum <= get_sum_of_rems<7, in_base<5>>(outer_48_bits_mask));
				__assume(b10_sum <= get_sum_of_rems<13, in_base<10>>(outer_48_bits_mask));

				lookup_4_ptr = (*b3m7_lookup)[b3_sum % 7].data() + byte_offset;
				lookup_5_ptr = (*b4m13_lookup)[b4_sum % 13].data() + byte_offset;
				lookup_6_ptr = (*b5m7_lookup)[b5_sum % 7].data() + byte_offset;
				lookup_7_ptr = (*b10m13_lookup)[b10_sum % 13].data() + byte_offset;
			}
		}
	}

	constexpr size_t last_pass = 2;

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
		const size_t elements_to_rollover = (pow_2_16 - ((number & detail::bits_1_16_mask) >> 1)) % pow_2_16; // map 65,536 -> 0

		uint64_t mask = *(uint64_t*)in;
		const uint64_t bit_patterns_mask = (*(uint64_t*)lookup_1_ptr &
											*(uint64_t*)lookup_2_ptr &
											*(uint64_t*)lookup_3_ptr &
											*(uint64_t*)lookup_4_ptr &
											*(uint64_t*)lookup_5_ptr &
											*(uint64_t*)lookup_6_ptr &
											*(uint64_t*)lookup_7_ptr);

		in += sizeof(uint64_t);
		lookup_1_ptr += sizeof(uint64_t);
		lookup_2_ptr += sizeof(uint64_t);
		lookup_3_ptr += sizeof(uint64_t);
		lookup_4_ptr += sizeof(uint64_t);
		lookup_5_ptr += sizeof(uint64_t);
		lookup_6_ptr += sizeof(uint64_t);
		lookup_7_ptr += sizeof(uint64_t);

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

			const uint64_t new_mask = (*(uint64_t*)lookup_1_ptr &
									   *(uint64_t*)lookup_2_ptr &
									   *(uint64_t*)lookup_3_ptr &
									   *(uint64_t*)lookup_4_ptr &
									   *(uint64_t*)lookup_5_ptr &
									   *(uint64_t*)lookup_6_ptr &
									   *(uint64_t*)lookup_7_ptr) << elements_to_rollover;

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
	size_t merge_bitmasks(size_t number, sieve_container& sieve)
	{
		constexpr static uint256_t static_pc_shuf_lookup{ .m256i_u8{
			0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
			0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 } };
		constexpr static uint256_t static_nybble_mask{ .m256i_u64{
			0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F } };

		size_t sieve_popcount = 0;

		const uint8_t* const aligned_end = sieve.data() + ((sieve.size() / 256) * sizeof(uint256_t));
		uint8_t* out = sieve.data();

		// while out != aligned_end:
		// - iterate by 4 blocks until we hit aligned_end, or a rollover of bits 1-16, whichever happens first
		// - if we broke before aligned_end, we hit a rollover
		//   - handle the next 4 blocks seperately, then continue
		// handle cleanup

		const uint8_t* in{};
		// pointers into our large bitmask lookups
		const uint8_t* lookup_1_ptr{};
		const uint8_t* lookup_2_ptr{};
		const uint8_t* lookup_3_ptr{};
		const uint8_t* lookup_4_ptr{};
		const uint8_t* lookup_5_ptr{};
		const uint8_t* lookup_6_ptr{};
		const uint8_t* lookup_7_ptr{};

		set_lookup_ptrs<pass>(out, 0, number, in,
							  lookup_1_ptr,
							  lookup_2_ptr,
							  lookup_3_ptr,
							  lookup_4_ptr,
							  lookup_5_ptr,
							  lookup_6_ptr,
							  lookup_7_ptr);

		for (; out != aligned_end; )
		{
			// The hot loop handles 256 elements (256 bits) per iteration.
			// How many times can we do this before reaching the lookups' ends?
			// Calculate this up front so the hot loop can run using a simpler condition.

			const size_t elements_to_rollover = pow_2_16 - ((number & detail::bits_1_16_mask) >> 1);

			constexpr size_t bytes_per_step = sizeof(uint256_t);
			constexpr size_t elements_per_step = bytes_per_step * 8;

			const size_t n_steps = std::min(elements_to_rollover / elements_per_step,
											(aligned_end - out) / bytes_per_step);

			uint256_t nybble_mask{};
			uint256_t pc_shuf_lookup{};
			uint256_t pc{};
			if constexpr (pass == last_pass)
			{
				nybble_mask = _mm256_loadu_si256(&static_nybble_mask);
				pc_shuf_lookup = _mm256_loadu_si256(&static_pc_shuf_lookup);
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
				const uint256_t data_5 = _mm256_loadu_si256((uint256_t*)(lookup_5_ptr + offset));
				const uint256_t data_6 = _mm256_loadu_si256((uint256_t*)(lookup_6_ptr + offset));
				const uint256_t data_7 = _mm256_loadu_si256((uint256_t*)(lookup_7_ptr + offset));

				if (offset % (2 * bytes_per_step) == 0) // every other iteration
				{
					constexpr size_t prefetch_offset = 2ull * 64;

					_mm_prefetch((char*)in + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_1_ptr + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_2_ptr + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_3_ptr + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_4_ptr + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_5_ptr + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_6_ptr + offset + prefetch_offset, _MM_HINT_T0);
					_mm_prefetch((char*)lookup_7_ptr + offset + prefetch_offset, _MM_HINT_T0);
				}

				const uint256_t m1 = _mm256_and_si256(data_0, data_1);
				const uint256_t m2 = _mm256_and_si256(data_2, data_3);
				const uint256_t m3 = _mm256_and_si256(data_4, data_5);
				const uint256_t m4 = _mm256_and_si256(data_6, data_7);

				uint256_t merged_data = _mm256_and_si256(_mm256_and_si256(m1, m2),
														 _mm256_and_si256(m3, m4));

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
			}

			out += n_steps * bytes_per_step;
			in += n_steps * bytes_per_step;
			lookup_1_ptr += n_steps * bytes_per_step;
			lookup_2_ptr += n_steps * bytes_per_step;
			lookup_3_ptr += n_steps * bytes_per_step;
			lookup_4_ptr += n_steps * bytes_per_step;
			lookup_5_ptr += n_steps * bytes_per_step;
			lookup_6_ptr += n_steps * bytes_per_step;
			lookup_7_ptr += n_steps * bytes_per_step;

			number += n_steps * elements_per_step * 2; // * 2 because the sieve only contains odd numbers

			if constexpr (pass == last_pass)
			{
				sieve_popcount += pc.m256i_u64[0];
				sieve_popcount += pc.m256i_u64[1];
				sieve_popcount += pc.m256i_u64[2];
				sieve_popcount += pc.m256i_u64[3];
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
		constexpr size_t leftover_elements = sieve.size() % 256;

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



	// 2x the expected number of candidates from the sieve passes
	constexpr size_t candidates_capacity = [] {
		double cleared = 0.0;
		for (size_t i = 1; i < small_primes_lookup.size(); ++i)
			cleared += (1.0 - cleared) * (1.0 / small_primes_lookup[i]);
		return 2 * size_t((1.0 - cleared) * sieve_container::size() * prime_sieve::steps);
	}();
	static alignas(64) std::array<uint64_t, candidates_capacity> candidates_storage;

	// buffer candidates for full primality testing until we have 64
	constexpr size_t pt_buffer_capacity = 64;
	static alignas(64) std::array<uint64_t, pt_buffer_capacity> pt_buffer;
	static size_t pt_buffer_size = 0;

	mbp::find_multibase_primes::find_multibase_primes()
	{
		gmp_rand.seed(mpir_ui{ 0xdeadbeef });

		count_passes(std::cout << "(counting passes)\n");
		count_passes(a = ps15 = b = c = d = e = f = g = bldt = 0);
		count_passes(bidt = b2 = b3 = b4 = b5 = passes = pc_hash = 0);
	}

	void mbp::find_multibase_primes::run()
	{
		constexpr size_t loop_size = 2ull * sieve_container::size() * prime_sieve::steps;

		size_t number = benchmark_mode ? bm_start : load_from_results();

		// Round starting number down to the nearest odd multiple of a product of primes
		number -= prime_sieve::product_of_static_sieve_primes; // n -= k
		number -= number % (2 * prime_sieve::product_of_static_sieve_primes); // n -= n % 2k
		number += prime_sieve::product_of_static_sieve_primes; // n += k

		// Align sieve on a multiple of 8, plus 1
		while ((number & 0b1111) != 0b0001)
			number -= (2 * prime_sieve::product_of_static_sieve_primes);

		prime_sieve::set_up_sieve_offsets_cache(number);

		size_t next_div_test_reorder = number + div_test::reorder_interval;

		// Start the clock after setup
		const auto start = util::current_time_in_ms();

		// (condition should optimize out)
		while (benchmark_mode ? number < bm_stop : true)
		{
			// The upper 32 bits of a 64 bit integer only change every 4 billion ints.
			// Detect iterations where the upper bits can not change, and allow
			// functions to optimize based on this.
			if (util::upper_32_bits_match(number, number + loop_size))
			{
				main_loop<true>(number);
			}
			else
			{
				main_loop<false>(number);
			}

			number += loop_size;



			if (next_div_test_reorder <= number)
			{
				full_div_tests.update_div_test_order();
				next_div_test_reorder += div_test::reorder_interval;

			#if analyze_div_tests
				full_div_tests.print_div_tests();
				//run_div_test_analysis(number);
			#endif
			}

			count_passes(++passes);
		} // end main loop



	#if analyze_div_tests
		  // Run one final step before exiting
		  // run_div_test_analysis();
		full_div_tests.print_div_tests();
	#endif

		std::cout << "Finished. " << util::current_time_in_ms() - start << " ms elapsed\n";

		count_passes(std::cout << passes << " main loop iters\n");

		log_pass_counts("Passed static sieve and\n"\
						"  bit pattern filters: ", a, (bm_size / 2));
		log_pass_counts("Passed sieve:          ", ps15, a);
		log_pass_counts("Passed 4-rem tests, p2:", b, ps15);
		log_pass_counts("Passed 5-rem tests:    ", c, b);
		log_pass_counts("Passed 10-rem tests:   ", d, c);
		log_pass_counts("Passed 8-rem tests:    ", e, d);
		log_pass_counts("Passed 12-rem tests:   ", f, e);
		log_pass_counts("Passed 16-rem tests:   ", g, f);
		log_pass_counts("P. branchless divtests:", bldt, g);
		log_pass_counts("P. branching divtests: ", bidt, bldt);
		log_pass_counts("Passed b2 BPSW test:   ", b2, bidt);
		log_pass_counts("Passed b3 prime test:  ", b3, b2);
		log_pass_counts("Passed b4 prime test:  ", b4, b3);
		log_pass_counts("Passed b5 prime test:  ", b5, b4);

		count_passes(std::cout << "\nhash of pass counts: " <<
					 std::hex << pc_hash << std::dec << '\n');
	}

	template<bool on_fast_path>
	void mbp::find_multibase_primes::main_loop(const uint64_t number)
	{
		uint64_t* const candidates = candidates_storage.data();
		uint64_t* candidates_end = candidates;

		// Merge static sieve, popcount, gcd, and div test bitmasks
		for (size_t i = 0; i < prime_sieve::steps; ++i)
		{
			const uint64_t sieve_start = number + (i * sieve_container::size() * 2);
			merge_bitmasks<1>(sieve_start, (*sieves)[i]);
		}
		for (size_t i = 0; i < prime_sieve::steps; ++i)
		{
			const uint64_t sieve_start = number + (i * sieve_container::size() * 2);
			sieve_popcounts[i] = merge_bitmasks<2>(sieve_start, (*sieves)[i]);
			count_passes(a += sieve_popcounts[i]);
		}

		// Sieve until one of our density thresholds is reached
		for (size_t i = 0; i < prime_sieve::steps; ++i)
		{
			const uint64_t sieve_start = number + (i * sieve_container::size() * 2);
			prime_sieve::partial_sieve(sieve_start, (*sieves)[i], sieve_popcounts[i]);
			count_passes(ps15 += (*sieves)[i].count_bits());
		}

		// Convert 1-bit candidates to 64-bit candidates
		for (size_t i = 0; i < prime_sieve::steps; ++i)
		{
			const uint64_t sieve_start = number + (i * sieve_container::size() * 2);
			candidates_end = prime_sieve::gather_sieve_results(candidates_end, (*sieves)[i], sieve_start);
		}

		// Perform some div tests separately when a specialized implementation is faster

		// base 12 mod 29, and 6 mod 37 (4 remainders, part 2)
		candidates_end = two_div_tests_with_four_rems<12, 29, 6, 37, on_fast_path>(candidates, candidates_end);
		count_passes(b += (candidates_end - candidates));

		// bases 3, 4, 5, and 9 mod 11 (5 remainders)
		candidates_end = div_tests_with_five_rems<on_fast_path>(candidates, candidates_end);
		count_passes(c += (candidates_end - candidates));

		// bases 6, 7, and 8 mod 11 (10 remainders)
		candidates_end = div_tests_with_10_rems<on_fast_path>(candidates, candidates_end);
		count_passes(d += (candidates_end - candidates));

		// bases 8 and 9 mod 17 (8 remainders)
		candidates_end = div_tests_with_8_rems<on_fast_path>(candidates, candidates_end);
		count_passes(e += (candidates_end - candidates));

		// bases 6, 7, and 11 mod 13 (12 remainders)
		candidates_end = div_tests_with_12_rems<on_fast_path>(candidates, candidates_end);
		count_passes(f += (candidates_end - candidates));

		// bases 3, 5, 6, 7, 10, 11 and 12 mod 17 (16 remainders)
		candidates_end = div_tests_with_16_rems<on_fast_path>(candidates, candidates_end);
		count_passes(g += (candidates_end - candidates));

		// bases 8 and 12 mod 19 (6 remainders)
		//candidates_end = two_div_tests_with_six_rems<on_fast_path>(candidates, candidates_end);
		//count_passes(dt6b += (candidates_end - candidates));

		// bases 4, 5, 6, and 9 mod 19 (9 remainders)
		//candidates_end = div_tests_with_nine_rems<on_fast_path>(candidates, candidates_end);
		//count_passes(dt9 += (candidates_end - candidates));

		// bases 7 and 11 mod 19 (3 remainders)
		//candidates_end = two_div_tests_with_three_rems<on_fast_path>(candidates, candidates_end);
		//count_passes(dt3b += (candidates_end - candidates));



		// Check for small prime factors across all bases
		candidates_end = full_div_tests.branchless_div_tests<on_fast_path>(candidates, candidates_end, div_test::n_of_branchless_tests);
		count_passes(bldt += (candidates_end - candidates));

		candidates_end = full_div_tests.branching_div_tests<on_fast_path>(candidates, candidates_end, div_test::n_of_branchless_tests);
		count_passes(bidt += (candidates_end - candidates));



		// Collect prime candidates until we have filled our buffer, then do full primality tests starting with base 2
		for (const auto* ptr = candidates; ptr != candidates_end; ++ptr)
		{
			pt_buffer[pt_buffer_size++] = *ptr;

			if (pt_buffer_size == pt_buffer_capacity)
			{
				full_primality_tests(pt_buffer.data(), pt_buffer.data() + pt_buffer.size());
				pt_buffer_size = 0;
			}
		}

	}



	void print_config()
	{
		std::stringstream ss{};

		if constexpr (benchmark_mode)
		{
			ss << "Benchmarking from ";
			if constexpr (bm_start == p11) ss << "p11";
			else if constexpr (bm_start == p12) ss << "p12";
			else ss << bm_start;
			ss << ", size: " << bm_size / 1'000'000'000 << " B\n";
		}

		using namespace prime_sieve;
		ss << "Static sieve size: " << static_sieve_size
			<< ", primes: 3-" << static_sieve_primes.back() << '\n';
		ss << "Sieve limit: " << largest_sieve_prime
			<< ", vector/scalar thresholds: " << vector_density_threshold << ", " << scalar_density_threshold
			<< ", steps: " << steps
			<< ", candidate capacity: " << candidates_capacity << '\n';

		//ss << div_test::div_tests.size() << " div tests ("
		//	<< div_test::n_of_branchless_tests << "branchless + "
		//	<< div_test::div_tests.size() - div_test::n_of_branchless_tests << " branching), "
		//	<< div_test::n_of_primes << " prime factors, "
		//	<< "bases 2-" << div_test::up_to_base
		//	<< ", partial reorder every " << div_test::reorder_interval / 1'000'000'000 << " B\n";

		ss << "SPRP rounds: " << prime_test::n_random_bases << ", td limit: " << largest_sieve_prime << '\n';

	#define stringify(macro) #macro

		std::cout << ss.str() << std::endl;
	}
}

template void mbp::find_multibase_primes::main_loop<true>(const uint64_t);
template void mbp::find_multibase_primes::main_loop<false>(const uint64_t);
