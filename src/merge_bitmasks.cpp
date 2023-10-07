
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



	// alternate lite version (one lookup per pass)

	template<size_t pass = 0>
	__forceinline void set_lookup_ptr(const uint8_t* out_ptr,
									  const size_t sieve_offset,
									  const uint64_t next_number,
									  const uint8_t*& in_ptr,
									  const uint8_t*& lookup_ptr)
	{
		using namespace div_test;
		using namespace div_test::detail;

		static_assert(pass >= 1 && pass <= 18);

		const uint64_t bits_1_16 = (next_number & bits_1_16_mask) >> 1;
		const size_t byte_offset = bits_1_16 / 8;

		if constexpr (pass != 1)
		{
			in_ptr = out_ptr; // wherever we're writing to in the sieve is also the location we have to read from
		}

		if constexpr (pass == 1) // static sieve + pc lookup
		{
			in_ptr = static_sieve.data() + sieve_offset; // the offset into the sieve is also the offset into the static sieve

			const size_t pc_idx = pc_lookup_idx(next_number);
			lookup_ptr = pc_lookup[pc_idx].data() + byte_offset;
		}
		else if constexpr (pass == 2) // gcd lookup
		{
			const size_t gcd_idx = gcd_lookup_idx(next_number);
			lookup_ptr = gcd_lookup[gcd_idx].data() + byte_offset;
		}
		else // passes 3-n use a base % prime bitmask
		{
			constexpr size_t base_mod_prime[18][2] = { { 3, 5 },
													  { 3, 7 },
													  { 3, 13 },
													  { 4, 7 },
													  { 4, 13 },
													  { 4, 17 },
													  { 5, 7 },
													  { 5, 13 },
													  { 8, 13 },
													  { 9, 13 },
													  { 10, 13 },
													  { 13, 17 },
													  { 3, 11 },
													  { 4, 11 },
													  { 5, 11 },
													  { 9, 11 } };

			constexpr size_t base = base_mod_prime[pass - 3][0];
			constexpr size_t prime = base_mod_prime[pass - 3][1];

			const size_t sum_of_rems = get_sum_of_rems<prime, in_base<base>>(next_number & outer_48_bits_mask);
			__assume(sum_of_rems <= get_sum_of_rems<prime, in_base<base>>(outer_48_bits_mask));

			bit_array<pow_2_16>* ptr{};
			if constexpr (pass == 3) { ptr = b3m5_lookup->data(); }
			if constexpr (pass == 4) { ptr = b3m7_lookup->data(); }
			if constexpr (pass == 5) { ptr = b3m13_lookup->data(); }
			if constexpr (pass == 6) { ptr = b4m7_lookup->data(); }
			if constexpr (pass == 7) { ptr = b4m13_lookup->data(); }
			if constexpr (pass == 8) { ptr = b4m17_lookup->data(); }
			if constexpr (pass == 9) { ptr = b5m7_lookup->data(); }
			if constexpr (pass == 10) { ptr = b5m13_lookup->data(); }
			if constexpr (pass == 11) { ptr = b8m13_lookup->data(); }
			if constexpr (pass == 12) { ptr = b9m13_lookup->data(); }
			if constexpr (pass == 13) { ptr = b10m13_lookup->data(); }
			if constexpr (pass == 14) { ptr = b13m17_lookup->data(); }
			if constexpr (pass == 15) { ptr = b3m11_lookup->data(); }
			if constexpr (pass == 16) { ptr = b4m11_lookup->data(); }
			if constexpr (pass == 17) { ptr = b5m11_lookup->data(); }
			if constexpr (pass == 18) { ptr = b9m11_lookup->data(); }

			lookup_ptr = ptr[sum_of_rems % prime].data() + byte_offset;
		}
	}

	using set_lookup_f = decltype(&set_lookup_ptr<0>);

	__forceinline void merge_one_block(uint8_t*& out,
									   const size_t sieve_offset,
									   uint64_t& number,
									   const uint8_t*& in,
									   const uint8_t*& lookup_ptr,
									   set_lookup_f set_lookup_f,
									   size_t& sieve_popcount)
	{
		const size_t elements_to_rollover = (pow_2_16 - ((number & bits_1_16_mask) >> 1)) % pow_2_16; // map 65,536 -> 0

		uint64_t mask = *(uint64_t*)in;
		uint64_t bit_patterns_mask = *(uint64_t*)lookup_ptr;

		in += sizeof(uint64_t);
		lookup_ptr += sizeof(uint64_t);

		if (elements_to_rollover >= 64) // copy another block using lookup data
		{
			mask &= bit_patterns_mask;
		}
		else // elements_to_rollover is 0 to 63
		{
			const uint8_t v{}; // we don't want set_lookup_ptr() to modify the input ptr
			auto dummy_ptr = &v;

			set_lookup_f(out, 0, number + (elements_to_rollover * 2), dummy_ptr, lookup_ptr);

			const uint64_t new_mask = (*(uint64_t*)lookup_ptr) << elements_to_rollover;

			const uint64_t select_from_old = (1ull << elements_to_rollover) - 1; // up to and including rollover
			const uint64_t select_from_new = ~select_from_old;

			mask &= ((bit_patterns_mask & select_from_old) |
					 (new_mask & select_from_new));

			// (re)set
			set_lookup_f(out, sieve_offset + sizeof(uint64_t), number + (64ull * 2), dummy_ptr, lookup_ptr);
		}

		*(uint64_t*)out = mask;
		out += sizeof(uint64_t);

		sieve_popcount += pop_count(mask);

		number += (64ull * 2);
	}

	size_t merge_bitmask(uint64_t number, sieve_container& sieve, set_lookup_f set_lookup_f,
						 const bool popcount_on_this_pass)
	{
		constexpr static uint256_t static_pc_shuf_lookup{.m256i_u8{
			0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
				0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 } };
		constexpr static uint256_t static_nybble_mask{.m256i_u64{
			0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F, 0x0F0F0F0F0F0F0F0F } };

		size_t sieve_popcount = 0;

		uint8_t* out = sieve.data();

		const uint8_t* in{}; // either the static sieve (pass 1) or the sieve (every other pass)
		const uint8_t* lookup_ptr{}; // pointer into a large bitmask lookup

		set_lookup_f(out, 0, number, in, lookup_ptr);

		constexpr size_t bytes_per_step = sizeof(uint256_t);

		const uint8_t* const aligned_end = out + ((sieve_container::size() / 256) * bytes_per_step);

		for (; out != aligned_end; )
		{
			constexpr size_t elements_per_step = bytes_per_step * 8;
			const size_t elements_to_rollover = pow_2_16 - ((number & bits_1_16_mask) >> 1);
			const size_t n_steps = std::min(elements_to_rollover / elements_per_step,
											(aligned_end - out) / bytes_per_step);

			const uint8_t* const in_a = in;
			const uint8_t* const in_b = lookup_ptr;

			if (popcount_on_this_pass)
			{
				uint256_t nybble_mask = _mm256_loadu_si256(&static_nybble_mask);
				uint256_t pc_shuf_lookup = _mm256_loadu_si256(&static_pc_shuf_lookup);
				uint256_t pc{};

				for (size_t offset = 0; offset < n_steps * bytes_per_step; offset += bytes_per_step)
				{
					const uint256_t sieve_data = _mm256_loadu_si256((uint256_t*)(/*in*/ in_a + offset));
					const uint256_t lookup_data = _mm256_loadu_si256((uint256_t*)(/*lookup_ptr*/ in_b + offset));
					const uint256_t merged_data = _mm256_and_si256(sieve_data, lookup_data);
					_mm256_storeu_si256((uint256_t*)(out + offset), merged_data);

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

				alignas(sizeof(uint256_t)) uint64_t buf[4]{};
				_mm256_storeu_si256((uint256_t*)buf, pc);

				sieve_popcount += buf[0];
				sieve_popcount += buf[1];
				sieve_popcount += buf[2];
				sieve_popcount += buf[3];
			}
			else
			{
				for (size_t offset = 0; offset < n_steps * bytes_per_step; offset += bytes_per_step)
				{
					const uint256_t sieve_data = _mm256_loadu_si256((uint256_t*)(/*in*/in_a + offset));
					const uint256_t lookup_data = _mm256_loadu_si256((uint256_t*)(/*lookup_ptr*/in_b + offset));
					const uint256_t merged_data = _mm256_and_si256(sieve_data, lookup_data);
					_mm256_storeu_si256((uint256_t*)(out + offset), merged_data);
				}
			}

			out += n_steps * bytes_per_step;
			in += n_steps * bytes_per_step;
			lookup_ptr += n_steps * bytes_per_step;

			number += n_steps * elements_per_step * 2;

			// if we are not at aligned_end, we are at a rollover
			if (out != aligned_end)
			{
				// handle 4 64-bit blocks
				for (size_t block = 0; block < 4; ++block)
				{
					merge_one_block(out, out - sieve.data(), number, in, lookup_ptr, set_lookup_f, sieve_popcount);
				}
			}
		} // end main copy/merge loop

		constexpr size_t leftover_elements = sieve_container::size() % 256;

		if constexpr (leftover_elements > 64ull * 0)
		{
			merge_one_block(out, out - sieve.data(), number, in, lookup_ptr, set_lookup_f, sieve_popcount);
		}
		if constexpr (leftover_elements > 64ull * 1)
		{
			merge_one_block(out, out - sieve.data(), number, in, lookup_ptr, set_lookup_f, sieve_popcount);
		}
		if constexpr (leftover_elements > 64ull * 2)
		{
			merge_one_block(out, out - sieve.data(), number, in, lookup_ptr, set_lookup_f, sieve_popcount);
		}
		if constexpr (leftover_elements > 64ull * 3)
		{
			merge_one_block(out, out - sieve.data(), number, in, lookup_ptr, set_lookup_f, sieve_popcount);
		}

		return sieve_popcount;
	}

	namespace detail
	{
		template<size_t n_of_masks, size_t mask = 1>
		__forceinline void merge_bitmasks_one_by_one_impl(const uint64_t number,
														  std::array<sieve_container, prime_sieve::steps>& sieves,
														  std::array<size_t, prime_sieve::steps>& sieve_popcounts)
		{
			for (size_t i = 0; i < prime_sieve::steps; ++i)
			{
				const uint64_t sieve_start = number + (i * sieve_container::size() * 2);
				//std::cout << "Calling merge_bitmask with mask=" << mask << " on sieve i=" << i << '\n';
				//std::cin.ignore();

				const size_t popcount = merge_bitmask(sieve_start, sieves[i], &set_lookup_ptr<mask>, mask == n_of_masks);
				if constexpr (mask == n_of_masks)
					sieve_popcounts[i] = popcount;
			}

			if constexpr (mask < n_of_masks)
				return merge_bitmasks_one_by_one_impl<n_of_masks, mask + 1>(number, sieves, sieve_popcounts);
		}
	}

	void merge_bitmasks_one_by_one(const uint64_t number,
								   std::array<sieve_container, prime_sieve::steps>& sieves,
								   std::array<size_t, prime_sieve::steps>& sieve_popcounts)
	{
		detail::merge_bitmasks_one_by_one_impl<18>(number, sieves, sieve_popcounts);
	}

}
