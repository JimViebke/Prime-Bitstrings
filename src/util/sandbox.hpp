#pragma once

#include <iostream>
#include <vector>

#pragma warning(push, 0)
#include "mpirxx.h"
#pragma warning(pop)

#include "math/franken_mpir.hpp"
#include "math/math.hpp"
#include "utility.hpp"

namespace mbp
{
	namespace detail
	{
		void find_p_n()
		{
			const auto start = current_time_in_ms();

			gmp_randclass r{ gmp_randinit_mt };
			r.seed(mpir_ui{ 0xdeadbeef });

			mpz_class mpz_number;

			size_t next_base = 4;
			const size_t max_base = 9;

			std::cout << "p2 = 10\n";
			std::cout << "p3 = 10\n";

			for (size_t number = 3; ; number += 2)
			{
				if (!franken::mpir_is_likely_prime_BPSW(number)) continue;

				for (size_t base = 3; ; ++base)
				{
					mpz_number = bin_to_base(number, base);
					if (!mpir_is_prime(mpz_number, r)) break;

					// if this is the next P
					if (base == next_base)
					{
						// print in base 10 to get the zeroes and ones
						std::cout << "p" << base << " = " << bin_to_base(number, 10) << '\n';
						++next_base;
						if (next_base > max_base) goto done;
					}
				}
			}

		done:
			std::cout << "Found in " << current_time_in_ms() - start << " ms.\n";
		}

		void calculate_static_sieve_sizes()
		{
			size_t sieve_size = 1;
			for (auto p : small_primes_lookup)
			{
				if (p == 2) continue;
				sieve_size *= p;
				std::cout << "Static sieve size would be size " << sieve_size << " using product of primes up to " << p << '\n';
				if (sieve_size > 1'000'000'000) return;
			}
		}

		void playing_with_sieves()
		{
			std::vector<size_t> numbers;
			for (size_t i = 0; i <= 100; ++i)
			{
				numbers.push_back(i);
			}
			numbers[1] = 0;

			// print
			for (auto n : numbers)
				if (n != 0)
					std::cout << n << ' ';
				else
					std::cout << ". ";
			std::cout << "\n\n";

			// zero the twos
			for (size_t i = 2 * 2; i < numbers.size(); i += 2)
				numbers[i] = 0;

			// print
			for (auto n : numbers)
				if (n != 0)
					std::cout << n << ' ';
				else
					std::cout << ". ";
			std::cout << "\n\n";

			// zero every other three
			for (size_t i = 3 * 3; i < numbers.size(); i += (2 * 3))
				numbers[i] = 0;

			// print
			for (auto n : numbers)
				if (n != 0)
					std::cout << n << ' ';
				else
					std::cout << ". ";
			std::cout << "\n\n";
		}

		void expiremental_SIMD_gather_sieve_results()
		{
			/*

			// Set some magic numbers
			std::cout << "Indices lookup:\n";
			const uint256_t sieve_gather_indices = _mm256_broadcastq_epi64(uint128_t{ .m128i_u8 = { 1, 2, 4, 7, 8, 11, 13, 14 } });
			print_vector256(sieve_gather_indices);

			std::cout << "block_offsets lookup:\n";
			const uint256_t block_offsets = { .m256i_u8 = {
								 0,  0,  0,  0,  0,  0,  0,  0,
								15, 15, 15, 15, 15, 15, 15, 15,
								30, 30, 30, 30, 30, 30, 30, 30,
								45, 45, 45, 45, 45, 45, 45, 45 } };
			print_vector256(block_offsets);

			const uint256_t ymm_number = _mm256_broadcastq_epi64(uint128_t{ .m128i_u64 = { number } });

			for (; sieve_ptr < sieve_end; sieve_ptr += 60, number += 120)
			{
				// Load 15 + 1 bytes from the sieve
				const uint128_t bytes_a = _mm_loadu_si128((const uint128_t*)sieve_ptr);
				const uint128_t bytes_b = _mm_loadu_si128((const uint128_t*)(sieve_ptr + 15));
				const uint128_t bytes_c = _mm_loadu_si128((const uint128_t*)(sieve_ptr + 30));
				const uint128_t bytes_d = _mm_loadu_si128((const uint128_t*)(sieve_ptr + 45));

				// Combine to two registers
				uint256_t ymm_ab = _mm256_inserti128_si256(_mm256_castsi128_si256(bytes_a), bytes_b, 1);
				uint256_t ymm_cd = _mm256_inserti128_si256(_mm256_castsi128_si256(bytes_c), bytes_d, 1);
				std::cout << "Aligned sieve bytes:\n";
				print_vector256(ymm_ab);
				print_vector256(ymm_cd);

				// Use the magic numbers to gather the bytes we want
				ymm_ab = _mm256_shuffle_epi8(ymm_ab, sieve_gather_indices);
				ymm_cd = _mm256_shuffle_epi8(ymm_cd, sieve_gather_indices);
				std::cout << "Gather relevant bytes into halves of each lane:\n";
				print_vector256(ymm_ab);
				print_vector256(ymm_cd);

				// Combine to one register
				uint256_t ymm_abcd = _mm256_blend_epi32(ymm_ab, ymm_cd, 0b00'11'00'11);
				std::cout << "Blend ab + cd:\n";
				print_vector256(ymm_abcd);

				// Use ymm_abcd to filter for the relevant block offsets and indicies
				const uint256_t filtered_block_offsets = _mm256_sign_epi8(block_offsets, ymm_abcd);
				const uint256_t filtered_indices = _mm256_sign_epi8(sieve_gather_indices, ymm_abcd);
				std::cout << "Filtered block offsets...:\n";
				print_vector256(filtered_block_offsets);
				std::cout << "...filtered indices:\n";
				print_vector256(filtered_indices);

				// add filters
				ymm_abcd = _mm256_add_epi8(filtered_block_offsets, filtered_indices);
				std::cout << "Offsets from sieve ptr:\n";
				print_vector256(ymm_abcd);

				// << 1 to multiply all offsets by 2
				ymm_abcd = _mm256_slli_epi64(ymm_abcd, 1);
				std::cout << "Offsets from N:\n";
				print_vector256(ymm_abcd);

				// move
				auto upperhalf = _mm256_permute4x64_epi64(ymm_abcd, 0b11'10); // out[0] = in[2], out[1] = in[3]
				auto ymm0 = _mm256_blend_epi32(_mm256_setzero_si256(), ymm_abcd, 0b0000'0001);
				auto ymm1 = _mm256_blend_epi32(_mm256_setzero_si256(), ymm_abcd, 0b0000'0010);
				auto ymm2 = _mm256_blend_epi32(_mm256_setzero_si256(), ymm_abcd, 0b0000'0100);
				auto ymm3 = _mm256_blend_epi32(_mm256_setzero_si256(), ymm_abcd, 0b0000'1000);
				auto ymm4 = _mm256_blend_epi32(_mm256_setzero_si256(), upperhalf, 0b0000'0001);
				auto ymm5 = _mm256_blend_epi32(_mm256_setzero_si256(), upperhalf, 0b0000'0010);
				auto ymm6 = _mm256_blend_epi32(_mm256_setzero_si256(), upperhalf, 0b0000'0100);
				auto ymm7 = _mm256_blend_epi32(_mm256_setzero_si256(), upperhalf, 0b0000'1000);
				std::cout << "Split sets of four byte offsets to eight registers:\n";
				print_vector256(ymm0);
				print_vector256(ymm1);
				print_vector256(ymm2);
				print_vector256(ymm3);
				print_vector256(ymm4);
				print_vector256(ymm5);
				print_vector256(ymm6);
				print_vector256(ymm7);

				// 0 already in place
				ymm1 = _mm256_srli_si256(ymm1, 1 * 4);
				ymm2 = _mm256_srli_si256(ymm2, 2 * 4);
				ymm3 = _mm256_srli_si256(ymm3, 3 * 4);
				// 4 already in place
				ymm5 = _mm256_srli_si256(ymm5, 1 * 4);
				ymm6 = _mm256_srli_si256(ymm6, 2 * 4);
				ymm7 = _mm256_srli_si256(ymm7, 3 * 4);
				std::cout << "Align sets of four byte offsets:\n";
				print_vector256(ymm0);
				print_vector256(ymm1);
				print_vector256(ymm2);
				print_vector256(ymm3);
				print_vector256(ymm4);
				print_vector256(ymm5);
				print_vector256(ymm6);
				print_vector256(ymm7);

				ymm0 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm0));
				ymm1 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm1));
				ymm2 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm2));
				ymm3 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm3));
				ymm4 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm4));
				ymm5 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm5));
				ymm6 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm6));
				ymm7 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm7));
				std::cout << "Zero-extend bytes into 64-bit positions:\n";
				print_vector256(ymm0);
				print_vector256(ymm1);
				print_vector256(ymm2);
				print_vector256(ymm3);
				print_vector256(ymm4);
				print_vector256(ymm5);
				print_vector256(ymm6);
				print_vector256(ymm7);

				ymm0 = _mm256_add_epi64(ymm0, ymm_number);
				ymm1 = _mm256_add_epi64(ymm1, ymm_number);
				ymm2 = _mm256_add_epi64(ymm2, ymm_number);
				ymm3 = _mm256_add_epi64(ymm3, ymm_number);
				ymm4 = _mm256_add_epi64(ymm4, ymm_number);
				ymm5 = _mm256_add_epi64(ymm5, ymm_number);
				ymm6 = _mm256_add_epi64(ymm6, ymm_number);
				ymm7 = _mm256_add_epi64(ymm7, ymm_number);
				std::cout << "Actually add the numbers:\n";
				print_vector256(ymm0);
				print_vector256(ymm1);
				print_vector256(ymm2);
				print_vector256(ymm3);
				print_vector256(ymm4);
				print_vector256(ymm5);
				print_vector256(ymm6);
				print_vector256(ymm7);

				for (size_t i = 0; i < 4; ++i)
				{
					auto val = ymm0.m256i_u64[i];
					*sieve_candidates = val;
					sieve_candidates += (val != number);
				}
				for (size_t i = 0; i < 4; ++i)
				{
					auto val = ymm1.m256i_u64[i];
					*sieve_candidates = val;
					sieve_candidates += (val != number);
				}
				for (size_t i = 0; i < 4; ++i)
				{
					auto val = ymm2.m256i_u64[i];
					*sieve_candidates = val;
					sieve_candidates += (val != number);
				}
				for (size_t i = 0; i < 4; ++i)
				{
					auto val = ymm3.m256i_u64[i];
					*sieve_candidates = val;
					sieve_candidates += (val != number);
				}
				for (size_t i = 0; i < 4; ++i)
				{
					auto val = ymm4.m256i_u64[i];
					*sieve_candidates = val;
					sieve_candidates += (val != number);
				}
				for (size_t i = 0; i < 4; ++i)
				{
					auto val = ymm5.m256i_u64[i];
					*sieve_candidates = val;
					sieve_candidates += (val != number);
				}
				for (size_t i = 0; i < 4; ++i)
				{
					auto val = ymm6.m256i_u64[i];
					*sieve_candidates = val;
					sieve_candidates += (val != number);
				}
				for (size_t i = 0; i < 4; ++i)
				{
					auto val = ymm7.m256i_u64[i];
					*sieve_candidates = val;
					sieve_candidates += (val != number);
				}

				std::cout << "-----------\n";
				std::cin.ignore();
			}

			// After the above loop, we still need cleanup at the end - the sieve is not a multiple of 15*4
			*/
		}
	}
}
