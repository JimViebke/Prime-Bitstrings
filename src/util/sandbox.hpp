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
			const auto start = util::current_time_in_ms();

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
			std::cout << "Found in " << util::current_time_in_ms() - start << " ms.\n";
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

		void expiremental_SIMD_gather_sieve_results(const uint8_t* const sieve_ptr, const size_t number)
		{
			// Set some magic numbers
			std::cout << "Indices lookup:\n";
			const uint256_t sieve_gather_indices = _mm256_broadcastq_epi64(uint128_t{ .m128i_u8 = { 1, 2, 4, 7, 8, 11, 13, 14 } });
			util::print_vector_as<uint8_t>(sieve_gather_indices);

			std::cout << "block_offsets lookup:\n";
			const uint256_t block_offsets = { .m256i_u8 = {
				0, 0, 0, 0, 0, 0, 0, 0,
				15, 15, 15, 15, 15, 15, 15, 15,
				30, 30, 30, 30, 30, 30, 30, 30,
				45, 45, 45, 45, 45, 45, 45, 45 } };
			util::print_vector_as<uint8_t>(block_offsets);

			const uint256_t ymm_number = _mm256_broadcastq_epi64(uint128_t{ .m128i_u64 = { number } });

			//for (; sieve_ptr < sieve_end; sieve_ptr += 60, number += 120)
			//{
				// Load 15 + 1 bytes from the sieve
			const uint128_t bytes_a = _mm_loadu_si128((const uint128_t*)sieve_ptr);
			const uint128_t bytes_b = _mm_loadu_si128((const uint128_t*)(sieve_ptr + 15));
			const uint128_t bytes_c = _mm_loadu_si128((const uint128_t*)(sieve_ptr + 30));
			const uint128_t bytes_d = _mm_loadu_si128((const uint128_t*)(sieve_ptr + 45));

			// Combine to two registers
			uint256_t ymm_ab = _mm256_inserti128_si256(_mm256_castsi128_si256(bytes_a), bytes_b, 1);
			uint256_t ymm_cd = _mm256_inserti128_si256(_mm256_castsi128_si256(bytes_c), bytes_d, 1);
			std::cout << "Aligned sieve bytes:\n";
			util::print_vector_as<uint8_t>(ymm_ab);
			util::print_vector_as<uint8_t>(ymm_cd);

			// Use the magic numbers to gather the bytes we want
			ymm_ab = _mm256_shuffle_epi8(ymm_ab, sieve_gather_indices);
			ymm_cd = _mm256_shuffle_epi8(ymm_cd, sieve_gather_indices);
			std::cout << "Gather relevant bytes into halves of each lane:\n";
			util::print_vector_as<uint8_t>(ymm_ab);
			util::print_vector_as<uint8_t>(ymm_cd);

			// Combine to one register
			uint256_t ymm_abcd = _mm256_blend_epi32(ymm_ab, ymm_cd, 0b00'11'00'11);
			std::cout << "Blend ab + cd:\n";
			util::print_vector_as<uint8_t>(ymm_abcd);

			// Use ymm_abcd to filter for the relevant block offsets and indicies
			const uint256_t filtered_block_offsets = _mm256_sign_epi8(block_offsets, ymm_abcd);
			const uint256_t filtered_indices = _mm256_sign_epi8(sieve_gather_indices, ymm_abcd);
			std::cout << "Filtered block offsets...:\n";
			util::print_vector_as<uint8_t>(filtered_block_offsets);
			std::cout << "...filtered indices:\n";
			util::print_vector_as<uint8_t>(filtered_indices);

			// add filters
			ymm_abcd = _mm256_add_epi8(filtered_block_offsets, filtered_indices);
			std::cout << "Offsets from sieve ptr:\n";
			util::print_vector_as<uint8_t>(ymm_abcd);

			// << 1 to multiply all offsets by 2
			ymm_abcd = _mm256_slli_epi64(ymm_abcd, 1);
			std::cout << "Offsets from N:\n";
			util::print_vector_as<uint8_t>(ymm_abcd);

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
			util::print_vector_as<uint8_t>(ymm0);
			util::print_vector_as<uint8_t>(ymm1);
			util::print_vector_as<uint8_t>(ymm2);
			util::print_vector_as<uint8_t>(ymm3);
			util::print_vector_as<uint8_t>(ymm4);
			util::print_vector_as<uint8_t>(ymm5);
			util::print_vector_as<uint8_t>(ymm6);
			util::print_vector_as<uint8_t>(ymm7);

			// 0 already in place
			ymm1 = _mm256_srli_si256(ymm1, 1 * 4);
			ymm2 = _mm256_srli_si256(ymm2, 2 * 4);
			ymm3 = _mm256_srli_si256(ymm3, 3 * 4);
			// 4 already in place
			ymm5 = _mm256_srli_si256(ymm5, 1 * 4);
			ymm6 = _mm256_srli_si256(ymm6, 2 * 4);
			ymm7 = _mm256_srli_si256(ymm7, 3 * 4);
			std::cout << "Align sets of four byte offsets:\n";
			util::print_vector_as<uint8_t>(ymm0);
			util::print_vector_as<uint8_t>(ymm1);
			util::print_vector_as<uint8_t>(ymm2);
			util::print_vector_as<uint8_t>(ymm3);
			util::print_vector_as<uint8_t>(ymm4);
			util::print_vector_as<uint8_t>(ymm5);
			util::print_vector_as<uint8_t>(ymm6);
			util::print_vector_as<uint8_t>(ymm7);

			ymm0 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm0));
			ymm1 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm1));
			ymm2 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm2));
			ymm3 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm3));
			ymm4 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm4));
			ymm5 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm5));
			ymm6 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm6));
			ymm7 = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(ymm7));
			std::cout << "Zero-extend bytes into 64-bit positions:\n";
			util::print_vector_as<uint8_t>(ymm0);
			util::print_vector_as<uint8_t>(ymm1);
			util::print_vector_as<uint8_t>(ymm2);
			util::print_vector_as<uint8_t>(ymm3);
			util::print_vector_as<uint8_t>(ymm4);
			util::print_vector_as<uint8_t>(ymm5);
			util::print_vector_as<uint8_t>(ymm6);
			util::print_vector_as<uint8_t>(ymm7);

			ymm0 = _mm256_add_epi64(ymm0, ymm_number);
			ymm1 = _mm256_add_epi64(ymm1, ymm_number);
			ymm2 = _mm256_add_epi64(ymm2, ymm_number);
			ymm3 = _mm256_add_epi64(ymm3, ymm_number);
			ymm4 = _mm256_add_epi64(ymm4, ymm_number);
			ymm5 = _mm256_add_epi64(ymm5, ymm_number);
			ymm6 = _mm256_add_epi64(ymm6, ymm_number);
			ymm7 = _mm256_add_epi64(ymm7, ymm_number);
			std::cout << "Actually add the numbers:\n";
			util::print_vector_as<uint64_t>(ymm0);
			util::print_vector_as<uint64_t>(ymm1);
			util::print_vector_as<uint64_t>(ymm2);
			util::print_vector_as<uint64_t>(ymm3);
			util::print_vector_as<uint64_t>(ymm4);
			util::print_vector_as<uint64_t>(ymm5);
			util::print_vector_as<uint64_t>(ymm6);
			util::print_vector_as<uint64_t>(ymm7);

			//for (size_t i = 0; i < 4; ++i)
			//{
			//	auto val = ymm0.m256i_u64[i];
			//	*sieve_candidates = val;
			//	sieve_candidates += (val != number);
			//}
			//for (size_t i = 0; i < 4; ++i)
			//{
			//	auto val = ymm1.m256i_u64[i];
			//	*sieve_candidates = val;
			//	sieve_candidates += (val != number);
			//}
			//for (size_t i = 0; i < 4; ++i)
			//{
			//	auto val = ymm2.m256i_u64[i];
			//	*sieve_candidates = val;
			//	sieve_candidates += (val != number);
			//}
			//for (size_t i = 0; i < 4; ++i)
			//{
			//	auto val = ymm3.m256i_u64[i];
			//	*sieve_candidates = val;
			//	sieve_candidates += (val != number);
			//}
			//for (size_t i = 0; i < 4; ++i)
			//{
			//	auto val = ymm4.m256i_u64[i];
			//	*sieve_candidates = val;
			//	sieve_candidates += (val != number);
			//}
			//for (size_t i = 0; i < 4; ++i)
			//{
			//	auto val = ymm5.m256i_u64[i];
			//	*sieve_candidates = val;
			//	sieve_candidates += (val != number);
			//}
			//for (size_t i = 0; i < 4; ++i)
			//{
			//	auto val = ymm6.m256i_u64[i];
			//	*sieve_candidates = val;
			//	sieve_candidates += (val != number);
			//}
			//for (size_t i = 0; i < 4; ++i)
			//{
			//	auto val = ymm7.m256i_u64[i];
			//	*sieve_candidates = val;
			//	sieve_candidates += (val != number);
			//}

			std::cout << "-----------\n";
			std::cin.ignore();
			//}

			// After the above loop, we still need cleanup at the end - the sieve is not a multiple of 15*4
		}

		void simd_popcount_64(size_t* output, size_t* input, const size_t* const end)
		{
			static constexpr uint256_t static_nibble_popcounts = []() consteval {
				uint256_t pcs{ .m256i_u8 = { 0 } };
				for (size_t i = 0; i < 16; ++i)
					pcs.m256i_u8[i] = pcs.m256i_u8[i + 16] = uint8_t(std::popcount(i));
				return pcs;
			}();

			const uint256_t nibble_popcounts = _mm256_loadu_si256((uint256_t*)&static_nibble_popcounts);
			const uint256_t nibble_mask = _mm256_set1_epi8(0b0000'1111);

			// load one iteration ahead
			uint256_t ymm_next = _mm256_loadu_si256((uint256_t*)input);

			for (; input < end; )
			{
				// shift high four bits right by four, so both registers have
				// byte values in range 0-F inclusive
				uint256_t ymm0 = _mm256_and_si256(ymm_next, nibble_mask);
				uint256_t ymm1 = _mm256_and_si256(_mm256_srli_epi64(ymm_next, 4), nibble_mask);

				// load one iteration ahead
				input += 4;
				ymm_next = _mm256_loadu_si256((uint256_t*)input);

				// replace all bytes (values 0-F) with their popcount (values 0-4)
				ymm0 = _mm256_shuffle_epi8(ymm0, nibble_popcounts);
				ymm1 = _mm256_shuffle_epi8(ymm1, nibble_popcounts);

				// add low and high together
				ymm0 = _mm256_add_epi8(ymm0, ymm1);

				// sum blocks of eight bytes
				ymm0 = _mm256_sad_epu8(ymm0, _mm256_setzero_si256());

				// extract four popcounts
				const auto pc_a = _mm256_extract_epi8(ymm0, 0);
				const auto pc_b = _mm256_extract_epi8(ymm0, 8);
				// extract the high half
				uint128_t xmm1 = _mm256_extracti128_si256(ymm0, 1);
				const auto pc_c = _mm_extract_epi8(xmm1, 0);
				const auto pc_d = _mm_extract_epi8(xmm1, 8);

				*(output + 0) = pc_a;
				*(output + 1) = pc_b;
				*(output + 2) = pc_c;
				*(output + 3) = pc_d;
				output += 4;

				// After calculating the four popcounts, we still need the original four
				// candidates (four loads, four cmovs, four stores)
			}

		}

		void test_simd_popcount()
		{

			std::array<size_t, 100> input{};
			std::array<size_t, 100> output{};
			output.fill(0);

			size_t some_val = 0;
			for (size_t i = 0; i < 100; ++i)
			{
				std::cin >> some_val;
				input[i] = some_val;

				if (some_val == 7) break;
			}

			simd_popcount_64(output.data(), input.data(), input.data() + 100);

			for (size_t i = 0; i < 100; ++i)
			{
				std::cout << output[i] << ' ';
			}
		}

		void denser_sieves()
		{
			std::cout << "Testing denser sieves...\n";

			std::array<uint8_t, 100> sieve{};
			sieve.fill(true);
			const size_t start = 11; // the starting point, the 0th value of the sieve

			constexpr std::array primes = { 3, 5, 7 };

			for (auto p : primes)
			{
				// Find out how far it is to the next multiple of p.
				size_t n = p - (start % p); // n = 3 - (11 % 3)
				// n = 1

				// odd multiples of 3: 99, 105

				// Start is always odd. Therefore, if n is odd, start + n is pointing to the next even multiple of p.
				// -- increase by p
				if (n % 2 == 1)
					n += p;
				// before: start + n == 11 + 1 == 12 (even multiple)
				// after:  start + n == 11 + 4 == 15 (next odd multiple of p)

				// We now have the distance to the next odd multiple of p.
				// Divide by 2 to store the *index* of the next odd multiple of p.
				std::cout << "Clear multiples of " << p << ":";
				for (size_t i = n / 2; i < sieve.size(); i += p)
				{
					sieve[i] = false;
					std::cout << ' ' << start + 2 * i;
				}
				std::cout << '\n';
			}

			std::cout << "\nDone. Prime candidates: ";
			for (size_t i = 0; i < sieve.size(); ++i)
				if (sieve[i])
					std::cout << start + 2 * i << ' ';
			std::cout << '\n';



			// for each prime N, calculate the gaps between multiples of N

			std::cout << "\n\n\n";

			std::vector<size_t> packed_number_line;
			for (size_t i = 3; packed_number_line.size() < 1'000; ++i)
				if (i % 2 != 0 && i % 3 != 0 && i % 5 != 0)
					packed_number_line.push_back(i);

			std::cout << "Packed number line: (no 2s, 3s, or 5s)\n";
			for (size_t i = 0; i < packed_number_line.size(); ++i)
			{
				if (i % 8 == 0) std::cout << '\n';

				std::cout << std::setw(5) << packed_number_line[i];
			}
			std::cout << '\n';

			auto print_gaps_for = [packed_number_line](int prime_factor) {
				size_t counter = 0;
				bool found_first = false;
				for (auto n : packed_number_line)
				{
					++counter;
					if (n % prime_factor == 0)
					{
						if (!found_first) // don't print before the first multiple
						{
							found_first = true;
							counter = 0;
							continue;
						}
						std::cout << ' ' << counter;
						counter = 0;
					}
				}
				std::cout << '\n';
			};

			std::cout << "Gaps for 7:";
			print_gaps_for(7);

			std::cout << "Gaps for 11:";
			print_gaps_for(11);

			std::cout << "Gaps for 13:";
			print_gaps_for(13);

			std::cout << "Gaps for 17:";
			print_gaps_for(17);

			std::cout << "Gaps for 19:";
			print_gaps_for(19);

			std::cout << "Gaps for 23:";
			print_gaps_for(23);

			std::cout << "Gaps for 29:";
			print_gaps_for(29);

			std::cout << "Gaps for 31:";
			print_gaps_for(31);

			std::cout << "Gaps for 37:";
			print_gaps_for(37);



		} // denser_sieves()

	} // namespace detail

} // namespace mbp
