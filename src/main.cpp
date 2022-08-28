
#define VCL_NAMESPACE vcl
#include "../lib/vcl/vectorclass.h"

#include <bitset>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#include "config.hpp"
#include "io/io.hpp"
#include "math/franken_mpir.hpp"
#include "math/math.hpp"
#include "trial_division/multibase_div_tests.hpp"
#include "util/find.hpp"
#include "util/sandbox.hpp"
#include "util/types.hpp"
#include "util/utility.hpp"


namespace mbp
{
	constexpr size_t static_sieve_size = std::accumulate(static_sieve_primes.begin(),
														 static_sieve_primes.end(), size_t(1), std::multiplies());

	using sieve_t = uint8_t;
	using sieve_container = std::vector<sieve_t>;
	// using sieve_container = std::array<sieve_t, static_sieve_size>;
	const sieve_container generate_static_sieve()
	{
		sieve_container sieve(static_sieve_size, true);
		// sieve_container sieve{}; std::fill(begin(sieve), end(sieve), true);

		// for each prime, mark off all multiples
		for (const auto p : static_sieve_primes)
			for (size_t i = 0; i < sieve.size(); i += p)
				sieve[i] = false;

		return sieve;
	}
	const sieve_container static_sieve = generate_static_sieve();

	using sieve_offset_t = narrowest_uint_for_val<static_sieve_size>;
	std::vector<sieve_offset_t> sieve_offsets_cache(small_primes_lookup.size());

	void set_up_sieve_offsets_cache(const size_t start)
	{
		// Start with the first prime not in the static sieve.
		for (size_t i = static_sieve_primes.size() + 1; i < small_primes_lookup.size(); ++i)
		{
			const sieve_prime_t p = small_primes_lookup[i];

			// Find out how far it is to the next multiple of p.
			sieve_prime_t n = p - (start % p);

			// Start is always odd. Therefore, if n is odd, it is pointing to the next even multiple of p.
			// -- increase by p
			if (n % 2 == 1)
				n += p;
			// However, if n is even, it is pointing to the next odd multiple of p.
			// -- do nothing

			// We now have the distance to the next odd multiple of p.
			// Divide by 2 to store the *index* of the next odd multiple of p.
			sieve_offsets_cache[i] = n / 2;
		}
	}

	void partial_sieve(sieve_container& sieve)
	{
		static_assert(static_sieve_size > sieve_primes_cap);

		sieve_t* begin = sieve.data();
		const sieve_t* end = begin + sieve.size();

		sieve_prime_t* cache_ptr = sieve_offsets_cache.data() + static_sieve_primes.size() + 1;

		// Start with the first prime not in the static sieve
		for (const sieve_prime_t* prime_ptr = small_primes_lookup.data() + static_sieve_primes.size() + 1;
			 prime_ptr < small_primes_lookup.data() + small_primes_lookup.size();
			 ++prime_ptr, ++cache_ptr)
		{
			// Get the next prime
			const size_t p = *prime_ptr;

			// Get the index of the next odd multiple of p
			sieve_t* j = begin + *cache_ptr;

			// Mark false each (implicitly odd) multiple of p
			do
			{
				*j = false;
				j += p;
			} while (j < end);

			// Update the cache for the next sieving
			*cache_ptr = sieve_prime_t(j - end);
		}
	}



	consteval auto generate_bitmask_lookup()
	{
		// +1 so bitmasks[n_of_rems] is always safe
		std::array<size_t, div_test::max_remainders + 1> bitmasks = { 0 };

		for (size_t i = 1; i < bitmasks.size(); ++i)
		{
			size_t bitmask = 0;
			for (size_t j = 0; j < 64; j += i)
			{
				bitmask <<= i;
				bitmask |= 1;
			}
			bitmasks[i] = bitmask;
		}

		return bitmasks;
	}
	constexpr std::array<size_t, div_test::max_remainders + 1> bitmask_lookup = generate_bitmask_lookup();

	// takes N^2 memory, even though we only need (N^2) / 2
	using popcount_t = uint16_t;
	static std::array<mbp::aligned64<popcount_t, 64>, 64> popcounts{};

	__forceinline bool has_small_divisor(const size_t number)
	{
		using namespace div_test;

		if (recursive_is_divisible_by<5, in_base<3>>(number)) return true;

		if (recursive_is_divisible_by<7, in_base<3>>(number)) return true;
		if (recursive_is_divisible_by<7, in_base<4>>(number)) return true;
		if (recursive_is_divisible_by<7, in_base<5>>(number)) return true;

	#if analyze_div_tests
		bool found_div = false;
	#endif

		bool which_way_boss = div_tests.begin()->is_first_with_n_remainders;
		n_of_remainders_t n_of_rems_boss = div_tests.begin()->n_of_remainders;

		for (div_test_const auto& div_test : div_tests)
		{
			size_t rem = 0;

			//const size_t n_of_rems = div_test.n_of_remainders;
			//__assume(n_of_rems > 0);
			//__assume(n_of_rems <= max_remainders);

			const auto& my_rems = div_test.remainders;

			__assume(n_of_rems_boss > 0);
			__assume(n_of_rems_boss <= max_remainders);


			if (which_way_boss)
			{
				which_way_boss = *((&div_test.is_first_with_n_remainders) + sizeof(div_test));

				const size_t my_bitmask = bitmask_lookup[n_of_rems_boss];
				auto& my_pcs = popcounts[n_of_rems_boss];

				// for switch (n), run cases n through 1, where the index is n-1 through 0
				constexpr size_t start = __LINE__ + 10;
			#define IDX(n) ((max_remainders - (n - start)) - 1)
			#define CASE(n) [[fallthrough]]; case(IDX(n) + 1): \
				{ \
				const auto pc = pop_count(number & (my_bitmask << IDX(n))); \
					my_pcs[IDX(n)] = popcount_t(pc); \
					rem += pc * my_rems[IDX(n)]; \
				}
				switch (n_of_rems_boss) // handle cases N through 1
				{
					CASE(__LINE__); // case (max)
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__); // case (1)
					static_assert(start + max_remainders == __LINE__);
					break;
				default:
					__assume(false);
				}
			#undef CASE
			#undef IDX
			}
			else
			{
				which_way_boss = *((&div_test.is_first_with_n_remainders) + sizeof(div_test));

				const auto& my_pcs = popcounts[n_of_rems_boss];

				// for switch (n), run cases n through 1, where the index is n-1 through 0
				constexpr size_t start = __LINE__ + 5;
			#define IDX(n) ((max_remainders - (n - start)) - 1)
			#define CASE(n) [[fallthrough]]; case(IDX(n) + 1): rem += size_t(my_pcs[IDX(n)]) * my_rems[IDX(n)];
				switch (n_of_rems_boss)
				{
					CASE(__LINE__); // case (max)
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__);
					CASE(__LINE__); // case (1)
					static_assert(start + max_remainders == __LINE__);
					break;
				default:
					__assume(false);
				}
			#undef CASE
			#undef IDX
			}

			n_of_rems_boss = *((&div_test.n_of_remainders) + sizeof(div_test));

			if (has_small_prime_factor(rem, div_test.prime_idx))
			{
			#if analyze_div_tests
				div_test.hits++;
				found_div = true;
				return true;
			#else
				return true;
			#endif
			}
		}

	#if analyze_div_tests
		return found_div;
	#else
		return false;
	#endif
	}

	void print_div_tests()
	{
	#if analyze_div_tests
		using namespace div_test;

		//std::sort(div_tests.begin(), div_tests.end(), [] (const auto& a, const auto& b)
		//		  {
		//			  //return a.hits < b.hits;
		//			  //return a.n_of_remainders < b.n_of_remainders;

		//			  if (a.n_of_remainders == b.n_of_remainders)
		//				  return a.base < b.base;
		//			  else
		//				  return a.n_of_remainders < b.n_of_remainders;

		//			  //if (a.prime_idx == b.prime_idx)
		//				 // return a.base < b.base;
		//			  //else
		//				 // return a.prime_idx < b.prime_idx;

		//			  //return
		//				 // a.hits / a.n_of_remainders <
		//				 // b.hits / b.n_of_remainders;
		//		  });

		auto w = std::setw;
		for (const auto& dt : div_tests)
		{
			std::cout << "   base " << std::setfill(' ') << w(2) << size_t(dt.base) << " % " << w(3) << size_t(small_primes_lookup[dt.prime_idx]) << ":  ";
			if (dt.hits == 0)
			{
				std::cout << "       -       ";
			}
			else
			{
				std::cout << w(8) << dt.hits << " hits  ";
			}

			std::cout << w(2) << size_t(dt.n_of_remainders) << " remainders: 1";
			for (size_t j = 1; j < dt.n_of_remainders; ++j)
			{
				std::cout << ' ' << w(3) << dt.remainders[j];
				if (j == 20)
				{
					std::cout << " ...";
					break;
				}
			}
			std::cout << '\n';
		}
	#endif
	}

	void run_div_test_analysis()
	{
	#if analyze_div_tests
		using namespace div_test;

		auto div_test_pred = [](auto a, auto b) { return a.hits < b.hits; };

		if (std::is_sorted(div_tests.rbegin(), div_tests.rend(), div_test_pred))
		{
			std::cout << "Div tests have not changed frequency ordering\n";
		}
		else
		{
			const auto copy = div_tests;

			// sort descending
			std::sort(div_tests.rbegin(), div_tests.rend(), div_test_pred);

			size_t moved = 0;
			for (size_t i = 0; i < div_tests.size(); ++i)
				if (div_tests[i].base != copy[i].base || div_tests[i].prime_idx != copy[i].prime_idx)
					moved++;

			// clear "firsts"
			for (auto& dt : div_tests)
			{
				dt.is_first_with_n_remainders = false;
			}

			// set firsts
			for (size_t i = 0; i <= div_test::max_remainders; ++i)
			{
				for (auto& div_test : div_tests)
				{
					if (div_test.n_of_remainders == i)
					{
						div_test.is_first_with_n_remainders = true;
						break;
					}
				}
			}

			print_div_tests();

			static std::stringstream ss;
			ss << ' ' << moved;

			std::cout << moved << " div tests changed position\n";
			std::cout << '(' << ss.str() << ")\n";
		}
	#endif
	}



	void find_multibase_primes()
	{
		gmp_randclass r{ gmp_randinit_mt };
		r.seed(mpir_ui{ 0xdeadbeef });

		size_t number = benchmark_mode ? bm_start : load_from_results();
		mpz_class mpz_number = 0ull; // it's a surprise tool that will help us later

		sieve_container sieve = static_sieve;
		if ((size_t)sieve.data() % 8 != 0) { std::cout << "unaligned!"; std::cin.ignore(); }

		// Round starting number down to the nearest odd multiple of the sieve sieze
		number -= static_sieve.size(); // n -= k
		number -= number % (2 * static_sieve.size()); // n -= n % 2k
		number += static_sieve.size(); // n += k

		set_up_sieve_offsets_cache(number);

		constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();
		constexpr size_t gcd_lookup = build_gcd_lookup();

	#if analyze_div_tests
		const size_t div_test_log_interval = 1'000'000;
		size_t next_div_test_checkpoint = div_test_log_interval;
	#endif

		constexpr size_t scratch_size = [] {
			double cleared = 0.0;
			for (size_t i = 1; i < small_primes_lookup.size(); ++i)
				cleared += ((1.0 - cleared) * (1.0 / small_primes_lookup[i]));
			return size_t((1.0 - cleared) * double(static_sieve_size) * 1.1);
		}();
		static size_t scratch[scratch_size]{};

		// Start the clock after setup
		const auto start = current_time_in_ms();



		count_passes(size_t a, b, c, d);
		count_passes(a = b = c = d = 0);

		// (condition optimizes out when not benchmarking)
		while (benchmark_mode ? number < bm_stop : true)
		{
			// Perform additional sieving on the static sieve
			sieve = static_sieve;
			partial_sieve(sieve);

			const size_t number_before_loop = number;



			// 1. Collect candidates that have not been marked composite by the sieve.

			// Safe to move these higher? Can v1 = const_v2 ever move v1?
			const char* const sieve_begin = (const char*)sieve.data();
			const char* const sieve_end = sieve_begin + sieve.size();
			const char* const sieve_end_rounded = sieve_end - (sieve.size() % 8);

			const char* sieve_ptr = sieve_begin;

			size_t* sieve_candidates = scratch;
			for (; sieve_ptr < sieve_end_rounded; sieve_ptr += 8, number += 16)
			{
				*sieve_candidates = number;
				sieve_candidates += *sieve_ptr;

				*sieve_candidates = number + 2;
				sieve_candidates += *(sieve_ptr + 1);

				*sieve_candidates = number + 4;
				sieve_candidates += *(sieve_ptr + 2);

				*sieve_candidates = number + 6;
				sieve_candidates += *(sieve_ptr + 3);

				*sieve_candidates = number + 8;
				sieve_candidates += *(sieve_ptr + 4);

				*sieve_candidates = number + 10;
				sieve_candidates += *(sieve_ptr + 5);

				*sieve_candidates = number + 12;
				sieve_candidates += *(sieve_ptr + 6);

				*sieve_candidates = number + 14;
				sieve_candidates += *(sieve_ptr + 7);
			}

			// handle last few elements
			for (; sieve_ptr < sieve_end; ++sieve_ptr, number += 2)
			{
				*sieve_candidates = number;
				sieve_candidates += *sieve_ptr;
			}

			count_passes(a += (sieve_candidates - scratch)); // How many ints passed the sieve?



			// 2. Collect candidates that have a prime number of bits set

			const size_t* const candidates_end = sieve_candidates; // already points one past the end
			const size_t* const candidates_end_rounded = candidates_end - ((candidates_end - scratch) % 4);
			const size_t* candidate_ptr = scratch;

			size_t* passed_pc_test_ptr = scratch;
			for (; candidate_ptr < candidates_end_rounded; candidate_ptr += 4)
			{
				const size_t c0 = *candidate_ptr;
				const size_t c1 = *(candidate_ptr + 1);
				const size_t c2 = *(candidate_ptr + 2);
				const size_t c3 = *(candidate_ptr + 3);

				const size_t pc_c0 = pop_count(c0);
				const size_t pc_c1 = pop_count(c1);
				const size_t pc_c2 = pop_count(c2);
				const size_t pc_c3 = pop_count(c3);

				*passed_pc_test_ptr = c0;
				passed_pc_test_ptr += bool(tiny_primes_lookup & (1ull << pc_c0));

				*passed_pc_test_ptr = c1;
				passed_pc_test_ptr += bool(tiny_primes_lookup & (1ull << pc_c1));

				*passed_pc_test_ptr = c2;
				passed_pc_test_ptr += bool(tiny_primes_lookup & (1ull << pc_c2));

				*passed_pc_test_ptr = c3;
				passed_pc_test_ptr += bool(tiny_primes_lookup & (1ull << pc_c3));
			}

			// handle last few elements
			for (; candidate_ptr < candidates_end; ++candidate_ptr)
			{
				const size_t c0 = *candidate_ptr;
				*passed_pc_test_ptr = c0;
				passed_pc_test_ptr += bool(tiny_primes_lookup & (1ull << pop_count(c0)));
			}

			count_passes(b += (passed_pc_test_ptr - scratch));



			// 3. Collect candidates with an alternating bitsum that shares a GCD of 1 with a product of primes

			const size_t* const gcd_candidates_end = passed_pc_test_ptr;
			const size_t* const gcd_candidates_end_rounded = gcd_candidates_end - ((gcd_candidates_end - scratch) % 4);
			candidate_ptr = scratch; // reset

			size_t* passed_gcd_test_ptr = scratch;
			for (; candidate_ptr < gcd_candidates_end_rounded; candidate_ptr += 4)
			{
				const size_t n0 = *candidate_ptr;
				const size_t n0_pca = pop_count(n0 & 0xAAAAAAAAAAAAAAAA);
				const size_t n0_pcb = pop_count(n0 & 0x5555555555555555);
				*passed_gcd_test_ptr = n0;
				if (gcd_lookup & (1ull << (n0_pca + 32 - n0_pcb))) ++passed_gcd_test_ptr;

				const size_t n1 = *(candidate_ptr + 1);
				const size_t n1_pca = pop_count(n1 & 0xAAAAAAAAAAAAAAAA);
				const size_t n1_pcb = pop_count(n1 & 0x5555555555555555);
				*passed_gcd_test_ptr = n1;
				if (gcd_lookup & (1ull << (n1_pca + 32 - n1_pcb))) ++passed_gcd_test_ptr;

				const size_t n2 = *(candidate_ptr + 2);
				const size_t n2_pca = pop_count(n2 & 0xAAAAAAAAAAAAAAAA);
				const size_t n2_pcb = pop_count(n2 & 0x5555555555555555);
				*passed_gcd_test_ptr = n2;
				if (gcd_lookup & (1ull << (n2_pca + 32 - n2_pcb))) ++passed_gcd_test_ptr;

				const size_t n3 = *(candidate_ptr + 3);
				const size_t n3_pca = pop_count(n3 & 0xAAAAAAAAAAAAAAAA);
				const size_t n3_pcb = pop_count(n3 & 0x5555555555555555);
				*passed_gcd_test_ptr = n3;
				if (gcd_lookup & (1ull << (n3_pca + 32 - n3_pcb))) ++passed_gcd_test_ptr;
			}

			// handle last few elements
			for (; candidate_ptr < gcd_candidates_end; ++candidate_ptr)
			{
				const size_t n0 = *candidate_ptr;
				const auto n0_pca = pop_count(n0 & 0xAAAAAAAAAAAAAAAA);
				const auto n0_pcb = pop_count(n0 & 0x5555555555555555);

				*passed_gcd_test_ptr = n0;
				passed_gcd_test_ptr += bool(gcd_lookup & (1ull << (n0_pca + 32 - n0_pcb)));
			}



			for (size_t* n = scratch; n < passed_gcd_test_ptr; ++n)
			{
				number = *n;

				count_passes(++c);

				// Bail if n has a small prime factor in any base
				if (has_small_divisor(number)) continue;

				count_passes(++d);



				// Do full primality tests, starting with base 2
				if (!franken::mpir_is_likely_prime_BPSW(number)) continue;

				// convert uint64_t to char array of ['0', '1'...] for MPIR
				char bin_str[64 + 1];
				auto result = std::to_chars(&bin_str[0], &bin_str[64], number, 2);
				*result.ptr = '\0';

				mpz_number.set_str(bin_str, 3);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 4);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 5);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 6);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 7);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 8);
				if (!franken::mpir_is_prime(mpz_number, r, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 9);
				if (!mpir_is_prime(mpz_number, r)) { log_result(number, 8); continue; }

				mpz_number.set_str(bin_str, 10);
				if (!mpir_is_prime(mpz_number, r)) { log_result(number, 9); continue; }

				mpz_number.set_str(bin_str, 11);
				if (!mpir_is_prime(mpz_number, r)) { log_result(number, 10); continue; }

				mpz_number.set_str(bin_str, 12);
				if (!mpir_is_prime(mpz_number, r)) { log_result(number, 11); continue; }

				mpz_number.set_str(bin_str, 13);
				if (!mpir_is_prime(mpz_number, r)) { log_result(number, 12); continue; }

			} // end hot loop

			number = number_before_loop + 2 * sieve.size();



		#if analyze_div_tests
			for (const auto& dt : div_test::div_tests)
			{
				if (dt.hits >= next_div_test_checkpoint)
				{
					run_div_test_analysis();
					next_div_test_checkpoint += div_test_log_interval;
					break;
				}
			}
		#endif

		} // end outer loop

		std::cout << "Finished. " << current_time_in_ms() - start << " ms elapsed\n";

		count_passes(std::cout << "Passed sieve:      " << std::setw(9) << a << '\n');
		count_passes(std::cout << "Passed prime test: " << std::setw(9) << b << '\n');
		count_passes(std::cout << "Passed GCD test:   " << std::setw(9) << c << '\n');
		count_passes(std::cout << "Passed div tests:  " << std::setw(9) << d << '\n');
	}

} // namespace mbp



int main()
{
	mbp::print_preamble();

	mbp::find_multibase_primes();
}
