
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
#include "multibase_div_tests.hpp"
#include "util/sandbox.hpp"
#include "util/types.hpp"
#include "util/utility.hpp"


namespace mbp
{
	constexpr size_t static_sieve_size = std::accumulate(static_sieve_primes.begin(),
														 static_sieve_primes.end(), size_t(1), std::multiplies());

	using sieve_t = uint8_t;
	const std::vector<sieve_t> generate_static_sieve()
	{
		std::vector<sieve_t> sieve(static_sieve_size, true);

		// for each prime, mark off all multiples
		for (const auto p : static_sieve_primes)
			for (size_t i = 0; i < sieve.size(); i += p)
				sieve[i] = false;

		return sieve;
	}

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

	void partial_sieve(std::vector<sieve_t>& sieve)
	{
		static_assert(static_sieve_size > sieve_primes_cap);

		const sieve_offset_t sieve_size = sieve_offset_t(sieve.size());

		// Start with the first prime not in the static sieve
		for (size_t i = static_sieve_primes.size() + 1; i < small_primes_lookup.size(); ++i)
		{
			const sieve_prime_t p = small_primes_lookup[i];

			// Get the index of the next odd multiple of p
			sieve_offset_t j = sieve_offsets_cache[i];

			// Mark false each (implicitly odd) multiple of p
			for (; j < sieve_size; j += p)
			{
				sieve[j] = false;
			}

			// Update the cache for the next sieving
			sieve_offsets_cache[i] = j - static_sieve_size;
		}
	}



	consteval auto generate_bitmask_lookup()
	{
		std::array<size_t, 64> bitmasks = { 0 };

		for (size_t i = 1; i < 64; ++i)
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
	constexpr std::array<size_t, 64> bitmask_lookup = generate_bitmask_lookup();

	using div_tests_t = std::array<div_test::div_test_t, div_test::div_tests_size>;
	static div_test_constexpr div_tests_t div_tests = div_test::generate_div_tests(); // intellisense false positive

	// takes N^2 memory, even though we only need (N^2) / 2
	static std::array<mbp::aligned64, 64> popcounts{};

	__forceinline bool has_small_divisor_cached(const size_t number)
	{
		using namespace div_test;

		if (recursive_is_divisible_by<5, in_base<3>>(number)) return true;

		if (recursive_is_divisible_by<7, in_base<3>>(number)) return true;
		if (recursive_is_divisible_by<7, in_base<4>>(number)) return true;
		if (recursive_is_divisible_by<7, in_base<5>>(number)) return true;

	#if analyze_div_tests
		bool found_div = false;
	#endif

		for (div_test_const auto& div_test : div_tests)
		{
			size_t rem = 0;

			const size_t n_of_rems = div_test.n_of_remainders;
			__assume(n_of_rems > 0);
			__assume(n_of_rems <= max_remainders);

			const auto& my_rems = div_test.remainders;

			if (div_test.is_first_with_n_remainders)
			{
				const size_t my_bitmask = bitmask_lookup[n_of_rems];
				auto& my_pcs = popcounts[n_of_rems];

				// for switch (n), run cases n through 1, where the index is n-1 through 0
				constexpr size_t start = __LINE__ + 10;
			#define IDX(n) ((max_remainders - (n - start)) - 1)
			#define CASE(n) [[fallthrough]]; case(IDX(n) + 1): \
				{ \
				const auto pc = pop_count(number & (my_bitmask << IDX(n))); \
					my_pcs[IDX(n)] = uint8_t(pc); \
					rem += pc * my_rems[IDX(n)]; \
				}
				switch (n_of_rems) // handle cases N through 1
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
				const auto& my_pcs = popcounts[n_of_rems];

				// for switch (n), run cases n through 1, where the index is n-1 through 0
				constexpr size_t start = __LINE__ + 5;
			#define IDX(n) ((max_remainders - (n - start)) - 1)
			#define CASE(n) [[fallthrough]]; case(IDX(n) + 1): rem += size_t(my_pcs[IDX(n)]) * my_rems[IDX(n)];
				switch (n_of_rems)
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
					CASE(__LINE__); // case (1)
					static_assert(start + max_remainders == __LINE__);
					break;
				default:
					__assume(false);
				}
			#undef CASE
			#undef IDX
			}

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

	void print_div_test_analysis()
	{
	#if analyze_div_tests
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
			std::cout << "   base " << std::setfill(' ') << w(2) << size_t(dt.base) << "^n % " << w(3) << size_t(small_primes_lookup[dt.prime_idx]) << ":  ";
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



	void find_multibase_primes()
	{
		gmp_random::r.seed(mpir_ui(0xdeadbeef));

		size_t number = benchmark_mode ? bm_start : load_from_results();
		mpz_class mpz_number = 0ull; // it's a surprise tool that will help us later

		const std::vector<sieve_t> static_sieve = generate_static_sieve();
		std::vector<sieve_t> sieve = static_sieve;

		/* The number must start on an odd multiple of the sieve size. To round N to the nearest odd multiple of K:
		 * n -= k;
		 * n -= n % 2k;
		 * n += k; */
		number -= static_sieve.size();
		number -= number % (2 * static_sieve.size());
		number += static_sieve.size();

		set_up_sieve_offsets_cache(number);

		constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();
		constexpr size_t gcd_lookup = build_gcd_lookup();

		// Start the clock after setup
		const auto start = current_time_in_ms();

		// (condition optimizes out when not benchmarking)
		while (benchmark_mode ? number < bm_stop : true)
		{
			// Perform additional sieving on the static sieve
			sieve = static_sieve;
			partial_sieve(sieve);

			for (size_t i = 0; i < sieve.size(); ++i, number += 2)
			{
				// Bail if this number is already known to have a small prime factor
				if (!sieve[i]) continue;

				// Bail if n does not have a prime number of bits set.
				if ((tiny_primes_lookup & (1ull << pop_count(number))) == 0) continue;

				// Bail if gcd(abs(alternating sums), 1155) is not equal to one.
				const auto pca = pop_count(number & 0xAAAAAAAAAAAAAAAA);
				const auto pcb = pop_count(number & 0x5555555555555555);
				if ((gcd_lookup & (1ull << abs(pca - pcb))) == 0) continue;

				// Run cheap trial division tests across multiple bases
				if (has_small_divisor_cached(number)) continue;



				// Do full primality tests, starting with base 2
				if (!franken::mpir_is_likely_prime_BPSW(number)) continue;

				// convert uint64_t to char array of ['0', '1'...] for MPIR
				char bin_str[64 + 1];
				auto result = std::to_chars(&bin_str[0], &bin_str[64], number, 2);
				*result.ptr = '\0';

				mpz_number.set_str(bin_str, 3);
				if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 4);
				if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 5);
				if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 6);
				if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 7);
				if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 8);
				if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

				mpz_number.set_str(bin_str, 9);
				if (!mpir_is_prime(mpz_number)) { log_result(number, 8); continue; }

				mpz_number.set_str(bin_str, 10);
				if (!mpir_is_prime(mpz_number)) { log_result(number, 9); continue; }

				mpz_number.set_str(bin_str, 11);
				if (!mpir_is_prime(mpz_number)) { log_result(number, 10); continue; }

				mpz_number.set_str(bin_str, 12);
				if (!mpir_is_prime(mpz_number)) { log_result(number, 11); continue; }

				mpz_number.set_str(bin_str, 13);
				if (!mpir_is_prime(mpz_number)) { log_result(number, 12); continue; }

			} // end hot loop

		#if analyze_div_tests
			for (const auto& dt : div_tests)
				if (dt.hits >= 1'000'000)
					goto div_test_summaries;
		#endif
		} // end sieve loop

	#if analyze_div_tests
		div_test_summaries :
		print_div_test_analysis();
	#endif

		std::cout << "Finished. " << current_time_in_ms() - start << " ms elapsed\n";
	}



	void find_multibase_primes_permute()
	{
		gmp_random::r.seed(mpir_ui(0xdeadbeef));

		size_t number = benchmark_mode ? bm_start : load_from_results();
		size_t bits = number >> 1; // The last bit is always 1; everything else is permuted
		mpz_class mpz_number = 0ull; // it's a surprise tool that will help us later

		constexpr size_t gcd_lookup = build_gcd_lookup();

		// Start the clock after setup
		const auto start = current_time_in_ms();

		// Search for primes across all permutations of p bits, where p is:
		// { 13, 17, 19, 23, 29, 31, ...? }

		// (condition optimizes out when not benchmarking)
		while (benchmark_mode ? number < bm_stop : true)
		{
			// The last bit is always 1; everything else is permuted
			lex_permute(bits);
			number = (bits << 1) | 0b1;

			// Bail if gcd(abs(alternating sums), 1155) is not equal to one.
			const auto pca = pop_count(number & 0xAAAAAAAAAAAAAAAA);
			const auto pcb = pop_count(number & 0x5555555555555555);
			if ((gcd_lookup & (1ull << abs(pca - pcb))) == 0) continue;

			// Run cheap trial division tests in base 2
			// if (b2_has_small_divisor<32>(number)) continue;

			// Run cheap trial division tests across multiple bases
			if (has_small_divisor_cached(number)) continue;

			// Run cheap trial division tests in base 2
			// if (b2_has_small_divisor<100, 32 + 1>(number)) continue;



			// Do full primality tests, starting with base 2
			if (!franken::mpir_is_likely_prime_BPSW(number)) continue;

			// convert uint64_t to char array of ['0', '1'...] for MPIR
			char bin_str[64 + 1];
			auto result = std::to_chars(&bin_str[0], &bin_str[64], number, 2);
			*result.ptr = '\0';

			mpz_number.set_str(bin_str, 3);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

			mpz_number.set_str(bin_str, 4);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

			mpz_number.set_str(bin_str, 5);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

			mpz_number.set_str(bin_str, 6);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

			mpz_number.set_str(bin_str, 7);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

			mpz_number.set_str(bin_str, 8);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes - 1)) continue;

			mpz_number.set_str(bin_str, 9);
			if (!mpir_is_prime(mpz_number)) { log_result(number, 8); continue; }

			mpz_number.set_str(bin_str, 10);
			if (!mpir_is_prime(mpz_number)) { log_result(number, 9); continue; }

			mpz_number.set_str(bin_str, 11);
			if (!mpir_is_prime(mpz_number)) { log_result(number, 10); continue; }

			mpz_number.set_str(bin_str, 12);
			if (!mpir_is_prime(mpz_number)) { log_result(number, 11); continue; }

			mpz_number.set_str(bin_str, 13);
			if (!mpir_is_prime(mpz_number)) { log_result(number, 12); continue; }
		}

	#if analyze_div_tests
		print_div_test_analysis();
	#endif

		std::cout << "Finished. " << current_time_in_ms() - start << " ms elapsed\n";
	}

} // namespace mbp



int main()
{
	mbp::find_multibase_primes();
}
