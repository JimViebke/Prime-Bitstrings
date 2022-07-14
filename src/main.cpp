
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <charconv>

#include "utility.hpp"
#include "math.hpp"
#include "multibase_div_tests.hpp"
#include "config.hpp"
#include "sandbox.hpp"
#pragma warning(push, 0)
#include "franken_mpir.hpp"
#pragma warning(pop)



size_t load_from_results()
{
	// Load the largest (not necessarily last) number from the "results" log.

	std::ifstream ifs(mbp::results_path);
	size_t number = 0;

	std::string str;
	while (std::getline(ifs, str))
	{
		if (str.empty()) continue;

		const auto start = str.find('(');
		const auto end = str.find(')');

		if (start == std::string::npos ||
			end == std::string::npos ||
			start > end) continue;

		str = str.substr(start + 1, (end - start) - 1);

		size_t v = 0;
		auto r = std::from_chars(str.c_str(), str.c_str() + str.size(), v);

		if (r.ec != std::errc())
		{
			std::cout << "Read bad value: " << v << ", error: " << size_t(r.ec) << std::endl;
			continue;
		}

		if (v > number)
			number = v;
	}

	std::cout << "Loaded starting point from " << mbp::results_path << ": " << number << std::endl;

	return number;
}

void log_time()
{
	time_t timestamp = time(0);
	tm now;
	localtime_s(&now, &timestamp);
	std::cout << std::setfill('0') << std::setw(2) << ((now.tm_hour % 12 == 0) ? 12 : now.tm_hour % 12) << ':' << std::setw(2) << now.tm_min << '\t';
}

// suppress "unreachable" warning while in benchmark mode
#pragma warning(push)
#pragma warning(disable: 4702)
void log_result(const mpz_class& n, size_t up_to_base)
{
	log_time();

	std::stringstream ss;
	ss << bin_to_base(n, 10) << " is a p" << up_to_base << " (" << n << ")\n";

	std::cout << ss.str();

	if constexpr (mbp::benchmark_mode) return;

	std::ofstream ofs(mbp::results_path, std::ofstream::app);
	ofs << ss.str();
}
#pragma warning(pop)

using sieve_t = uint8_t;
const std::vector<sieve_t> generate_static_sieve()
{
	std::vector<sieve_t> sieve(mbp::static_sieve_size, true);

	// for each prime, mark off all multiples
	for (const auto p : mbp::static_sieve_primes)
		for (size_t i = 0; i < sieve.size(); i += p)
			sieve[i] = false;

	return sieve;
}

std::vector<size_t> sieve_offsets_cache(small_primes_lookup.size());

void set_up_sieve_offsets_cache(const size_t start)
{
	// Start with the first prime not in the static sieve.
	for (size_t i = mbp::static_sieve_primes.size() + 1; i < small_primes_lookup.size(); ++i)
	{
		const size_t p = small_primes_lookup[i];

		// Find out how far it is to the next multiple of p.
		size_t n = p - (start % p);

		// Start is always odd. Therefore, if n is odd, it is pointing to the next even multiple of p.
		// -- increase by p
		// However, if n is even, it is pointing to the next odd multiple of p.
		// -- do nothing

		// if (n % 2 == 1) // branchful
		//	n += p;
		n += p & -(n % 2 == 1); // branchless

		// We now have the distance to the next odd multiple of p.
		// Divide by 2 to store the *index* of the next odd multiple of p.
		sieve_offsets_cache[i] = n / 2;
	}
}

void partial_sieve(std::vector<sieve_t>& sieve)
{
	static_assert(mbp::static_sieve_size > mbp::sieve_primes_cap);

	// Start with the first prime not in the static sieve
	for (size_t i = mbp::static_sieve_primes.size() + 1; i < small_primes_lookup.size(); ++i)
	{
		const size_t p = small_primes_lookup[i];

		// Get the index of the next odd multiple of p
		size_t j = sieve_offsets_cache[i];

		// Mark false each (implicitly odd) multiple of p
		for (; j < sieve.size(); j += p)
		{
			sieve[j] = false;
		}

		// Update the cache for the next sieving
		sieve_offsets_cache[i] = j - mbp::static_sieve_size;
	}
}

__forceinline bool has_small_divisor(const size_t number,
									 const std::vector<std::vector<uint8_t>>& remainders,
									 const std::array<size_t, mbp::div_test::mod_remainders_size>& bitmasks)
{
	using namespace mbp;
	using namespace div_test;

	if (is_divisible_by<5, in_base<3>>(number)) return true;

	if (is_divisible_by<7, in_base<3>>(number)) return true;
	if (is_divisible_by<7, in_base<4>>(number)) return true;
	if (is_divisible_by<7, in_base<5>>(number)) return true;

	// div_test::n_of_bases * [the number of primes to skip, starting from 3]
	for (size_t i = 0; i < remainders.size(); ++i)
	{
		// skip two expensive tests
		if (remainders[i].size() == 64) continue;

		// Do this here, because first shift is always 0 and first rem is always 1
		size_t rem = pop_count(number & (bitmasks[i]));

		for (size_t l = 1; l < remainders[i].size(); ++l)
		{
			rem += pop_count(number & (bitmasks[i] << l)) * remainders[i][l];
		}

		// see if the sum of remainders is evenly divisible by a given prime
		if (div_test::detail::has_small_prime_factor(rem, (i / div_test::n_of_bases) + 1 + div_test::n_of_primes_with_hardcoded_divtests))
		{
#if 0
			static std::vector<uint8_t> active_indexes(remainders.size(), false);
			if (!active_indexes[i]) // if this isn't known to be an active index
			{
				active_indexes[i] = true; // note it as active

				std::cout << "\n\n\nA number was divisible by " << small_primes_lookup[(i / div_test::n_of_bases) + 1]
					<< " in base " << (i % div_test::n_of_bases) + 3 << " (" << remainders[i].size() << " remainders)\n\n";
				//bool any = false;
				//for (size_t j = 0; j < active_indexes.size(); ++j)
				//{
				//	if (!active_indexes[j])
				//	{
				//		any = true;
				//		std::cout << "No numbers found divisible by " << small_primes_lookup[(j / div_test::n_of_bases) + 1]
				//			<< " in base " << (j % div_test::n_of_bases) + 3 << " (idx = " << j << ")\n";
				//	}
				//}
				//if (!any) std::cout << "All values in remainders lookup used\n";

				mbp::detail::print_active_mod_remainders(active_indexes);
			}
#endif
			return true;
		}
	}

	return false;
}

#pragma warning(push)
#pragma warning(disable: 26450)
constexpr std::array<size_t, 64> bitmask_lookup = []
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
} ();
#pragma warning(pop)

__forceinline bool has_small_divisor(const size_t number,
									 const std::array<mbp::div_test::div_test_t, mbp::div_test::div_tests_size>& div_tests)
{
	using namespace mbp;
	using namespace mbp::div_test;

	if (is_divisible_by<5, in_base<3>>(number)) return true;

	if (is_divisible_by<7, in_base<3>>(number)) return true;
	if (is_divisible_by<7, in_base<4>>(number)) return true;
	if (is_divisible_by<7, in_base<5>>(number)) return true;

	for (const div_test_t& div_test : div_tests)
	{
		// Perform the first popcount here, because first shift is always 0 and first rem is always 1
		size_t rem = pop_count(number & bitmask_lookup[div_test.n_of_remainders]);

		for (size_t i = 1; i < div_test.n_of_remainders; ++i)
		{
			rem += pop_count(number & (bitmask_lookup[div_test.n_of_remainders] << i)) * div_test.remainders[i];
		}

		if (div_test::detail::has_small_prime_factor(rem, div_test.prime_idx))
		{
			// std::cout << number << " is divisible by " << small_primes_lookup[div_test.prime_idx] << " in base " << size_t(div_test.base) << '\n';
			return true;
		}
	}

	return false;
}

void find_multibase_primes()
{
	using namespace mbp;

	gmp_random::r.seed(mpir_ui(0xdeadbeef));

	size_t number = mbp::benchmark_mode ? mbp::bm_start : load_from_results();
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
	constexpr size_t gcd_1155_lookup = build_gcd_1155_lookup();

	static const std::array<div_test::div_test_t, div_test::div_tests_size> div_tests = div_test::generate_div_tests();

	// Start the clock after setup
	const auto start = current_time_in_ms();

	// (condition optimizes out when not benchmarking)
	while (mbp::benchmark_mode ? number < bm_stop : true)
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
			if ((gcd_1155_lookup & (1ull << abs(pca - pcb))) == 0) continue;

			// Run cheap trial division tests across multiple bases
			if (has_small_divisor(number, div_tests)) continue;



			// Do full primality tests, starting with base 2
			if (!franken::mpir_is_likely_prime_BPSW(number)) continue;

			// convert uint64_t to char array of ['0', '1'...] for MPIR
			char bin_str[64 + 1];
			auto result = std::to_chars(&bin_str[0], &bin_str[64], number, 2);
			*result.ptr = '\0';

			mpz_number.set_str(bin_str, 3);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes)) continue;

			mpz_number.set_str(bin_str, 4);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes)) continue;

			mpz_number.set_str(bin_str, 5);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes)) continue;

			mpz_number.set_str(bin_str, 6);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes)) continue;

			mpz_number.set_str(bin_str, 7);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes)) continue;

			mpz_number.set_str(bin_str, 8);
			if (!franken::mpir_is_prime(mpz_number, div_test::n_of_primes)) continue;

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
	}

	std::cout << "Finished. " << current_time_in_ms() - start << " ms elapsed\n";
}


int main()
{
	find_multibase_primes();
}
