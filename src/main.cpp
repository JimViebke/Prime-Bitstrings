
#include <iostream>
#include <vector>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <sstream>
#include <charconv>
#include <bitset>

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

void log_result(const mpz_class& n, size_t up_to_base)
{
	log_time();

	std::stringstream ss;
	ss << bin_to_base(n, 10) << " is a p" << up_to_base << " (" << n << ")\n";

	std::cout << ss.str();

	std::ofstream ofs(mbp::results_path, std::ofstream::app);
	ofs << ss.str();
}

void partial_sieve(const size_t start, std::vector<uint8_t>& sieve)
{
	constexpr static size_t idx = mbp::static_sieve_primes.size() + 1;

	for (size_t i = idx; i < small_primes_lookup.size(); ++i)
	{
		const size_t p = small_primes_lookup[i];

		// Find out how far past we are from the previous multiple of p.
		// ie 3 == 10 % 7
		// This will always be <= p.
		const size_t offset_from_last_pn = (start % p) / 2;
		// Divide by 2 because the sieve only represents odd numbers

		// Now mark false each multiple of p, starting from [0 + p - distance from previous p],
		// where 0 is the n we started with.
		for (size_t j = p - offset_from_last_pn; j < sieve.size(); j += p)
		{
			sieve[j] = false;
		}
	}
}

const std::vector<uint8_t> generate_static_sieve()
{
	std::vector<uint8_t> sieve(mbp::static_sieve_size, true);

	// for each prime, mark off all multiples
	for (const auto p : mbp::static_sieve_primes)
		for (auto i = p; i < sieve.size(); i += p)
			sieve[i] = false;

	return sieve;
}

inline bool has_small_divisor(const size_t number,
							  const std::vector<std::vector<uint8_t>>& remainders,
							  const std::array<size_t, mbp::div_test::mod_remainders_size>& bitmasks)
{
	using namespace mbp;

	for (size_t i = 0; i < remainders.size(); ++i)
	{
		// skip six expensive tests, four of which are always false
		if (remainders[i].size() == 64) continue;

		// Do this here, because first shift is always 0 and first rem is always 1
		size_t rem = pop_count(number & (bitmasks[i]));

		for (size_t l = 1; l < remainders[i].size(); ++l)
		{
			rem += pop_count(number & (bitmasks[i] << l)) * remainders[i][l];
		}

		// see if the sum of remainders is evenly divisible by a given prime
		if (div_test::has_small_prime_factor(rem, (i / div_test::n_of_bases) + 1)) return true;
	}

	return false;
}

void find_multibase_primes()
{
	using namespace mbp;

	gmp_random::r.seed(mpir_ui(0xdeadbeef));

	size_t number = mbp::benchmark_mode ? mbp::bm_start : load_from_results();
	mpz_class mpz_number = 0ull; // it's a surprise tool that will help us later

	const std::vector<uint8_t> static_sieve = generate_static_sieve();
	std::vector<uint8_t> sieve;

	/* The number must start on an odd multiple of the sieve size. To round N to the nearest odd multiple of K:
	 * n -= k;
	 * n -= n % 2k;
	 * n += k; */
	number -= static_sieve.size();
	number -= number % (2 * static_sieve.size());
	number += static_sieve.size();

	constexpr size_t tiny_primes_lookup = build_tiny_primes_lookup();
	constexpr size_t gcd_1155_lookup = build_gcd_1155_lookup();

	// Dimensions are [primes * bases][remainders]
	const std::vector<std::vector<uint8_t>> remainders = div_test::generate_mod_remainders();
	// Dimensions are [primes * bases]
	constexpr std::array<size_t, div_test::mod_remainders_size> bitmasks = div_test::generate_mod_remainder_bitmasks();

	// Don't start the clock until here
	const auto start = current_time_in_ms();

	// Condition optimizes out when not benchmarking
	while (mbp::benchmark_mode ? number < bm_stop : true)
	{
		// perform additional sieving on the static sieve
		sieve = static_sieve;
		partial_sieve(number, sieve);

		for (size_t i = 0; i < mbp::static_sieve_size; ++i, number += 2)
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
			if (has_small_divisor(number, remainders, bitmasks)) continue;

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
