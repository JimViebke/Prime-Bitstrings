// Prime Bitstrings.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <sstream>
#include <charconv>
#include <bitset>

//#include "pk_prime.hpp"

#include "utility.hpp"
#include "math.hpp"
#include "multibase_div_tests.hpp"
#include "config.hpp"
// #include "sandbox.hpp"

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

	// Only log large primes to file
	if (up_to_base >= mbp::smallest_base_to_log)
	{
		std::ofstream ofs(mbp::results_path, std::ofstream::app);
		ofs << ss.str();
	}
}

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

void partial_sieve(const size_t start, std::vector<uint8_t>& sieve)
{
	for (const size_t p : small_primes_lookup)
	{
		// only mark off the small primes not already here
		if (p <= mbp::static_sieve_primes.back()) continue;

		// Find out how far past we are from the previous multiple of p.
		// ie 3 == 10 % 7
		// This will always be <= p.
		size_t offset_from_last_pn = (start % p) / 2;
		// Divide by 2 because the sieve only represents odd numbers

		// Now mark false each multiple of p, starting from [0 + p - distance from previous p],
		// where 0 is the n we started with.
		for (size_t i = p - offset_from_last_pn; i < sieve.size(); i += p)
		{
			sieve[i] = false;
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
							  const std::vector<std::vector<std::vector<uint8_t>>>& remainders,
							  const std::vector<std::vector<size_t>>& bitmasks)
{
	using namespace mbp;

	for (size_t j = 0; j < div_test::n_of_primes; j += div_test::primes_per_round)
	{
		// for each base 3..8
		for (size_t base = 3; base <= div_test::up_to_base; ++base)
		{
			// for each small prime
			for (size_t k = j; k < j + div_test::primes_per_round; ++k)
			{
				if (k == 0) continue;

				// mask against bitmask[base][k] to collect remainders in each set of positions
				size_t rem = 0;
				for (size_t l = 0; l < remainders[base][k].size(); ++l)
				{
					rem += pop_count(number & (bitmasks[base][k] << l)) * remainders[base][k][l];
				}

				// see if the sum of remainders is evenly divisible by a given prime
				if (divides_evenly(rem, k))
				{
					std::stringstream ss;
					ss << number << " would be divisible by " << small_primes_lookup[k] << " in base " << base << '\n';
					std::cout << ss.str();
					return true;
				}
			}
		}
	}

	return false;
}

inline bool has_small_divisor_compact_edition(const size_t number,
											  const std::vector<std::vector<uint8_t>>& remainders,
											  const std::vector<std::vector<size_t>>& bitmasks)
{
	using namespace mbp;

	for (size_t j = 0; j < div_test::n_of_primes; j += div_test::primes_per_round)
	{
		// for each base 3..8
		for (size_t base = 3; base <= div_test::up_to_base; ++base)
		{
			// for each small prime
			for (size_t k = j; k < j + div_test::primes_per_round; ++k)
			{
				// mask against bitmask[base][k] to collect remainders in each set of positions
				size_t rem = 0;
				for (size_t l = 0; l < remainders[base * div_test::n_of_primes + k].size(); ++l)
				{
					rem += pop_count(number & (bitmasks[base][k] << l)) * remainders[base * div_test::n_of_primes + k][l];
				}

				// see if the sum of remainders is evenly divisible by a given prime
				if (divides_evenly(rem, k)) return true;
			}
		}
	}

	return false;
}

inline bool has_small_divisor_compacter_er_edition(const size_t number,
												   const std::vector<uint8_t>& remainders,
												   const std::vector<std::vector<size_t>>& bitmasks)
{
	using namespace mbp;

	for (size_t j = 0; j < div_test::n_of_primes; j += div_test::primes_per_round)
	{
		// for each base 3..8
		for (size_t base = 3; base <= div_test::up_to_base; ++base)
		{
			// for each small prime
			for (size_t k = j; k < j + div_test::primes_per_round; ++k)
			{
				// mask against bitmask[base][k] to collect remainders in each set of positions
				size_t rem = 0;
				const size_t next_index = base * div_test::n_of_primes * 64 + k * 64;
				for (size_t l = 0; remainders[next_index + l] != uint8_t(-1) && l < 64; ++l)
				{
					rem += pop_count(number & (bitmasks[base][k] << l)) * remainders[next_index + l];
				}

				// see if the sum of remainders is evenly divisible by a given prime
				if (divides_evenly(rem, k)) return true;
			}
		}
	}

	return false;
}

inline bool has_small_divisor_compacter_er_er_edition(const size_t number,
													  const std::vector<uint8_t>& remainders,
													  const std::vector<size_t>& bitmasks)
{
	using namespace mbp;

	for (size_t j = 0; j < div_test::n_of_primes; j += div_test::primes_per_round)
	{
		// for each base 3..8
		for (size_t base = 3; base <= div_test::up_to_base; ++base)
		{
			// for each small prime
			for (size_t k = j; k < j + div_test::primes_per_round; ++k)
			{
				// mask against bitmask[base][k] to collect remainders in each set of positions
				size_t rem = 0;
				const size_t next_index = base * div_test::n_of_primes * 64 + k * 64;
				for (size_t l = 0; remainders[next_index + l] != uint8_t(-1) && l < 64; ++l)
				{
					rem += pop_count(number & (bitmasks[base * div_test::n_of_primes + k] << l)) * remainders[next_index + l];
				}

				// see if the sum of remainders is evenly divisible by a given prime
				if (divides_evenly(rem, k)) return true;
			}
		}
	}

	return false;
}

inline bool has_small_divisor_another_edition(const size_t number,
											  const std::vector<std::vector<uint8_t>>& remainders,
											  const std::vector<size_t>& bitmasks)
{
	using namespace mbp;

	for (size_t j = 0; j < div_test::n_of_primes; j += div_test::primes_per_round)
	{
		// for each base 3..8
		for (size_t base = 3; base <= div_test::up_to_base; ++base)
		{
			// for each small prime
			for (size_t k = j; k < j + div_test::primes_per_round; ++k)
			{
				// mask against bitmask[base][k] to collect remainders in each set of positions
				size_t rem = 0;
				for (size_t l = 0; l < remainders[base * div_test::n_of_primes + k].size(); ++l)
				{
					rem += pop_count(number & (bitmasks[base * div_test::n_of_primes + k] << l)) * remainders[base * div_test::n_of_primes + k][l];
				}

				// see if the sum of remainders is evenly divisible by a given prime
				if (divides_evenly(rem, k)) return true;
			}
		}
	}

	return false;
}

inline bool has_small_divisor_v6(const size_t number,
								 const std::vector<std::vector<std::vector<uint8_t>>>& remainders,
								 const std::vector<std::vector<size_t>>& bitmasks)
{
	using namespace mbp;

	// for each small prime
	for (size_t i = 0; i < div_test::n_of_primes; ++i)
	{
		// for each base 3..8
		for (size_t base = 3; base <= div_test::up_to_base; ++base)
		{
			// mask against bitmask[i][base]
			size_t rem = 0;
			for (size_t l = 0; l < remainders[i][base].size(); ++l)
			{
				rem += pop_count(number & (bitmasks[i][base] << l)) * remainders[i][base][l];
			}

			// see if the sum of remainders is evenly divisible by a given prime
			if (divides_evenly(rem, i)) return true;
		}
	}

	return false;
}

inline bool has_small_divisor_v7(const size_t number,
								 const std::vector<std::vector<uint8_t>>& remainders,
								 const std::vector<size_t>& bitmasks)
{
	using namespace mbp;

	constexpr size_t n_of_bases = (div_test::up_to_base + 1) - 3;

	for (size_t i = 0; i < remainders.size(); ++i)
	{
		size_t rem = 0;
		for (size_t l = 0; l < remainders[i].size(); ++l)
		{
			rem += pop_count(number & (bitmasks[i] << l)) * remainders[i][l];
		}

		// see if the sum of remainders is evenly divisible by a given prime
		if (divides_evenly(rem, i / n_of_bases)) return true;
	}

	return false;
}

inline bool has_small_divisor_v8(const size_t number,
								 const std::vector<std::vector<uint8_t>>& remainders,
								 const std::vector<size_t>& bitmasks)
{
	using namespace mbp;

	size_t i = 0; // :(

	for (size_t j = 0; j < div_test::n_of_primes; j += div_test::primes_per_round)
	{
		// for each base 3..8
		for (size_t base = 3; base <= div_test::up_to_base; ++base)
		{
			// for each small prime
			for (size_t k = j; k < j + div_test::primes_per_round; ++k, ++i)// :(
			{
				if (k == 0) continue;

				// mask to collect remainders in each set of positions
				size_t rem = 0;
				for (size_t l = 0; l < remainders[i].size(); ++l)
				{
					rem += pop_count(number & (bitmasks[i] << l)) * remainders[i][l];
				}

				// see if the sum of remainders is evenly divisible by a given prime
				if (divides_evenly(rem, k)) return true;
			}
		}
	}

	return false;
}

/*__declspec(noinline)*/ inline bool has_small_divisor_v9(const size_t number,
														  const std::vector<std::vector<uint8_t>>& remainders,
														  const std::vector<size_t>& bitmasks)
{
	using namespace mbp;

	constexpr size_t n_of_bases = (div_test::up_to_base + 1) - 3;

	for (size_t i = 0; i < remainders.size(); ++i)
	{
		size_t rem = 0;
		for (size_t l = 0; l < remainders[i].size(); ++l)
		{
			rem += pop_count(number & (bitmasks[i] << l)) * remainders[i][l];
		}

		// see if the sum of remainders is evenly divisible by a given prime
		if (divides_evenly(rem, (i / n_of_bases) + 1)) return true;
	}

	return false;
}


void find_multibase_primes()
{
	using namespace mbp;

	gmp_random::r.seed(rand());

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

	// Dimensions are [base 3..n][primes][remainders]
	const std::vector<std::vector<std::vector<uint8_t>>> remainders = generate_remainders_for_bases();
	// Dimensions are [base 3..n][bitmasks for p]
	const std::vector<std::vector<size_t>> bitmasks = generate_mod_remainder_bitmasks(remainders);

	// hmmm
	// const std::vector<std::vector<uint8_t>> remainders_compact = generate_remainders_for_bases_compacter_version();

	// hmmmmmm
	// const std::vector<uint8_t> remainders_compacter_er = generate_remainders_for_bases_compacter_er_version();

	// hmmmmmmmmm                             -- I'm building this using the old (3D) version of the remainders vector,
	//                                        -- but for now that should be fine.
	// const std::vector<size_t> bitmasks_compact = generate_mod_remainder_bitmasks_compact_edition(remainders);

	const std::vector<std::vector<std::vector<uint8_t>>> remainders_v6 = generate_remainders_v6();
	const std::vector<std::vector<size_t>> bitmasks_v6 = generate_bitmasks_v6(remainders_v6);

	const std::vector<std::vector<uint8_t>> remainders_v7 = generate_remainders_v7();
	const std::vector<size_t> bitmasks_v7 = generate_bitmasks_v7(remainders_v7);

	// Note the use of the OG remainders, to generate v8 with the same order
	const std::vector<std::vector<uint8_t>> remainders_v8 = generate_remainders_v8(remainders);
	const std::vector<size_t> bitmasks_v8 = generate_bitmasks_v7(remainders_v8); // reuse v7 gen

	const std::vector<std::vector<uint8_t>> remainders_v9 = generate_remainders_v9();
	const std::vector<size_t> bitmasks_v9 = generate_bitmasks_v7(remainders_v9); // reuse v7 gen

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

			// Bail if gcd(abs(# of even bits - # of odd bits), 1155) is not equal to one.
			const int pca = (int)pop_count(number & 0xAAAAAAAAAAAAAAAA);
			const int pcb = (int)pop_count(number & 0x5555555555555555);
			if ((gcd_1155_lookup & (1ull << abs(pca - pcb))) == 0) continue;

			// Run cheap trial division tests across multiple bases
			// if (has_small_divisor(number, remainders, bitmasks)) continue;
			// if (has_small_divisor_compact_edition(number, remainders_compact, bitmasks)) continue;
			// if (has_small_divisor_compacter_er_edition(number, remainders_compacter_er, bitmasks)) continue;
			// if (has_small_divisor_compacter_er_er_edition(number, remainders_compacter_er, bitmasks_compact)) continue;
			// if (has_small_divisor_another_edition(number, remainders_compact, bitmasks_compact)) continue;
			// if (has_small_divisor_v6(number, remainders_v6, bitmasks_v6)) continue;
			// if (has_small_divisor_v7(number, remainders_v7, bitmasks_v7)) continue;
			// if (has_small_divisor_v8(number, remainders_v8, bitmasks_v8)) continue;
			if (has_small_divisor_v9(number, remainders_v9, bitmasks_v9)) continue;

			// bool v1 = has_small_divisor(number, remainders, bitmasks);
			// bool v2 = has_small_divisor_compact_edition(number, remainders_compact, bitmasks);
			// bool v3 = has_small_divisor_compacter_er_edition(number, remainders_compacter_er, bitmasks);
			// bool v4 = has_small_divisor_compacter_er_er_edition(number, remainders_compacter_er, bitmasks_compact);
			// bool v5 = has_small_divisor_another_edition(number, remainders_compact, bitmasks_compact);
			// bool v6 = has_small_divisor_v6(number, remainders_v6, bitmasks_v6);
			// bool v7 = has_small_divisor_v7(number, remainders_v7, bitmasks_v7);
			// bool v8 = has_small_divisor_v8(number, remainders_v8, bitmasks_v8);
			// bool v9 = has_small_divisor_v9(number, remainders_v9, bitmasks_v9);

			//if (v1 != v9)
			//	std::cout << "fffffff\n";

			//if (v1) continue;

			// Do full primality tests, starting with base 2
			if (!franken::mpir_is_likely_prime_BPSW(number)) continue;

			// convert uint64_t to char array of ['0', '1'...] for MPIR
			char bin_str[64 + 1];
			auto result = std::to_chars(&bin_str[0], &bin_str[64], number, 2);
			*result.ptr = '\0';

			mpz_number.set_str(bin_str, 3);
			if (!mpir_is_prime(mpz_number)) continue;

			mpz_number.set_str(bin_str, 4);
			if (!mpir_is_prime(mpz_number)) continue;

			mpz_number.set_str(bin_str, 5);
			if (!mpir_is_prime(mpz_number)) continue;

			mpz_number.set_str(bin_str, 6);
			if (!mpir_is_prime(mpz_number)) continue;

			mpz_number.set_str(bin_str, 7);
			if (!mpir_is_prime(mpz_number)) continue;

			mpz_number.set_str(bin_str, 8);
			if (!mpir_is_prime(mpz_number)) continue;

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
	//mbp::detail::cheaper_divtests();

	//mbp::detail::diminishing_returns();

	find_multibase_primes();
}
