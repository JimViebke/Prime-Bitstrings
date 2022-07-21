
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cassert>
#include <bitset>

#include "utility.hpp"
#include "math.hpp"
#include "multibase_div_tests.hpp"
#include "config.hpp"
#include "sandbox.hpp"
#pragma warning(push, 0)
#include "franken_mpir.hpp"
#pragma warning(pop)


namespace mbp
{
	static std::filesystem::path results_path;

	void set_up_results_path()
	{
		namespace fs = std::filesystem;
		fs::path dir = fs::current_path();

		// Search upward for Prime Bitstrings folder
		while (dir.filename() != "Prime Bitstrings")
		{
			if (dir.has_parent_path())
			{
				dir = dir.parent_path();
			}
			else
			{
				std::cout << "Could not find " << results_filename << " along " << fs::current_path() << '\n';
				exit(EXIT_FAILURE);
			}
		}

		results_path = dir;
		results_path.append(results_filename);
	}

	size_t load_from_results()
	{
		// Load the largest (not necessarily last) number from the "results" log.
		set_up_results_path();

		std::ifstream ifs(results_path);
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

			if (r.ec != std::errc{})
			{
				std::cout << "Read bad entry: " << v << ", error: " << size_t(r.ec) << '\n';
				continue;
			}

			if (v > number)
				number = v;
		}

		if (div_test::max_pn_bitwidth <= std::bit_width(number))
		{
			std::cout << "error: max_pn_bitwidth should be larger than std::bit_width(number)\n";
			exit(EXIT_FAILURE);
		}

		std::cout << "Loaded " << number << " from " << results_path.generic_string() << '\n';

		return number;
	}

	void log_time()
	{
		time_t timestamp = time(0);
		tm now;
		localtime_s(&now, &timestamp);
		std::cout << std::setfill(' ') << std::setw(2) << ((now.tm_hour % 12 == 0) ? 12 : now.tm_hour % 12) << ':' << std::setfill('0') << std::setw(2) << now.tm_min << '\t';
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

		if constexpr (benchmark_mode) return;

		std::ofstream ofs(results_path, std::ofstream::app);
		ofs << ss.str();
	}
#pragma warning(pop)



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

	std::vector<sieve_prime_t> sieve_offsets_cache(small_primes_lookup.size());

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
		static_assert(static_sieve_size > sieve_primes_cap);

		// Start with the first prime not in the static sieve
		for (size_t i = static_sieve_primes.size() + 1; i < small_primes_lookup.size(); ++i)
		{
			const sieve_prime_t p = small_primes_lookup[i];

			// Get the index of the next odd multiple of p
			sieve_prime_t j = sieve_offsets_cache[i];

			// Mark false each (implicitly odd) multiple of p
			for (; j < sieve.size(); j += p)
			{
				sieve[j] = false;
			}

			// Update the cache for the next sieving
			sieve_offsets_cache[i] = j - static_sieve_size;
		}
	}



#pragma warning(push)
#pragma warning(disable: 26450)
	constexpr std::array<size_t, 64> generate_bitmask_lookup()
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
#pragma warning(pop)

	namespace div_test
	{
#if analyze_div_tests
		using div_tests_t = std::array<div_test_t, div_tests_size>;
#else
		using div_tests_t = const std::array<div_test_t, div_tests_size>;
#endif
	}

	static div_test::div_tests_t div_tests = div_test::generate_div_tests();

	_declspec(noinline) /*__forceinline*/ bool has_small_divisor(const size_t number)
	{
		using namespace div_test;

		if (is_divisible_by<5, in_base<3>>(number)) return true;

		if (is_divisible_by<7, in_base<3>>(number)) return true;
		if (is_divisible_by<7, in_base<4>>(number)) return true;
		if (is_divisible_by<7, in_base<5>>(number)) return true;

#if analyze_div_tests
		bool found_div = false;

		for (auto& div_test : div_tests)
#else
		for (const auto& div_test : div_tests)
#endif
		{
			// Perform the first popcount here, because first shift is always 0 and first rem is always 1
			size_t rem = pop_count(number & bitmask_lookup[div_test.bitmask_idx]);

			for (size_t i = 1; i < div_test.n_of_remainders; ++i)
			{
				rem += pop_count(number & (bitmask_lookup[div_test.bitmask_idx] << i)) * div_test.remainders[i];
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

	__declspec(noinline) /*__forceinline*/ bool has_small_divisor_vectorized(const size_t number)
	{
		using namespace div_test;

		if (is_divisible_by<5, in_base<3>>(number)) return true;

		if (is_divisible_by<7, in_base<3>>(number)) return true;
		if (is_divisible_by<7, in_base<4>>(number)) return true;
		if (is_divisible_by<7, in_base<5>>(number)) return true;

		constexpr size_t avx2 = 256 / 8;
		constexpr size_t steps = avx2 / sizeof(uint64_t);

		alignas(avx2) const uint64_t num[steps] = { number, number, number, number };

#if analyze_div_tests
		bool found_div = false;

		for (auto& div_test : div_tests)
#else
		for (const auto& div_test : div_tests)
#endif
		{
			// Show the compiler we're on a multiple of 4
			const uint64_t iters = uint64_t(div_test.n_of_remainders >> 2) << 2;
			assert(iters % 4 == 0);
			const uint64_t mask = bitmask_lookup[div_test.bitmask_idx];

			alignas(avx2) uint32_t sums[steps * 2] = { 0 };

			for (size_t i = 0; i < iters; i += steps)
			{
				// copy each remainder into two adjacent 32-bit words
				alignas(avx2) uint64_t rems[steps];
				uint32_t* rems_p = (uint32_t*)rems;
				for (size_t k = 0; k < steps * 2; ++k)
					rems_p[k] = div_test.remainders[i + (k / 2)];

				// load the next four bitmasks
				alignas(avx2) uint64_t data[steps];
				uint32_t* data_p = (uint32_t*)data;
				for (size_t k = 0; k < steps; ++k)
					data[k] = mask << (i + k);

#pragma omp simd
				for (size_t k = 0; k < steps; ++k)
					data[k] &= num[k];

				for (size_t k = 0; k < steps; ++k)
					data[k] = (uint64_t)std::popcount(data[k]);

#pragma omp simd
				for (size_t k = 0; k < steps * 2; ++k)
				{
					// AVX512 has a popcount instruction, but not AVX2
					// (might be faster to do this using four x64 popcounts?)
					// data_p[d] = data_p[d] - ((data_p[d] >> 1) & 0x55555555);
					// data_p[d] = (data_p[d] & 0x33333333) + ((data_p[d] >> 2) & 0x33333333);
					// data_p[d] = ((data_p[d] + (data_p[d] >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;

					// multiply remainder counts by remainders
					data_p[k] *= rems_p[k];

					// save results
					sums[k] += data_p[k];
				}
			}

			size_t sum = std::accumulate(sums, sums + (steps * 2), size_t(0));

			if (has_small_prime_factor(sum, div_test.prime_idx))
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
		constexpr size_t gcd_1155_lookup = build_gcd_1155_lookup();

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
				if ((gcd_1155_lookup & (1ull << abs(pca - pcb))) == 0) continue;

				// Run cheap trial division tests across multiple bases
				bool a = has_small_divisor(number);
				bool b = has_small_divisor_vectorized(number);
				if (a != b)
					std::cout << "uhoh ";

				if (a) continue;

				// if (has_small_divisor(number)) continue;
				// if (has_small_divisor_vectorized(number)) continue;



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
				std::cout << "      -       ";
			}
			else
			{
				std::cout << w(7) << dt.hits << " hits  ";
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

		std::cout << "Finished. " << current_time_in_ms() - start << " ms elapsed\n";
	}

} // namespace mbp



int main()
{
	mbp::find_multibase_primes();
}
