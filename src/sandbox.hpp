#pragma once

#include <iostream>
#include <vector>
#include <bitset>
#include <map>

#pragma warning(push, 0)
#include "mpirxx.h"
#pragma warning(pop)

#include "math.hpp"
#include "utility.hpp"
#include "multibase_div_tests.hpp"
#pragma warning(push, 0)
#include "franken_mpir.hpp"
#pragma warning(pop)

namespace mbp
{
	namespace detail
	{
		void test_binary_to_decimal()
		{
			for (size_t i = 0; i < 50; ++i)
			{
				std::cout << i << '\t' << binary_to_base((size_t)i, 10) << '\n';
			}
		}

		void find_p_n()
		{
			const auto start = current_time_in_ms();

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
					if (!mpir_is_prime(mpz_number)) break;

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

		void mpir_testing()
		{
			gmp_randclass r{ gmp_randinit_mt };
			r.seed(rand());

			const mpz_class million{ 1000 * 1000 };
			const mpz_class billion{ million * 1000 };
			const mpz_class trillion{ billion * 1000 };
			const mpz_class quadrillion{ trillion * 1000 };
			const mpz_class quintillion{ quadrillion * 1000 };
			const mpz_class sextillion{ quintillion * 1000 };
			const mpz_class septillion{ sextillion * 1000 };
			const mpz_class octillion{ septillion * 1000 };
			const mpz_class nonillion{ octillion * 1000 };
			const mpz_class decillion{ nonillion * 1000 };

			const size_t target = 10'000;

			std::cout << mpz_class{ "1000000000000000035000061" } << std::endl;

			const mpz_class num{ "1000000000000000035000061" };
			std::cout << num << std::endl;

			size_t count = 0;
			auto start = current_time_in_ms();
			for (mpz_class i = num; count < target; ++i)
				if (mpz_likely_prime_p(i.get_mpz_t(), r.get_randstate_t(), 0))
					++count;

			std::cout << "A linear search found the next " << count << " primes in " << current_time_in_ms() - start << "ms\n";

			start = current_time_in_ms();
			count = 0;
			for (mpz_class i = num; count < target; mpz_next_prime_candidate(i.get_mpz_t(), i.get_mpz_t(), r.get_randstate_t()))
				if (mpz_likely_prime_p(i.get_mpz_t(), r.get_randstate_t(), 0))
					++count;

			std::cout << "A next-candidate search found the next " << count << " primes in " << current_time_in_ms() - start << "ms\n";
		}

		void mpir_testing_2()
		{
			gmp_randclass r{ gmp_randinit_mt };
			r.seed(rand());

			const mpz_class num{ "413" };
			std::cout << num;

			auto start = current_time_in_ms();

			if (mpz_likely_prime_p(num.get_mpz_t(), r.get_randstate_t(), 0))
			{
				std::cout << " is a prime number\n";
			}
			else
			{
				std::cout << " is not a prime number, but\n";

				mpz_class next_prime = num;

				do
				{
					mpz_next_prime_candidate(next_prime.get_mpz_t(), next_prime.get_mpz_t(), r.get_randstate_t());
				} while (!mpz_likely_prime_p(next_prime.get_mpz_t(), r.get_randstate_t(), 0));

				std::cout << next_prime << " is the next prime up." << std::endl;
			}

			std::cout << "Found in " << current_time_in_ms() - start << " ms\n";
		}

		void pk_testing()
		{
			// print n largest primes, counting down from 2^64
			for (size_t i = size_t(-1), n = 0; n < 1000; i--)
			{
				if (pk::is_prime(i))
				{
					n++;
					std::cout << i << '\n';
				}
			}
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
			std::cout << "(no suitable sieve size found)\n";
		}

		void print_mod_remainders()
		{
			using namespace mbp::div_test;

			const std::vector<std::vector<uint8_t>> remainders = generate_mod_remainders();

			auto it = remainders.begin();

			// std::multimap<size_t, std::string> ordered_checks;

			for (size_t i = 1; i < n_of_primes; ++i)
			{
				for (size_t b = 3; b <= up_to_base; ++b)
				{
					std::stringstream ss;
					ss << b << "^n % " << small_primes_lookup[i] << " = (" << it->size() << " place values)";
					std::cout << ss.str() << '\n';

					// ordered_checks.insert({ small_primes_lookup[i] * it->size(), ss.str() });

					++it;
				}
			}

			//for (const auto& kv : ordered_checks)
			//	std::cout << kv.first << '\t' << kv.second << '\n';
		}

		void print_active_mod_remainders(const std::vector<uint8_t>& active)
		{
			// Log the list of remainders for each [prime, base] divtests.
			// Only log remainders if the specific test is known to have eliminated at least one candidate.

			using namespace mbp::div_test;

			static std::vector<std::vector<uint8_t>> remainders = generate_mod_remainders();

			size_t r_idx = 0;

			// std::multimap<size_t, std::string> ordered_checks;

			for (size_t i = 1; i < n_of_primes; ++i)
			{
				for (size_t b = 3; b <= up_to_base; ++b)
				{
					if (active[r_idx])
					{
						std::cout << b << "^n % " << small_primes_lookup[i] << " has " << remainders[r_idx].size() << " values: ";

						if (remainders[r_idx].size() == 64)
						{
							std::cout << "(skipping)";
						}
						else
						{
							for (const auto k : remainders[r_idx])
								std::cout << size_t(k) << ' ';
						}

						std::cout << '\n';
					}

					// ordered_checks.insert({ small_primes_lookup[i] * it->size(), ss.str() });

					++r_idx;
				}
			}

			//for (const auto& kv : ordered_checks)
			//	std::cout << kv.first << '\t' << kv.second << '\n';
		}
	}
}
