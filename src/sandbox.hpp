#pragma once

#include <iostream>
#include <vector>

#pragma warning(push, 0)
#include "mpirxx.h"
#pragma warning(pop)

#include "math.hpp"
#include "utility.hpp"
#pragma warning(push, 0)
#include "franken_mpir.hpp"
#pragma warning(pop)

namespace mbp
{
	namespace detail
	{
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
	}
}
