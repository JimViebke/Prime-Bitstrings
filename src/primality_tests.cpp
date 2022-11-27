
#include <charconv>

#include "find_multibase_primes.hpp"

#include "config.hpp"
#include "io/io.hpp"
#include "math/franken_mpir.hpp"

void mbp::find_multibase_primes::full_primality_tests(const uint64_t* candidates,
													  const uint64_t* const candidates_end)
{
	for (; candidates < candidates_end; ++candidates)
	{
		const size_t candidate = *candidates;

		if (!franken::mpir_is_likely_prime_BPSW(candidate)) continue;
		count_passes(++b2);

		// convert uint64_t to char array of ['0', '1'...] for MPIR
		char bin_str[64 + 1];
		auto result = std::to_chars(&bin_str[0], &bin_str[64], candidate, 2);
		*result.ptr = '\0';

		mpz_number.set_str(bin_str, 3);
		if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;
		count_passes(++b3);

		mpz_number.set_str(bin_str, 4);
		if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;
		count_passes(++b4);

		mpz_number.set_str(bin_str, 5);
		if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;
		count_passes(++b5);

		mpz_number.set_str(bin_str, 6);
		if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;

		mpz_number.set_str(bin_str, 7);
		if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1)) continue;

		mpz_number.set_str(bin_str, 8);
		if (!franken::mpir_is_prime(mpz_number, gmp_rand, div_test::n_of_primes - 1))
		{
			if constexpr (!benchmark_mode) log_result(candidate, 7);
			continue;
		}

		mpz_number.set_str(bin_str, 9);
		if (!mpir_is_prime(mpz_number, gmp_rand))
		{
			if constexpr (!benchmark_mode) log_result(candidate, 8);
			continue;
		}

		mpz_number.set_str(bin_str, 10);
		if (!mpir_is_prime(mpz_number, gmp_rand)) { log_result(candidate, 9); continue; }

		mpz_number.set_str(bin_str, 11);
		if (!mpir_is_prime(mpz_number, gmp_rand)) { log_result(candidate, 10); continue; }

		mpz_number.set_str(bin_str, 12);
		if (!mpir_is_prime(mpz_number, gmp_rand)) { log_result(candidate, 11); continue; }

		mpz_number.set_str(bin_str, 13);
		if (!mpir_is_prime(mpz_number, gmp_rand)) { log_result(candidate, 12); continue; }

		log_result(candidate, 13);
	}
}
