#pragma once

#pragma warning(push, 0)
#include "mpirxx.h"
#pragma warning(pop)

#include <stdint.h>
#include <nmmintrin.h>

#include "config.hpp"
#pragma warning(push, 0)
#include "franken_fermat.hpp"
#pragma warning(pop)

const std::vector<size_t> small_primes_lookup = build_small_primes_lookup();

namespace bpsw_1_native
{
	int pow(int pow_a, unsigned int pow_b, int pow_c)
	{
		int result = 1;
		pow_a = pow_a % pow_c;
		while (pow_b > 0) {
			if (pow_b & 1)
				result = (result * pow_a) % pow_c;
			pow_b = pow_b >> 1;
			pow_a = (pow_a * pow_a) % pow_c;
		}
		return result;
	}
	bool MiillerTest(int MT_dt, int MT_num)
	{
		int MT_a = 2 + rand() % (MT_num - 4);
		int MT_x = pow(MT_a, MT_dt, MT_num);
		if (MT_x == 1 || MT_x == MT_num - 1)
			return true;
		while (MT_dt != MT_num - 1) {
			MT_x = (MT_x * MT_x) % MT_num;
			MT_dt *= 2;
			if (MT_x == 1)
				return false;
			if (MT_x == MT_num - 1)
				return true;
		}
		return false;
	}
	bool prime(int P_N, int P_K)
	{
		if (P_N <= 1 || P_N == 4)
			return false;
		if (P_N <= 3)
			return true;
		int P_D = P_N - 1;
		while (P_D % 2 == 0)
			P_D /= 2;
		for (int i = 0; i < P_K; i++)
			if (MiillerTest(P_D, P_N) == false)
				return false;
		return true;
	}
}

namespace gmp_random
{
	gmp_randclass r{ gmp_randinit_mt };
}

bool mpir_is_prime(const mpz_class& p, const size_t div = 0)
{
	return bool(mpz_likely_prime_p(p.get_mpz_t(), gmp_random::r.get_randstate_t(), div));
}

bool mpir_is_probable_prime(const mpz_class& p, const int prob, const size_t div)
{
	return bool(mpz_probable_prime_p(p.get_mpz_t(), gmp_random::r.get_randstate_t(), prob, div));
}

constexpr bool brute_force_is_prime(const size_t n)
{
	if (n % 2 == 0) return false;

	const size_t sqrt_n = size_t(sqrt(n));
	for (size_t i = 3; i <= sqrt_n; i += 2)
	{
		if (n % i == 0) return false;
	}

	return true;
}

size_t binary_to_base(size_t binary, const size_t base)
{
	constexpr size_t highest_bit = size_t(1) << 63;

	size_t num = 0;

	// for each bit
	for (auto i = 0; i < 64; binary <<= 1, ++i)
	{
		num *= base;

		// if the (i-1)th bit is set, +1
		if ((binary & highest_bit) == highest_bit)
		{
			num += 1;
		}
	}

	return num;
}

mpz_class bin_to_base(const mpz_class& binary, const int base)
{
	return mpz_class{ binary.get_str(2), base };
}

inline uint64_t pop_count(uint64_t n)
{
	return uint64_t(_mm_popcnt_u64(n));

	//n -= ((n >> 1) & 0x5555555555555555ull);
	//n = (n & 0x3333333333333333ull) + (n >> 2 & 0x3333333333333333ull);
	//return ((n + (n >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

constexpr uint64_t build_tiny_primes_lookup()
{
	// Generate a 64 bit lookup, where the prime-numbered bits are set high
	// 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 | 47 | 53 | 59 | 61

	uint64_t lookup = 0;

	// We don't care about tiny pNs, therefore leave out the smallest of the primes.
	//lookup |= (1ull << 2);
	//lookup |= (1ull << 3);
	//lookup |= (1ull << 5);
	//lookup |= (1ull << 7);
	//lookup |= (1ull << 11);
	lookup |= (1ull << 13);
	lookup |= (1ull << 17);
	lookup |= (1ull << 19);
	lookup |= (1ull << 23);
	lookup |= (1ull << 29);
	lookup |= (1ull << 31);
	lookup |= (1ull << 37);
	lookup |= (1ull << 41);
	lookup |= (1ull << 43);
	lookup |= (1ull << 47);
	lookup |= (1ull << 53);
	lookup |= (1ull << 59);
	lookup |= (1ull << 61);

	return lookup;
}

const std::vector<size_t> build_small_primes_lookup()
{
	std::vector<size_t> primes;

	for (size_t i = 2; i <= mbp::sieve_primes_cap; ++i)
		if (mpir_is_prime(i))
			primes.push_back(i);

	return primes;
}

constexpr size_t gcd(size_t a, size_t b)
{
	if (b == 0)
		return a;
	return gcd(b, a % b);
}

constexpr size_t build_gcd_1155_lookup()
{
	size_t lookup = 0;

	// set bit i high if and only if the GCD of (i, 1155) is 1
	for (size_t i = 0; i < 32; ++i)
	{
		lookup |= size_t(gcd(i, 1155u) == 1ull) << i;
	}

	return lookup;
}
