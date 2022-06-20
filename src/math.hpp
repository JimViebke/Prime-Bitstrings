#pragma once

#include <stdint.h>
#include <nmmintrin.h>

#include "mpirxx.h"
#include "pk_prime.hpp"

namespace bpsw_1_native
{
	int pow(int pow_a, unsigned int pow_b, int pow_c) {
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
	bool MiillerTest(int MT_dt, int MT_num) {
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
	bool prime(int P_N, int P_K) {
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
	return mpz_likely_prime_p(p.get_mpz_t(), gmp_random::r.get_randstate_t(), div);
}

bool mpir_is_probable_prime(const mpz_class& p, const int prob, const size_t div)
{
	return mpz_probable_prime_p(p.get_mpz_t(), gmp_random::r.get_randstate_t(), prob, div);
}

bool brute_force_is_prime(const size_t n)
{
	const size_t sqrt_n = size_t(sqrt(n));

	for (size_t i = 2; i <= sqrt_n; ++i)
	{
		if (n / i * i == n) return false;
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
	return _mm_popcnt_u64(n);

	//n -= ((n >> 1) & 0x5555555555555555ull);
	//n = (n & 0x3333333333333333ull) + (n >> 2 & 0x3333333333333333ull);
	//return ((n + (n >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

// Generate a 64 bit lookup, where the prime-numbered bits are set high
// 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 | 47 | 53 | 59 | 61
constexpr uint64_t build_tiny_primes_lookup()
{
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

constexpr size_t small_primes_cap = 1621; // 1000

std::vector<size_t> generate_small_primes()
{
	std::vector<size_t> primes;

	for (size_t i = 2; i <= small_primes_cap; ++i)
		if (mpir_is_prime(i))
			primes.push_back(i);

	return primes;
}

static const std::vector<size_t> small_primes_lookup = generate_small_primes();

constexpr uint16_t gcd(uint16_t a, uint16_t b)
{
	if (b == 0)
		return a;
	return gcd(b, a % b);
}

constexpr size_t build_gcd_1155_lookup()
{
	size_t lookup = 0;

	// set bit i high if and only if the GCD of (i, 1155) is 1
	for (uint16_t i = 0; i < 32; ++i)
	{
		lookup |= (gcd(i, 1155) == 1ull) << i;
	}

	return lookup;
}



// Test code only below this point.



void compare_implementations()
{
	size_t num = 282607273285049; // p11 - too large for 32-bit BPSW
	// size_t num = 113;
	std::cout << "Naive + native implementation:\n";
	std::cout << num << " is " << (mpir_is_prime(num) ? "prime" : "not prime") << std::endl;

	std::cout << "BPSW + native implementation:\n";
	std::cout << num << " is " << (bpsw_1_native::prime((int)num, 50) ? "prime" : "not prime") << std::endl;
}

void test_bpsw_1_native()
{
	size_t primes_found = 0;

	for (size_t num = 2; ; ++num)
	{
		if (bpsw_1_native::prime((int)num, 50))
		{
			std::cout << num << ' ';
			++primes_found;
		}
	}
}

void test_binary_to_decimal()
{
	for (size_t i = 0; i < 50; ++i)
	{
		std::cout << i << '\t' << binary_to_base((size_t)i, 10) << '\n';
	}
}

// This shouldn't work - converting large primes between bases should be larger than a size_t
void find_p2_8()
{
	auto start = current_time_in_ms();

	for (size_t binary = 3; ; binary += 2)
	{
		if (!mpir_is_prime(binary)) continue;

		const size_t b3 = binary_to_base(binary, 3);
		if (!mpir_is_prime(b3)) continue;

		const size_t b4 = binary_to_base(binary, 4);
		if (!mpir_is_prime(b4)) continue;

		const size_t b5 = binary_to_base(binary, 5);
		if (!mpir_is_prime(b5)) continue;

		const size_t b6 = binary_to_base(binary, 6);
		if (!mpir_is_prime(b6)) continue;

		const size_t b7 = binary_to_base(binary, 7);
		if (!mpir_is_prime(b7)) continue;

		const size_t b8 = binary_to_base(binary, 8);
		if (!mpir_is_prime(b8)) continue;

		const size_t b10 = binary_to_base(binary, 10);

		std::cout << b10 << " is prime in bases 2-8 (" <<
			binary << ", " <<
			b3 << ", " <<
			b4 << ", " <<
			b5 << ", " <<
			b6 << ", " <<
			b7 << ", " <<
			b8 << ")\n";
		break;
	}

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
	for (size_t i = -1, n = 0; n < 1000; i--)
	{
		if (pk::is_prime(i))
		{
			n++;
			std::cout << i << '\n';
		}
	}
}
