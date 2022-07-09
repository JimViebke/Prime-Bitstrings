#pragma once

#include <vector>

#include "pk_prime.hpp"
#include "math.hpp"

/*
To determine if a bitstring has a small prime divisor in base b, we can do better than converting the bitstring to base b, then calculating bistring % smallprime == 0.

Given any bitstring in the original base 2, the value of each digit is either 0 or base^(digit position). For example, in base 3, the bits can only represent 3^0, 3^1, 3^2... respectively.

Therefore, we can easily calculate the remainder of [place value] % [smallprime]. This could result in a 64-entry lookup, however, these remainders follow a short, repetitive pattern. Therefore, instead of comparing 64 bits with 64 different remainders, we only need to grab each unique remainder. Given k unique remainders, we can use a mask k times to select every k-th bit in the original number, perform a popcount, and multiple that popcount by the remainder. Adding these k remainders together gives us a small number, not greater than perhaps 100-200.

Then, cheaply determine if this sum is divisible by p. This determines if the bitstring would be divisible by p in base b.
*/

namespace mbp::div_test
{
	// Dimensions are [primes * bases][remainders], for testing base^place_value % prime == 0
	constexpr const std::vector<std::vector<uint8_t>> generate_mod_remainders()
	{
		std::vector<std::vector<uint8_t>> remainders;
		remainders.reserve(mod_remainders_size);

		for (size_t i = 1; i < n_of_primes; ++i) // for each small prime
		{
			for (size_t base = 3; base <= up_to_base; ++base) // for each base 3..n
			{
				if (small_primes_lookup[i] > (uint8_t)-1)
					std::cout << "N mod " << small_primes_lookup[i] << " doesn't fit uint8_t" << std::endl;

				std::vector<uint8_t> rems;
				for (size_t j = 0; j < 64; ++j) // calculate base^j MOD prime
				{
					uint8_t rem = uint8_t(pk::powMod(base, j, small_primes_lookup[i]));
					if (rem == 1 && j > 0) break; // break when the pattern repeats
					rems.push_back(rem);
				}
				remainders.push_back(rems);
			}
		}

		return remainders;
	}

	constexpr std::array<size_t, mod_remainders_size> generate_mod_remainder_bitmasks()
	{
		std::array<size_t, mod_remainders_size> bitmasks = {};

		const auto remainders = generate_mod_remainders();

		for (size_t i = 0; i < mod_remainders_size; ++i)
		{
			const size_t k = remainders[i].size();
			size_t bitmask = 1;

			for (size_t j = 0; j < 64 && k < 64; j += k)
			{
				bitmask <<= k;
				bitmask |= 1;
			}

			bitmasks[i] = bitmask;
		}

		return bitmasks;
	}

	namespace detail
	{
		template<size_t prime, size_t base>
		struct bitmask_for
		{
		private:
			static constexpr size_t f()
			{
				size_t i = 0;
				for (; i < 64; ++i) // calculate base^i MOD prime
				{
					uint8_t rem = uint8_t(pk::powMod(base, i, prime));
					if (rem == 1 && i > 0) break; // break when the pattern repeats
				}

				size_t bitmask = 0;
				for (size_t j = 0; j < 64 && i < 64; j += i)
				{
					bitmask <<= i;
					bitmask |= 1;
				}

				return bitmask;
			}
		public:
			static constexpr size_t val = f();
		};

		template<size_t a, size_t b, size_t m>
		struct pow_mod
		{
			static constexpr size_t rem = pk::powMod(a, b, m);
		};

		template<size_t bitmask>
		struct period_of
		{
		private:
			static constexpr size_t f()
			{
				// shift bitmask until we find a 1
				static_assert(bitmask & 1);
				for (size_t i = 1; i < 64; ++i)
					if ((bitmask >> i) & 1)
						return i;
			}
		public:
			static constexpr size_t val = f();
		};

		constexpr size_t largest_remainder()
		{
			// find the largest remainder
			size_t largest_remainder = 0;
			const auto remainders = generate_mod_remainders();
			for (const auto& a : remainders)
				for (const auto& b : a)
					if (b > largest_remainder)
						largest_remainder = b;
			return largest_remainder;
		}

		// The sum of remainders won't be larger than (largest prime <56) * (the largest remainder)
		constexpr size_t prime_factor_lookup_size = 53 * largest_remainder();

		std::vector<size_t> build_prime_factor_lookup()
		{
			std::vector<size_t> lookup;
			lookup.reserve(prime_factor_lookup_size);

			// for every possible summation of remainders
			for (size_t i = 0; i < prime_factor_lookup_size; ++i)
			{
				// for the first 64 primes
				size_t entry = 0;
				for (size_t p = 0; p < 64; ++p)
				{
					// calculate if that summation of remainders is divisible by that prime
					entry |= (size_t((i % small_primes_lookup[p]) == 0) << p);
				}

				lookup.push_back(entry);
			}

			return lookup;
		}

		// Replaces "n % prime[k] == 0" with "lookup[n] & (1 << k)"
		const std::vector<size_t> prime_factor_lookup = build_prime_factor_lookup();

		__forceinline bool has_small_prime_factor(const size_t n, const size_t prime_index)
		{
			return (prime_factor_lookup[n] & (1ull << prime_index)) != 0;
		}
	}

	// Suppress warnings about bitmasks having upper bits moved
#pragma warning(push)
#pragma warning(disable: 26450)

	__forceinline bool divisible_by_5_in_base_3(const size_t number)
	{
		using namespace detail;

		constexpr size_t b3_mask = bitmask_for<5, 3>::val;
		static_assert(period_of<b3_mask>::val == 4); // 3^n % 5 has 4 values
		size_t rem = 0;
		rem += pop_count(number & (b3_mask << 0)) * pow_mod<3, 0, 5>::rem;
		rem += pop_count(number & (b3_mask << 1)) * pow_mod<3, 1, 5>::rem;
		rem += pop_count(number & (b3_mask << 2)) * pow_mod<3, 2, 5>::rem;
		rem += pop_count(number & (b3_mask << 3)) * pow_mod<3, 3, 5>::rem;
		return detail::has_small_prime_factor(rem, 2); // idx 2; 5 is the 3rd prime
	}

	__forceinline bool divisible_by_5(const size_t number)
	{
		if (divisible_by_5_in_base_3(number)) return true;

		return false;
	}

	__forceinline bool divisible_by_7_in_base_3(const size_t number)
	{
		using namespace detail;
		constexpr size_t b3_mask = bitmask_for<7, 3>::val;
		static_assert(period_of<b3_mask>::val == 6);
		size_t rem = 0;
		rem += pop_count(number & (b3_mask << 0)) * pow_mod<3, 0, 7>::rem;
		rem += pop_count(number & (b3_mask << 1)) * pow_mod<3, 1, 7>::rem;
		rem += pop_count(number & (b3_mask << 2)) * pow_mod<3, 2, 7>::rem;
		rem += pop_count(number & (b3_mask << 3)) * pow_mod<3, 3, 7>::rem;
		rem += pop_count(number & (b3_mask << 4)) * pow_mod<3, 4, 7>::rem;
		rem += pop_count(number & (b3_mask << 5)) * pow_mod<3, 5, 7>::rem;
		return detail::has_small_prime_factor(rem, 3); // idx 3; 7 is the 4th prime
	}

	__forceinline bool divisible_by_7_in_base_4(const size_t number)
	{
		using namespace detail;
		constexpr size_t b4_mask = bitmask_for<7, 4>::val;
		static_assert(period_of<b4_mask>::val == 3);
		size_t rem = 0;
		rem += pop_count(number & (b4_mask << 0)) * pow_mod<4, 0, 7>::rem;
		rem += pop_count(number & (b4_mask << 1)) * pow_mod<4, 1, 7>::rem;
		rem += pop_count(number & (b4_mask << 2)) * pow_mod<4, 2, 7>::rem;
		return detail::has_small_prime_factor(rem, 3); // idx 3; 7 is the 4th prime
	}

	__forceinline bool divisible_by_7_in_base_5(const size_t number)
	{
		using namespace detail;
		constexpr size_t b5_mask = bitmask_for<7, 5>::val;
		static_assert(period_of<b5_mask>::val == 6);
		size_t rem = 0;
		rem += pop_count(number & (b5_mask << 0)) * pow_mod<5, 0, 7>::rem;
		rem += pop_count(number & (b5_mask << 1)) * pow_mod<5, 1, 7>::rem;
		rem += pop_count(number & (b5_mask << 2)) * pow_mod<5, 2, 7>::rem;
		rem += pop_count(number & (b5_mask << 3)) * pow_mod<5, 3, 7>::rem;
		rem += pop_count(number & (b5_mask << 4)) * pow_mod<5, 4, 7>::rem;
		rem += pop_count(number & (b5_mask << 5)) * pow_mod<5, 5, 7>::rem;
		return detail::has_small_prime_factor(rem, 3); // idx 3; 7 is the 4th prime
	}

	__forceinline bool divisible_by_7(const size_t number)
	{
		if (divisible_by_7_in_base_3(number)) return true;
		if (divisible_by_7_in_base_4(number)) return true;
		if (divisible_by_7_in_base_5(number)) return true;
		return false;
	}

#pragma warning(pop)

}
