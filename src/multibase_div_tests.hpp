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
	// Dimensions are [primes * bases][remainders]
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
		// Replaces "n % prime[k] == 0" with "lookup[n] & (1 << k)"
		const std::vector<size_t> build_divides_evenly_lookup()
		{
			// Given a % b, a must be smaller than (the largest prime <64) * (the largest remainder)

			// find the largest remainder
			size_t largest_remainder = 0;

			{
				const auto remainders = generate_mod_remainders();
				for (const auto& a : remainders)
					for (const auto& b : a)
						if (b > largest_remainder)
							largest_remainder = b;
			} // deallocate remainders copy

			std::vector<size_t> lookup;
			lookup.reserve(largest_remainder * 61); // 61 is the smallest prime <= 64

			// for every possible summation of remainders
			for (size_t i = 0; i < largest_remainder * 61; ++i)
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

		const std::vector<size_t> divides_evenly_lookup = build_divides_evenly_lookup();
	}

	inline bool divides_evenly(const size_t n, const size_t prime_index)
	{
		return (detail::divides_evenly_lookup[n] & (1ull << prime_index)) != 0;
	}
}
