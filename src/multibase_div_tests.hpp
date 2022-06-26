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

namespace detail
{
	// For numbers in base b, dimensions are [primes][residues]
	std::vector<std::vector<uint8_t>> generate_remainders_for_base(const size_t base, const size_t n_of_primes)
	{
		std::vector<std::vector<uint8_t>> remainders;

		for (size_t i = 0; i < n_of_primes; ++i) // for each small prime
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

		return remainders;
	}
}

// dimensions are [base 3..n][primes][residues]
std::vector<std::vector<std::vector<uint8_t>>> generate_remainders_for_bases(
	const size_t stop_base,
	const size_t n_of_primes)
{
	std::vector<std::vector<std::vector<uint8_t>>> remainders;
	remainders.push_back(decltype(remainders)::value_type()); // Dummy values for base 0
	remainders.push_back(decltype(remainders)::value_type()); // b1
	remainders.push_back(decltype(remainders)::value_type()); // b2

	for (size_t base = 3; base <= stop_base; ++base)
	{
		remainders.push_back(detail::generate_remainders_for_base(base, n_of_primes));
	}

	return remainders;
}

// Dimensions are [base 3..n][bitmasks for p]
// This will correctly handle the dummy elements 0..2
std::vector<std::vector<size_t>> generate_mod_remainder_bitmasks(const std::vector<std::vector<std::vector<uint8_t>>>& remainders)
{
	std::vector<std::vector<size_t>> bitmasks;

	for (auto& r : remainders)
	{
		std::vector<size_t> b;

		for (size_t i = 0; i < r.size(); ++i)
		{
			const size_t k = r[i].size();
			size_t bitmask = 0;

			for (size_t j = 0; j < 64; j += k)
			{
				bitmask <<= k;
				bitmask |= 1;
			}
			b.push_back(bitmask);
		}

		bitmasks.push_back(b);
	}

	return bitmasks;
}
