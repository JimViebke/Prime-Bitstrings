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
	// For numbers in base b, dimensions are [primes][remainders]
	const std::vector<std::vector<uint8_t>> generate_remainders_for_base(const size_t base, const size_t n_of_primes)
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

// dimensions are [base 3..n][primes][remainders] (v1)
const std::vector<std::vector<std::vector<uint8_t>>> generate_remainders_for_bases()
{
	std::vector<std::vector<std::vector<uint8_t>>> remainders;
	remainders.push_back(decltype(remainders)::value_type()); // Dummy values for base 0
	remainders.push_back(decltype(remainders)::value_type()); // b1
	remainders.push_back(decltype(remainders)::value_type()); // b2

	for (size_t base = 3; base <= mbp::div_test::up_to_base; ++base)
	{
		remainders.push_back(detail::generate_remainders_for_base(base, mbp::div_test::n_of_primes));
	}

	return remainders;
}

// dimensions are [base 3..n * primes][remainders] (v2, v5)
// instead of [base][prime][remainder]
//            [base*N_of_primes + prime][remainder]
const std::vector<std::vector<uint8_t>> generate_remainders_for_bases_compacter_version()
{
	std::vector<std::vector<uint8_t>> remainders;

	for (size_t i = 0; i < 3 * mbp::div_test::n_of_primes; ++i)
	{
		remainders.push_back(decltype(remainders)::value_type()); // dummy values for bases 0, 1, and 2
	}

	for (size_t base = 3; base <= mbp::div_test::up_to_base; ++base)
	{
		for (size_t i = 0; i < mbp::div_test::n_of_primes; ++i) // for each small prime
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

// dimensions are [base 3..n * primes][remainders] (v3, v4)
// instead of [base][prime][remainder]
// instead of [base*N_of_primes + prime][remainder]
//            [base*N_of_primes*64 + prime*64 + remainders]
const std::vector<uint8_t> generate_remainders_for_bases_compacter_er_version()
{
	std::vector<uint8_t> remainders;

	for (size_t i = 0; i < 3 * mbp::div_test::n_of_primes * 64; ++i)
	{
		remainders.push_back(decltype(remainders)::value_type()); // dummy values for bases 0, 1, and 2
	}

	for (size_t base = 3; base <= mbp::div_test::up_to_base; ++base)
	{
		for (size_t i = 0; i < mbp::div_test::n_of_primes; ++i) // for each small prime
		{
			if (small_primes_lookup[i] > (uint8_t)-1)
				std::cout << "N mod " << small_primes_lookup[i] << " doesn't fit uint8_t" << std::endl;

			size_t j = 0;
			for (; j < 64; ++j) // calculate base^j MOD prime
			{
				uint8_t rem = uint8_t(pk::powMod(base, j, small_primes_lookup[i]));
				if (rem == 1 && j > 0) break; // break when the pattern repeats
				remainders.push_back(rem);
			}

			// :(
			for (; j < 64; ++j)
				remainders.push_back(uint8_t(-1));
		}
	}

	return remainders;
}

// Dimensions are [base 3..n][bitmasks for p], with blank elements for bases 0..2 (v1, v2, 3)
const std::vector<std::vector<size_t>> generate_mod_remainder_bitmasks(const std::vector<std::vector<std::vector<uint8_t>>>& remainders)
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

// Dimensions are [base 3..n][bitmasks for p], with blank elements for bases 0..2 (v4, v5)
const std::vector<size_t> generate_mod_remainder_bitmasks_compact_edition(const std::vector<std::vector<std::vector<uint8_t>>>& remainders)
{
	std::vector<size_t> bitmasks;

	// 3 * N_of_primes to pad the lookup
	for (size_t i = 0; i < 3 * mbp::div_test::n_of_primes; ++i)
	{
		bitmasks.push_back(0);
	}

	for (auto& r : remainders)
	{
		for (size_t i = 0; i < r.size(); ++i)
		{
			const size_t k = r[i].size();
			size_t bitmask = 0;

			for (size_t j = 0; j < 64; j += k)
			{
				bitmask <<= k;
				bitmask |= 1;
			}
			bitmasks.push_back(bitmask);
		}
	}

	return bitmasks;

}

// Dimensions are [primes][base 3..n][remainders]
const std::vector<std::vector<std::vector<uint8_t>>> generate_remainders_v6()
{
	std::vector<std::vector<std::vector<uint8_t>>> remainders;

	for (size_t i = 0; i < mbp::div_test::n_of_primes; ++i)
	{
		decltype(remainders)::value_type r2;
		r2.push_back(decltype(r2)::value_type()); // blank values for bases 0-2
		r2.push_back(decltype(r2)::value_type());
		r2.push_back(decltype(r2)::value_type());

		for (size_t base = 3; base <= mbp::div_test::up_to_base; ++base)
		{
			decltype(r2)::value_type r3;

			for (size_t j = 0; j < 64; ++j)
			{
				uint8_t rem = uint8_t(pk::powMod(base, j, small_primes_lookup[i]));
				if (rem == 1 && j > 0) break; // break when the pattern repeats
				r3.push_back(rem);
			}
			r2.push_back(r3);
		}
		remainders.push_back(r2);
	}

	return remainders;
}

// Dimensions are [prime][base 3..n], with blank elements for bases 0..2
const std::vector<std::vector<size_t>> generate_bitmasks_v6(const std::vector<std::vector<std::vector<uint8_t>>>& remainders)
{
	std::vector<std::vector<size_t>> bitmasks;

	// for each prime
	for (auto& r : remainders)
	{
		std::vector<size_t> b;
		b.push_back(0); // empty values for bases 0..2
		b.push_back(0);
		b.push_back(0);

		// for each base
		for (size_t i = 3; i < r.size(); ++i)
		{
			// get the period of remainders
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

// Dimensions are [primes * bases][remainders]
const std::vector<std::vector<uint8_t>> generate_remainders_v7()
{
	std::vector<std::vector<uint8_t>> remainders;
	// if up_to_base==8, then (3..n values) == (6 values) == (n + 1 - 3 values)
	remainders.reserve(mbp::div_test::n_of_primes *
					   ((mbp::div_test::up_to_base + 1) - 3));

	for (size_t i = 0; i < mbp::div_test::n_of_primes; ++i) // for each small prime
	{
		for (size_t base = 3; base <= mbp::div_test::up_to_base; ++base) // for each base 3..n
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

// Dimensions are [prime * base 3..n], with no blank elements for padding
const std::vector<size_t> generate_bitmasks_v7(const std::vector<std::vector<uint8_t>>& remainders)
{
	std::vector<size_t> bitmasks;
	bitmasks.reserve(remainders.size());

	for (auto& rems : remainders)
	{
		const size_t k = rems.size();
		size_t bitmask = 0;

		for (size_t j = 0; j < 64; j += k)
		{
			bitmask <<= k;
			bitmask |= 1;
		}

		bitmasks.push_back(bitmask);
	}

	return bitmasks;
}

// Dimensions are [primes and bases in v1 order][remainders]
const std::vector<std::vector<uint8_t>> generate_remainders_v8(const decltype(generate_remainders_for_bases())& remainders_v1)
{
	using namespace mbp;

	std::vector<std::vector<uint8_t>> remainders;
	remainders.reserve(remainders.size() * remainders[0].size());

	// loops are ripped from the OG divtest
	for (size_t j = 0; j < div_test::n_of_primes; j += div_test::primes_per_round)
	{
		// for each base 3..8
		for (size_t base = 3; base <= div_test::up_to_base; ++base)
		{
			// for each small prime
			for (size_t k = j; k < j + div_test::primes_per_round; ++k)
			{
				remainders.push_back(remainders_v1[base][k]);
			}
		}
	}

	// sanity check
	if (remainders.size() != remainders.size() * remainders[0].size())
	{
		std::cout << "Something went not okay.\n";
	}

	return remainders;
}

// Dimensions are [primes * bases][remainders]
const std::vector<std::vector<uint8_t>> generate_remainders_v9()
{
	std::vector<std::vector<uint8_t>> remainders;
	// if up_to_base==8, then (3..n values) == (6 values) == (n + 1 - 3 values)
	remainders.reserve(mbp::div_test::n_of_primes - 1 *
					   ((mbp::div_test::up_to_base + 1) - 3));

	for (size_t i = 1; i < mbp::div_test::n_of_primes; ++i) // for each small prime
	{
		for (size_t base = 3; base <= mbp::div_test::up_to_base; ++base) // for each base 3..n
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


namespace detail
{
	// Replaces "n % prime[k] == 0" with "lookup[n] & (1 << k)"
	const std::vector<size_t> build_divides_evenly_lookup()
	{
		// Given a % b, a must be smaller than (the largest prime <64) * (the largest remainder)

		// find the largest remainder
		size_t largest_remainder = 0;
		const auto remainders = generate_remainders_for_bases();
		for (const auto& a : remainders)
			for (const auto& b : a)
				for (const auto& c : b)
					if (c > largest_remainder)
						largest_remainder = c;

		std::vector<size_t> lookup;
		lookup.reserve(largest_remainder * 61); // 61 is the smallest prime <= 64

		// for every possible summation of remainders
		for (size_t i = 0; i < largest_remainder * 61; ++i)
		{
			size_t entry = 0;
			// for the first 64 primes
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
