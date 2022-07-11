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
		// bitmask for all base^n mod prime
		template<size_t base, size_t prime>
		struct bitmask_for
		{
			static constexpr size_t val = []
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
			}();
		};

		// calculate a^b mod m
		template<size_t a, size_t b, size_t m>
		struct pow_mod
		{
			static constexpr size_t rem = pk::powMod(a, b, m);
		};

		template<size_t bitmask>
		struct period_of
		{
			static constexpr size_t val = []
			{
				for (size_t i = 1; i < 64; ++i)
					if ((bitmask >> i) & 1)
						return i;
			}();
		};

		template<size_t prime>
		struct get_prime_index
		{
			static constexpr size_t idx = []
			{
				for (size_t i = 0; i < 64; ++i)
					if (small_primes_lookup[i] == prime)
						return i;
			}();
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

	// 1/3 - hardcoded version:

	__forceinline bool divisible_by_5_in_base_3(const size_t number)
	{
		using namespace detail;

		constexpr size_t mask = bitmask_for<3, 5>::val;
		static_assert(period_of<mask>::val == 4); // 3^n % 5 has 4 values
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<3, 0, 5>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<3, 1, 5>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<3, 2, 5>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<3, 3, 5>::rem;
		return detail::has_small_prime_factor(rem, get_prime_index<5>::idx);
	}

	__forceinline bool divisible_by_7_in_base_3(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<3, 7>::val;
		static_assert(period_of<mask>::val == 6);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<3, 0, 7>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<3, 1, 7>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<3, 2, 7>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<3, 3, 7>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<3, 4, 7>::rem;
		rem += pop_count(number & (mask << 5)) * pow_mod<3, 5, 7>::rem;
		return detail::has_small_prime_factor(rem, get_prime_index<7>::idx);
	}

	__forceinline bool divisible_by_7_in_base_4(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<4, 7>::val;
		static_assert(period_of<mask>::val == 3);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<4, 0, 7>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<4, 1, 7>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<4, 2, 7>::rem;
		return detail::has_small_prime_factor(rem, get_prime_index<7>::idx);
	}

	__forceinline bool divisible_by_7_in_base_5(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<5, 7>::val;
		static_assert(period_of<mask>::val == 6);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<5, 0, 7>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<5, 1, 7>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<5, 2, 7>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<5, 3, 7>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<5, 4, 7>::rem;
		rem += pop_count(number & (mask << 5)) * pow_mod<5, 5, 7>::rem;
		return detail::has_small_prime_factor(rem, get_prime_index<7>::idx);
	}

	__forceinline bool divisible_by_11_in_base_3(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<3, 11>::val;
		static_assert(period_of<mask>::val == 5);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<3, 0, 11>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<3, 1, 11>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<3, 2, 11>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<3, 3, 11>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<3, 4, 11>::rem;
		return detail::has_small_prime_factor(rem, 4); // idx 4; 11 is the 5th prime
	}

	__forceinline bool divisible_by_11_in_base_4(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<4, 11>::val;
		static_assert(period_of<mask>::val == 5);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<4, 0, 11>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<4, 1, 11>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<4, 2, 11>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<4, 3, 11>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<4, 4, 11>::rem;
		return detail::has_small_prime_factor(rem, 4);
	}

	__forceinline bool divisible_by_11_in_base_5(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<5, 11>::val;
		static_assert(period_of<mask>::val == 5);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<5, 0, 11>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<5, 1, 11>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<5, 2, 11>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<5, 3, 11>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<5, 4, 11>::rem;
		return detail::has_small_prime_factor(rem, 4);
	}

	__forceinline bool divisible_by_11_in_base_6(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<6, 11>::val;
		static_assert(period_of<mask>::val == 10);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<6, 0, 11>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<6, 1, 11>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<6, 2, 11>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<6, 3, 11>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<6, 4, 11>::rem;
		rem += pop_count(number & (mask << 5)) * pow_mod<6, 5, 11>::rem;
		rem += pop_count(number & (mask << 6)) * pow_mod<6, 6, 11>::rem;
		rem += pop_count(number & (mask << 7)) * pow_mod<6, 7, 11>::rem;
		rem += pop_count(number & (mask << 8)) * pow_mod<6, 8, 11>::rem;
		rem += pop_count(number & (mask << 9)) * pow_mod<6, 9, 11>::rem;
		return detail::has_small_prime_factor(rem, 4);
	}

	__forceinline bool divisible_by_11_in_base_7(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<7, 11>::val;
		static_assert(period_of<mask>::val == 10);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<7, 0, 11>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<7, 1, 11>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<7, 2, 11>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<7, 3, 11>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<7, 4, 11>::rem;
		rem += pop_count(number & (mask << 5)) * pow_mod<7, 5, 11>::rem;
		rem += pop_count(number & (mask << 6)) * pow_mod<7, 6, 11>::rem;
		rem += pop_count(number & (mask << 7)) * pow_mod<7, 7, 11>::rem;
		rem += pop_count(number & (mask << 8)) * pow_mod<7, 8, 11>::rem;
		rem += pop_count(number & (mask << 9)) * pow_mod<7, 9, 11>::rem;
		return detail::has_small_prime_factor(rem, 4);
	}

	__forceinline bool divisible_by_11_in_base_8(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<8, 11>::val;
		static_assert(period_of<mask>::val == 10);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<8, 0, 11>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<8, 1, 11>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<8, 2, 11>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<8, 3, 11>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<8, 4, 11>::rem;
		rem += pop_count(number & (mask << 5)) * pow_mod<8, 5, 11>::rem;
		rem += pop_count(number & (mask << 6)) * pow_mod<8, 6, 11>::rem;
		rem += pop_count(number & (mask << 7)) * pow_mod<8, 7, 11>::rem;
		rem += pop_count(number & (mask << 8)) * pow_mod<8, 8, 11>::rem;
		rem += pop_count(number & (mask << 9)) * pow_mod<8, 9, 11>::rem;
		return detail::has_small_prime_factor(rem, 4);
	}

	// 2/3 - templated version:

	template<size_t base>
	struct in_base
	{
		static constexpr size_t val = base;
		static constexpr bool is(size_t b) { return b == base; }
	};

	template<size_t divisor, typename base_t, const size_t base = base_t::val> requires (divisor == 5 && base == 3)
		__forceinline bool is_divisible_by(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<base, divisor>::val;
		static_assert(period_of<mask>::val == 4);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<base, 0, divisor>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<base, 1, divisor>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<base, 2, divisor>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<base, 3, divisor>::rem;

		return detail::has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

	template<size_t divisor, typename base_t, const size_t base = base_t::val> requires (divisor == 7 && base == 3)
		__forceinline bool is_divisible_by(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<base, divisor>::val;
		static_assert(period_of<mask>::val == 6);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<base, 0, divisor>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<base, 1, divisor>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<base, 2, divisor>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<base, 3, divisor>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<base, 4, divisor>::rem;
		rem += pop_count(number & (mask << 5)) * pow_mod<base, 5, divisor>::rem;

		return detail::has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

	template<size_t divisor, typename base_t, const size_t base = base_t::val> requires (divisor == 7 && base == 4)
		__forceinline bool is_divisible_by(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<base, divisor>::val;
		static_assert(period_of<mask>::val == 3);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<base, 0, divisor>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<base, 1, divisor>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<base, 2, divisor>::rem;

		return detail::has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

	template<size_t divisor, typename base_t, const size_t base = base_t::val> requires (divisor == 7 && base == 5)
		__forceinline bool is_divisible_by(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<base, divisor>::val;
		static_assert(period_of<mask>::val == 6);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<base, 0, divisor>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<base, 1, divisor>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<base, 2, divisor>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<base, 3, divisor>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<base, 4, divisor>::rem;
		rem += pop_count(number & (mask << 5)) * pow_mod<base, 5, divisor>::rem;

		return detail::has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

	template<size_t divisor, typename base_t, const size_t base = base_t::val> requires (divisor == 11 && base == 8)
		__forceinline bool is_divisible_by(const size_t number)
	{
		using namespace detail;
		constexpr size_t mask = bitmask_for<base, divisor>::val;
		static_assert(period_of<mask>::val == 10);
		size_t rem = 0;
		rem += pop_count(number & (mask << 0)) * pow_mod<base, 0, divisor>::rem;
		rem += pop_count(number & (mask << 1)) * pow_mod<base, 1, divisor>::rem;
		rem += pop_count(number & (mask << 2)) * pow_mod<base, 2, divisor>::rem;
		rem += pop_count(number & (mask << 3)) * pow_mod<base, 3, divisor>::rem;
		rem += pop_count(number & (mask << 4)) * pow_mod<base, 4, divisor>::rem;
		rem += pop_count(number & (mask << 5)) * pow_mod<base, 5, divisor>::rem;
		rem += pop_count(number & (mask << 6)) * pow_mod<base, 6, divisor>::rem;
		rem += pop_count(number & (mask << 7)) * pow_mod<base, 7, divisor>::rem;
		rem += pop_count(number & (mask << 8)) * pow_mod<base, 8, divisor>::rem;
		rem += pop_count(number & (mask << 9)) * pow_mod<base, 9, divisor>::rem;

		return detail::has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

	// 3/3 - recursive templated version:

	template<size_t divisor, size_t base, size_t mask, size_t place_value>
	__forceinline void recursive_is_divisible_by(size_t& rem, const size_t number)
	{
		rem += pop_count(number & (mask << place_value)) * detail::pow_mod<base, place_value, divisor>::rem;
		if constexpr (place_value > 0)
			recursive_is_divisible_by<divisor, base, mask, place_value - 1>(rem, number);
		return;
	}

	template<size_t divisor, typename base_t, const size_t base = base_t::val>
	__forceinline bool recursive_is_divisible_by(const size_t number)
	{
		static_assert(base >= 3 && base <= up_to_base);
		constexpr size_t bitmask = detail::bitmask_for<base, divisor>::val;
		size_t rem = 0;
		recursive_is_divisible_by<divisor, base, bitmask, detail::period_of<bitmask>::val - 1>(rem, number);
		return detail::has_small_prime_factor(rem, detail::get_prime_index<divisor>::idx);
	}

#pragma warning(pop)

}
