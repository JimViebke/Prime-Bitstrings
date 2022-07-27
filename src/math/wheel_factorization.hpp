#pragma once

#include "config.hpp"
#include "math.hpp"

namespace mbp
{
	/*
	An odd number N is definitely composite if N >= wheel_size and N % wheel_size is
	- a composite number, or
	- one of the primes used to construct the wheel
	*/

	constexpr size_t wheel_size = std::accumulate(wheel_primes.begin(), wheel_primes.end(),
												  size_t(1), std::multiplies<>());

	namespace detail
	{
		consteval std::array<uint8_t, wheel_size> build_wheel()
		{
			std::array<uint8_t, wheel_size> wheel = { uint8_t(false) };
			// Mark off each multiple of each prime
			for (const auto prime : wheel_primes)
				for (size_t i = 0; i < wheel_size; i += prime)
					wheel[i] = true;
			return wheel;
		}

		// True if N % wheel_size is composite, for N >= wheel_size
		constexpr std::array<uint8_t, wheel_size> wheel = build_wheel();
	}

	// True if n has a prime factor <= wheel_primes.back()
	__forceinline bool is_composite(const size_t n)
	{
		return detail::wheel[n % wheel_size];
	}
}
