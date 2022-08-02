#pragma once

#include "config.hpp"
#include "utility.hpp"

namespace mbp
{
	using sieve_prime_t = narrowest_uint_for_val<sieve_primes_cap>;

	class alignas(64) aligned64
	{
	public:
		uint8_t& operator[](size_t i) { return data[i]; }
		const uint8_t& operator[](size_t i) const { return data[i]; }
	private:
		std::array<uint8_t, 64> data{};
	};

}
