#pragma once

#include "config.hpp"
#include "utility.hpp"

namespace mbp
{
	using sieve_prime_t = util::narrowest_uint_for_val<sieve_primes_cap>;

	using uint128_t = __m128i;
	using uint256_t = __m256i;

	template<typename T, size_t N>
	class alignas(64) aligned64
	{
	public:
		T& operator[](size_t i) { return data[i]; }
		const T& operator[](size_t i) const { return data[i]; }
	private:
		std::array<T, N> data{};
	};

}
