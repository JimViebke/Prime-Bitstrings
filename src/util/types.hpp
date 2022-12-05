#pragma once

#include "../config.hpp"
#include "utility.hpp"
#include "bit_array.hpp"

namespace mbp
{
	using sieve_prime_t = util::narrowest_uint_for_val<sieve_primes_cap>;

	using sieve_container = mbp::bit_array<static_sieve_size>;
}
