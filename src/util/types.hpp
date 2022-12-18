#pragma once

#include "../config.hpp"
#include "utility.hpp"
#include "bit_array.hpp"

namespace mbp
{
	using sieve_prime_t = util::narrowest_uint_for_val<prime_sieve::largest_sieve_prime>;

	using sieve_container = mbp::bit_array<prime_sieve::static_sieve_size>;
}
