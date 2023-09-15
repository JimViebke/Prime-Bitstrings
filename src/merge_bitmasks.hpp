#pragma once

#include <memory>

#include "config.hpp"
#include "math/math.hpp"
#include "trial_division/multibase_div_tests.hpp"

namespace mbp
{
	template<size_t pass>
	size_t merge_bitmasks(size_t number, sieve_container& sieve);

	void merge_bitmasks_one_by_one(uint64_t number,
								   std::array<sieve_container, prime_sieve::steps>& sieves,
								   std::array<size_t, prime_sieve::steps>& sieve_popcounts);
}
