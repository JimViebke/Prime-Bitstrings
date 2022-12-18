#pragma once

#include "util/types.hpp"

namespace mbp::prime_sieve
{
	const sieve_container generate_static_sieve();

	void set_up_sieve_offsets_cache(const size_t start);

	void partial_sieve(const uint64_t number,
					   sieve_container& sieve
					   count_passes(, size_t& ps15));

	uint64_t* gather_sieve_results(uint64_t* candidates,
								   const sieve_container& sieve,
								   const uint64_t number);

}
