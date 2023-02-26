#pragma once

#include "util/types.hpp"

namespace mbp::prime_sieve
{
	const sieve_container generate_static_sieve();

	void set_up_sieve_offsets_cache(const size_t start);

	void partial_sieve(sieve_container& sieve,
					   const size_t sieve_popcount);

	uint64_t* gather_sieve_results(uint64_t* candidates,
								   const sieve_container& sieve,
								   const uint64_t number);

	namespace detail
	{
		void verify_sieve_offset_cache(const uint64_t start);
	}
}
