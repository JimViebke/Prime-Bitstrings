#pragma once

#include "../config.hpp"
#include "../math/math.hpp"
#include "../util/types.hpp"
#include "../util/utility.hpp"

namespace mbp::div_test
{
	using base_t = narrowest_uint_for_val<up_to_base>;
	using prime_idx_t = narrowest_uint_for_val<div_test::n_of_primes>;
	using n_of_remainders_t = narrowest_uint_for_val<div_test::max_remainders>;
	using remainder_t = narrowest_uint_for_val<mbp::small_primes_lookup[div_test::n_of_primes - 1]>;

	class div_test_t
	{
	public:
	#if analyze_div_tests
		bool used = false;
		uint32_t hits = 0;
		base_t base = 0;
	#endif

		prime_idx_t prime_idx = 0;
		n_of_remainders_t n_of_remainders = 0; // is also the index of the req'd bitmask
		bool is_first_with_n_remainders = false;

		std::array<remainder_t, 64> remainders{ 0 };
	};

	namespace detail
	{
		class uncompressed_div_test_t
		{
		public:
			uint32_t hits = 0;
			base_t base = 0;

			prime_idx_t prime_idx = 0;
			n_of_remainders_t n_of_remainders = 0;
			std::array<remainder_t, 64> remainders{ 0 };
		};
	}

}
