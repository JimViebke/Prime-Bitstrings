#pragma once

#include "config.hpp"

#include "math/math.hpp"

namespace mbp
{
	class find_multibase_primes
	{
	public:
		find_multibase_primes();

		void run();

	private:

		template<bool on_fast_path>
		void main_loop(const size_t number);

		void full_primality_tests(const uint64_t* candidates_begin,
								  const uint64_t* const candidates_end);

		gmp_randclass gmp_rand{ gmp_randinit_mt };
		mpz_class mpz_number = 0ull;

		count_passes(size_t a, ps15, b, c, d, e, f, g, h, bldt, bidt, b2, b3, b4, b5, passes, pc_hash = 0);
	};

	void print_config();
}
