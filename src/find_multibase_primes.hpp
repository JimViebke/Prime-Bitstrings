#pragma once

#include "config.hpp"
#include "math/math.hpp"
#include "trial_division/multibase_div_tests.hpp"

namespace mbp
{
	class find_multibase_primes
	{
	public:
		find_multibase_primes();

		void run();

	private:

		template<bool on_fast_path>
		inline_toggle uint64_t* div_tests(uint64_t* candidates_end);

		void full_primality_tests(const uint64_t* candidates,
								  const uint64_t* const candidates_end);

		div_test::full_div_tests full_div_tests{};

		gmp_randclass gmp_rand{ gmp_randinit_mt };
		mpz_class mpz_number{};

		count_passes(size_t a, ps15, b, c, d, e, f, bldt, bidt, b2, b3, b4, b5, passes, pc_hash);
	};

	void print_config();
}
