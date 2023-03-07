#pragma once

#include <memory>

#include "config.hpp"
#include "math/math.hpp"
#include "trial_division/multibase_div_tests.hpp"

namespace mbp
{
	template<size_t pass>
	size_t merge_bitmasks(size_t number, sieve_container& sieve);
}
