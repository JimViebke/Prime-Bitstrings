#pragma once

#include <type_traits>
#include <cmath>

namespace mbp::franken_boost
{
	namespace detail
	{
		template <typename Real>
		inline constexpr Real sqrt_impl_2(Real x, Real s, Real s2)
		{
			return !(s < s2) ? s2 : sqrt_impl_2(x, (x / s + s) / 2, s);
		}

		template <typename Real>
		inline constexpr Real sqrt_impl_1(Real x, Real s)
		{
			return sqrt_impl_2(x, (x / s + s) / 2, s);
		}

		template <typename Real>
		inline constexpr Real sqrt_impl(Real x)
		{
			return sqrt_impl_1(x, x > 1 ? x : Real(1));
		}
	}

	template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
	inline constexpr double sqrt(Z x)
	{
		return detail::sqrt_impl<double>(static_cast<double>(x));
	}
}
