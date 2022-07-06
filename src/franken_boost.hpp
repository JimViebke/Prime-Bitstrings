#pragma once

#include <type_traits>
#include <cmath>

namespace franken_boost
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

		template <typename T>
		inline constexpr bool isinf(T x)
		{
			if (BOOST_MATH_IS_CONSTANT_EVALUATED(x))
			{
				return x == std::numeric_limits<T>::infinity() || -x == std::numeric_limits<T>::infinity();
			}
			else
			{
				using std::isinf;

				if constexpr (!std::is_integral_v<T>)
				{
					return isinf(x);
				}
				else
				{
					return isinf(static_cast<double>(x));
				}
			}
		}

		template <typename T>
		inline constexpr bool isnan(T x)
		{
			if (BOOST_MATH_IS_CONSTANT_EVALUATED(x))
			{
				return x != x;
			}
			else
			{
				using std::isnan;

				if constexpr (!std::is_integral_v<T>)
				{
					return isnan(x);
				}
				else
				{
					return isnan(static_cast<double>(x));
				}
			}
		}
	}

	template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
	inline constexpr Real sqrt(Real x)
	{
		if (BOOST_MATH_IS_CONSTANT_EVALUATED(x))
		{
			return detail::isnan(x) ? std::numeric_limits<Real>::quiet_NaN() :
				detail::isinf(x) ? std::numeric_limits<Real>::infinity() :
				detail::sqrt_impl<Real>(x);
		}
		else
		{
			using std::sqrt;
			return sqrt(x);
		}
	}

	template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
	inline constexpr double sqrt(Z x)
	{
		return detail::sqrt_impl<double>(static_cast<double>(x));
	}
}
