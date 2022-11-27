#pragma once

#include "../config.hpp"
#include "utility.hpp"

namespace mbp
{
	using sieve_prime_t = util::narrowest_uint_for_val<sieve_primes_cap>;

	// Bitset with cache line alignment.
	// Supress warning C4324: 'structure was padded due to alignment'
#pragma warning(push)
#pragma warning(disable: 4324)
	template<size_t n_elements>
	class alignas(64) bit_array
	{
	public:
		bit_array() { clear_all(); }

		__forceinline void clear_bit(const size_t idx)
		{			
			_data[idx / 8] &= ~(1 << (idx % 8));
		}

		__forceinline void set_bit(const size_t idx)
		{
			_data[idx / 8] |= (1 << (idx % 8));
		}

		void set_all()
		{
			for (auto& el : _data)
				el = uint8_t(-1);

			constexpr size_t extra_bits = n_elements % 64;
			if constexpr (extra_bits > 0)
			{
				uint64_t* ptr = (uint64_t*)_data.data();
				ptr += (_data.size() / 8) - 1;

				*ptr >>= (64 - extra_bits);
			}
		}

		void clear_all()
		{
			for (auto& el : _data)
				el = uint8_t(0);
		}

		__forceinline uint64_t load_u64_chunk(const size_t idx) const
		{
			return *(uint64_t*)(&_data[idx / 8]);
		}

		__forceinline uint8_t* data()
		{
			return _data.data();
		}

		__forceinline const uint8_t* data() const
		{
			return _data.data();
		}

		constexpr static size_t size() { return n_elements; }

		size_t count_bits() const
		{
			const uint64_t* in = (const uint64_t*)_data.data();
			const uint64_t* const end = in + (_data.size() / 8);

			size_t pc = 0;

			for (; in != end; ++in)
			{
				pc += _mm_popcnt_u64(*in);
			}

			return pc;
		}

	private:
		std::array<uint8_t, ((n_elements / 64) + 1) * 8> _data{};
	};
#pragma warning(pop)

	using sieve_container = mbp::bit_array<static_sieve_size>;
}
