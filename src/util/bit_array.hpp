#pragma once

// Bitset with cache line alignment

namespace mbp
{
	// Suppress warning C4324: 'structure was padded due to alignment'
#pragma warning(push)
#pragma warning(disable: 4324)
	template<size_t n_elements>
	class alignas(64) bit_array
	{
	public:
		constexpr bit_array()
		{
			static_assert(size_in_bytes() % sizeof(uint64_t) == 0);

			clear_all();
		}

		constexpr __forceinline void clear_bit(const size_t idx)
		{
			_data[idx / 8] &= ~(1 << (idx % 8));
		}

		constexpr __forceinline void set_bit(const size_t idx)
		{
			_data[idx / 8] |= (1 << (idx % 8));
		}

		constexpr void set_all()
		{
			if (std::is_constant_evaluated())
			{
				for (auto& el : _data)
					el = uint8_t(-1);
			}
			else
			{
				uint64_t* in = (uint64_t*)_data.data();
				const uint64_t* const end = (uint64_t*)(_data.data() + _data.size());

				for (; in != end; ++in)
					*in = uint64_t(-1);
			}
		}

		constexpr void clear_all()
		{
			if (std::is_constant_evaluated())
			{
				for (auto& el : _data)
					el = 0;
			}
			else
			{
				uint64_t* in = (uint64_t*)&_data[0];
				const uint64_t* const end = (uint64_t*)(_data.data() + _data.size());

				for (; in != end; ++in)
					*in = 0;
			}
		}

		__forceinline uint64_t load_u64_chunk(const size_t idx) const
		{
			return *(uint64_t*)(&_data[idx / 8]);
		}

		__forceinline uint8_t* data() { return _data.data(); }

		__forceinline const uint8_t* data() const { return _data.data(); }

		constexpr static size_t size() { return n_elements; }

		size_t count_bits() const
		{
			const uint64_t* in = (uint64_t*)_data.data();
			const uint64_t* const end = (uint64_t*)(_data.data() + _data.size());

			size_t pc = 0;

			for (; in != end - 1; ++in)
				pc += _mm_popcnt_u64(*in);

			constexpr size_t extra_bits = n_elements % 64;
			if constexpr (extra_bits > 0)
			{
				// shift extra bits away so we get the right count
				pc += _mm_popcnt_u64(*in << (64 - extra_bits));
			}

			return pc;
		}

	private:
		constexpr static size_t size_in_bits()
		{
			return n_elements + 63 - ((n_elements + 63) % 64);
		}

		constexpr static size_t size_in_bytes()
		{
			return size_in_bits() / 8;
		}

		std::array<uint8_t, size_in_bytes()> _data{};
	};
#pragma warning(pop)
}
