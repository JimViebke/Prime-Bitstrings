#pragma once

// Bitset with 32-byte / 256-bit alignment

#include <algorithm>

#include "simd.hpp"

namespace mbp
{
	// Suppress warning C4324: 'structure was padded due to alignment'
#pragma warning(push)
#pragma warning(disable: 4324)
	template<size_t n_elements>
	class alignas(32) bit_array
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

		__forceinline size_t count_bits() const
		{
			constexpr static uint256_t static_nybble_lookup{ .m256i_u8{
				0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
				0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 } };

			const uint256_t nybble_mask = _mm256_set1_epi8(0b0000'1111);
			const uint256_t nybble_pc_lookup = _mm256_loadu_si256(&static_nybble_lookup);

			const uint256_t* const rounded_end = ((uint256_t*)_data.data()) + (n_elements / 256);

			uint256_t ymm_pc{};
			for (const uint256_t* in = (uint256_t*)_data.data(); in != rounded_end; )
			{
				const uint256_t* const stride_end = in + std::min(rounded_end - in,
																  15ll); // uint8::max / 8

				uint256_t inner_pc{};
				for (; in < stride_end; ++in)
				{
					const uint256_t ymm0 = _mm256_loadu_si256(in);

					const uint256_t nybbles_lo = _mm256_and_si256(ymm0, nybble_mask);
					const uint256_t nybbles_hi = _mm256_and_si256(_mm256_srli_epi64(ymm0, 4), nybble_mask);

					inner_pc = _mm256_add_epi8(inner_pc, _mm256_shuffle_epi8(nybble_pc_lookup, nybbles_lo));
					inner_pc = _mm256_add_epi8(inner_pc, _mm256_shuffle_epi8(nybble_pc_lookup, nybbles_hi));
				}

				ymm_pc = _mm256_add_epi64(ymm_pc, _mm256_sad_epu8(inner_pc, _mm256_setzero_si256()));
			}

			// add high and low lanes
			ymm_pc = _mm256_add_epi64(ymm_pc, _mm256_permute2x128_si256(ymm_pc, ymm_pc, 1));
			// add high and low elements
			ymm_pc = _mm256_add_epi64(ymm_pc, _mm256_srli_si256(ymm_pc, 8));
			// extract
			size_t pc = ymm_pc.m256i_u64[0];

			const uint64_t* in = (const uint64_t*)rounded_end;

			constexpr size_t extra_bytes = (n_elements % 256) / 8;
			if constexpr (extra_bytes >= 8 * 1) pc += _mm_popcnt_u64(*in++);
			if constexpr (extra_bytes >= 8 * 2) pc += _mm_popcnt_u64(*in++);
			if constexpr (extra_bytes >= 8 * 3) pc += _mm_popcnt_u64(*in++);

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
