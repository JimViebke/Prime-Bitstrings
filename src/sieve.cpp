
#include "math/math.hpp"
#include "sieve.hpp"

namespace mbp
{
	namespace detail
	{
		constexpr uint64_t n_sieve_chunks = (sieve_container::size() / 64) + 1;
		constexpr size_t max_pc = 27; // 27 is the highest pc of 64-bit chunks in the static sieve

		using chunk_count_t = util::narrowest_uint_for_val<n_sieve_chunks>;
		static std::array<chunk_count_t, max_pc + 1> n_chunks_with_pc{};
		// [popcount][chunk, index]
		static std::array<std::array<uint64_t, n_sieve_chunks * 2>, max_pc + 1> sorted_chunks{};

		template<size_t n_bits>
		__forceinline void extract_candidates(uint64_t& chunk,
											  const uint64_t idx,
											  uint64_t*& candidates)
		{
			// read the next bit in the chunk
			const size_t tzcnt = _tzcnt_u64(chunk);

			// reset the bit we just read
			if constexpr (n_bits > 1)
				chunk = _blsr_u64(chunk);

			// store the chunk's index plus the bit's index
			*candidates++ = idx + tzcnt;

			if constexpr (n_bits - 1 > 0)
				extract_candidates<n_bits - 1>(chunk, idx, candidates);
		}

		template<size_t popcount>
		__forceinline void extract_candidates_with_popcount(uint64_t*& candidates)
		{
			for (size_t j = 0, n_chunks = n_chunks_with_pc[popcount]; j < n_chunks; ++j)
			{
				uint64_t chunk = sorted_chunks[popcount][j * 2];
				uint64_t index = sorted_chunks[popcount][j * 2 + 1] * 64;

				// generate instructions to extract n candidates
				extract_candidates<popcount>(chunk, index, candidates);
			}
		}
	}

	// sort 64-bit chunks of the sieve by popcount, then use custom loops that perform the exact number of required reads
	uint64_t* gather_sieve_results(uint64_t* candidates,
								   const sieve_container& sieve,
								   const uint64_t number)
	{
		using namespace detail;

		const uint64_t* sieve_data = (uint64_t*)sieve.data();

		// autovectorizes to three large writes
		for (auto& pc : n_chunks_with_pc)
			pc = 0;

		for (size_t i = 0; i < n_sieve_chunks; ++i)
		{
			const uint64_t chunk = sieve_data[i];

			const size_t pc = pop_count(chunk);
			const size_t idx = n_chunks_with_pc[pc]++;

			sorted_chunks[pc][idx * 2] = chunk;
			sorted_chunks[pc][idx * 2 + 1] = i;
		}

		uint64_t* const candidates_start = candidates;

		extract_candidates_with_popcount<1>(candidates);
		extract_candidates_with_popcount<2>(candidates);
		extract_candidates_with_popcount<3>(candidates);
		extract_candidates_with_popcount<4>(candidates);
		extract_candidates_with_popcount<5>(candidates);
		extract_candidates_with_popcount<6>(candidates);
		extract_candidates_with_popcount<7>(candidates);
		extract_candidates_with_popcount<8>(candidates);
		extract_candidates_with_popcount<9>(candidates);
		extract_candidates_with_popcount<10>(candidates);
		extract_candidates_with_popcount<11>(candidates);
		extract_candidates_with_popcount<12>(candidates);
		extract_candidates_with_popcount<13>(candidates);
		extract_candidates_with_popcount<14>(candidates);
		extract_candidates_with_popcount<15>(candidates);
		extract_candidates_with_popcount<16>(candidates);
		extract_candidates_with_popcount<17>(candidates);
		extract_candidates_with_popcount<18>(candidates);
		extract_candidates_with_popcount<19>(candidates);
		extract_candidates_with_popcount<20>(candidates);
		extract_candidates_with_popcount<21>(candidates);
		extract_candidates_with_popcount<22>(candidates);
		extract_candidates_with_popcount<23>(candidates);
		extract_candidates_with_popcount<24>(candidates);
		extract_candidates_with_popcount<25>(candidates);
		extract_candidates_with_popcount<26>(candidates);
		extract_candidates_with_popcount<27>(candidates);

		// autovectorizes
		for (uint64_t* ptr = candidates_start; ptr < candidates; ++ptr)
		{
			uint64_t candidate = *ptr;
			candidate *= 2;
			candidate += number;
			*ptr = candidate;
		}

		return candidates;
	}

}
