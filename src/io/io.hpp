#pragma once

#include <vector>

namespace mbp
{
	void print_preamble();

	void set_up_results_path();
	size_t load_from_results();

	void log_time();
	void log_result(const size_t n, const size_t up_to_base);



	// New threadsafe mechanisms. Not necessarily compatible with the above.

	namespace io
	{
		constexpr size_t num_threads = 1; // should be passed by command line instead
		constexpr size_t block_size = 1'000'000'000;

		struct multibase_prime
		{
			size_t bitstring;
			size_t up_to_base;
		};

		struct block
		{
			const size_t start;
			std::vector<multibase_prime> multibase_primes;

			bool operator<(const auto& other) const { return start < other.start; }
		};

		block get_block();

		block log_block_and_get_next_block(const block& completed);
	}

}
