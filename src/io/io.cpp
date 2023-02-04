#include "io.hpp"

#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <set>

#include "../config.hpp"
#include "../math/math.hpp"

void mbp::print_preamble()
{
	std::cout << "Built " __DATE__ " " __TIME__ << ((mbp::benchmark_mode) ? "... " : "\n");
}

namespace mbp
{
	static std::filesystem::path results_path;
}

void mbp::set_up_results_path()
{
	namespace fs = std::filesystem;
	fs::path dir = fs::current_path();

	// Search upward for Prime Bitstrings folder
	while (dir.filename() != "Prime Bitstrings")
	{
		if (dir.has_parent_path())
		{
			dir = dir.parent_path();
		}
		else
		{
			std::cout << "Could not find " << results_filename << " along " << fs::current_path() << '\n';
			exit(EXIT_FAILURE);
		}
	}

	results_path = dir;
	results_path.append(results_filename);
}

size_t mbp::load_from_results()
{
	// Load the largest (not necessarily last) number from the results log.
	set_up_results_path();

	std::ifstream ifs(results_path);
	size_t number = 0;

	std::string str;
	while (std::getline(ifs, str))
	{
		if (str.empty()) continue;

		const auto start = str.find('(');
		const auto end = str.find(')');

		if (start == std::string::npos ||
			end == std::string::npos ||
			start > end) continue;

		str = str.substr(start + 1, (end - start) - 1);

		size_t v = 0;
		auto r = std::from_chars(str.c_str(), str.c_str() + str.size(), v);

		if (r.ec != std::errc{})
		{
			std::cout << "Read bad entry: " << v << ", error: " << size_t(r.ec) << '\n';
			continue;
		}

		if (v > number)
			number = v;
	}

	std::cout << "Loaded " << number << " from " << results_path.generic_string() << '\n';

	//number -= (number % mbp::block_size);
	//std::cout << "Assigning search blocks from " << number << '\n';

	return number;
}

void mbp::log_time()
{
	time_t timestamp = time(0);
	tm now;
	localtime_s(&now, &timestamp);
	std::cout << std::setw(2) << ((now.tm_hour % 12 == 0) ? 12 : now.tm_hour % 12) << ':' << std::setfill('0') << std::setw(2) << now.tm_min << std::setfill(' ') << '\t';
}

void mbp::log_result(const uint64_t n, const size_t up_to_base)
{
	log_time();

	std::stringstream ss;
	ss << std::bitset<64>(n).to_string().c_str() + std::countl_zero(n) << " is a p" << up_to_base;
	if (up_to_base < 10) ss << ' '; // padding space for alignment
	ss << " (" << n << ')';

	static auto last_perf_time = util::current_time_in_ms();
	static uint64_t last_n = 0;
	static uint64_t results_hash = 0xdeadbeef;

	const auto perf_time = util::current_time_in_ms();

	// Don't log to file or print perf info for the first result
	if (last_n != 0)
	{
		// Don't log in benchmark mode
		if constexpr (!benchmark_mode)
		{
			std::ofstream ofs(results_path, std::ofstream::app);
			ofs << ss.str() << '\n';
		}

		// Continue populating the stringstream

		results_hash = util::hash(results_hash ^ n);
		ss << "    " << std::hex << std::setfill('0') << std::setw(16) << results_hash << std::dec;

		const double elapsed_seconds = double(perf_time - last_perf_time) / 1'000.0;
		const double ints_per_second = double(n - last_n) / elapsed_seconds;
		const double billions_of_ints_per_second = ints_per_second / 1'000'000'000.0;
		ss << "    " << std::setfill(' ') << std::setw(4) << std::setprecision(1) << std::fixed << std::right
			<< billions_of_ints_per_second << " B ints/second";
	}

	ss << '\n';
	std::cout << ss.str();

	last_perf_time = perf_time;
	last_n = n;
}



namespace mpb::io
{
	using namespace mbp::io;

	namespace detail
	{
		struct
		{
			std::mutex mutex;

			size_t next_block_to_assign;
			std::set<size_t> assigned_blocks;

			// Don't log results unless all earlier results have been logged
			size_t next_block_to_log;
			std::set<block> completed_blocks;
		} state;

		block get_next_block()
		{
			// Construct an object representing the next block of work
			block next{ .start = state.next_block_to_assign };

			// Increment the next block to start
			state.next_block_to_assign += block_size;

			// Note that the block has been assigned
			state.assigned_blocks.insert(next.start);

			return next;
		}
	}

	block get_block()
	{
		std::lock_guard<std::mutex> lock(detail::state.mutex);
		return detail::get_next_block();
	}

	block log_block_and_get_next_block(const block& completed)
	{
		std::lock_guard<std::mutex> lock(detail::state.mutex);

		using namespace detail;

		{
			auto it = state.assigned_blocks.find(completed.start);
			if (it == state.assigned_blocks.cend())
			{
				std::cout << "A block was completed, but it was not known to have been assigned. This should never happen." << std::endl;
				std::cin.ignore(); // Should never happen. Pause execution.
			}
			state.assigned_blocks.erase(it);
		}

		// For the sake of simpler code, always append the completed block to the queue.
		state.completed_blocks.insert(completed);

		// Log as many completed blocks as we can
		for (auto it = state.completed_blocks.begin();
			 it != state.completed_blocks.end() && it->start == state.next_block_to_log;
			 it = state.completed_blocks.begin())
		{
			for (const auto& mbp : it->multibase_primes)
				mbp::log_result(mbp.bitstring, mbp.up_to_base);

			state.next_block_to_log += block_size;

			state.completed_blocks.erase(it);
		}

		// Return the next block
		return detail::get_next_block();
	}

}

