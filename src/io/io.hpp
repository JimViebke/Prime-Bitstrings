#pragma once

#include <filesystem>
#include <iostream>
#include <fstream>

#include "../math/math.hpp"
#include "../config.hpp"

namespace mbp
{
	static std::filesystem::path results_path;

	void set_up_results_path()
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

	size_t load_from_results()
	{
		// Load the largest (not necessarily last) number from the "results" log.
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

		if (div_test::max_pn_bitwidth <= std::bit_width(number))
		{
			std::cout << "error: max_pn_bitwidth should be larger than std::bit_width(number)\n";
			exit(EXIT_FAILURE);
		}

		std::cout << "Loaded " << number << " from " << results_path.generic_string() << '\n';

		return number;
	}

	void log_time()
	{
		time_t timestamp = time(0);
		tm now;
		localtime_s(&now, &timestamp);
		std::cout << std::setfill(' ') << std::setw(2) << ((now.tm_hour % 12 == 0) ? 12 : now.tm_hour % 12) << ':' << std::setfill('0') << std::setw(2) << now.tm_min << '\t';
	}

	void log_result(const mpz_class& n, size_t up_to_base)
	{
		log_time();

		std::stringstream ss;
		ss << bin_to_base(n, 10) << " is a p" << up_to_base << " (" << n << ")\n";

		std::cout << ss.str();

		// Don't write to file in benchmark mode
		if constexpr (!benchmark_mode)
		{
			std::ofstream ofs(results_path, std::ofstream::app);
			ofs << ss.str();
		}
	}

}
