#include "utility.hpp"

#include <bitset>
#include <chrono>
#include <iostream>

namespace mbp::util
{
	long long current_time_in_ms()
	{
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	}

	void print_as_bits(const uint64_t n)
	{
		std::stringstream ss;
		ss << std::bitset<8>(uint8_t(n >> (8 * 7))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 6))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 5))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 4))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 3))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 2))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 1))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 0)));
		std::cout << ss.str();
	}
	void print_as_bits(const uint32_t n)
	{
		std::stringstream ss;
		ss << std::bitset<8>(uint8_t(n >> (8 * 3))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 2))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 1))) << ' ';
		ss << std::bitset<8>(uint8_t(n >> (8 * 0)));
		std::cout << ss.str();
	}

	// via https://stackoverflow.com/a/12996028/2924233
	uint64_t hash(uint64_t x)
	{
		x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
		x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
		x = x ^ (x >> 31);
		return x;
	}
}
