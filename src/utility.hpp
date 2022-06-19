#pragma once

#include <chrono>

#include "direct.h"

inline auto current_time_in_ms()
{
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

inline auto current_time_in_us()
{
	return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

void create_folder(const std::string & path)
{
	_mkdir(path.c_str());
}
