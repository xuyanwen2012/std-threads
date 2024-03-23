#pragma once

#include <string>
#include <unordered_map>

[[nodiscard]] std::unordered_map<std::string, size_t> read_config(
	const std::string& filename);
