#include <fstream>
#include <iostream>

#include "json.hpp"
#include <unordered_map>

using json = nlohmann::json;

std::unordered_map<std::string, size_t> read_config(
	const std::string& filename)
{
	std::ifstream config_file("config.json");
	if (!config_file.is_open())
	{
		//throw std::runtime_error("Could not open config file" + filename);

		std::cout << "Could not open config file " << filename << '\n';
		std::cout << "!!!!!!!!!!!!!!!!!!Returning default configuration!!!!!!!!!!!!!!!!!!!!!!!\n";
		std::cout << "!!!!!!!!!!!!!!!!!!Returning default configuration!!!!!!!!!!!!!!!!!!!!!!!\n";
		std::cout << "!!!!!!!!!!!!!!!!!!Returning default configuration!!!!!!!!!!!!!!!!!!!!!!!\n";

		// return a default configuration
		return {
			{"radix_tree", 1}, {"unique", 1}, {"sorting", 1}, {"morton", 1}
		};
	}

	// Parse JSON
	json config;
	config_file >> config;

	// Hashtable to store configuration
	std::unordered_map<std::string, size_t> config_map;

	// Iterate over stages
	for (auto& stage : config["stages"])
	{
		std::string stage_name = stage["name"];
		int num_threads = stage["num_threads"];

		// Store stage configuration in hashtable
		config_map[stage_name] = num_threads;
	}

	for (const auto& stage : config_map)
	{
		std::cout << "Stage: " << stage.first << '\n';
		std::cout << "  Num Threads: " << stage.second << '\n';
		std::cout << '\n';
	}

	return config_map;
}
