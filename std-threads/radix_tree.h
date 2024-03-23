#pragma once

#include <cassert>
#include <stdexcept>

// I am using only pointers because this gives me a unified front end for both CPU/and GPU
struct radix_tree
{
#define UNINITIALIZED (-1)

	const size_t capacity;

	// ------------------------
	// Essential Data
	// ------------------------
	int n_brt_nodes = UNINITIALIZED;

	uint8_t* u_prefix_n;
	bool* u_has_leaf_left;
	bool* u_has_leaf_right;
	int* u_left_child;
	int* u_parent;

	// ------------------------
	// Constructors
	// ------------------------

	radix_tree() = delete;

	explicit radix_tree(size_t n_to_allocate);

	radix_tree(const radix_tree&) = delete;
	radix_tree& operator=(const radix_tree&) = delete;
	radix_tree(radix_tree&&) = delete;
	radix_tree& operator=(radix_tree&&) = delete;

	~radix_tree();

	// ------------------------
	// Getter/Setters
	// ------------------------

	void set_n_nodes(const size_t n_nodes)
	{
		assert(n_nodes < capacity);
		n_brt_nodes = static_cast<int>(n_nodes);
	}

	[[nodiscard]] int n_nodes() const
	{
		if (n_brt_nodes == UNINITIALIZED)
			throw std::runtime_error("BRT nodes unset!!!");
		return n_brt_nodes;
	}
};
