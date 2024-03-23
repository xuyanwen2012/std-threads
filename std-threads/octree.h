#pragma once

#include <stdexcept>
#include <glm/glm.hpp>

struct octree
{
#define UNINITIALIZED (-1)

	const size_t capacity;

	// ------------------------
	// Essential Data
	// ------------------------

	int n_oct_nodes = UNINITIALIZED;

	// [Outputs]
	int (*u_children)[8];
	glm::vec4* u_corner;
	float* u_cell_size;
	int* u_child_node_mask;
	int* u_child_leaf_mask;

	// ------------------------
	// Constructors
	// ------------------------

	octree() = delete;

	explicit octree(size_t capacity);

	octree(const octree&) = delete;
	octree& operator=(const octree&) = delete;
	octree(octree&&) = delete;
	octree& operator=(octree&&) = delete;

	~octree();

	// ------------------------
	// Getter/Setters
	// ------------------------
	void set_n_nodes(const size_t n_nodes)
	{
		assert(n_nodes < capacity);
		n_oct_nodes = static_cast<int>(n_nodes);
	}

	[[nodiscard]] int n_nodes() const
	{
		if (n_oct_nodes == UNINITIALIZED)
			throw std::runtime_error("BRT nodes unset!!!");
		return n_oct_nodes;
	}
};
