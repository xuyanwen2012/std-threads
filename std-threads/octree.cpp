#include "octree.h"

octree::octree(const size_t capacity) : capacity(capacity)
{
	u_children = new int[capacity][8];
	u_corner = new glm::vec4[capacity];
	u_cell_size = new float[capacity];
	u_child_node_mask = new int[capacity];
	u_child_leaf_mask = new int[capacity];
}

octree::~octree()
{
	delete[] u_children;
	delete[] u_corner;
	delete[] u_cell_size;
	delete[] u_child_node_mask;
	delete[] u_child_leaf_mask;
}
