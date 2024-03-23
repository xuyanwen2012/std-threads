#include <thread>
#include <iostream>
#include <random>
#include <algorithm>
#include <mutex>
#include <condition_variable>
#include <unordered_map>


#include "morton.h"
#include "oct_v2.h"

#include "BS_thread_pool.hpp"
#include "BS_thread_pool_utils.hpp"
#include "octree.h"
#include "radix_tree.h"
#include "read_config.hpp"

#if defined(__linux__) || defined(__APPLE__)
#include <pthread.h> // assume on linux

void stick_this_thread_to_core(
	std::thread::native_handle_type thread_native_handle, size_t core_id)
{
	cpu_set_t cpu_set;
	CPU_ZERO(&cpu_set);
	CPU_SET(core_id, &cpu_set);

	auto ret =
		pthread_setaffinity_np(thread_native_handle, sizeof(cpu_set_t), &cpu_set);
	if (ret != 0)
	{
		std::cerr << "Error calling pthread_setaffinity_np: " << ret << '\n';
	}
}
#endif

class barrier
{
public:
	explicit barrier(const size_t count) : thread_count_(count), count_(count), generation_(0)
	{
	}

	void wait()
	{
		std::unique_lock lock(mutex_);
		int gen = generation_;
		if (--count_ == 0)
		{
			generation_++;
			count_ = thread_count_;
			cond_.notify_all();
		}
		else
		{
			cond_.wait(lock, [this, gen] { return gen != generation_; });
		}
	}

private:
	std::mutex mutex_;
	std::condition_variable cond_;
	size_t thread_count_;
	size_t count_;
	int generation_;
};


// ---------------------------------------------------------------------
// Global variables
// ---------------------------------------------------------------------

BS::thread_pool pool;
BS::synced_stream sync_out;

void monitor_tasks()
{
	sync_out.println(pool.get_tasks_total(), " tasks total, ",
	                 pool.get_tasks_running(), " tasks running, ",
	                 pool.get_tasks_queued(), " tasks queued.");
}

#ifdef NDEBUG
#define DEBUG_PRINT(...) ((void)0)
#else
#define DEBUG_PRINT(...) sync_out.println(__VA_ARGS__)
#endif


// ---------------------------------------------------------------------
// Morton encoding (1->1 relation)
// ---------------------------------------------------------------------

[[nodiscard]]
BS::multi_future<void>
dispatch_morton_code(
	const size_t desired_n_threads,
	const std::vector<glm::vec4>& u_points,
	std::vector<morton_t>& u_morton,
	const float min_coord,
	const float range)
{
	const auto n = static_cast<int>(u_points.size());

	return pool.submit_blocks(
		0, n,
		[&u_points, &u_morton, min_coord, range](const int start, const int end)
		{
			DEBUG_PRINT("[tid ", std::this_thread::get_id(), "] started. (morton)");

			std::transform(
				u_points.begin() + start, u_points.begin() + end,
				u_morton.begin() + start,
				[min_coord, range](const glm::vec4& xyz)
				{
					return shared::xyz_to_morton32(xyz, min_coord, range);
				});

			DEBUG_PRINT("[tid ", std::this_thread::get_id(), "] ended. (morton)");
		},
		desired_n_threads);
}


// ---------------------------------------------------------------------
// Radix Sort (challenging)
// ---------------------------------------------------------------------

constexpr int BASE_BITS = 8;
constexpr int BASE = (1 << BASE_BITS); // 256
constexpr int MASK = (BASE - 1); // 0xFF

constexpr int DIGITS(const unsigned int v, const int shift)
{
	return (v >> shift) & MASK;
}

struct
{
	std::mutex mtx;
	int bucket[BASE] = {}; // shared among threads
	std::condition_variable cv;
	size_t current_thread = 0;
} sort;

void k_binning_pass(
	const size_t tid,
	barrier& barrier,
	const morton_t* u_sort_begin,
	const morton_t* u_sort_end,
	std::vector<morton_t>& u_sort_alt, // output
	const int shift
)
{
	DEBUG_PRINT("[tid ", tid, "] started. (Binning, shift=", shift, ")");

	int local_bucket[BASE] = {};

	// compute histogram (local)
	std::for_each(u_sort_begin, u_sort_end, [shift, &local_bucket](const morton_t& code)
	{
		++local_bucket[DIGITS(code, shift)];
	});

	std::unique_lock lck(sort.mtx);

	// update to shared bucket
	for (auto i = 0; i < BASE; ++i)
	{
		sort.bucket[i] += local_bucket[i];
	}

	lck.unlock();

	barrier.wait();


	if (tid == 0)
	{
		std::partial_sum(std::begin(sort.bucket), std::end(sort.bucket),
		                 std::begin(sort.bucket));
	}

	barrier.wait();

	lck.lock();
	sort.cv.wait(lck, [&] { return tid == sort.current_thread; });

	// update the local_bucket from the shared bucket
	for (auto i = 0; i < BASE; i++)
	{
		sort.bucket[i] -= local_bucket[i];
		local_bucket[i] = sort.bucket[i];
	}

	--sort.current_thread;
	sort.cv.notify_all();

	lck.unlock();

	std::for_each(u_sort_begin, u_sort_end,
	              [shift, &local_bucket, &u_sort_alt](const morton_t& code)
	              {
		              u_sort_alt[local_bucket[DIGITS(code, shift)]++] = code;
	              });

	DEBUG_PRINT("[tid ", tid, "] ended. (Binning, shift=", shift, ")");
}

// I know they have implemented this in the library, but it was private, and I cant modify the library because
// I want to make sure for all other project can use the library as is. 
//
template <typename T>
class [[nodiscard]] my_blocks
{
public:
	/**
	 * @brief Construct a `blocks` object with the given specifications.
	 *
	 * @param first_index_ The first index in the range.
	 * @param index_after_last_ The index after the last index in the range.
	 * @param num_blocks_ The desired number of blocks to divide the range into.
	 */
	my_blocks(const T first_index_, const T index_after_last_,
	          const size_t num_blocks_)
		: first_index(first_index_),
		  index_after_last(index_after_last_),
		  num_blocks(num_blocks_)
	{
		if (index_after_last > first_index)
		{
			const size_t total_size =
				static_cast<size_t>(index_after_last - first_index);
			if (num_blocks > total_size) num_blocks = total_size;
			block_size = total_size / num_blocks;
			remainder = total_size % num_blocks;
			if (block_size == 0)
			{
				block_size = 1;
				num_blocks = (total_size > 1) ? total_size : 1;
			}
		}
		else
		{
			num_blocks = 0;
		}
	}

	/**
	 * @brief Get the first index of a block.
	 *
	 * @param block The block number.
	 * @return The first index.
	 */
	[[nodiscard]] T start(const size_t block) const
	{
		return first_index + static_cast<T>(block * block_size) +
			static_cast<T>(block < remainder ? block : remainder);
	}

	/**
	 * @brief Get the index after the last index of a block.
	 *
	 * @param block The block number.
	 * @return The index after the last index.
	 */
	[[nodiscard]] T end(const size_t block) const
	{
		return (block == num_blocks - 1) ? index_after_last : start(block + 1);
	}

	/**
	 * @brief Get the number of blocks. Note that this may be different than the
	 * desired number of blocks that was passed to the constructor.
	 *
	 * @return The number of blocks.
	 */
	[[nodiscard]] size_t get_num_blocks() const { return num_blocks; }

private:
	/**
	 * @brief The size of each block (except possibly the last block).
	 */
	size_t block_size = 0;

	/**
	 * @brief The first index in the range.
	 */
	T first_index = 0;

	/**
	 * @brief The index after the last index in the range.
	 */
	T index_after_last = 0;

	/**
	 * @brief The number of blocks.
	 */
	size_t num_blocks = 0;

	/**
	 * @brief The remainder obtained after dividing the total size by the number
	 * of blocks.
	 */
	size_t remainder = 0;
}; // class blocks


[[nodiscard]]
BS::multi_future<void>
dispatch_binning_pass(
	const size_t n_threads,
	barrier& barrier,
	const std::vector<morton_t>& u_sort,
	std::vector<morton_t>& u_sort_alt,
	const int shift)
{
	constexpr auto first_index = 0;
	const auto index_after_last = static_cast<int>(u_sort.size());

	const my_blocks blks(first_index, index_after_last, n_threads);

	BS::multi_future<void> future;
	future.reserve(blks.get_num_blocks());

	std::fill_n(sort.bucket, BASE, 0);
	sort.current_thread = n_threads - 1;

	// I could have used the simpler API, but I need the 'blk' index for my kernel

	for (size_t blk = 0; blk < blks.get_num_blocks(); ++blk)
	{
		future.push_back(
			pool.submit_task([start = blks.start(blk), end = blks.end(blk), blk, &barrier, &u_sort, &u_sort_alt, shift]
			{
				k_binning_pass(static_cast<int>(blk), barrier, u_sort.data() + start, u_sort.data() + end, u_sort_alt,
				               shift);
			}));
	}

	return future;
}

// ---------------------------------------------------------------------
// Unique (for CPU, should only use a single thread), we have small problem size
// 
// And this future will return the number of unique elements
// 
// ---------------------------------------------------------------------

[[nodiscard]]
std::future<int>
dispatch_unique(
	const std::vector<morton_t>& u_sort,
	std::vector<morton_t>& u_sort_unique)
{
	return pool.submit_task([&]
	{
		DEBUG_PRINT("[tid ", std::this_thread::get_id(), "] started. (unique)");

		const auto last =
			std::unique_copy(u_sort.begin(), u_sort.end(), u_sort_unique.begin());

		const auto n_unique = std::distance(u_sort_unique.begin(), last);

		DEBUG_PRINT("[tid ", std::this_thread::get_id(), "] ended. (unqiue = ", n_unique, ")");

		return static_cast<int>(n_unique);
	});
}


// ---------------------------------------------------------------------
// Radix Tree (1->1 relation)
// ---------------------------------------------------------------------


#if defined(__GNUC__) || defined(__clang__)
#define CLZ(x) __builtin_clz(x)
#elif defined(_MSC_VER)
#include <intrin.h>
#define CLZ(x) _lzcnt_u32(x)
#else
#error "CLZ not supported on this platform"
#endif

unsigned int ceil_div_u32(const unsigned int a, const unsigned int b)
{
	assert(b != 0);
	return (a + b - 1) / b;
}

uint8_t delta_u32(const unsigned int a, const unsigned int b)
{
	[[maybe_unused]] constexpr unsigned int bit1_mask =
		static_cast<unsigned int>(1) << (sizeof(a) * 8 - 1);
	assert((a & bit1_mask) == 0);
	assert((b & bit1_mask) == 0);
	return static_cast<uint8_t>(CLZ(a ^ b) - 1);
}

int log2_ceil_u32(const unsigned int x)
{
	// Counting from LSB to MSB, number of bits before last '1'
	// This is floor(log(x))
	const auto n_lower_bits = ((8 * sizeof(x)) - CLZ(x) - 1);

	// Add 1 if 2^n_lower_bits is less than x
	//     (i.e. we rounded down because x was not a power of 2)
	return static_cast<int>(n_lower_bits + ((1 << n_lower_bits) < x));
}


inline void process_radix_tree_i(const int i, const int n /*n_brt_nodes*/,
                                 const morton_t* codes,
                                 const radix_tree* out_brt)
{
	// 'i' is the iterator within a chunk
	// 'codes' is the base address of the whole data, for each chunk, we need to
	// use the offset 'out_brt' is the base address of the whole data, for each
	// chunk, we need to use the offset

	const auto code_i = codes[i];

	const auto prefix_n = out_brt->u_prefix_n;
	const auto has_leaf_left = out_brt->u_has_leaf_left;
	const auto has_leaf_right = out_brt->u_has_leaf_right;
	const auto left_child = out_brt->u_left_child;
	const auto parent = out_brt->u_parent;

	// Determine direction of the range (+1 or -1)
	int d;
	if (i == 0)
	{
		d = 1;
	}
	else
	{
		const auto delta_diff_right = delta_u32(code_i, codes[i + 1]);
		const auto delta_diff_left = delta_u32(code_i, codes[i - 1]);
		const auto direction_difference = delta_diff_right - delta_diff_left;
		d = (direction_difference > 0) - (direction_difference < 0);
	}

	// Compute upper bound for the length of the range

	auto l = 0;
	if (i == 0)
	{
		// First node is root, covering whole tree
		l = n - 1;
	}
	else
	{
		const auto delta_min = delta_u32(code_i, codes[i - d]);
		auto l_max = 2;
		// Cast to ptrdiff_t so in case the result is negative (since d is +/- 1),
		// we can catch it and not index out of bounds
		while (i + static_cast<std::ptrdiff_t>(l_max) * d >= 0 &&
			i + l_max * d <= n &&
			delta_u32(code_i, codes[i + l_max * d]) > delta_min)
		{
			l_max *= 2;
		}
		const auto l_cutoff = (d == -1) ? i : n - i;
		int t;
		int divisor;
		// Find the other end using binary search
		for (t = l_max / 2, divisor = 2; t >= 1;
		     divisor *= 2, t = l_max / divisor)
		{
			if (l + t <= l_cutoff &&
				delta_u32(code_i, codes[i + (l + t) * d]) > delta_min)
			{
				l += t;
			}
		}
	}

	const auto j = i + l * d;

	// Find the split position using binary search
	const auto delta_node = delta_u32(codes[i], codes[j]);
	prefix_n[i] = delta_node;
	auto s = 0;
	const auto max_divisor = 1 << log2_ceil_u32(l);
	auto divisor = 2;
	const auto s_cutoff = (d == -1) ? i - 1 : n - i - 1;
	for (auto t = ceil_div_u32(l, 2); divisor <= max_divisor;
	     divisor <<= 1, t = ceil_div_u32(l, divisor))
	{
		if (s + t <= s_cutoff &&
			delta_u32(code_i, codes[i + (s + t) * d]) > delta_node)
		{
			s += t;
		}
	}

	// Split position
	const auto gamma = i + s * d + std::min(d, 0);
	left_child[i] = gamma;
	has_leaf_left[i] = (std::min(i, j) == gamma);
	has_leaf_right[i] = (std::max(i, j) == gamma + 1);
	// Set parents of left and right children, if they aren't leaves
	// can't set this node as parent of its leaves, because the
	// leaf also represents an internal node with a differnent parent
	if (!has_leaf_left[i])
	{
		parent[gamma] = i;
	}
	if (!has_leaf_right[i])
	{
		parent[gamma + 1] = i;
	}
}

[[nodiscard]]
BS::multi_future<void>
dispatch_build_radix_tree(
	const size_t n_desired_threads,
	const std::vector<morton_t>& sorted_morton,
	radix_tree* radix_tree)
{
	return
		pool.submit_blocks(
			0, radix_tree->n_nodes(),
			[&](const int start, const int end)
			{
				for (auto i = start; i < end; ++i)
				{
					process_radix_tree_i(i, radix_tree->n_nodes(), sorted_morton.data(),
					                     radix_tree);
				}
			},
			n_desired_threads);
}


// ---------------------------------------------------------------------
// Edge Count (1->1 relation)
// ---------------------------------------------------------------------
#include "defines.h"

namespace shared
{
	H_D_I void process_edge_count_i(const int i, const uint8_t* prefix_n,
	                                const int* parents, int* edge_count)
	{
		const auto my_depth = prefix_n[i] / 3;
		const auto parent_depth = prefix_n[parents[i]] / 3;
		edge_count[i] = my_depth - parent_depth;
	}
} // namespace shared

[[nodiscard]]
BS::multi_future<void>
dispatch_edge_count(
	const size_t n_desired_threads,
	const std::unique_ptr<radix_tree>& radix_tree,
	std::vector<int>& edge_count)
{
	return
		pool.submit_blocks(
			0, radix_tree->n_nodes(),
			[&](const int start, const int end)
			{
				for (auto i = start; i < end; ++i)
				{
					shared::process_edge_count_i(i, radix_tree->u_prefix_n, radix_tree->u_parent,
					                             edge_count.data());
				}
			}
			, n_desired_threads);
}

// ---------------------------------------------------------------------
// Edge Offset (prefix sum, for CPU, should use a single thread)
// this will return the number of octree nodes
// ---------------------------------------------------------------------

[[nodiscard]]
std::future<int>
dispatch_edge_offset(
	const int n_brt_nodes,
	const std::vector<int>& edge_count,
	std::vector<int>& edge_offset)
{
	// should not use ".end()", we allocated more than actual data
	// we need to use our own "n_brt_nodes"
	return pool.submit_task([&]
	{
		std::partial_sum(edge_count.begin(), edge_count.begin() + n_brt_nodes,
		                 edge_offset.begin());
		return edge_offset[n_brt_nodes - 1];
	});
}


// ---------------------------------------------------------------------
// Octree (1->1 relation, but has a lot of input)
// ---------------------------------------------------------------------

[[nodiscard]]
BS::multi_future<void>
dispatch_build_octree(
	const size_t n_desired_threads,
	const std::vector<int>& edge_count,
	const std::vector<int>& edge_offset,
	const std::vector<morton_t>& sorted_unique_morton,
	const std::unique_ptr<radix_tree>& radix_tree,
	const std::unique_ptr<octree>& octree, // output
	const float min_coord,
	const float range
)
{
	return pool.submit_blocks(
		0, radix_tree->n_nodes(),
		[&](const int start, const int end)
		{
			for (auto i = start; i < end; ++i)
			{
				shared::v2::process_oct_node(i, octree->u_children, octree->u_corner,
				                             octree->u_cell_size, octree->u_child_node_mask,
				                             edge_offset.data(), edge_count.data(),
				                             sorted_unique_morton.data(), radix_tree->u_prefix_n,
				                             radix_tree->u_parent, min_coord, range);
			}
		},
		n_desired_threads);
}

int main(const int argc, const char* argv[])
{
	auto test_n_threads = 1u;
	if (argc > 1)
	{
		test_n_threads = std::stoi(argv[1]);
	}

	const auto max_cores = std::thread::hardware_concurrency();


	std::cout << "Max CPU Cores: " << max_cores << '\n';

	if (test_n_threads > max_cores)
	{
		std::cerr << "Number of threads exceeds max threads.\n";
		return 1;
	}

	const auto config = read_config("config.json");

	const auto morton_n_thread = config.at("morton");
	const auto sort_n_thread = config.at("sorting");
	//constexpr auto unique_n_thread = 1;
	const auto brt_n_thread = config.at("radix_tree");
	constexpr auto edge_count_n_thread = 1;
	//constexpr auto edge_offset_n_thread = 1;
	constexpr auto octree_n_thread = 1;


	// make sure they add up to the total number of threads

	//if (const auto n_threads = morton_n_thread + sort_n_thread +
	//		unique_n_thread + brt_n_thread;
	//	n_threads > max_cores)
	//{
	//	std::cerr << "Number of threads does not match the sum of the stages.\n";
	//	return 1;
	//}

	// ---------------------------------------------------------------------
	// Initialization
	// ---------------------------------------------------------------------


	//constexpr auto capacity = 1 << 20; // ~1M
	constexpr auto n = 1920 * 1080; // ~2M

	constexpr auto min_coord = 0.0f;
	constexpr auto range = 1024.0f;
	constexpr auto educated_guess = 0.55; // from my experiment, #octree nodes < 55%, so allocate less 

	std::vector<glm::vec4> u_points(n);
	std::vector<morton_t> u_morton(n);
	std::vector<morton_t> u_morton_alt(n);
	auto brt = std::make_unique<radix_tree>(n);
	std::vector<int> u_edge_count(n);
	std::vector<int> u_edge_offset(n);
	auto oct = std::make_unique<octree>(n * educated_guess);

	std::mt19937 gen(114514); // NOLINT(cert-msc51-cpp)
	std::uniform_real_distribution dis(min_coord, min_coord + range);
	std::generate(u_points.begin(), u_points.end(), [&dis, &gen]
	{
		return glm::vec4(dis(gen), dis(gen), dis(gen), 1.0f);
	});


	BS::timer timer;

	barrier sort_barrier(sort_n_thread);


	// ---------------------------------------------------------------------
	// Main
	// ---------------------------------------------------------------------

	timer.start();

	dispatch_morton_code(morton_n_thread, u_points, u_morton, min_coord, range).wait();

	const auto morton_timestamp = timer.current_ms();

	dispatch_binning_pass(sort_n_thread, sort_barrier, u_morton, u_morton_alt, 0).wait();
	dispatch_binning_pass(sort_n_thread, sort_barrier, u_morton_alt, u_morton, 8).wait();
	dispatch_binning_pass(sort_n_thread, sort_barrier, u_morton, u_morton_alt, 16).wait();
	dispatch_binning_pass(sort_n_thread, sort_barrier, u_morton_alt, u_morton, 24).wait();

	const auto sort_timestamp = timer.current_ms();

	auto unique_future =
		dispatch_unique(u_morton, u_morton_alt);
	auto n_unique = unique_future.get();
	brt->set_n_nodes(n_unique - 1);


	const auto unique_timestamp = timer.current_ms();

	dispatch_build_radix_tree(brt_n_thread, u_morton_alt, brt.get()).wait();

	const auto brt_timestamp = timer.current_ms();

	dispatch_edge_count(edge_count_n_thread, brt, u_edge_count).wait();

	auto offset_future =
		dispatch_edge_offset(brt->n_nodes(), u_edge_count, u_edge_offset);
	auto n_octree_nodes = offset_future.get();
	oct->set_n_nodes(n_octree_nodes);

	const auto edge_offset_timestamp = timer.current_ms();

	dispatch_build_octree(octree_n_thread, u_edge_count, u_edge_offset, u_morton, brt, oct, min_coord, range).wait();

	timer.stop();

	// ---------------------------------------------------------------------
	// Validation
	// ---------------------------------------------------------------------

	std::cout << "==========================\n";
	std::cout << " Total Time spent: " << timer.ms() << " ms\n";
	std::cout << " n_unique = " << n_unique << '\n';
	std::cout << " n_brt_nodes = " << brt->n_nodes() << '\n';
	std::cout << " n_octree_nodes = " << n_octree_nodes << " (" << static_cast<double>(n_octree_nodes) / n * 100.0 <<
		"%)\n";
	std::cout << "--------------------------\n";
	std::cout << " Morton: " << morton_timestamp << " ms\n";
	std::cout << " Sort: " << sort_timestamp - morton_timestamp << " ms\n";
	std::cout << " Unique: " << unique_timestamp - sort_timestamp << " ms\n";
	std::cout << " BRT: " << brt_timestamp - unique_timestamp << " ms\n";
	std::cout << " Edge Count && Offset: " << edge_offset_timestamp - brt_timestamp << " ms\n";
	std::cout << " Octree: " << timer.current_ms() - edge_offset_timestamp << " ms\n";
	std::cout << "==========================\n";

	const auto is_sorted = std::is_sorted(u_morton.begin(), u_morton.end());
	std::cout << "Is sorted: " << std::boolalpha
		<< is_sorted << '\n';

	// peek 32 brt nodes

	//for (auto i = 0; i < 32; ++i)
	//{
	//	std::cout << "Node " << i << '\n';
	//	std::cout << "  Prefix: " << static_cast<int>(brt->u_prefix_n[i]) << '\n';
	//	std::cout << "  Has Left Leaf: " << std::boolalpha << brt->u_has_leaf_left[i] << '\n';
	//	std::cout << "  Has Right Leaf: " << std::boolalpha << brt->u_has_leaf_right[i] << '\n';
	//	std::cout << "  Left Child: " << brt->u_left_child[i] << '\n';
	//	std::cout << "  Parent: " << brt->u_parent[i] << '\n';
	//	std::cout << '\n';
	//}

	return 0;
}
