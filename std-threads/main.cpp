#include <thread>
#include <iostream>
#include <random>
#include <algorithm>
#include <mutex>
#include <condition_variable>

//#include <pthread.h> // assume on linux

#include "morton.h"

#include "BS_thread_pool.hpp"
#include "BS_thread_pool_utils.hpp"

//
//void stick_this_thread_to_core(
//	std::thread::native_handle_type thread_native_handle, size_t core_id)
//{
//	cpu_set_t cpu_set;
//	CPU_ZERO(&cpu_set);
//	CPU_SET(core_id, &cpu_set);
//
//	auto ret =
//		pthread_setaffinity_np(thread_native_handle, sizeof(cpu_set_t), &cpu_set);
//	if (ret != 0)
//	{
//		std::cerr << "Error calling pthread_setaffinity_np: " << ret << '\n';
//	}
//}


class barrier
{
public:
	explicit barrier(const int count) : thread_count_(count), count_(count), generation_(0)
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
	int thread_count_;
	int count_;
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
#define DEBUG_PRINT(...) sync_out.println(__va_args__)
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
			std::transform(
				u_points.begin() + start, u_points.begin() + end,
				u_morton.begin() + start,
				[min_coord, range](const glm::vec4& xyz)
				{
					return shared::xyz_to_morton32(xyz, min_coord, range);
				});
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
	int current_thread = 0;
} sort;

void k_binning_pass(
	const int tid,
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
//void
dispatch_binning_pass(
	const int n_threads,
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


static std::unordered_map<std::string, size_t> config = {
	{"morton", 1},
	{"sort", 2},
	{"unique", 1},
	{"radix_tree", 0}
};


int main(const int argc, const char* argv[])
{
	auto n_threads = 1u;
	if (argc > 1)
	{
		n_threads = std::stoi(argv[1]);
	}

	const auto max_threads = std::thread::hardware_concurrency();

	std::cout << "Number of threads: " << n_threads << '\n';
	std::cout << "Max threads: " << max_threads << '\n';

	if (n_threads > max_threads)
	{
		std::cerr << "Number of threads exceeds max threads.\n";
		return 1;
	}

	//constexpr auto n = 1 << 20; // ~1M
	constexpr auto n = 1920 * 1080; // ~2M

	constexpr auto min_coord = 0.0f;
	constexpr auto range = 1024.0f;

	std::vector<glm::vec4> u_points(n);

	std::vector<morton_t> u_morton(n);
	std::vector<morton_t> u_morton_alt(n);

	std::mt19937 gen(114514); // NOLINT(cert-msc51-cpp)
	std::uniform_real_distribution dis(min_coord, min_coord + range);
	std::generate(u_points.begin(), u_points.end(), [&dis, &gen]
	{
		return glm::vec4(dis(gen), dis(gen), dis(gen), 1.0f);
	});


	BS::timer timer;

	barrier sort_barrier(n_threads);

	timer.start();

	dispatch_morton_code(n_threads, u_points, u_morton, min_coord, range).wait();

	auto morton_timestamp = timer.current_ms();

	dispatch_binning_pass(n_threads, sort_barrier, u_morton, u_morton_alt, 0).wait();
	dispatch_binning_pass(n_threads, sort_barrier, u_morton_alt, u_morton, 8).wait();
	dispatch_binning_pass(n_threads, sort_barrier, u_morton, u_morton_alt, 16).wait();
	dispatch_binning_pass(n_threads, sort_barrier, u_morton_alt, u_morton, 24).wait();

	auto sort_timestamp = timer.current_ms() - morton_timestamp;

	timer.stop();

	std::cout << "==========================\n";
	std::cout << " Total Time spent: " << timer.ms() << " ms\n";
	std::cout << "--------------------------\n";
	std::cout << " Morton: " << morton_timestamp << " ms\n";
	std::cout << " Sort: " << sort_timestamp << " ms\n";
	std::cout << "==========================\n";

	const auto is_sorted = std::is_sorted(u_morton.begin(), u_morton.end());
	std::cout << "Is sorted: " << std::boolalpha
		<< is_sorted << '\n';

	return 0;
}
