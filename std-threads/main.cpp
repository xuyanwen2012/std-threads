#include <thread>
#include <iostream>
#include <random>
#include <algorithm>
#include <mutex>
#include <condition_variable>
#include <pthread.h> // assume on linux

#include "morton.h"

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

template <typename NumeratorT, typename DenominatorT>
constexpr NumeratorT divide_and_round_up(const NumeratorT n,
                                         const DenominatorT d)
{
	static_assert(
		std::is_integral_v<NumeratorT> && std::is_integral_v<DenominatorT>,
		"DivideAndRoundUp is only intended for integral types.");
	// Static cast to undo integral promotion.
	return static_cast<NumeratorT>(n / d + (n % d != 0 ? 1 : 0));
}

void stick_this_thread_to_core(
	std::thread::native_handle_type thread_native_handle,
	int core_id)
{
	cpu_set_t cpu_set;
	CPU_ZERO(&cpu_set);
	CPU_SET(core_id, &cpu_set);

	auto ret = pthread_setaffinity_np(thread_native_handle, sizeof(cpu_set_t), &cpu_set);
	if (ret != 0)
	{
		std::cerr << "Error calling pthread_setaffinity_np: " << ret << '\n';
	}
}

// ---------------------------------------------------------------------
// Morton encoding (1->1 relation)
// ---------------------------------------------------------------------


// 1 -> 1 relation, so we can use iterator
void k_morton_code(
	const std::vector<glm::vec4>::iterator input_begin,
	const std::vector<glm::vec4>::iterator input_end,
	const std::vector<morton_t>::iterator morton_begin,
	const float min_coord,
	const float range) noexcept
{
	std::transform(input_begin, input_end, morton_begin,
	               [min_coord, range](const glm::vec4& xyz)
	               {
		               return shared::xyz_to_morton32(xyz, min_coord, range);
	               });
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
	std::vector<morton_t>& u_sort_alt,
	const int shift
)
{
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
}

void dispatch_binning_pass(const int n_threads,
                           const std::vector<morton_t>& u_sort,
                           std::vector<morton_t>& u_sort_alt,
                           const int shift)
{
	barrier barrier(n_threads);

	const auto chunk_size = divide_and_round_up(u_sort.size(), n_threads);

	sort.current_thread = n_threads - 1;
	std::fill_n(sort.bucket, BASE, 0);

	std::vector<std::thread> threads;
	threads.reserve(n_threads);
	for (auto i = 0; i < n_threads; i++)
	{
		threads.emplace_back(k_binning_pass,
		                     i,
		                     std::ref(barrier),
		                     u_sort.data() + i * chunk_size,
		                     (i == n_threads - 1)
			                     ? u_sort.data() + u_sort.size()
			                     : u_sort.data() + (i + 1) * chunk_size,
		                     std::ref(u_sort_alt),
		                     shift);
	}

	for (auto& t : threads)
	{
		t.join();
	}
}

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
	//constexpr auto n = 667;
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


	std::vector<std::thread> threads;
	const auto chunk_size = divide_and_round_up(n, n_threads);
	threads.reserve(n_threads);
	for (auto i = 0; i < n_threads; ++i)
	{
		threads.emplace_back(
			k_morton_code, u_points.begin() + i * chunk_size,
			(i == n_threads - 1)
				? u_points.end()
				: u_points.begin() + (i + 1) * chunk_size,
			u_morton.begin() + i * chunk_size, min_coord, range);
		const auto handle = threads.back().native_handle();
		stick_this_thread_to_core(handle, i);
	}

	for (auto& t : threads)
	{
		t.join();
	}

	std::vector u_sort_ref{u_morton};

	auto start = std::chrono::high_resolution_clock::now();

	dispatch_binning_pass(n_threads, u_morton, u_morton_alt, 0);
	dispatch_binning_pass(n_threads, u_morton_alt, u_morton, 8);
	dispatch_binning_pass(n_threads, u_morton, u_morton_alt, 16);
	dispatch_binning_pass(n_threads, u_morton_alt, u_morton, 24);

	auto end = std::chrono::high_resolution_clock::now();

	const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Time: " << duration << " ms\n";

	const auto is_sorted = std::is_sorted(u_morton.begin(), u_morton.end());
	std::cout << "Is sorted: " << std::boolalpha
		<< is_sorted << '\n';

	start = std::chrono::high_resolution_clock::now();

	std::sort(u_sort_ref.begin(), u_sort_ref.end());

	end = std::chrono::high_resolution_clock::now();

	std::cout << "std::sort time: "
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
		<< " ms\n";

	{
		const auto is_equal = std::equal(u_morton.begin(), u_morton.end(), u_sort_ref.begin());
		std::cout << "Is equal: " << std::boolalpha << is_equal << '\n';
	}

	return 0;
}
