/*
** Developed by Henning Lindemann
 */

#include <atomic>
#include <cassert>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <random>
#include <thread>
#include <unordered_set>
#include <vector>
#include "utilities/blocked_bloom_filter.hpp"
#include "utilities/one_hashing_blocked_bloom_filter.hpp"
#include "utilities/bloom_filter.hpp"
#include "utilities/runtime_class.hpp"

class StdSet
{
 private:
  std::unordered_set<uint64_t> set;

 public:
  void insert(uint64_t value) { set.insert(value); }

  bool contains(uint64_t value)
  {
    auto search = set.find(value);
    return search != set.end();
  }

  size_t size_in_bytes() { return 0; }

  size_t num_hash_functions_get() { return 1; }
};

std::vector<uint64_t> generate_testdata(uint64_t size, uint64_t seed)
{
  std::vector<uint64_t> test_data;
  std::mt19937 rng(seed);
  std::uniform_int_distribution<uint64_t> dist;

  test_data.reserve(size);
  for (uint64_t i = 0; i < size; i++)
  {
    test_data.push_back(dist(rng));
  }

  return test_data;
}

template <class AMD>
void benchmark(const char *name, AMD &amd, std::vector<uint64_t> &insert_data,
               std::vector<uint64_t> &test_data)
{
  RunTimeClass rt_insert{};
  for (const uint64_t &v : insert_data)
  {
    amd.insert(v);
  }
  const size_t micro_insert = rt_insert.elapsed();

  RunTimeClass rt_check_random{};
  size_t false_positives = 0;
  for (const uint64_t &v : test_data)
  {
    if (amd.contains(v))
    {
      false_positives++;
    }
  }
  const size_t micro_check_random = rt_check_random.elapsed();

  RunTimeClass rt_check_inserted{};
  for ([[maybe_unused]] const uint64_t &v : insert_data)
  {
    assert(amd.contains(v));
  }
  const size_t micro_check_inserted = rt_check_inserted.elapsed();
  const float false_positive_percent =
      100.0 * (float) false_positives / ((float) test_data.size());
  const uint64_t size_in_bytes = amd.size_in_bytes();
  const float bits_per_element = 8.0 * (float) size_in_bytes /
                                      ((float) insert_data.size());
  printf("%35s\t 1\t%7zu\t%7zu\t%7zu\t%10zu"
         "\t%10.5f\t%11" PRIu64 "\t%7.2f\t%2zu\n",
         name,
         micro_insert / size_t(1000),
         micro_check_random / size_t(1000),
         micro_check_inserted / size_t(1000),
         false_positives,
         false_positive_percent,
         size_in_bytes,
         bits_per_element,
         amd.num_hash_functions_get());
}

template <class AMD>
void benchmark_threaded(const char *name, AMD &amd,
                        std::vector<uint64_t> &insert_data,
                        std::vector<uint64_t> &test_data, uint32_t num_threads)
{
  if (num_threads == 1)
  {
    benchmark(name, amd, insert_data, test_data);
    return;
  }

  std::vector<std::thread> threads;
  RunTimeClass rt_insert{};
  threads.reserve(num_threads);
  for (uint32_t t = 0; t < num_threads; t++)
  {
    threads.push_back(std::thread([t, &amd, &insert_data, &num_threads]() {
      const uint64_t size = insert_data.size();
      const uint64_t start = (size / num_threads) * t;
      const uint64_t end = (t == num_threads - 1)
                             ? size
                             : (size / num_threads) * (t + 1);
      for (uint64_t i = start; i < end; i++)
      {
        amd.insert(insert_data[i]);
      }
    }));
  }

  for (auto &t : threads)
  {
    t.join();
  }
  const size_t micro_insert = rt_insert.elapsed();

  threads.clear();

  std::atomic<uint64_t> false_positives = 0;
  RunTimeClass rt_check_random{};
  for (uint32_t t = 0; t < num_threads; t++)
  {
    threads.push_back(
        std::thread([t, &amd, &test_data, &num_threads, &false_positives]() {
          uint64_t local_false_positives = 0;
          const uint64_t size = test_data.size();
          const uint64_t start = (size / num_threads) * t;
          const uint64_t end =
              (t == num_threads - 1) ? size : (size / num_threads) * (t + 1);
          for (uint64_t i = start; i < end; i++)
          {
            if (amd.contains(test_data[i]))
            {
              local_false_positives++;
            }
          }
          false_positives.fetch_add(local_false_positives,
                                    std::memory_order_relaxed);
        }));
  }

  for (auto &t : threads)
  {
    t.join();
  }
  const size_t micro_check_random = rt_check_random.elapsed();

  threads.clear();

  RunTimeClass rt_check_inserted{};
  for (uint32_t t = 0; t < num_threads; t++)
  {
    threads.push_back(std::thread([t, &amd, &insert_data, &num_threads]() {
      const uint64_t size = insert_data.size();
      const uint64_t start = (size / num_threads) * t;
      const uint64_t end =
          (t == num_threads - 1) ? size : (size / num_threads) * (t + 1);
      for (uint64_t i = start; i < end; i++)
      {
        assert(amd.contains(insert_data[i]));
      }
    }));
  }

  for (auto &t : threads)
  {
    t.join();
  }
  const size_t micro_check_inserted = rt_check_inserted.elapsed();
  const float false_positive_percent =
      100.0 * (float) false_positives.load() / ((float) test_data.size());
  const uint64_t size_in_bytes = amd.size_in_bytes();
  const float bits_per_element = 8.0 * (float) size_in_bytes /
                                      ((float) insert_data.size());
  printf("%35s\t%2u\t%7zu\t%7zu\t%7zu\t%10" PRIu64
         "\t%10.5f\t%11" PRIu64 "\t%7.2f\t%2zu\n",
         name,
         num_threads,
         micro_insert / size_t(1000),
         micro_check_random / size_t(1000),
         micro_check_inserted / size_t(1000),
         false_positives.load(),
         false_positive_percent,
         size_in_bytes,
         bits_per_element,
         amd.num_hash_functions_get());
}

int main(int /*argc*/, char* /**argv*/[])
{
  const uint64_t insert_size = 1 << 24;
  const uint64_t test_size = insert_size;
#if 1
  const std::vector<double> error{0.005, 0.0005, 0.00005};
  const std::vector<uint32_t> threads{1, 2, 3, 4, 5, 6, 8, 10, 14};
  const std::vector<size_t> multipliers{8, 12, 16, 24};
  const std::vector<size_t> num_hash_functions{3, 5, 9};
#else
  const std::vector<double> error{0.005, 0.0005, 0.00005};
  const std::vector<uint32_t> threads{6};
  const std::vector<size_t> multipliers{8, 12, 16, 24};
  const std::vector<size_t> num_hash_functions{3, 5, 9};
#endif

  std::vector<uint64_t> insert_data;
  std::vector<uint64_t> test_data;

  insert_data = generate_testdata(insert_size, 0);
  test_data = generate_testdata(test_size, 1);

  // Disable buffering output
  setvbuf(stdout, nullptr, _IONBF, 0);
  printf("Number of insert elements: %zu\n", insert_data.size());
  printf("Number of test elements: %zu\n", test_data.size());

  printf("name\tthreads\tinsert time(ms)\tcheck random(ms)\tcheck insert(ms)"
         "\tfalse positives\tfalse positive rate(%%)\tsize(byte)\tbits per "
         "element\tnum hash functions\n");
  {
    StdSet std;
    benchmark("std::unordered_set", std, insert_data, test_data);
  }

  for (auto n : num_hash_functions)
  {
    for (auto m : multipliers)
    {
      BloomFilter<false> amd(insert_data.size() * m, n);
      benchmark("bloom filter", amd, insert_data, test_data);
    }
  }

  for (const double e : error)
  {
    char str[64];
    snprintf(str, 64, "bloom filte by error %.1e", e);
    BloomFilter<false> amd(e, insert_data.size());
    benchmark(str, amd, insert_data, test_data);
  }

  for (auto n : num_hash_functions)
  {
    for (auto m : multipliers)
    {
      BlockedBloomFilter<false> amd(
          insert_data.size() * m / BlockedBloomFilterBlockSize, n);
      benchmark("blocked bloom filter", amd, insert_data, test_data);
    }
  }

  for (auto m : multipliers)
  {
    OneHashingBlockedBloomFilter<false> amd(
        insert_data.size() * m / OneHashingBlockedBloomFilterBlockSize);
    benchmark("one hashing blocked bloom filter", amd, insert_data, test_data);
  }

  for (auto n : num_hash_functions)
  {
    for (auto m : multipliers)
    {
      for (auto t : threads)
      {
        BloomFilter<true> amd(insert_data.size() * m, n);
        benchmark_threaded("ts bloom filter", amd, insert_data, test_data, t);
      }
    }
  }

  for (const double e : error)
  {
    for (auto t : threads)
    {
      char str[64];
      snprintf(str, 64, "ts bloom filte by error %.1e", e);
      BloomFilter<true> amd(e, insert_data.size());
      benchmark_threaded(str, amd, insert_data, test_data, t);
    }
  }

  for (auto n : num_hash_functions)
  {
    for (auto m : multipliers)
    {
      for (auto t : threads)
      {
        BlockedBloomFilter<true> amd(
            insert_data.size() * m / BlockedBloomFilterBlockSize, n);
        benchmark_threaded("ts blocked bloom filter", amd, insert_data,
                           test_data, t);
      }
    }
  }

  for (auto m : multipliers)
  {
    for (auto t : threads)
    {
      OneHashingBlockedBloomFilter<true> amd(
          insert_data.size() * m / OneHashingBlockedBloomFilterBlockSize);
      benchmark_threaded("ts one hashing blocked bloom filter", amd,
                         insert_data, test_data, t);
    }
  }
}
