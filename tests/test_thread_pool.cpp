#include "doctest.h"

#include <atomic>
#include <chrono>
#include <stdexcept>
#include <thread>
#include <vector>

#include "covidx/detail/thread_pool.h"

TEST_CASE("Thread pool executes all tasks") {
  covidx::detail::ThreadPool pool(4);
  std::atomic<int> counter{0};
  for (int i = 0; i < 100; ++i) {
    pool.submit(
        [&counter] { counter.fetch_add(1, std::memory_order_relaxed); });
  }
  pool.wait();
  CHECK(counter.load() == 100);
  pool.stop();
}

TEST_CASE("Thread pool rethrows the first task exception") {
  covidx::detail::ThreadPool pool(2);
  pool.submit([] { throw std::runtime_error("boom"); });
  CHECK_THROWS_WITH_AS(pool.wait(), "boom", std::runtime_error);
}

TEST_CASE("Thread pool aborts and drops queued tasks") {
  covidx::detail::ThreadPool pool(1);
  std::atomic<int> ran{0};
  pool.submit([&ran] {
    ran.fetch_add(1, std::memory_order_relaxed);
    throw std::runtime_error("abort");
  });
  for (int i = 0; i < 50; ++i) {
    pool.submit([&ran] { ran.fetch_add(1, std::memory_order_relaxed); });
  }
  CHECK_THROWS_AS(pool.wait(), std::runtime_error);
  // The first task may have completed; the remaining queued tasks must not run.
  CHECK(ran.load() == 1);
}

TEST_CASE("Thread pool wait can be used repeatedly") {
  covidx::detail::ThreadPool pool(2);
  std::atomic<int> counter{0};
  for (int round = 0; round < 5; ++round) {
    for (int i = 0; i < 10; ++i) {
      pool.submit(
          [&counter] { counter.fetch_add(1, std::memory_order_relaxed); });
    }
    pool.wait();
  }
  CHECK(counter.load() == 50);
  pool.stop();
}
