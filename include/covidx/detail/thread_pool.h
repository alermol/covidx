#pragma once
//
// Fixed-size worker pool with a FIFO task queue.
//
// Adapted from gcrank (MIT), pinned at
// 704eee87135ff221ec3c4d60068c39c48bd0bfaa.
//
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace covidx::detail {

// A pool of at least one worker. Exceptions thrown by a task are captured and
// rethrown from wait()/stop() on the thread that owns the pool. Once a task
// throws, the pool aborts: queued tasks are dropped and submit() becomes a
// no-op.
class ThreadPool {
public:
  explicit ThreadPool(unsigned num_threads);
  ~ThreadPool();

  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator=(const ThreadPool&) = delete;

  [[nodiscard]] unsigned size() const {
    return static_cast<unsigned>(workers_.size());
  }

  // Enqueues a task; ignored once the pool is stopping or aborted.
  void submit(std::function<void()> task);

  // Blocks until all submitted tasks finish (or the pool aborts) and rethrows
  // the first task exception, if any.
  void wait();

  // Signals workers to drain and exit, joins them, and rethrows a task
  // exception, if any. Called implicitly by the destructor.
  void stop();

private:
  void worker_loop();

  std::mutex mutex_;
  std::condition_variable cv_task_;
  std::condition_variable cv_done_;
  std::deque<std::function<void()>> tasks_;
  std::vector<std::thread> workers_;
  std::size_t active_ = 0;
  bool stopping_ = false;
  bool abort_ = false;
  std::exception_ptr error_;
};

} // namespace covidx::detail
