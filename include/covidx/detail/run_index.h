#pragma once
//
// Owning in-memory run index with invariant validation.
//
#include <cstdint>
#include <vector>

#include "covidx/detail/bucket_table.h"
#include "covidx/detail/run.h"

namespace covidx::detail {

class RunIndex {
public:
  RunIndex() = default;

  void reset(uint64_t length);
  void set_runs(std::vector<Run> runs);
  void validate() const;

  [[nodiscard]] RunIndexView view() const;
  [[nodiscard]] uint64_t covered_count(uint64_t x, uint64_t y) const;
  [[nodiscard]] uint64_t total_covered() const {
    return total_covered_;
  }
  [[nodiscard]] uint64_t runs() const {
    return start_.size();
  }
  [[nodiscard]] uint32_t w() const {
    return w_;
  }
  [[nodiscard]] uint64_t length() const {
    return length_;
  }
  [[nodiscard]] const std::vector<Run>& run_list() const {
    return runs_;
  }

private:
  uint64_t length_ = 0;
  uint32_t w_ = 20;
  uint64_t total_covered_ = 0;
  std::vector<Run> runs_;
  std::vector<uint32_t> start_;
  std::vector<uint32_t> end_;
  std::vector<uint64_t> cum_;
  std::vector<uint32_t> first_;
};

} // namespace covidx::detail
