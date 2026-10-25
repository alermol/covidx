#pragma once
//
// BED3+ parsing shared by covidx-build and covidx-query.
//
#include <cstdint>
#include <istream>
#include <string>
#include <vector>

namespace covidx::detail {

struct BedRecord {
  std::string contig;
  uint64_t start = 0;
  uint64_t end = 0;
  std::size_t line = 0;
};

// Parses BED3+ from a stream. Blank, '#', 'track', and 'browser' lines are
// skipped. Extra columns are ignored. Malformed lines throw InputError with
// the origin and 1-based line number.
[[nodiscard]] std::vector<BedRecord> read_bed_stream(std::istream& in,
                                                     const std::string& origin);

// Opens a BED path, or stdin when path is "-".
[[nodiscard]] std::vector<BedRecord> read_bed(const std::string& path);

} // namespace covidx::detail
