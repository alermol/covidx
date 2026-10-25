#pragma once
//
// Public format constants and read-only metadata types.
//
#include <cstdint>
#include <string>
#include <vector>

#include "covidx/error.h"

namespace covidx {

// Version of the on-disk format implemented by this library.
inline constexpr const char* kFormatVersion = "1.0.0";
// Version of the builder implementation (independent of the format version).
inline constexpr const char* kBuilderVersion = "0.1.0";

// Unit file magic, exactly eight bytes with no terminator.
inline constexpr char kUnitMagic[8] = {'C', 'O', 'V', 'I', 'D', 'X', 'U', '1'};

inline constexpr uint32_t kMetaRecordSize = 48;
inline constexpr uint32_t kPackAlignment = 8;
inline constexpr uint64_t kMaxContigLength = 2147483647ULL; // INT32_MAX

// Layout of units built by covidx-build.
enum class Layout { PerSample, Pooled, Both };

// How positions covered by M/=/X are treated when the read base is ambiguous.
enum class AmbiguousBasePolicy { Exclude, Cover };

// Fixed group of samples, built as one union unit.
struct Group {
  std::string name;
  std::vector<std::string> samples;
};

// Half-open interval [start, end) on a named contig, 0-based.
struct Region {
  std::string contig;
  uint64_t start = 0;
  uint64_t end = 0;
};

// Load-time validation level. Structural is the default.
enum class Validation { None, Structural, Full };

// Read-only contig metadata (indexed contigs only, in header order).
struct ContigInfo {
  uint32_t id = 0;
  std::string name;
  uint64_t length = 0;
};

// Read-only unit metadata.
struct UnitInfo {
  uint32_t id = 0;
  std::string name;
  enum class Kind : uint8_t {
    Sample = 0,
    Group = 1,
    Pooled = 2
  } kind = Kind::Sample;
  std::vector<std::string> members;
  std::string file; // relative POSIX path
  uint64_t file_size = 0;
  uint32_t file_crc32 = 0;
  uint64_t total_runs = 0;
  uint64_t total_covered = 0;
  uint64_t contigs_with_runs = 0;
};

// Size and CRC32 of one file listed in the manifest.
struct FileInfo {
  uint64_t size = 0;
  uint32_t crc32 = 0;
};

// Parses a layout name.
[[nodiscard]] Layout parse_layout(const std::string& value);
// Returns the canonical layout name.
[[nodiscard]] const char* layout_name(Layout layout);
// Parses an ambiguous-base policy name.
[[nodiscard]] AmbiguousBasePolicy
parse_ambiguous_base_policy(const std::string& value);
// Returns the canonical policy name.
[[nodiscard]] const char*
ambiguous_base_policy_name(AmbiguousBasePolicy policy);

} // namespace covidx
