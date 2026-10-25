#pragma once
//
// The canonical manifest.json schema.
//
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "covidx/format.h"

namespace covidx {

struct SampleInfo {
  std::string name;
  std::vector<std::string> read_groups;
};

struct SourceInfo {
  std::string format = "bam";
  std::string header_so = "coordinate";
  std::string index_path;
  std::string path;
  uint64_t records_accepted = 0;
  uint64_t records_examined = 0;
  std::optional<std::string> sha256;
  uint64_t size = 0;
};

struct PolicyInfo {
  std::string ambiguous_read_bases = "cover";
  bool count_deletions = true;
  uint16_t effective_exclude_flags = 0x0F04;
  uint16_t exclude_flags = 0x0704;
  bool include_supplementary = false;
  uint8_t min_mapq = 0;
  std::string reference_base_identity = "ignored";
  std::string unmapped_policy = "skip";
};

struct ScopeInfo {
  std::string contig;
  uint64_t start = 0;
  uint64_t end = 0;
};

struct IndexInfo {
  uint64_t contig_count = 0;
  uint32_t meta_record_size = kMetaRecordSize;
  uint32_t pack_alignment = kPackAlignment;
  bool partial_scope = false;
  uint64_t total_runs = 0;
  uint64_t unit_count = 0;
};

class Manifest {
public:
  std::string format = "covidx";
  std::string format_version = kFormatVersion;
  std::string builder_version = kBuilderVersion;
  std::string created;
  std::string endianness = "little";
  std::string layout = "per-sample";
  SourceInfo source;
  PolicyInfo policy;
  std::vector<ContigInfo> contigs;
  std::vector<ScopeInfo> scope;
  std::vector<SampleInfo> samples;
  std::vector<UnitInfo> units;
  std::map<std::string, FileInfo> files;
  IndexInfo index;
  std::string semantics_covered = "reference positions covered by at least one "
                                  "accepted alignment span under policy";
  std::string semantics_features =
      "accepted alignments overlapping the window, each counted once over its "
      "full reference span";
  std::vector<std::string> warnings;

  [[nodiscard]] std::string to_json() const;
  [[nodiscard]] static Manifest from_json(std::string_view text,
                                          std::string_view origin);
  void save(const std::string& path) const;
  [[nodiscard]] static Manifest load(const std::string& path);
};

// ISO-8601 UTC timestamp from SOURCE_DATE_EPOCH, or the current time when the
// variable is unset. Throws UsageError for a malformed value.
[[nodiscard]] std::string created_from_env();

} // namespace covidx
