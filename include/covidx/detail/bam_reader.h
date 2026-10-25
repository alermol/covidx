#pragma once
//
// Validated BAM header/sample parsing and per-task region iteration.
//
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/sam.h>

namespace covidx::detail {

struct ContigRecord {
  uint32_t id = 0;
  std::string name;
  uint64_t length = 0;
};

struct SampleRecord {
  std::string name;
  std::vector<std::string> read_groups;
};

struct BamHeaderInfo {
  std::string path;
  std::string header_so;
  std::string index_path;
  uint64_t file_size = 0;
  int64_t file_mtime = 0;
  int64_t index_mtime = 0;
  std::vector<ContigRecord> contigs;
  std::vector<SampleRecord> samples;
  std::unordered_map<std::string, uint32_t> rg_to_sample;
  std::unordered_map<std::string, uint32_t> sample_by_name;
};

// Opens and validates a local BAM header and index. Throws InputError for
// remote paths, non-BAM input, missing/unsorted @HD, invalid @SQ/@RG, or a
// missing index when require_index is true.
[[nodiscard]] BamHeaderInfo read_bam_header(const std::string& path,
                                            bool require_index = true);

// Resolves the sample of a read from its RG aux tag. Throws InputError when
// the tag is missing/unknown and unknown_sample is empty.
[[nodiscard]] uint32_t sample_id_for_read(const bam1_t& read,
                                          const BamHeaderInfo& header,
                                          const std::string& unknown_sample);

struct SharedBam;

// Iterates all reads overlapping one contig's merged scope intervals. Each
// read is returned once even when it overlaps several intervals, and positions
// are non-decreasing. Readers are not thread-safe; use one reader per worker.
// Immutable headers/indexes are cached per worker thread, while every reader
// owns its own BAM handle and iterator.
class BamRegionReader {
public:
  BamRegionReader(std::string path, int32_t tid,
                  std::vector<std::pair<uint64_t, uint64_t>> intervals);
  ~BamRegionReader();

  BamRegionReader(const BamRegionReader&) = delete;
  BamRegionReader& operator=(const BamRegionReader&) = delete;

  // Returns false at end of input. Throws InputError on truncated/corrupt BAM
  // or non-coordinate input.
  [[nodiscard]] bool next(bam1_t& record);

private:
  void open_next_interval();
  void prune_spanners(uint64_t interval_start);
  [[nodiscard]] std::string read_key(const bam1_t& record, int64_t end) const;

  std::string path_;
  int32_t tid_ = -1;
  std::vector<std::pair<uint64_t, uint64_t>> intervals_;
  std::size_t interval_index_ = 0;
  samFile* fp_ = nullptr;
  std::shared_ptr<SharedBam> shared_;
  sam_hdr_t* hdr_ = nullptr;
  hts_idx_t* idx_ = nullptr;
  hts_itr_t* itr_ = nullptr;
  uint64_t current_interval_start_ = 0;
  uint64_t current_interval_end_ = 0;
  int64_t last_pos_ = -1;
  std::unordered_map<std::string, uint64_t> spanning_;
};

} // namespace covidx::detail
