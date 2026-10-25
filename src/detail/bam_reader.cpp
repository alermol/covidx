#include "covidx/detail/bam_reader.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <sys/stat.h>

#include "covidx/detail/text.h"
#include "covidx/error.h"
#include "covidx/format.h"

namespace covidx::detail {

// Immutable BAM header/index shared by successive readers on one worker
// thread. Each reader still owns its own samFile and iterator, so readers can
// coexist and interleave safely; only the read-only metadata is shared.
struct SharedBam {
  sam_hdr_t* hdr = nullptr;
  hts_idx_t* idx = nullptr;

  ~SharedBam() {
    if (idx != nullptr)
      hts_idx_destroy(idx);
    if (hdr != nullptr)
      sam_hdr_destroy(hdr);
  }
};

namespace {

namespace fs = std::filesystem;

[[nodiscard]] bool has_remote_scheme(const std::string& path) {
  return path.find("://") != std::string::npos;
}

[[nodiscard]] std::string index_candidate(const std::string& path) {
  std::vector<std::string> candidates;
  if (path.size() > 4 && path.ends_with(".bam")) {
    const std::string base = path.substr(0, path.size() - 4);
    candidates = {path + ".bai", base + ".bai", path + ".csi", base + ".csi"};
  } else {
    candidates = {path + ".bai", path + ".csi", path + ".bam.bai",
                  path + ".bam.csi"};
  }
  for (const std::string& candidate : candidates) {
    std::error_code ec;
    if (fs::exists(candidate, ec))
      return candidate;
  }
  return path + ".bai";
}

// A cache key that changes when the BAM or its index is replaced. The size and
// mtime checks are the same ones the build uses to detect a changing source.
[[nodiscard]] std::string bam_cache_key(const std::string& path) {
  struct stat bam_stat = {};
  if (::stat(path.c_str(), &bam_stat) != 0)
    throw InputError("cannot stat BAM '" + path + "'");
  std::string key = path;
  key.push_back('\x1f');
  key += std::to_string(bam_stat.st_size);
  key.push_back('\x1f');
  key += std::to_string(bam_stat.st_mtime);
  const std::string index_path = index_candidate(path);
  struct stat index_stat = {};
  if (::stat(index_path.c_str(), &index_stat) == 0) {
    key.push_back('\x1f');
    key += index_path;
    key.push_back('\x1f');
    key += std::to_string(index_stat.st_size);
    key.push_back('\x1f');
    key += std::to_string(index_stat.st_mtime);
  }
  return key;
}

using BamCache = std::unordered_map<std::string, std::shared_ptr<SharedBam>>;

// At most a few BAMs are open per worker in practice. Clearing the map drops
// only cache references; readers that still hold a shared_ptr keep it alive.
[[nodiscard]] BamCache& bam_cache() {
  thread_local BamCache cache;
  return cache;
}

[[nodiscard]] std::string take_kstring(kstring_t& ks) {
  std::string out;
  if (ks.s != nullptr) {
    out.assign(ks.s, ks.l);
    std::free(ks.s);
    ks.s = nullptr;
    ks.l = 0;
    ks.m = 0;
  }
  return out;
}

void validate_header_name(std::string_view value, const char* kind) {
  if (value.empty())
    throw InputError(std::string(kind) + " must not be empty");
  if (!is_valid_utf8(value))
    throw InputError(std::string(kind) + " is not valid UTF-8");
  if (value.find('\0') != std::string_view::npos ||
      value.find('\t') != std::string_view::npos ||
      value.find('\r') != std::string_view::npos ||
      value.find('\n') != std::string_view::npos) {
    throw InputError(std::string(kind) + " contains a control character");
  }
}

} // namespace

BamHeaderInfo read_bam_header(const std::string& path, bool require_index) {
  if (has_remote_scheme(path)) {
    throw InputError("remote BAM URLs are not supported: '" + path + "'");
  }
  struct stat st = {};
  if (::stat(path.c_str(), &st) != 0) {
    throw InputError("cannot stat BAM '" + path + "'");
  }
  if (!S_ISREG(st.st_mode)) {
    throw InputError("BAM path is not a regular file: '" + path + "'");
  }

  samFile* fp = sam_open(path.c_str(), "r");
  if (fp == nullptr)
    throw InputError("cannot open BAM '" + path + "'");
  const htsFormat* format = hts_get_format(fp);
  if (format == nullptr || format->format != bam) {
    sam_close(fp);
    throw InputError("input is not BAM (SAM/CRAM are not supported): '" + path +
                     "'");
  }
  sam_hdr_t* hdr = sam_hdr_read(fp);
  if (hdr == nullptr) {
    sam_close(fp);
    throw InputError("cannot read BAM header: '" + path + "'");
  }

  BamHeaderInfo info;
  info.path = path;
  info.file_size = static_cast<uint64_t>(st.st_size);
  info.file_mtime = static_cast<int64_t>(st.st_mtime);
  info.index_path = index_candidate(path);

  hts_idx_t* idx = nullptr;
  if (require_index) {
    idx = sam_index_load(fp, path.c_str());
    if (idx == nullptr) {
      sam_hdr_destroy(hdr);
      sam_close(fp);
      throw InputError("BAM index (.bai/.csi) is missing or unreadable for '" +
                       path + "'; run samtools index");
    }
  }

  struct stat idx_st = {};
  if (::stat(info.index_path.c_str(), &idx_st) == 0) {
    info.index_mtime = static_cast<int64_t>(idx_st.st_mtime);
  }

  try {
    kstring_t ks{0, 0, nullptr};
    const int so_rc = sam_hdr_find_tag_hd(hdr, "SO", &ks);
    if (so_rc != 0) {
      throw InputError("BAM header must contain @HD SO:coordinate: '" + path +
                       "'");
    }
    info.header_so = take_kstring(ks);
    if (info.header_so != "coordinate") {
      throw InputError("BAM is not coordinate sorted (SO=" + info.header_so +
                       "); run samtools sort");
    }

    const int nref = sam_hdr_nref(hdr);
    if (nref < 0)
      throw InputError("BAM header has an invalid reference count");
    std::unordered_set<std::string> contig_names;
    info.contigs.reserve(static_cast<std::size_t>(nref));
    for (int tid = 0; tid < nref; ++tid) {
      const char* name = sam_hdr_tid2name(hdr, tid);
      if (name == nullptr)
        throw InputError("BAM header has a null contig name");
      validate_header_name(name, "contig name");
      if (!contig_names.insert(name).second) {
        throw InputError("duplicate @SQ name '" + std::string(name) +
                         "' in BAM header");
      }
      const hts_pos_t length = sam_hdr_tid2len(hdr, tid);
      if (length <= 0) {
        throw InputError("contig '" + std::string(name) +
                         "' has non-positive length " + std::to_string(length));
      }
      if (static_cast<uint64_t>(length) > kMaxContigLength) {
        throw InputError("contig '" + std::string(name) + "' length " +
                         std::to_string(length) +
                         " exceeds BAM alignment coordinate limit INT32_MAX");
      }
      info.contigs.push_back(ContigRecord{static_cast<uint32_t>(tid), name,
                                          static_cast<uint64_t>(length)});
    }

    // Parse @RG lines from the raw header text. htslib de-duplicates its
    // lookup table, but the format contract requires rejecting duplicate IDs.
    std::unordered_map<std::string, std::vector<std::string>> sample_to_rgs;
    std::unordered_set<std::string> rg_ids;
    const char* header_text = sam_hdr_str(hdr);
    if (header_text == nullptr)
      throw InputError("BAM header text is unavailable");
    std::string_view header_view(header_text);
    std::size_t line_start = 0;
    int rg_line = 0;
    while (line_start <= header_view.size()) {
      const std::size_t line_end = header_view.find('\n', line_start);
      const std::string_view line =
          header_view.substr(line_start, line_end == std::string_view::npos
                                             ? std::string_view::npos
                                             : line_end - line_start);
      line_start = line_end == std::string_view::npos ? header_view.size() + 1
                                                      : line_end + 1;
      if (!line.starts_with("@RG\t") && line != "@RG")
        continue;
      std::string rg_id;
      std::string sample;
      std::size_t field_start = 3;
      if (line.size() > 3 && line[3] == '\t')
        ++field_start;
      for (;;) {
        const std::size_t tab = line.find('\t', field_start);
        const std::string_view field = line.substr(
            field_start, tab == std::string_view::npos ? std::string_view::npos
                                                       : tab - field_start);
        if (field.starts_with("ID:"))
          rg_id = std::string(field.substr(3));
        else if (field.starts_with("SM:"))
          sample = std::string(field.substr(3));
        if (tab == std::string_view::npos)
          break;
        field_start = tab + 1;
      }
      if (rg_id.empty()) {
        throw InputError("BAM @RG line " + std::to_string(rg_line) +
                         " has no ID");
      }
      validate_header_name(rg_id, "read-group ID");
      if (!rg_ids.insert(rg_id).second) {
        throw InputError("duplicate @RG ID '" + rg_id + "' in BAM header");
      }
      if (sample.empty()) {
        throw InputError("BAM @RG '" + rg_id + "' has no non-empty SM tag");
      }
      try {
        validate_unit_name(sample, "sample", false);
      } catch (const UsageError& e) {
        throw InputError(std::string("invalid sample name in BAM header: ") +
                         e.what());
      }
      sample_to_rgs[sample].push_back(rg_id);
      ++rg_line;
    }

    std::vector<std::string> sample_names;
    sample_names.reserve(sample_to_rgs.size());
    for (const auto& entry : sample_to_rgs)
      sample_names.push_back(entry.first);
    std::sort(sample_names.begin(), sample_names.end());
    info.samples.reserve(sample_names.size());
    for (const std::string& name : sample_names) {
      std::vector<std::string>& rgs = sample_to_rgs[name];
      std::sort(rgs.begin(), rgs.end());
      const uint32_t sample_id = static_cast<uint32_t>(info.samples.size());
      info.sample_by_name.emplace(name, sample_id);
      for (const std::string& rg : rgs)
        info.rg_to_sample.emplace(rg, sample_id);
      info.samples.push_back(SampleRecord{name, std::move(rgs)});
    }
  } catch (...) {
    if (idx != nullptr)
      hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    sam_close(fp);
    throw;
  }

  if (idx != nullptr)
    hts_idx_destroy(idx);
  sam_hdr_destroy(hdr);
  sam_close(fp);
  return info;
}

uint32_t sample_id_for_read(const bam1_t& read, const BamHeaderInfo& header,
                            const std::string& unknown_sample) {
  const uint8_t* aux = bam_aux_get(&read, "RG");
  if (aux == nullptr) {
    if (unknown_sample.empty()) {
      throw InputError(
          "read '" + std::string(bam_get_qname(&read)) +
          "' has no RG tag; use --unknown-sample NAME or build a pooled index");
    }
    const auto it = header.sample_by_name.find(unknown_sample);
    if (it == header.sample_by_name.end()) {
      throw InputError("unknown sample '" + unknown_sample +
                       "' is not declared");
    }
    return it->second;
  }
  const char* rg = bam_aux2Z(aux);
  if (rg == nullptr) {
    throw InputError("read '" + std::string(bam_get_qname(&read)) +
                     "' has a malformed RG aux tag");
  }
  const auto it = header.rg_to_sample.find(rg);
  if (it != header.rg_to_sample.end())
    return it->second;
  if (unknown_sample.empty()) {
    throw InputError("read '" + std::string(bam_get_qname(&read)) +
                     "' refers to unknown RG '" + std::string(rg) +
                     "'; use --unknown-sample NAME");
  }
  const auto sample = header.sample_by_name.find(unknown_sample);
  if (sample == header.sample_by_name.end()) {
    throw InputError("unknown sample '" + unknown_sample + "' is not declared");
  }
  return sample->second;
}

BamRegionReader::BamRegionReader(
    std::string path, int32_t tid,
    std::vector<std::pair<uint64_t, uint64_t>> intervals)
    : path_(std::move(path)), tid_(tid), intervals_(std::move(intervals)) {
  if (intervals_.empty())
    throw UsageError("region reader requires at least one interval");
  for (std::size_t i = 0; i < intervals_.size(); ++i) {
    const auto& [start, end] = intervals_[i];
    if (start >= end)
      throw UsageError("region reader interval is empty or reversed");
    if (i > 0 && intervals_[i - 1].second >= start) {
      throw UsageError("region reader intervals must be sorted and disjoint");
    }
  }

  const std::string cache_key = bam_cache_key(path_);
  BamCache& cache = bam_cache();
  const auto cached = cache.find(cache_key);
  if (cached != cache.end()) {
    shared_ = cached->second;
  } else {
    fp_ = sam_open(path_.c_str(), "r");
    if (fp_ == nullptr)
      throw InputError("cannot open BAM '" + path_ + "'");
    auto shared = std::make_shared<SharedBam>();
    shared->hdr = sam_hdr_read(fp_);
    if (shared->hdr == nullptr) {
      sam_close(fp_);
      fp_ = nullptr;
      throw InputError("cannot read BAM header: '" + path_ + "'");
    }
    shared->idx = sam_index_load(fp_, path_.c_str());
    if (shared->idx == nullptr) {
      sam_close(fp_);
      fp_ = nullptr;
      throw InputError("BAM index is missing or unreadable for '" + path_ +
                       "'");
    }
    if (cache.size() >= 8)
      cache.clear();
    cache.emplace(cache_key, shared);
    shared_ = std::move(shared);
  }
  // A cached reader did not read the header through its own fp, so open the
  // BAM only when it was loaded fresh above.
  if (fp_ == nullptr) {
    fp_ = sam_open(path_.c_str(), "r");
    if (fp_ == nullptr)
      throw InputError("cannot open BAM '" + path_ + "'");
  }
  hdr_ = shared_->hdr;
  idx_ = shared_->idx;
  if (tid_ < 0 || tid_ >= sam_hdr_nref(hdr_)) {
    sam_close(fp_);
    fp_ = nullptr;
    hdr_ = nullptr;
    idx_ = nullptr;
    throw InputError("region reader has an invalid contig ID");
  }
}

BamRegionReader::~BamRegionReader() {
  if (itr_ != nullptr)
    hts_itr_destroy(itr_);
  // hdr_/idx_ are owned by shared_, which releases the last reference when
  // the reader is destroyed.
  if (fp_ != nullptr)
    sam_close(fp_);
}

void BamRegionReader::prune_spanners(uint64_t interval_start) {
  for (auto it = spanning_.begin(); it != spanning_.end();) {
    if (it->second <= interval_start) {
      it = spanning_.erase(it);
    } else {
      ++it;
    }
  }
}

std::string BamRegionReader::read_key(const bam1_t& record, int64_t end) const {
  std::string key = bam_get_qname(&record);
  key.push_back('\x1f');
  key += std::to_string(record.core.tid);
  key.push_back('\x1f');
  key += std::to_string(record.core.pos);
  key.push_back('\x1f');
  key += std::to_string(record.core.flag);
  key.push_back('\x1f');
  key += std::to_string(record.core.l_qseq);
  key.push_back('\x1f');
  key += std::to_string(end);
  // Two alignments with the same name/position/flag/span can still have
  // different CIGARs; include the operation list so deduplication cannot
  // collapse distinct alignments.
  const uint32_t* cigar = bam_get_cigar(&record);
  key.push_back('\x1f');
  key += std::to_string(record.core.n_cigar);
  for (uint32_t i = 0; i < record.core.n_cigar; ++i) {
    key.push_back('\x1f');
    key += std::to_string(cigar[i]);
  }
  return key;
}

void BamRegionReader::open_next_interval() {
  if (itr_ != nullptr) {
    hts_itr_destroy(itr_);
    itr_ = nullptr;
  }
  if (interval_index_ >= intervals_.size())
    return;
  current_interval_start_ = intervals_[interval_index_].first;
  current_interval_end_ = intervals_[interval_index_].second;
  // A new iterator re-emits every read overlapping the interval, including
  // reads whose POS precedes the interval start. Those records restart the
  // coordinate order, so the sortedness watermark is per interval, not
  // global. Spanning reads are still deduplicated by spanning_.
  last_pos_ = -1;
  prune_spanners(current_interval_start_);
  itr_ = sam_itr_queryi(idx_, tid_,
                        static_cast<hts_pos_t>(current_interval_start_),
                        static_cast<hts_pos_t>(current_interval_end_));
  if (itr_ == nullptr) {
    throw InputError("cannot create BAM iterator for '" + path_ +
                     "' contig ID " + std::to_string(tid_));
  }
}

bool BamRegionReader::next(bam1_t& record) {
  for (;;) {
    if (itr_ == nullptr) {
      if (interval_index_ >= intervals_.size())
        return false;
      open_next_interval();
    }
    const int rc = sam_itr_next(fp_, itr_, &record);
    if (rc == -1) {
      hts_itr_destroy(itr_);
      itr_ = nullptr;
      ++interval_index_;
      continue;
    }
    if (rc < -1) {
      throw InputError("truncated or corrupt BAM while reading '" + path_ +
                       "'");
    }
    if (record.core.tid != tid_)
      continue;
    if (record.core.pos < last_pos_) {
      throw InputError("BAM is not coordinate sorted at contig ID " +
                       std::to_string(tid_) + ", position " +
                       std::to_string(record.core.pos));
    }
    last_pos_ = record.core.pos;
    const int64_t end = bam_endpos(&record);
    const std::string key = read_key(record, end);
    const auto seen = spanning_.find(key);
    if (seen != spanning_.end()) {
      if (static_cast<uint64_t>(end) <= current_interval_end_)
        spanning_.erase(seen);
      continue;
    }
    if (static_cast<uint64_t>(end) > current_interval_end_) {
      spanning_.emplace(key, static_cast<uint64_t>(end));
    }
    return true;
  }
}

} // namespace covidx::detail
