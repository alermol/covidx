#include "covidx/detail/spool.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <random>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <unistd.h>

#include "covidx/detail/build.h"
#include "covidx/detail/cigar.h"
#include "covidx/detail/run_sweep.h"
#include "covidx/detail/scope.h"
#include "covidx/detail/thread_pool.h"
#include "covidx/error.h"

namespace covidx::detail {

namespace {

namespace fs = std::filesystem;

constexpr std::size_t kSpoolBufferBytes = 1u << 16;
constexpr std::size_t kRecordBytes = 16;
constexpr std::size_t kMaxOpenSpools = 256;

void encode_u64(uint8_t* out, uint64_t value) {
  for (int i = 0; i < 8; ++i)
    out[i] = static_cast<uint8_t>(value >> (8 * i));
}

uint64_t decode_u64(const uint8_t* in) {
  uint64_t value = 0;
  for (int i = 0; i < 8; ++i)
    value |= static_cast<uint64_t>(in[i]) << (8 * i);
  return value;
}

[[nodiscard]] bool passes_pre_filter(const bam1_t& read,
                                     const ReadFilter& filter) {
  if ((read.core.flag & BAM_FUNMAP) != 0 || read.core.tid < 0 ||
      read.core.pos < 0)
    return false;
  if ((read.core.flag & filter.effective_exclude) != 0)
    return false;
  if (read.core.qual < filter.min_mapq)
    return false;
  return true;
}

struct BatchStats {
  BuildStats stats;
  std::mutex mutex;

  void add(const BuildStats& other) {
    std::lock_guard<std::mutex> lock(mutex);
    stats.records_examined += other.records_examined;
    stats.records_unmapped += other.records_unmapped;
    stats.records_flag_filtered += other.records_flag_filtered;
    stats.records_mapq_filtered += other.records_mapq_filtered;
    stats.records_no_blocks += other.records_no_blocks;
    stats.records_accepted += other.records_accepted;
  }
};

[[nodiscard]] uint64_t spool_records(const fs::path& path) {
  std::error_code ec;
  const uint64_t size = fs::file_size(path, ec);
  if (ec) {
    throw IoError("cannot stat spool '" + path.string() + "': " + ec.message());
  }
  if (size % kRecordBytes != 0) {
    throw FormatError("spool '" + path.string() + "' has a partial record");
  }
  return size / kRecordBytes;
}

} // namespace

RunSpoolWriter::RunSpoolWriter(std::string path) : path_(std::move(path)) {
  const fs::path parent = fs::path(path_).parent_path();
  if (!parent.empty()) {
    std::error_code ec;
    fs::create_directories(parent, ec);
    if (ec)
      throw OutputError("cannot create spool directory '" + parent.string() +
                        "'");
  }
  file_ = std::fopen(path_.c_str(), "wb");
  if (file_ == nullptr) {
    throw OutputError("cannot create spool '" + path_ +
                      "': " + std::strerror(errno));
  }
  try {
    buffer_.reserve(kSpoolBufferBytes + kRecordBytes);
  } catch (...) {
    std::fclose(file_);
    file_ = nullptr;
    throw;
  }
}

RunSpoolWriter::~RunSpoolWriter() {
  try {
    close();
  } catch (...) {
  }
  // close() can throw before it reaches fclose (for example while flushing a
  // failing stream); the destructor must still release the handle.
  if (file_ != nullptr) {
    std::fclose(file_);
    file_ = nullptr;
  }
}

void RunSpoolWriter::write(const Run& run) {
  const std::size_t offset = buffer_.size();
  buffer_.resize(offset + kRecordBytes);
  encode_u64(buffer_.data() + offset, run.start);
  encode_u64(buffer_.data() + offset + 8, run.end);
  ++records_;
  if (buffer_.size() >= kSpoolBufferBytes)
    flush();
}

void RunSpoolWriter::flush() {
  if (buffer_.empty())
    return;
  const std::size_t written =
      std::fwrite(buffer_.data(), 1, buffer_.size(), file_);
  if (written != buffer_.size()) {
    throw OutputError("write error on spool '" + path_ + "'");
  }
  buffer_.clear();
}

void RunSpoolWriter::close() {
  if (file_ == nullptr)
    return;
  flush();
  if (std::fflush(file_) != 0) {
    std::fclose(file_);
    file_ = nullptr;
    throw OutputError("flush error on spool '" + path_ + "'");
  }
  if (::fsync(::fileno(file_)) != 0) {
    std::fclose(file_);
    file_ = nullptr;
    throw OutputError("fsync error on spool '" + path_ + "'");
  }
  if (std::fclose(file_) != 0) {
    file_ = nullptr;
    throw OutputError("close error on spool '" + path_ + "'");
  }
  file_ = nullptr;
}

RunSpoolReader::RunSpoolReader(std::string path) : path_(std::move(path)) {
  file_ = std::fopen(path_.c_str(), "rb");
  if (file_ == nullptr) {
    if (errno == ENOENT)
      return;
    throw IoError("cannot open spool '" + path_ + "': " + std::strerror(errno));
  }
}

RunSpoolReader::~RunSpoolReader() {
  if (file_ != nullptr)
    std::fclose(file_);
}

bool RunSpoolReader::next(Run& run) {
  if (file_ == nullptr)
    return false;
  uint8_t data[kRecordBytes];
  const std::size_t got = std::fread(data, 1, kRecordBytes, file_);
  if (got == 0) {
    if (std::ferror(file_) != 0)
      throw IoError("read error on spool '" + path_ + "'");
    return false;
  }
  if (got != kRecordBytes)
    throw FormatError("truncated spool '" + path_ + "'");
  run.start = decode_u64(data);
  run.end = decode_u64(data + 8);
  if (run.start >= run.end)
    throw FormatError("invalid run in spool '" + path_ + "'");
  if (have_previous_ && run.start <= previous_.end) {
    throw FormatError("spool '" + path_ + "' is not sorted and canonical");
  }
  previous_ = run;
  have_previous_ = true;
  return true;
}

std::vector<Run> RunSpoolReader::read_all() {
  std::vector<Run> runs;
  Run run;
  while (next(run))
    runs.push_back(run);
  return runs;
}

void merge_run_spools(std::vector<std::unique_ptr<RunSpoolReader>>& readers,
                      RunSpoolWriter& writer) {
  struct Node {
    Run run;
    std::size_t index = 0;
  };
  struct Compare {
    bool operator()(const Node& a, const Node& b) const {
      if (a.run.start != b.run.start)
        return a.run.start > b.run.start;
      if (a.index != b.index)
        return a.index > b.index;
      return a.run.end > b.run.end;
    }
  };
  std::priority_queue<Node, std::vector<Node>, Compare> heap;
  for (std::size_t i = 0; i < readers.size(); ++i) {
    Run run;
    if (readers[i] != nullptr && readers[i]->next(run))
      heap.push(Node{run, i});
  }
  bool active = false;
  uint64_t cs = 0;
  uint64_t ce = 0;
  while (!heap.empty()) {
    const Node node = heap.top();
    heap.pop();
    const Run run = node.run;
    if (!active) {
      active = true;
      cs = run.start;
      ce = run.end;
    } else if (run.start <= ce) {
      ce = std::max(ce, run.end);
    } else {
      writer.write(Run{cs, ce});
      cs = run.start;
      ce = run.end;
    }
    Run next_run;
    if (readers[node.index]->next(next_run))
      heap.push(Node{next_run, node.index});
  }
  if (active)
    writer.write(Run{cs, ce});
}

BuildWorkspace::BuildWorkspace(std::string parent, bool keep) : keep_(keep) {
  fs::path parent_path;
  if (parent.empty()) {
    parent_path = fs::temp_directory_path();
  } else {
    parent_path = fs::path(std::move(parent));
  }
  std::error_code ec;
  fs::create_directories(parent_path, ec);
  if (ec)
    throw OutputError("cannot create work parent '" + parent_path.string() +
                      "'");

  static std::atomic<unsigned> counter{0};
  std::random_device random;
  bool created = false;
  for (int attempt = 0; attempt < 100; ++attempt) {
    const auto now =
        std::chrono::steady_clock::now().time_since_epoch().count();
    const std::string name = "covidx-work-" + std::to_string(::getpid()) + "-" +
                             std::to_string(counter.fetch_add(1)) + "-" +
                             std::to_string(now) + "-" +
                             std::to_string(random());
    root_ = parent_path / name;
    ec.clear();
    if (fs::create_directory(root_, ec)) {
      created = true;
      break;
    }
    if (ec && attempt == 99) {
      throw OutputError("cannot create work directory under '" +
                        parent_path.string() + "'");
    }
  }
  if (!created) {
    throw OutputError("cannot create a unique work directory under '" +
                      parent_path.string() + "'");
  }
  try {
    fs::permissions(root_, fs::perms::owner_all, fs::perm_options::replace, ec);
    fs::create_directories(data_dir(), ec);
    if (ec)
      throw OutputError("cannot create data directory in work directory");
    fs::create_directories(spool_dir(), ec);
    if (ec)
      throw OutputError("cannot create spool directory in work directory");
    root_string_ = root_.string();
  } catch (...) {
    // A throwing constructor does not run ~BuildWorkspace; remove the
    // directory this constructor already created.
    std::error_code cleanup;
    fs::remove_all(root_, cleanup);
    throw;
  }
}

BuildWorkspace::~BuildWorkspace() {
  if (keep_)
    return;
  std::error_code ec;
  fs::remove_all(root_, ec);
}

fs::path BuildWorkspace::sample_spool(uint32_t sample_id,
                                      uint32_t contig_id) const {
  char sample[32];
  char contig[32];
  std::snprintf(sample, sizeof(sample), "s%06u", sample_id);
  std::snprintf(contig, sizeof(contig), "c%06u", contig_id);
  return spool_dir() / sample / (std::string(contig) + ".bin");
}

fs::path BuildWorkspace::unit_spool(uint32_t unit_id,
                                    uint32_t contig_id) const {
  char unit[32];
  char contig[32];
  std::snprintf(unit, sizeof(unit), "u%06u", unit_id);
  std::snprintf(contig, sizeof(contig), "c%06u", contig_id);
  return spool_dir() / unit / (std::string(contig) + ".bin");
}

void BuildWorkspace::remove_spools() {
  std::error_code ec;
  fs::remove_all(spool_dir(), ec);
}

fs::path spool_path_for_unit(const BuildWorkspace& workspace,
                             const UnitPlan& unit, uint32_t contig_id) {
  if (unit.kind == UnitInfo::Kind::Sample && !unit.sample_ids.empty()) {
    return workspace.sample_spool(unit.sample_ids.front(), contig_id);
  }
  return workspace.unit_spool(unit.id, contig_id);
}

namespace {

BuildStats
build_contig_spools(const BamHeaderInfo& header, const BuildOptions& options,
                    const std::vector<UnitPlan>& units, uint32_t bam_tid,
                    uint32_t contig_id,
                    const std::vector<std::pair<uint64_t, uint64_t>>& intervals,
                    BuildWorkspace& workspace, std::size_t max_open_spools) {
  BuildStats stats;
  if (bam_tid >= header.contigs.size())
    throw UsageError("contig ID is out of range");
  const ContigRecord& contig = header.contigs[bam_tid];

  ReadFilter filter;
  filter.min_mapq = options.min_mapq;
  filter.effective_exclude = effective_exclude_flags(options);
  filter.extract =
      ExtractPolicy{options.count_deletions, options.ambiguous_read_bases};
  filter.contig_length = contig.length;

  struct BamDeleter {
    void operator()(bam1_t* record) const {
      bam_destroy1(record);
    }
  };

  if (options.layout == Layout::Pooled) {
    if (units.size() != 1)
      throw UsageError("pooled layout must have exactly one unit");
    RunSpoolWriter writer(
        workspace.unit_spool(units[0].id, contig_id).string());
    RunSweep sweep([&writer](Run run) { writer.write(run); });
    BamRegionReader reader(header.path, static_cast<int32_t>(bam_tid),
                           intervals);
    std::unique_ptr<bam1_t, BamDeleter> record(bam_init1());
    std::vector<Block> blocks;
    while (reader.next(*record)) {
      ++stats.records_examined;
      if (!passes_pre_filter(*record, filter)) {
        if ((record->core.flag & BAM_FUNMAP) != 0 || record->core.tid < 0 ||
            record->core.pos < 0) {
          ++stats.records_unmapped;
        } else if ((record->core.flag & filter.effective_exclude) != 0) {
          ++stats.records_flag_filtered;
        } else {
          ++stats.records_mapq_filtered;
        }
        continue;
      }
      if (!filter.accept(*record, blocks)) {
        ++stats.records_no_blocks;
        continue;
      }
      const std::vector<Block> clipped = clip_blocks(blocks, intervals);
      if (clipped.empty()) {
        ++stats.records_no_blocks;
        continue;
      }
      ++stats.records_accepted;
      sweep.advance(static_cast<uint64_t>(record->core.pos));
      for (const Block& block : clipped)
        sweep.add(block);
    }
    sweep.finish();
    writer.close();
    return stats;
  }

  const std::size_t sample_count = header.samples.size();
  const std::size_t batch_size =
      std::min<std::size_t>(std::max<std::size_t>(1, max_open_spools),
                            std::max<std::size_t>(1, sample_count));
  const std::size_t batches =
      sample_count == 0 ? 1 : (sample_count + batch_size - 1) / batch_size;
  if (batches > 1 && options.verbose && options.log) {
    options.log(LogLevel::Info, "contig '" + contig.name + "' uses " +
                                    std::to_string(batches) +
                                    " sample batches");
  }

  for (std::size_t batch = 0; batch < batches; ++batch) {
    const std::size_t begin = sample_count == 0 ? 0 : batch * batch_size;
    const std::size_t end =
        sample_count == 0 ? 0 : std::min(sample_count, begin + batch_size);
    std::vector<std::unique_ptr<RunSpoolWriter>> writers;
    std::vector<std::unique_ptr<RunSweep>> sweeps;
    writers.reserve(end - begin);
    sweeps.reserve(end - begin);
    for (std::size_t sample = begin; sample < end; ++sample) {
      writers.push_back(std::make_unique<RunSpoolWriter>(
          workspace.sample_spool(static_cast<uint32_t>(sample), contig_id)
              .string()));
      sweeps.push_back(std::make_unique<RunSweep>(
          [writer = writers.back().get()](Run run) { writer->write(run); }));
    }

    BamRegionReader reader(header.path, static_cast<int32_t>(bam_tid),
                           intervals);
    std::unique_ptr<bam1_t, BamDeleter> record(bam_init1());
    std::vector<Block> blocks;
    while (reader.next(*record)) {
      if (batch == 0) {
        ++stats.records_examined;
        if (!passes_pre_filter(*record, filter)) {
          if ((record->core.flag & BAM_FUNMAP) != 0 || record->core.tid < 0 ||
              record->core.pos < 0) {
            ++stats.records_unmapped;
          } else if ((record->core.flag & filter.effective_exclude) != 0) {
            ++stats.records_flag_filtered;
          } else {
            ++stats.records_mapq_filtered;
          }
          continue;
        }
      } else {
        if (!passes_pre_filter(*record, filter))
          continue;
      }

      if (!filter.accept(*record, blocks)) {
        if (batch == 0)
          ++stats.records_no_blocks;
        continue;
      }
      const uint32_t sample =
          sample_id_for_read(*record, header, options.unknown_sample);
      if (sample < begin || sample >= end)
        continue;
      const std::vector<Block> clipped = clip_blocks(blocks, intervals);
      if (clipped.empty()) {
        ++stats.records_no_blocks;
        continue;
      }
      ++stats.records_accepted;
      const std::size_t local = sample - begin;
      sweeps[local]->advance(static_cast<uint64_t>(record->core.pos));
      for (const Block& block : clipped)
        sweeps[local]->add(block);
    }
    for (auto& sweep : sweeps)
      sweep->finish();
    for (auto& writer : writers)
      writer->close();
  }

  return stats;
}

} // namespace

BuildStats build_run_spools(
    const BamHeaderInfo& header, const BuildOptions& options,
    const std::vector<UnitPlan>& units, const std::vector<uint32_t>& bam_tids,
    const std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& intervals,
    BuildWorkspace& workspace, std::size_t max_open_spools_override) {
  if (bam_tids.size() != intervals.size()) {
    throw UsageError("contig ID and scope interval lists disagree");
  }
  BuildStats total;
  total.units = units.size();
  total.contigs = bam_tids.size();

  unsigned threads = options.threads;
  if (threads == 0) {
    const unsigned hardware = std::thread::hardware_concurrency();
    threads = std::min<unsigned>(8, hardware == 0 ? 1 : hardware);
  }
  threads = std::max<unsigned>(1, std::min<unsigned>(threads, 1024));
  const std::size_t max_open_spools =
      max_open_spools_override != 0
          ? max_open_spools_override
          : std::max<std::size_t>(1, kMaxOpenSpools / threads);

  std::vector<std::size_t> order(bam_tids.size());
  for (std::size_t i = 0; i < order.size(); ++i)
    order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    uint64_t len_a = 0;
    uint64_t len_b = 0;
    for (const auto& interval : intervals[a])
      len_a += interval.second - interval.first;
    for (const auto& interval : intervals[b])
      len_b += interval.second - interval.first;
    if (len_a != len_b)
      return len_a > len_b;
    return a < b;
  });

  BatchStats combined;
  ThreadPool pool(threads);
  for (std::size_t task : order) {
    pool.submit([&, task] {
      const BuildStats local = build_contig_spools(
          header, options, units, bam_tids[task], static_cast<uint32_t>(task),
          intervals[task], workspace, max_open_spools);
      combined.add(local);
    });
  }
  pool.wait();

  // Merge group and pooled units from member sample spools.
  for (const UnitPlan& unit : units) {
    if (unit.kind == UnitInfo::Kind::Sample || unit.direct_from_reads)
      continue;
    for (std::size_t contig = 0; contig < bam_tids.size(); ++contig) {
      RunSpoolWriter writer(
          workspace.unit_spool(unit.id, static_cast<uint32_t>(contig))
              .string());
      std::vector<std::unique_ptr<RunSpoolReader>> readers;
      readers.reserve(unit.sample_ids.size());
      for (uint32_t sample : unit.sample_ids) {
        readers.push_back(std::make_unique<RunSpoolReader>(
            workspace.sample_spool(sample, static_cast<uint32_t>(contig))
                .string()));
      }
      merge_run_spools(readers, writer);
      writer.close();
    }
  }

  total.records_examined = combined.stats.records_examined;
  total.records_unmapped = combined.stats.records_unmapped;
  total.records_flag_filtered = combined.stats.records_flag_filtered;
  total.records_mapq_filtered = combined.stats.records_mapq_filtered;
  total.records_no_blocks = combined.stats.records_no_blocks;
  total.records_accepted = combined.stats.records_accepted;

  for (const UnitPlan& unit : units) {
    for (std::size_t contig = 0; contig < bam_tids.size(); ++contig) {
      total.runs += spool_records(
          spool_path_for_unit(workspace, unit, static_cast<uint32_t>(contig)));
    }
  }
  return total;
}

} // namespace covidx::detail
