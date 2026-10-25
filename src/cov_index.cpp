#include "covidx/cov_index.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include "covidx/detail/bam_reader.h"
#include "covidx/detail/bucket_table.h"
#include "covidx/detail/build.h"
#include "covidx/detail/crc32.h"
#include "covidx/detail/merge.h"
#include "covidx/detail/scope.h"
#include "covidx/detail/sha256.h"
#include "covidx/detail/spool.h"
#include "covidx/detail/text.h"
#include "covidx/detail/thread_pool.h"
#include "covidx/detail/unit_file.h"
#include "covidx/manifest.h"

namespace covidx {

namespace {

namespace fs = std::filesystem;
using detail::BamHeaderInfo;
using detail::BuildWorkspace;
using detail::ContigRecord;
using detail::Scope;
using detail::UnitFileReader;
using detail::UnitPlan;

[[nodiscard]] std::string unit_filename(uint32_t id) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "u%06u.bin", id);
  return buffer;
}

[[nodiscard]] std::vector<std::vector<std::string>>
read_tsv(const fs::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw IoError("cannot open TSV '" + path.string() + "'");
  std::vector<std::vector<std::string>> rows;
  std::string line;
  while (std::getline(in, line)) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
      const std::size_t tab = line.find('\t', start);
      fields.push_back(line.substr(
          start, tab == std::string::npos ? std::string::npos : tab - start));
      if (tab == std::string::npos)
        break;
      start = tab + 1;
    }
    rows.push_back(std::move(fields));
  }
  if (in.bad())
    throw IoError("read error on TSV '" + path.string() + "'");
  return rows;
}

void write_text_file(const fs::path& path, const std::string& contents) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out)
    throw OutputError("cannot create '" + path.string() + "'");
  out.write(contents.data(), static_cast<std::streamsize>(contents.size()));
  if (!out)
    throw OutputError("write error on '" + path.string() + "'");
  out.close();
  if (!out)
    throw OutputError("close error on '" + path.string() + "'");
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd >= 0) {
    (void)::fsync(fd);
    ::close(fd);
  }
}

void fsync_directory(const fs::path& path) {
  const int fd = ::open(path.c_str(), O_RDONLY | O_DIRECTORY);
  if (fd < 0)
    throw OutputError("cannot open directory for fsync: '" + path.string() +
                      "'");
  if (::fsync(fd) != 0) {
    ::close(fd);
    throw OutputError("fsync failed on directory '" + path.string() + "'");
  }
  ::close(fd);
}

void copy_or_link(const fs::path& source, const fs::path& target) {
  std::error_code ec;
  fs::create_directories(target.parent_path(), ec);
  if (ec)
    throw OutputError("cannot create directory '" +
                      target.parent_path().string() + "'");
  fs::create_hard_link(source, target, ec);
  if (!ec)
    return;
  ec.clear();
  fs::copy_file(source, target, fs::copy_options::overwrite_existing, ec);
  if (ec) {
    throw OutputError("cannot copy '" + source.string() + "' to '" +
                      target.string() + "': " + ec.message());
  }
}

[[nodiscard]] std::vector<std::string>
expected_contig_rows(const std::vector<ContigInfo>& contigs) {
  std::vector<std::string> rows;
  for (const ContigInfo& contig : contigs) {
    rows.push_back(std::to_string(contig.id) + "\t" + contig.name + "\t" +
                   std::to_string(contig.length));
  }
  return rows;
}

[[nodiscard]] std::vector<std::string>
expected_sample_rows(const std::vector<detail::SampleRecord>& samples) {
  std::vector<std::string> rows;
  for (const auto& sample : samples) {
    if (sample.read_groups.empty()) {
      rows.push_back(sample.name + "\t-");
    } else {
      for (const std::string& rg : sample.read_groups)
        rows.push_back(sample.name + "\t" + rg);
    }
  }
  return rows;
}

[[nodiscard]] std::string kind_string(UnitInfo::Kind kind) {
  switch (kind) {
  case UnitInfo::Kind::Sample:
    return "sample";
  case UnitInfo::Kind::Group:
    return "group";
  case UnitInfo::Kind::Pooled:
    return "pooled";
  }
  return "sample";
}

[[nodiscard]] std::vector<std::string>
expected_unit_rows(const std::vector<UnitInfo>& units) {
  std::vector<std::string> rows;
  for (const UnitInfo& unit : units) {
    const std::string prefix = std::to_string(unit.id) + "\t" +
                               kind_string(unit.kind) + "\t" + unit.name + "\t";
    if (unit.kind == UnitInfo::Kind::Pooled) {
      rows.push_back(prefix + "-");
    } else {
      std::vector<std::string> members = unit.members;
      std::sort(members.begin(), members.end());
      for (const std::string& member : members)
        rows.push_back(prefix + member);
    }
  }
  return rows;
}

void compare_rows(const fs::path& path,
                  const std::vector<std::string>& expected) {
  const auto actual = read_tsv(path);
  if (actual.size() != expected.size()) {
    throw FormatError("TSV '" + path.string() + "' has " +
                      std::to_string(actual.size()) + " rows, expected " +
                      std::to_string(expected.size()));
  }
  for (std::size_t i = 0; i < expected.size(); ++i) {
    std::string joined;
    for (std::size_t j = 0; j < actual[i].size(); ++j) {
      if (j != 0)
        joined.push_back('\t');
      joined += actual[i][j];
    }
    if (joined != expected[i]) {
      throw FormatError("TSV '" + path.string() + "' row " +
                        std::to_string(i + 1) + " disagrees with the manifest");
    }
  }
}

[[nodiscard]] bool safe_relative_path(const std::string& path) {
  if (path.empty() || path.front() == '/' ||
      path.find('\0') != std::string::npos)
    return false;
  const fs::path parsed(path);
  for (const auto& part : parsed) {
    if (part == "..")
      return false;
  }
  return true;
}

} // namespace

struct CovIndex::Impl {
  Manifest manifest;
  std::vector<ContigInfo> contigs;
  std::vector<UnitInfo> units;
  std::unordered_map<std::string, uint32_t> unit_by_name;
  std::unordered_map<std::string, uint32_t> contig_by_name;
  std::unordered_map<std::string, std::vector<Region>> scope_by_contig;
  std::vector<std::unique_ptr<UnitFileReader>> readers;
  std::unique_ptr<BuildWorkspace> workspace;
  BuildStats stats;
  bool built = false;
  bool loaded = false;

  void clear() {
    readers.clear();
    manifest = Manifest{};
    contigs.clear();
    units.clear();
    unit_by_name.clear();
    contig_by_name.clear();
    scope_by_contig.clear();
    workspace.reset();
    stats = BuildStats{};
    built = false;
    loaded = false;
  }

  void rebuild_maps() {
    unit_by_name.clear();
    contig_by_name.clear();
    scope_by_contig.clear();
    for (const ContigInfo& contig : contigs) {
      contig_by_name.emplace(contig.name, contig.id);
    }
    for (const UnitInfo& unit : units) {
      unit_by_name.emplace(unit.name, unit.id);
    }
    for (const ScopeInfo& interval : manifest.scope) {
      scope_by_contig[interval.contig].push_back(
          Region{interval.contig, interval.start, interval.end});
    }
  }

  [[nodiscard]] const UnitInfo& require_unit(std::string_view name) const {
    const auto it = unit_by_name.find(std::string(name));
    if (it == unit_by_name.end()) {
      std::string available;
      for (const UnitInfo& unit : units) {
        if (!available.empty())
          available += ", ";
        available += unit.name;
      }
      throw QueryError("unknown unit '" + std::string(name) +
                       "' (available: " + available + ")");
    }
    return units.at(it->second);
  }

  [[nodiscard]] const ContigInfo& require_contig(std::string_view name) const {
    const auto it = contig_by_name.find(std::string(name));
    if (it == contig_by_name.end()) {
      throw QueryError("unknown contig '" + std::string(name) + "'");
    }
    return contigs.at(it->second);
  }

  [[nodiscard]] const std::vector<Region>&
  require_scope(std::string_view contig) const {
    const auto it = scope_by_contig.find(std::string(contig));
    if (it == scope_by_contig.end())
      throw QueryError("unknown contig '" + std::string(contig) + "'");
    return it->second;
  }

  void validate_window(std::string_view unit_name, std::string_view contig_name,
                       uint64_t x, uint64_t y) const {
    const UnitInfo& unit = require_unit(unit_name);
    const ContigInfo& contig = require_contig(contig_name);
    if (x >= y)
      throw QueryError("empty or reversed window");
    if (y > contig.length)
      throw QueryError("window exceeds contig length");
    if (!detail::window_in_scope(require_scope(contig_name), x, y)) {
      throw QueryError("window is outside the indexed scope");
    }
    (void)unit;
  }

  [[nodiscard]] uint64_t count(std::string_view unit, std::string_view contig,
                               uint64_t x, uint64_t y) const {
    validate_window(unit, contig, x, y);
    const uint32_t unit_id = unit_by_name.at(std::string(unit));
    const uint32_t contig_id = contig_by_name.at(std::string(contig));
    const detail::RunIndexView view = readers.at(unit_id)->view(contig_id);
    return detail::covered_count(view, x, y);
  }
};

CovIndex::CovIndex() : impl_(std::make_unique<Impl>()) {}
CovIndex::~CovIndex() = default;
CovIndex::CovIndex(CovIndex&&) noexcept = default;
CovIndex& CovIndex::operator=(CovIndex&&) noexcept = default;

void CovIndex::build(const std::string& bam_path, const BuildOptions& options) {
  if (!impl_)
    impl_ = std::make_unique<Impl>();
  impl_->clear();
  BamHeaderInfo header = detail::read_bam_header(bam_path, true);
  if (options.require_fresh_index && header.index_mtime > 0 &&
      header.index_mtime < header.file_mtime) {
    throw InputError(
        "BAM index is older than the BAM; rebuild it with samtools index");
  }
  std::vector<std::string> warnings;
  if (header.index_mtime > 0 && header.index_mtime < header.file_mtime) {
    warnings.push_back("BAM index is older than the BAM; rebuild recommended");
  }

  if (options.unknown_sample == "pooled") {
    throw UsageError("'pooled' cannot be used as the unknown-sample name");
  }
  bool synthetic_unknown = false;
  if (!options.unknown_sample.empty() && options.layout != Layout::Pooled) {
    detail::validate_unit_name(options.unknown_sample, "sample", false);
    if (!header.sample_by_name.contains(options.unknown_sample)) {
      header.samples.push_back(
          detail::SampleRecord{options.unknown_sample, {}});
      synthetic_unknown = true;
    }
  }
  // Rebuild the sample tables in sorted name order, preserving RG mappings.
  {
    std::vector<detail::SampleRecord> samples = header.samples;
    std::sort(samples.begin(), samples.end(),
              [](const auto& a, const auto& b) { return a.name < b.name; });
    std::unordered_map<std::string, uint32_t> by_name;
    std::unordered_map<std::string, uint32_t> rg_map;
    for (std::size_t i = 0; i < samples.size(); ++i) {
      by_name.emplace(samples[i].name, static_cast<uint32_t>(i));
      for (const std::string& rg : samples[i].read_groups) {
        rg_map.emplace(rg, static_cast<uint32_t>(i));
      }
    }
    header.samples = std::move(samples);
    header.sample_by_name = std::move(by_name);
    header.rg_to_sample = std::move(rg_map);
  }
  if (synthetic_unknown) {
    warnings.push_back("reads without or with an undeclared RG were assigned "
                       "to unknown sample '" +
                       options.unknown_sample + "'");
  }
  if (options.log) {
    for (const std::string& warning : warnings)
      options.log(LogLevel::Warning, warning);
  }

  const Scope scope = detail::make_scope(header, options);
  const std::vector<UnitPlan> units = detail::make_units(header, options);

  std::vector<uint32_t> bam_tids;
  std::vector<std::vector<std::pair<uint64_t, uint64_t>>> intervals;
  std::vector<ContigInfo> contigs;
  for (const ContigRecord& record : header.contigs) {
    const auto per_contig = detail::scope_intervals(scope, record.name);
    if (per_contig.empty())
      continue;
    const uint32_t id = static_cast<uint32_t>(contigs.size());
    contigs.push_back(ContigInfo{id, record.name, record.length});
    bam_tids.push_back(record.id);
    intervals.push_back(per_contig);
  }

  std::optional<std::string> sha256;
  if (options.checksum_source)
    sha256 = detail::sha256_file_hex(bam_path);

  auto workspace =
      std::make_unique<BuildWorkspace>(options.work_dir, options.keep_work);
  if (options.keep_work && options.log) {
    options.log(LogLevel::Info,
                "keeping work directory '" + workspace->path_string() + "'");
  }
  BuildStats stats;
  try {
    stats = detail::build_run_spools(header, options, units, bam_tids,
                                     intervals, *workspace);

    std::vector<UnitInfo> unit_infos(units.size());
    std::mutex unit_mutex;
    std::vector<uint64_t> lengths;
    lengths.reserve(contigs.size());
    for (const ContigInfo& contig : contigs)
      lengths.push_back(contig.length);

    unsigned threads = options.threads;
    if (threads == 0) {
      const unsigned hardware = std::thread::hardware_concurrency();
      threads = std::min<unsigned>(8, hardware == 0 ? 1 : hardware);
    }
    threads = std::max<unsigned>(1, std::min<unsigned>(threads, 1024));
    detail::ThreadPool pool(threads);
    for (const UnitPlan& plan : units) {
      pool.submit([&, plan] {
        const std::string filename = unit_filename(plan.id);
        const fs::path path = workspace->data_dir() / filename;
        detail::UnitFileWriter writer(path.string(), plan.id,
                                      static_cast<uint32_t>(contigs.size()));
        UnitInfo info;
        info.id = plan.id;
        info.name = plan.name;
        info.kind = plan.kind;
        info.members = plan.members;
        info.file = "data/" + filename;
        for (std::size_t contig_id = 0; contig_id < contigs.size();
             ++contig_id) {
          detail::RunSpoolReader reader(
              detail::spool_path_for_unit(*workspace, plan,
                                          static_cast<uint32_t>(contig_id))
                  .string());
          const std::vector<detail::Run> runs = reader.read_all();
          const uint64_t length = contigs[contig_id].length;
          detail::UnitContigMeta meta;
          if (runs.empty()) {
            writer.append_contig(meta, {}, {}, {}, {});
            continue;
          }
          meta.runs = runs.size();
          meta.w = detail::choose_w(length, runs.size());
          std::vector<uint32_t> start;
          std::vector<uint32_t> end;
          std::vector<uint64_t> cum(runs.size() + 1, 0);
          start.reserve(runs.size());
          end.reserve(runs.size());
          for (std::size_t i = 0; i < runs.size(); ++i) {
            start.push_back(static_cast<uint32_t>(runs[i].start));
            end.push_back(static_cast<uint32_t>(runs[i].end));
            cum[i + 1] = cum[i] + (runs[i].end - runs[i].start);
          }
          const std::vector<uint32_t> first =
              detail::build_first(runs, length, meta.w);
          uint64_t max_population = 0;
          for (std::size_t b = 1; b < first.size(); ++b) {
            max_population =
                std::max<uint64_t>(max_population, first[b] - first[b - 1]);
          }
          meta.total_covered = cum.back();
          meta.max_bucket_population = max_population;
          info.total_runs += runs.size();
          info.total_covered += cum.back();
          ++info.contigs_with_runs;
          writer.append_contig(meta, start, end, cum, first);
        }
        info.file_crc32 = writer.close_and_crc();
        std::error_code size_error;
        info.file_size = fs::file_size(path, size_error);
        if (size_error) {
          throw OutputError("cannot stat unit file '" + path.string() +
                            "': " + size_error.message());
        }
        std::lock_guard<std::mutex> lock(unit_mutex);
        unit_infos[plan.id] = std::move(info);
      });
    }
    pool.wait();

    struct stat after = {};
    if (::stat(bam_path.c_str(), &after) != 0 ||
        static_cast<uint64_t>(after.st_size) != header.file_size ||
        static_cast<int64_t>(after.st_mtime) != header.file_mtime) {
      throw InputError("BAM changed during the build: '" + bam_path + "'");
    }

    // Spools have served their purpose once the final per-unit files exist.
    // Reclaim the disk space before the manifest and save() copy phases.
    if (!options.keep_work)
      workspace->remove_spools();

    // TSVs.
    const fs::path root = workspace->root();
    std::vector<std::string> contig_rows = expected_contig_rows(contigs);
    std::vector<std::string> sample_rows = expected_sample_rows(header.samples);
    std::vector<std::string> unit_rows = expected_unit_rows(unit_infos);
    const auto join_rows = [](const std::vector<std::string>& rows) {
      std::string text;
      for (const std::string& row : rows) {
        text += row;
        text.push_back('\n');
      }
      return text;
    };
    write_text_file(root / "contigs.tsv", join_rows(contig_rows));
    write_text_file(root / "samples.tsv", join_rows(sample_rows));
    write_text_file(root / "units.tsv", join_rows(unit_rows));
    detail::write_crc_sidecar(
        (root / "contigs.tsv").string() + ".crc32",
        detail::crc32_file((root / "contigs.tsv").string()));
    detail::write_crc_sidecar(
        (root / "samples.tsv").string() + ".crc32",
        detail::crc32_file((root / "samples.tsv").string()));
    detail::write_crc_sidecar(
        (root / "units.tsv").string() + ".crc32",
        detail::crc32_file((root / "units.tsv").string()));

    Manifest manifest;
    manifest.created = created_from_env();
    manifest.layout = layout_name(options.layout);
    manifest.contigs = contigs;
    for (const ContigInfo& contig : contigs) {
      for (const Region& region : scope.intervals) {
        if (region.contig == contig.name) {
          manifest.scope.push_back(
              ScopeInfo{region.contig, region.start, region.end});
        }
      }
    }
    for (const auto& sample : header.samples) {
      manifest.samples.push_back(SampleInfo{sample.name, sample.read_groups});
    }
    manifest.units = std::move(unit_infos);
    manifest.source.format = "bam";
    manifest.source.header_so = header.header_so;
    manifest.source.index_path = header.index_path;
    manifest.source.path = bam_path;
    manifest.source.records_accepted = stats.records_accepted;
    manifest.source.records_examined = stats.records_examined;
    manifest.source.sha256 = std::move(sha256);
    manifest.source.size = header.file_size;
    manifest.policy.ambiguous_read_bases =
        ambiguous_base_policy_name(options.ambiguous_read_bases);
    manifest.policy.count_deletions = options.count_deletions;
    manifest.policy.effective_exclude_flags =
        detail::effective_exclude_flags(options);
    manifest.policy.exclude_flags = options.exclude_flags;
    manifest.policy.include_supplementary = options.include_supplementary;
    manifest.policy.min_mapq = options.min_mapq;
    manifest.index.contig_count = contigs.size();
    manifest.index.unit_count = manifest.units.size();
    manifest.index.partial_scope = scope.partial;
    uint64_t total_runs = 0;
    for (const UnitInfo& unit : manifest.units)
      total_runs += unit.total_runs;
    manifest.index.total_runs = total_runs;
    manifest.warnings = warnings;
    for (const UnitInfo& unit : manifest.units) {
      manifest.files.emplace(unit.file,
                             FileInfo{unit.file_size, unit.file_crc32});
    }
    for (const char* name : {"contigs.tsv", "samples.tsv", "units.tsv"}) {
      const fs::path path = root / name;
      std::error_code size_error;
      const uint64_t size = fs::file_size(path, size_error);
      if (size_error) {
        throw OutputError("cannot stat '" + path.string() +
                          "': " + size_error.message());
      }
      manifest.files.emplace(name,
                             FileInfo{size, detail::crc32_file(path.string())});
    }
    manifest.save((root / "manifest.json").string());

    // Open readers for immediate queries.
    std::vector<std::unique_ptr<UnitFileReader>> readers;
    readers.reserve(manifest.units.size());
    for (const UnitInfo& unit : manifest.units) {
      auto reader = std::make_unique<UnitFileReader>(
          (root / unit.file).string(), lengths, true);
      if (reader->header().unit_id != unit.id) {
        throw FormatError("unit file ID disagrees with the manifest");
      }
      readers.push_back(std::move(reader));
    }

    impl_->manifest = std::move(manifest);
    impl_->contigs = std::move(contigs);
    impl_->units = impl_->manifest.units;
    impl_->readers = std::move(readers);
    impl_->workspace = std::move(workspace);
    impl_->stats = stats;
    impl_->stats.units = impl_->units.size();
    impl_->stats.contigs = impl_->contigs.size();
    impl_->stats.runs = impl_->manifest.index.total_runs;
    impl_->built = true;
    impl_->rebuild_maps();
  } catch (...) {
    if (!options.keep_work)
      workspace.reset();
    throw;
  }
}

void CovIndex::save(const std::string& out_dir, const std::string& name,
                    bool overwrite) const {
  if (!impl_ || !impl_->built || impl_->loaded || !impl_->workspace) {
    throw UsageError("save() requires a built index; call build() first");
  }
  if (!detail::valid_index_name(name)) {
    throw UsageError("invalid index name '" + name + "'");
  }
  const fs::path out_path(out_dir);
  std::error_code ec;
  fs::create_directories(out_path, ec);
  if (ec)
    throw OutputError("cannot create output directory '" + out_dir + "'");
  const fs::path final_path = out_path / (name + ".covidx");
  if (fs::is_symlink(final_path, ec)) {
    throw OutputError("refusing to replace symlink '" + final_path.string() +
                      "'");
  }
  if (fs::exists(final_path, ec) && !overwrite) {
    throw OutputError("index already exists: '" + final_path.string() + "'");
  }

  static std::atomic<unsigned> counter{0};
  std::random_device random;
  const std::string suffix = std::to_string(::getpid()) + "-" +
                             std::to_string(counter.fetch_add(1)) + "-" +
                             std::to_string(random());
  const fs::path temp_path = out_path / ("." + name + ".covidx.tmp." + suffix);
  fs::remove_all(temp_path, ec);
  if (!fs::create_directory(temp_path, ec)) {
    throw OutputError("cannot create temporary directory '" +
                      temp_path.string() + "'");
  }
  fs::permissions(temp_path, fs::perms::owner_all, fs::perm_options::replace,
                  ec);

  try {
    const fs::path root = impl_->workspace->root();
    for (const auto& [relative, info] : impl_->manifest.files) {
      (void)info;
      if (!safe_relative_path(relative))
        throw OutputError("unsafe manifest path");
      const fs::path source = root / relative;
      const fs::path target = temp_path / relative;
      copy_or_link(source, target);
      const fs::path sidecar = source.string() + ".crc32";
      if (fs::exists(sidecar)) {
        copy_or_link(sidecar, target.string() + ".crc32");
      }
    }
    copy_or_link(root / "manifest.json", temp_path / "manifest.json");

    const Manifest verify =
        Manifest::load((temp_path / "manifest.json").string());
    for (const auto& [relative, info] : verify.files) {
      const fs::path path = temp_path / relative;
      if (!fs::exists(path) || fs::file_size(path) != info.size) {
        throw OutputError("verification failed for '" + relative + "'");
      }
      if (detail::crc32_file(path.string()) != info.crc32) {
        throw OutputError("CRC verification failed for '" + relative + "'");
      }
    }
    // fsync every regular file and directory.
    for (const auto& entry : fs::recursive_directory_iterator(temp_path)) {
      if (entry.is_regular_file()) {
        const int fd = ::open(entry.path().c_str(), O_RDONLY);
        if (fd < 0)
          throw OutputError("cannot open '" + entry.path().string() +
                            "' for fsync");
        if (::fsync(fd) != 0) {
          ::close(fd);
          throw OutputError("fsync failed on '" + entry.path().string() + "'");
        }
        ::close(fd);
      } else if (entry.is_directory()) {
        fsync_directory(entry.path());
      }
    }
    fsync_directory(temp_path);

    if (fs::exists(final_path, ec)) {
      const fs::path backup = out_path / ("." + name + ".covidx.old." + suffix);
      fs::remove_all(backup, ec);
      fs::rename(final_path, backup, ec);
      if (ec)
        throw OutputError("cannot move existing index aside: " + ec.message());
      fs::rename(temp_path, final_path, ec);
      if (ec) {
        std::error_code restore;
        fs::rename(backup, final_path, restore);
        throw OutputError("cannot publish index: " + ec.message());
      }
      fsync_directory(out_path);
      fs::remove_all(backup, ec);
    } else {
      fs::rename(temp_path, final_path, ec);
      if (ec)
        throw OutputError("cannot publish index: " + ec.message());
      fsync_directory(out_path);
    }
  } catch (...) {
    std::error_code cleanup;
    fs::remove_all(temp_path, cleanup);
    throw;
  }
}

void CovIndex::load(const std::string& dir, const OpenOptions& options) {
  if (!impl_)
    impl_ = std::make_unique<Impl>();
  impl_->clear();
  const fs::path root(dir);
  Manifest manifest = Manifest::load((root / "manifest.json").string());
  if (manifest.format != "covidx")
    throw FormatError("not a covidx index directory");

  for (const auto& [relative, info] : manifest.files) {
    if (!safe_relative_path(relative)) {
      throw FormatError("unsafe path in manifest files: '" + relative + "'");
    }
    const fs::path path = root / relative;
    std::error_code ec;
    if (!fs::exists(path, ec)) {
      throw FormatError("index file is missing: '" + relative + "'");
    }
    if (fs::file_size(path) != info.size) {
      throw FormatError("index file size disagrees with the manifest: '" +
                        relative + "'");
    }
    if (options.validate_crc &&
        detail::crc32_file(path.string()) != info.crc32) {
      throw FormatError("CRC mismatch for '" + relative + "'");
    }
    const fs::path sidecar = path.string() + ".crc32";
    if (fs::exists(sidecar)) {
      if (detail::read_crc_sidecar(sidecar.string()) != info.crc32) {
        throw FormatError("CRC sidecar disagrees with the manifest for '" +
                          relative + "'");
      }
    } else if (options.log) {
      options.log(LogLevel::Warning,
                  "CRC sidecar missing for '" + relative + "'");
    }
  }

  compare_rows(root / "contigs.tsv", expected_contig_rows(manifest.contigs));
  std::vector<detail::SampleRecord> samples;
  samples.reserve(manifest.samples.size());
  for (const SampleInfo& sample : manifest.samples) {
    samples.push_back(detail::SampleRecord{sample.name, sample.read_groups});
  }
  compare_rows(root / "samples.tsv", expected_sample_rows(samples));
  compare_rows(root / "units.tsv", expected_unit_rows(manifest.units));

  std::vector<uint64_t> lengths;
  lengths.reserve(manifest.contigs.size());
  for (const ContigInfo& contig : manifest.contigs)
    lengths.push_back(contig.length);
  std::unordered_map<std::string, std::vector<Region>> scope_by_contig;
  for (const ScopeInfo& interval : manifest.scope) {
    scope_by_contig[interval.contig].push_back(
        Region{interval.contig, interval.start, interval.end});
  }
  std::vector<std::unique_ptr<UnitFileReader>> readers;
  readers.reserve(manifest.units.size());
  for (const UnitInfo& unit : manifest.units) {
    auto reader = std::make_unique<UnitFileReader>(
        (root / unit.file).string(), lengths, options.mmap, options.validation);
    if (reader->header().unit_id != unit.id) {
      throw FormatError("unit file ID disagrees with the manifest");
    }
    if (options.mmap && !reader->used_mmap() && options.log) {
      options.log(LogLevel::Warning,
                  "mmap unavailable for '" + unit.file + "'; using heap");
    }
    if (options.validate_crc && reader->file_crc32() != unit.file_crc32) {
      throw FormatError("unit file CRC disagrees with the manifest: '" +
                        unit.file + "'");
    }
    readers.push_back(std::move(reader));
  }

  // Runs must describe positions inside the declared scope. Structural
  // validation samples the same indices as the unit-file checks; Full scans
  // every run. This catches a unit file whose CRC and internal arrays are
  // self-consistent but whose coordinates escape the scope.
  if (options.validation != Validation::None) {
    for (std::size_t unit_id = 0; unit_id < readers.size(); ++unit_id) {
      const UnitInfo& unit = manifest.units[unit_id];
      uint64_t unit_runs = 0;
      uint64_t unit_covered = 0;
      uint64_t unit_contigs = 0;
      for (uint32_t contig_id = 0; contig_id < manifest.contigs.size();
           ++contig_id) {
        const detail::RunIndexView view = readers[unit_id]->view(contig_id);
        unit_runs += view.runs;
        unit_covered += view.total_covered;
        if (view.runs != 0)
          ++unit_contigs;
        if (view.runs == 0)
          continue;
        const auto scope_it =
            scope_by_contig.find(manifest.contigs[contig_id].name);
        if (scope_it == scope_by_contig.end()) {
          throw FormatError("manifest scope is missing indexed contig '" +
                            manifest.contigs[contig_id].name + "'");
        }
        const auto check_run = [&](uint64_t index) {
          const uint64_t start = view.start[static_cast<std::size_t>(index)];
          const uint64_t end = view.end[static_cast<std::size_t>(index)];
          if (!detail::window_in_scope(scope_it->second, start, end)) {
            throw FormatError("unit '" + unit.name + "' run " +
                              std::to_string(index) + " on contig '" +
                              manifest.contigs[contig_id].name +
                              "' is outside the indexed scope");
          }
        };
        if (options.validation == Validation::Full) {
          for (uint64_t i = 0; i < view.runs; ++i)
            check_run(i);
        } else {
          check_run(0);
          check_run(view.runs / 2);
          check_run(view.runs - 1);
        }
      }
      if (unit_runs != unit.total_runs || unit_covered != unit.total_covered ||
          unit_contigs != unit.contigs_with_runs) {
        throw FormatError("unit '" + unit.name +
                          "' totals disagree with its unit file and metadata");
      }
    }
  }

  impl_->manifest = std::move(manifest);
  impl_->contigs = impl_->manifest.contigs;
  impl_->units = impl_->manifest.units;
  impl_->readers = std::move(readers);
  impl_->stats = BuildStats{};
  impl_->stats.records_examined = impl_->manifest.source.records_examined;
  impl_->stats.records_accepted = impl_->manifest.source.records_accepted;
  impl_->stats.units = impl_->units.size();
  impl_->stats.contigs = impl_->contigs.size();
  impl_->stats.runs = impl_->manifest.index.total_runs;
  impl_->loaded = true;
  impl_->rebuild_maps();
}

uint64_t CovIndex::covered_count(std::string_view unit, std::string_view contig,
                                 uint64_t x, uint64_t y) const {
  if (!impl_ || (!impl_->built && !impl_->loaded)) {
    throw QueryError("index is not built or loaded");
  }
  return impl_->count(unit, contig, x, y);
}

double CovIndex::covered_fraction(std::string_view unit,
                                  std::string_view contig, uint64_t x,
                                  uint64_t y) const {
  if (!impl_ || (!impl_->built && !impl_->loaded)) {
    throw QueryError("index is not built or loaded");
  }
  const uint64_t count = impl_->count(unit, contig, x, y);
  return static_cast<double>(count) / static_cast<double>(y - x);
}

uint64_t CovIndex::total_covered(std::string_view unit,
                                 std::string_view contig) const {
  if (!impl_ || (!impl_->built && !impl_->loaded)) {
    throw QueryError("index is not built or loaded");
  }
  (void)impl_->require_unit(unit);
  (void)impl_->require_contig(contig);
  const uint32_t unit_id = impl_->unit_by_name.at(std::string(unit));
  const uint32_t contig_id = impl_->contig_by_name.at(std::string(contig));
  return impl_->readers.at(unit_id)->view(contig_id).total_covered;
}

const std::vector<UnitInfo>& CovIndex::units() const {
  if (!impl_ || (!impl_->built && !impl_->loaded)) {
    throw QueryError("index is not built or loaded");
  }
  return impl_->units;
}

const std::vector<ContigInfo>& CovIndex::contigs() const {
  if (!impl_ || (!impl_->built && !impl_->loaded)) {
    throw QueryError("index is not built or loaded");
  }
  return impl_->contigs;
}

const UnitInfo* CovIndex::unit(std::string_view name) const {
  if (!impl_ || (!impl_->built && !impl_->loaded))
    return nullptr;
  const auto it = impl_->unit_by_name.find(std::string(name));
  return it == impl_->unit_by_name.end() ? nullptr
                                         : &impl_->units.at(it->second);
}

const ContigInfo* CovIndex::contig(std::string_view name) const {
  if (!impl_ || (!impl_->built && !impl_->loaded))
    return nullptr;
  const auto it = impl_->contig_by_name.find(std::string(name));
  return it == impl_->contig_by_name.end() ? nullptr
                                           : &impl_->contigs.at(it->second);
}

const std::vector<Region>& CovIndex::scope(std::string_view contig) const {
  if (!impl_ || (!impl_->built && !impl_->loaded)) {
    throw QueryError("index is not built or loaded");
  }
  return impl_->require_scope(contig);
}

const Manifest& CovIndex::manifest() const {
  if (!impl_ || (!impl_->built && !impl_->loaded)) {
    throw QueryError("index is not built or loaded");
  }
  return impl_->manifest;
}

const BuildStats& CovIndex::build_stats() const {
  static const BuildStats empty;
  if (!impl_)
    return empty;
  return impl_->stats;
}

} // namespace covidx
