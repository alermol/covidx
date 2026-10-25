#include "covidx/detail/scope.h"

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "covidx/detail/build.h"
#include "covidx/detail/run_sweep.h"
#include "covidx/error.h"

namespace covidx::detail {

namespace {

[[nodiscard]] std::unordered_map<std::string, uint32_t>
contig_ids(const BamHeaderInfo& header) {
  std::unordered_map<std::string, uint32_t> ids;
  ids.reserve(header.contigs.size());
  for (const ContigRecord& contig : header.contigs)
    ids.emplace(contig.name, contig.id);
  return ids;
}

} // namespace

Scope make_scope(const BamHeaderInfo& header, const BuildOptions& options) {
  if (!options.contigs.empty() && !options.regions.empty()) {
    throw UsageError(
        "--contigs/--contigs-file and --regions are mutually exclusive");
  }
  const auto ids = contig_ids(header);
  Scope scope;

  if (!options.contigs.empty()) {
    std::unordered_set<std::string> selected;
    for (const std::string& name : options.contigs) {
      if (name.empty())
        throw UsageError("--contigs contains an empty contig name");
      if (ids.find(name) == ids.end()) {
        throw UsageError("unknown contig '" + name + "' in --contigs");
      }
      selected.insert(name);
    }
    for (const ContigRecord& contig : header.contigs) {
      if (selected.contains(contig.name)) {
        scope.intervals.push_back(Region{contig.name, 0, contig.length});
      }
    }
    scope.partial = false;
    return scope;
  }

  if (!options.regions.empty()) {
    std::vector<Region> input = options.regions;
    for (const Region& region : input) {
      if (region.contig.empty())
        throw InputError("scope region has an empty contig name");
      const auto it = ids.find(region.contig);
      if (it == ids.end()) {
        throw InputError("unknown contig '" + region.contig + "' in --regions");
      }
      const uint64_t length = header.contigs[it->second].length;
      if (region.start >= region.end) {
        throw InputError("zero-length or reversed --regions interval on '" +
                         region.contig + "': [" + std::to_string(region.start) +
                         ", " + std::to_string(region.end) + ")");
      }
      if (region.end > length) {
        throw InputError("--regions interval on '" + region.contig +
                         "' ends at " + std::to_string(region.end) +
                         " beyond contig length " + std::to_string(length));
      }
    }
    std::sort(input.begin(), input.end(),
              [&ids](const Region& a, const Region& b) {
                const uint32_t ta = ids.at(a.contig);
                const uint32_t tb = ids.at(b.contig);
                if (ta != tb)
                  return ta < tb;
                if (a.start != b.start)
                  return a.start < b.start;
                return a.end < b.end;
              });

    for (const Region& region : input) {
      if (!scope.intervals.empty() &&
          scope.intervals.back().contig == region.contig &&
          region.start <= scope.intervals.back().end) {
        scope.intervals.back().end =
            std::max(scope.intervals.back().end, region.end);
      } else {
        scope.intervals.push_back(region);
      }
    }

    std::unordered_map<std::string, std::vector<std::pair<uint64_t, uint64_t>>>
        by_contig;
    for (const Region& region : scope.intervals) {
      by_contig[region.contig].push_back({region.start, region.end});
    }
    for (const ContigRecord& contig : header.contigs) {
      const auto it = by_contig.find(contig.name);
      if (it == by_contig.end())
        continue;
      const auto& intervals = it->second;
      if (intervals.size() != 1 || intervals.front().first != 0 ||
          intervals.front().second != contig.length) {
        scope.partial = true;
        break;
      }
    }
    return scope;
  }

  scope.intervals.reserve(header.contigs.size());
  for (const ContigRecord& contig : header.contigs) {
    scope.intervals.push_back(Region{contig.name, 0, contig.length});
  }
  scope.partial = false;
  return scope;
}

std::vector<std::pair<uint64_t, uint64_t>>
scope_intervals(const Scope& scope, std::string_view contig) {
  std::vector<std::pair<uint64_t, uint64_t>> out;
  for (const Region& region : scope.intervals) {
    if (region.contig == contig)
      out.push_back({region.start, region.end});
  }
  return out;
}

bool window_in_scope(const std::vector<Region>& intervals, uint64_t x,
                     uint64_t y) {
  if (x >= y)
    return false;
  auto it = std::upper_bound(intervals.begin(), intervals.end(), x,
                             [](uint64_t value, const Region& region) {
                               return value < region.start;
                             });
  if (it == intervals.begin())
    return false;
  --it;
  return it->start <= x && y <= it->end;
}

std::vector<Block>
clip_blocks(const std::vector<Block>& blocks,
            const std::vector<std::pair<uint64_t, uint64_t>>& intervals) {
  std::vector<Block> out;
  std::size_t interval_index = 0;
  for (const Block& block : blocks) {
    while (interval_index < intervals.size() &&
           intervals[interval_index].second <= block.start) {
      ++interval_index;
    }
    std::size_t index = interval_index;
    while (index < intervals.size() && intervals[index].first < block.end) {
      const uint64_t start = std::max(block.start, intervals[index].first);
      const uint64_t end = std::min(block.end, intervals[index].second);
      if (start < end) {
        if (!out.empty() && out.back().end == start) {
          out.back().end = end;
        } else {
          out.push_back(Block{start, end});
        }
      }
      ++index;
    }
  }
  return out;
}

void build_single_unit_scoped(
    const BamHeaderInfo& header, const BuildOptions& options, uint32_t tid,
    const std::vector<std::pair<uint64_t, uint64_t>>& intervals,
    RunIndex& out) {
  if (tid >= header.contigs.size())
    throw UsageError("contig ID is out of range");
  if (intervals.empty())
    throw UsageError("scoped build requires at least one interval");
  const ContigRecord& contig = header.contigs[tid];
  ReadFilter filter;
  filter.min_mapq = options.min_mapq;
  filter.effective_exclude = effective_exclude_flags(options);
  filter.extract =
      ExtractPolicy{options.count_deletions, options.ambiguous_read_bases};
  filter.contig_length = contig.length;

  std::vector<Run> runs;
  RunSweep sweep([&runs](Run run) { runs.push_back(run); });
  BamRegionReader reader(header.path, static_cast<int32_t>(tid), intervals);
  struct BamDeleter {
    void operator()(bam1_t* record) const {
      bam_destroy1(record);
    }
  };
  std::unique_ptr<bam1_t, BamDeleter> record(bam_init1());
  std::vector<Block> blocks;
  while (reader.next(*record)) {
    if (!filter.accept(*record, blocks))
      continue;
    sweep.advance(static_cast<uint64_t>(record->core.pos));
    const std::vector<Block> clipped = clip_blocks(blocks, intervals);
    for (const Block& block : clipped)
      sweep.add(block);
  }
  sweep.finish();

  out.reset(contig.length);
  out.set_runs(std::move(runs));
}

} // namespace covidx::detail
