#include "covidx/detail/merge.h"

#include <algorithm>
#include <memory>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "covidx/detail/build.h"
#include "covidx/detail/cigar.h"
#include "covidx/detail/run_sweep.h"
#include "covidx/detail/scope.h"
#include "covidx/detail/text.h"
#include "covidx/error.h"

namespace covidx::detail {

namespace {

struct HeapNode {
  Run run;
  std::size_t list = 0;
};

struct HeapCompare {
  bool operator()(const HeapNode& a, const HeapNode& b) const {
    if (a.run.start != b.run.start)
      return a.run.start > b.run.start;
    if (a.list != b.list)
      return a.list > b.list;
    return a.run.end > b.run.end;
  }
};

void append_clipped(const std::vector<Block>& blocks,
                    const std::vector<std::pair<uint64_t, uint64_t>>& intervals,
                    RunSweep& sweep) {
  const std::vector<Block> clipped = clip_blocks(blocks, intervals);
  for (const Block& block : clipped)
    sweep.add(block);
}

} // namespace

std::vector<UnitPlan> make_units(const BamHeaderInfo& header,
                                 const BuildOptions& options) {
  if (options.layout == Layout::Pooled && !options.groups.empty()) {
    throw UsageError(
        "--group requires a layout with sample demultiplexing (not pooled)");
  }

  std::unordered_map<std::string, uint32_t> sample_ids;
  for (std::size_t i = 0; i < header.samples.size(); ++i) {
    sample_ids.emplace(header.samples[i].name, static_cast<uint32_t>(i));
  }

  std::vector<UnitPlan> units;
  if (options.layout == Layout::Pooled) {
    UnitPlan pooled;
    pooled.id = 0;
    pooled.name = "pooled";
    pooled.kind = UnitInfo::Kind::Pooled;
    pooled.direct_from_reads = true;
    units.push_back(std::move(pooled));
    return units;
  }

  units.reserve(header.samples.size() + options.groups.size() + 1);
  for (std::size_t i = 0; i < header.samples.size(); ++i) {
    UnitPlan unit;
    unit.id = static_cast<uint32_t>(units.size());
    unit.name = header.samples[i].name;
    unit.kind = UnitInfo::Kind::Sample;
    unit.members = {header.samples[i].name};
    unit.sample_ids = {static_cast<uint32_t>(i)};
    units.push_back(std::move(unit));
  }

  std::set<std::string> group_names;
  std::vector<Group> groups = options.groups;
  for (const Group& group : groups) {
    validate_unit_name(group.name, "group", false);
    if (sample_ids.contains(group.name)) {
      throw UsageError("group name '" + group.name +
                       "' collides with a sample name");
    }
    if (!group_names.insert(group.name).second) {
      throw UsageError("duplicate group name '" + group.name + "'");
    }
    if (group.samples.empty()) {
      throw UsageError("group '" + group.name + "' has no members");
    }
  }
  std::sort(groups.begin(), groups.end(),
            [](const Group& a, const Group& b) { return a.name < b.name; });

  for (const Group& group : groups) {
    std::set<std::string> members;
    for (const std::string& member : group.samples) {
      if (!sample_ids.contains(member)) {
        throw UsageError("group '" + group.name +
                         "' refers to unknown sample '" + member + "'");
      }
      if (!members.insert(member).second) {
        throw UsageError("group '" + group.name + "' lists sample '" + member +
                         "' twice");
      }
    }
    UnitPlan unit;
    unit.id = static_cast<uint32_t>(units.size());
    unit.name = group.name;
    unit.kind = UnitInfo::Kind::Group;
    unit.members.assign(members.begin(), members.end());
    for (const std::string& member : unit.members) {
      unit.sample_ids.push_back(sample_ids.at(member));
    }
    units.push_back(std::move(unit));
  }

  if (options.layout == Layout::Both) {
    UnitPlan pooled;
    pooled.id = static_cast<uint32_t>(units.size());
    pooled.name = "pooled";
    pooled.kind = UnitInfo::Kind::Pooled;
    for (std::size_t i = 0; i < header.samples.size(); ++i) {
      pooled.members.push_back(header.samples[i].name);
      pooled.sample_ids.push_back(static_cast<uint32_t>(i));
    }
    units.push_back(std::move(pooled));
  }
  return units;
}

std::vector<Run>
merge_run_lists(const std::vector<const std::vector<Run>*>& lists) {
  std::priority_queue<HeapNode, std::vector<HeapNode>, HeapCompare> heap;
  std::vector<std::size_t> positions(lists.size(), 0);
  for (std::size_t i = 0; i < lists.size(); ++i) {
    if (lists[i] != nullptr && !lists[i]->empty()) {
      heap.push(HeapNode{lists[i]->front(), i});
    }
  }

  std::vector<Run> out;
  bool active = false;
  uint64_t cs = 0;
  uint64_t ce = 0;
  while (!heap.empty()) {
    const HeapNode node = heap.top();
    heap.pop();
    const Run run = node.run;
    if (!active) {
      active = true;
      cs = run.start;
      ce = run.end;
    } else if (run.start <= ce) {
      ce = std::max(ce, run.end);
    } else {
      out.push_back(Run{cs, ce});
      cs = run.start;
      ce = run.end;
    }
    ++positions[node.list];
    if (positions[node.list] < lists[node.list]->size()) {
      heap.push(HeapNode{(*lists[node.list])[positions[node.list]], node.list});
    }
  }
  if (active)
    out.push_back(Run{cs, ce});
  return out;
}

void build_contig_all_units(
    const BamHeaderInfo& header, const BuildOptions& options,
    const std::vector<UnitPlan>& units, uint32_t tid,
    const std::vector<std::pair<uint64_t, uint64_t>>& intervals,
    std::vector<std::vector<Run>>& out) {
  if (tid >= header.contigs.size())
    throw UsageError("contig ID is out of range");
  if (intervals.empty())
    throw UsageError("contig build requires a non-empty scope");
  const ContigRecord& contig = header.contigs[tid];

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
  std::unique_ptr<bam1_t, BamDeleter> record(bam_init1());
  std::vector<Block> blocks;
  BamRegionReader reader(header.path, static_cast<int32_t>(tid), intervals);

  out.assign(units.size(), {});

  const bool pooled_only = options.layout == Layout::Pooled;
  if (pooled_only) {
    RunSweep sweep([&out](Run run) { out[0].push_back(run); });
    while (reader.next(*record)) {
      if (!filter.accept(*record, blocks))
        continue;
      sweep.advance(static_cast<uint64_t>(record->core.pos));
      append_clipped(blocks, intervals, sweep);
    }
    sweep.finish();
    return;
  }

  std::vector<std::vector<Run>> sample_runs(header.samples.size());
  std::vector<std::unique_ptr<RunSweep>> sweeps;
  sweeps.reserve(header.samples.size());
  for (std::size_t i = 0; i < header.samples.size(); ++i) {
    sweeps.push_back(std::make_unique<RunSweep>(
        [&sample_runs, i](Run run) { sample_runs[i].push_back(run); }));
  }

  while (reader.next(*record)) {
    if (!filter.accept(*record, blocks))
      continue;
    const uint32_t sample =
        sample_id_for_read(*record, header, options.unknown_sample);
    if (sample >= sweeps.size())
      throw InputError("sample ID is out of range");
    sweeps[sample]->advance(static_cast<uint64_t>(record->core.pos));
    append_clipped(blocks, intervals, *sweeps[sample]);
  }
  for (auto& sweep : sweeps)
    sweep->finish();

  for (const UnitPlan& unit : units) {
    if (unit.kind == UnitInfo::Kind::Sample) {
      out[unit.id] = sample_runs.at(unit.sample_ids.at(0));
    } else {
      std::vector<const std::vector<Run>*> lists;
      lists.reserve(unit.sample_ids.size());
      for (uint32_t sample_id : unit.sample_ids) {
        if (sample_id >= sample_runs.size())
          throw InputError("unit sample ID is out of range");
        lists.push_back(&sample_runs[sample_id]);
      }
      out[unit.id] = merge_run_lists(lists);
    }
  }
}

} // namespace covidx::detail
