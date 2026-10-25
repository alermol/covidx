#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "common.h"
#include "covidx/cov_index.h"
#include "covidx/detail/bed.h"
#include "covidx/detail/text.h"
#include "covidx/format.h"

namespace fs = std::filesystem;

namespace {

void print_help() {
  std::printf(
      "Usage:\n"
      "  covidx-build --bam <in.bam> --out <dir> [options]\n"
      "\n"
      "Required:\n"
      "  --bam <path>                         Coordinate-sorted BAM with "
      ".bai/.csi\n"
      "  --out <dir>                          Directory that receives "
      "<name>.covidx\n"
      "\n"
      "General:\n"
      "  --name <name>                        Index name (default: BAM "
      "basename)\n"
      "  --layout <per-sample|pooled|both>    Units to build (default: "
      "per-sample)\n"
      "  --group <NAME=S1,S2,...>             Fixed group; repeatable\n"
      "  --groups-file <path>                 TSV group file; repeatable\n"
      "  --unknown-sample <name>              Assign missing/unknown RG reads "
      "here\n"
      "  --min-mapq <0..255>                  MAPQ cutoff (default: 0)\n"
      "  --exclude-flags <int>                0x704 default; decimal or 0x "
      "hex\n"
      "  --include-supplementary              Include 0x800 alignments\n"
      "  --count-deletions                    Count D as covered (default)\n"
      "  --no-count-deletions                 Treat D as an uncovered gap\n"
      "  --ambiguous-read-bases <cover|exclude>  Read-base policy (default: "
      "cover)\n"
      "  --contigs <name[,name...]>           Index only these contigs, full "
      "length\n"
      "  --contigs-file <path>                One contig name per line; "
      "repeatable\n"
      "  --regions <bed>                      Index only these 0-based "
      "half-open intervals\n"
      "  --threads <1..1024>                  Worker threads (default: min(8, "
      "CPUs))\n"
      "  --work-dir <path>                    Parent for the owned work "
      "directory\n"
      "  --keep-work                          Keep spool/work files for "
      "debugging\n"
      "  --checksum-source                    Compute the BAM SHA-256 in a "
      "separate pass\n"
      "  --require-fresh-index                Error if the BAM index is older "
      "than the BAM\n"
      "  --overwrite                          Atomically replace an existing "
      "index\n"
      "  --verbose                            Progress and timing on stderr\n"
      "  --version, --help                    Print and exit 0\n"
      "\n"
      "Examples:\n"
      "  covidx-build --bam sample.bam --out indexes\n"
      "  covidx-build --bam cohort.bam --out indexes --layout both \\\n"
      "    --group trio=NA12878,NA12891,NA12892\n"
      "  covidx-build --bam cohort.bam --out indexes --contigs chr1,chr2\n"
      "  covidx-build --bam cohort.bam --out indexes --regions targets.bed\n");
}

std::vector<covidx::Group> parse_groups(const cli::Args& args) {
  std::vector<covidx::Group> groups;
  for (const std::string& text : args.all("group")) {
    const std::size_t equals = text.find('=');
    if (equals == std::string::npos || equals == 0 ||
        equals + 1 >= text.size()) {
      throw covidx::UsageError("--group expects NAME=S1,S2,...; got '" + text +
                               "'");
    }
    covidx::Group group;
    group.name = text.substr(0, equals);
    const std::vector<std::string> samples =
        cli::split_csv(text.substr(equals + 1));
    for (const std::string& sample : samples) {
      if (sample.empty())
        throw covidx::UsageError("--group contains an empty sample name");
      group.samples.push_back(sample);
    }
    groups.push_back(std::move(group));
  }
  for (const std::string& path : args.all("groups-file")) {
    const std::vector<std::string> lines = cli::read_lines(path, "groups file");
    for (const std::string& line : lines) {
      if (line.starts_with('#'))
        continue;
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
      if (fields.size() < 2) {
        throw covidx::UsageError("groups file '" + path +
                                 "' needs group<TAB>sample<TAB>... rows");
      }
      covidx::Group group;
      group.name = fields[0];
      for (std::size_t i = 1; i < fields.size(); ++i) {
        if (fields[i].empty()) {
          throw covidx::UsageError("groups file '" + path +
                                   "' has an empty sample name");
        }
        group.samples.push_back(fields[i]);
      }
      groups.push_back(std::move(group));
    }
  }
  return groups;
}

std::vector<std::string> parse_contigs(const cli::Args& args) {
  std::vector<std::string> contigs;
  if (args.has("contigs")) {
    for (const std::string& name : cli::split_csv(args.get("contigs"))) {
      if (name.empty())
        throw covidx::UsageError("--contigs contains an empty name");
      contigs.push_back(name);
    }
  }
  for (const std::string& path : args.all("contigs-file")) {
    for (const std::string& name : cli::read_lines(path, "contigs file")) {
      if (name.starts_with('#'))
        continue;
      contigs.push_back(name);
    }
  }
  if ((args.has("contigs") || args.has("contigs-file")) && contigs.empty()) {
    throw covidx::UsageError(
        "--contigs/--contigs-file produced an empty contig list");
  }
  return contigs;
}

int run(int argc, char** argv) {
  cli::Args args(
      argc, argv,
      {"bam", "out", "name", "layout", "group", "groups-file", "unknown-sample",
       "min-mapq", "exclude-flags", "ambiguous-read-bases", "contigs",
       "contigs-file", "regions", "threads", "work-dir"},
      {"include-supplementary", "count-deletions", "no-count-deletions",
       "keep-work", "checksum-source", "require-fresh-index", "overwrite",
       "verbose", "version", "help"},
      {"group", "groups-file", "contigs-file"});

  if (args.flag("help")) {
    print_help();
    return 0;
  }
  if (args.flag("version")) {
    std::printf("covidx-build %s (format %s)\n", covidx::kBuilderVersion,
                covidx::kFormatVersion);
    return 0;
  }
  if (!args.has("bam"))
    throw covidx::UsageError("--bam is required");
  if (!args.has("out"))
    throw covidx::UsageError("--out is required");

  const std::string bam = args.get("bam");
  const std::string out = args.get("out");
  if (bam.empty())
    throw covidx::UsageError("--bam must not be empty");
  if (out.empty())
    throw covidx::UsageError("--out must not be empty");
  std::string name = args.get("name");
  if (name.empty()) {
    name = fs::path(bam).filename().string();
    if (name.ends_with(".bam"))
      name.resize(name.size() - 4);
  }
  if (!covidx::detail::valid_index_name(name)) {
    throw covidx::UsageError("invalid --name '" + name +
                             "'; use [A-Za-z0-9._-]+ up to 128 bytes");
  }

  covidx::BuildOptions options;
  options.layout = covidx::parse_layout(args.get("layout", "per-sample"));
  options.groups = parse_groups(args);
  options.unknown_sample = args.get("unknown-sample");
  if (args.has("min-mapq")) {
    const uint64_t value = cli::parse_u64(args.get("min-mapq"), "--min-mapq");
    if (value > 255)
      throw covidx::UsageError("--min-mapq must be in 0..255");
    options.min_mapq = static_cast<uint8_t>(value);
  }
  if (args.has("exclude-flags")) {
    options.exclude_flags = static_cast<uint16_t>(
        cli::parse_flags(args.get("exclude-flags"), "--exclude-flags"));
  }
  options.include_supplementary = args.flag("include-supplementary");
  if (args.flag("count-deletions") && args.flag("no-count-deletions")) {
    throw covidx::UsageError(
        "--count-deletions and --no-count-deletions are mutually exclusive");
  }
  options.count_deletions = !args.flag("no-count-deletions");
  options.ambiguous_read_bases = covidx::parse_ambiguous_base_policy(
      args.get("ambiguous-read-bases", "cover"));
  options.contigs = parse_contigs(args);
  if (args.has("regions")) {
    if (!options.contigs.empty()) {
      throw covidx::UsageError(
          "--contigs/--contigs-file and --regions are mutually exclusive");
    }
    const auto records = covidx::detail::read_bed(args.get("regions"));
    if (records.empty())
      throw covidx::UsageError("--regions contains no intervals");
    for (const auto& record : records) {
      options.regions.push_back(
          covidx::Region{record.contig, record.start, record.end});
    }
  }
  if (args.has("threads")) {
    const uint64_t value = cli::parse_u64(args.get("threads"), "--threads");
    if (value < 1 || value > 1024)
      throw covidx::UsageError("--threads must be in 1..1024");
    options.threads = static_cast<unsigned>(value);
  } else {
    options.threads = cli::default_threads();
  }
  options.work_dir =
      args.get("work-dir", (fs::path(out) / ".covidx-work").string());
  options.keep_work = args.flag("keep-work");
  options.checksum_source = args.flag("checksum-source");
  options.require_fresh_index = args.flag("require-fresh-index");
  options.verbose = args.flag("verbose");

  const bool verbose = options.verbose;
  options.log = [](covidx::LogLevel level, std::string_view message) {
    std::fprintf(stderr, "covidx-build: %s: %.*s\n",
                 level == covidx::LogLevel::Warning ? "warning" : "info",
                 static_cast<int>(message.size()), message.data());
  };

  const auto started = std::chrono::steady_clock::now();
  covidx::CovIndex index;
  index.build(bam, options);
  index.save(out, name, args.flag("overwrite"));
  const auto elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started);
  if (verbose) {
    const covidx::BuildStats& stats = index.build_stats();
    std::fprintf(
        stderr,
        "covidx-build: info: contigs=%llu units=%llu runs=%llu features=%llu "
        "accepted=%llu elapsed=%.3fs out=%s\n",
        static_cast<unsigned long long>(stats.contigs),
        static_cast<unsigned long long>(stats.units),
        static_cast<unsigned long long>(stats.runs),
        static_cast<unsigned long long>(stats.features),
        static_cast<unsigned long long>(stats.records_accepted),
        elapsed.count(), (fs::path(out) / (name + ".covidx")).c_str());
  }
  return 0;
}

} // namespace

int main(int argc, char** argv) {
  try {
    return run(argc, argv);
  } catch (const std::exception& error) {
    std::fprintf(stderr, "covidx-build: %s: %s\n",
                 cli::exception_category(error).c_str(), error.what());
    return cli::exception_exit_code(error);
  }
}
