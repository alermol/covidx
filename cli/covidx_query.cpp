#include <algorithm>
#include <charconv>
#include <clocale>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

#include "common.h"
#include "covidx/cov_index.h"
#include "covidx/detail/bed.h"
#include "covidx/format.h"
#include "covidx/manifest.h"

namespace fs = std::filesystem;

namespace {

void print_help() {
  std::printf(
      "Usage:\n"
      "  covidx-query --index <dir> [--unit <name>] \\\n"
      "    (--bed <path|-> | --region <chrom:a-b> | --window <N> [--step <M>]) "
      "[options]\n"
      "\n"
      "Required:\n"
      "  --index <dir>                    Path to <name>.covidx\n"
      "\n"
      "Query:\n"
      "  --unit <name>                    Unit to query (required if not "
      "unambiguous)\n"
      "  --bed <path|->                   BED3+ input, 0-based half-open; '-' "
      "= stdin\n"
      "  --region <chrom:a-b>             Region query; 1-based inclusive by "
      "default\n"
      "  --one-based                      Explicit 1-based inclusive for "
      "--region\n"
      "  --zero-based                     Interpret --region as 0-based "
      "half-open\n"
      "  --window <N>                     Tile size; N > 0\n"
      "  --step <M>                       Tile step; M > 0 (default: N)\n"
      "  --contig <name>                  Restrict tiling to one contig\n"
      "  --no-header                      Omit the TSV header\n"
      "  --output <path|->                Output file or stdout (default: "
      "stdout)\n"
      "\n"
      "Loading:\n"
      "  --mmap / --no-mmap               Zero-copy load when possible "
      "(default: mmap)\n"
      "  --validate-crc / --no-validate-crc  Full CRC pass on load (default: "
      "validate)\n"
      "  --validation <none|structural|full> Structural validation level\n"
      "\n"
      "Listing:\n"
      "  --list-units                     Print unit names/kinds and exit 0\n"
      "  --list-contigs                   Print contig names/lengths and exit "
      "0\n"
      "  --version, --help                Print and exit 0\n"
      "\n"
      "Examples:\n"
      "  covidx-query --index indexes/sample.covidx --unit S1 --region "
      "chr1:1000-2000\n"
      "  covidx-query --index indexes/sample.covidx --window 10000 --step "
      "10000\n"
      "  covidx-query --index indexes/sample.covidx --bed targets.bed\n");
}

[[nodiscard]] std::string kind_name(covidx::UnitInfo::Kind kind) {
  switch (kind) {
  case covidx::UnitInfo::Kind::Sample:
    return "sample";
  case covidx::UnitInfo::Kind::Group:
    return "group";
  case covidx::UnitInfo::Kind::Pooled:
    return "pooled";
  }
  return "sample";
}

[[nodiscard]] bool window_inside_scope(const std::vector<covidx::Region>& scope,
                                       uint64_t x, uint64_t y) {
  if (x >= y)
    return false;
  auto it = std::upper_bound(scope.begin(), scope.end(), x,
                             [](uint64_t value, const covidx::Region& region) {
                               return value < region.start;
                             });
  if (it == scope.begin())
    return false;
  --it;
  return it->start <= x && y <= it->end;
}

[[nodiscard]] uint64_t parse_coordinate(const std::string& text,
                                        const std::string& what) {
  if (text.empty())
    throw covidx::UsageError("empty " + what + " in --region");
  uint64_t value = 0;
  const auto result =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    throw covidx::UsageError("invalid " + what + " in --region: '" + text +
                             "'");
  }
  return value;
}

struct Window {
  std::string contig;
  uint64_t start = 0;
  uint64_t end = 0;
};

class Output {
public:
  explicit Output(const std::string& path) {
    if (path == "-" || path.empty()) {
      stream_ = &std::cout;
    } else {
      file_ = std::make_unique<std::ofstream>(path, std::ios::binary |
                                                        std::ios::trunc);
      if (!*file_)
        throw covidx::OutputError("cannot open output file '" + path + "'");
      stream_ = file_.get();
    }
  }

  ~Output() = default;

  std::ostream& stream() {
    return *stream_;
  }

  void finish() {
    stream_->flush();
    if (!*stream_)
      throw covidx::OutputError("write error on query output");
  }

private:
  std::unique_ptr<std::ofstream> file_;
  std::ostream* stream_ = nullptr;
};

[[nodiscard]] std::string format_fraction(double value) {
  char buffer[64];
  const int written = std::snprintf(buffer, sizeof(buffer), "%.6f", value);
  if (written < 0)
    throw covidx::OutputError("failed to format covered_fraction");
  return std::string(buffer, static_cast<std::size_t>(written));
}

void write_row(std::ostream& out, const covidx::CovIndex& index,
               const std::string& unit, const Window& window) {
  const uint64_t covered =
      index.covered_count(unit, window.contig, window.start, window.end);
  const uint64_t length = window.end - window.start;
  const uint64_t uncovered = length - covered;
  const double fraction =
      static_cast<double>(covered) / static_cast<double>(length);
  out << window.contig << '\t' << window.start << '\t' << window.end << '\t'
      << unit << '\t' << covered << '\t' << uncovered << '\t'
      << format_fraction(fraction) << '\n';
  if (!out)
    throw covidx::OutputError("write error on query output");
}

[[nodiscard]] std::vector<Window> tile(uint64_t start, uint64_t end,
                                       uint64_t window, uint64_t step) {
  std::vector<Window> windows;
  if (window == 0 || step == 0)
    return windows;
  uint64_t position = start;
  while (position < end) {
    const uint64_t remaining = end - position;
    // Both additions are guarded by the remaining interval length so a
    // user-supplied window or step near UINT64_MAX cannot wrap.
    const uint64_t stop = window >= remaining ? end : position + window;
    windows.push_back(Window{"", position, stop});
    if (step >= remaining)
      break;
    position += step;
  }
  return windows;
}

struct ParsedRegion {
  std::string contig;
  uint64_t start = 0;
  uint64_t end = 0;
};

[[nodiscard]] ParsedRegion parse_region(const covidx::CovIndex& index,
                                        const std::string& text,
                                        bool zero_based) {
  const std::size_t colon = text.rfind(':');
  auto is_known_contig = [&](const std::string& name) {
    return index.contig(name) != nullptr;
  };
  if (colon == std::string::npos) {
    const covidx::ContigInfo* contig = index.contig(text);
    if (contig == nullptr)
      throw covidx::QueryError("unknown contig '" + text + "'");
    return ParsedRegion{text, 0, contig->length};
  }
  const std::string chrom = text.substr(0, colon);
  const std::string suffix = text.substr(colon + 1);
  const std::size_t dash = suffix.find('-');
  bool coordinate_like =
      dash != std::string::npos && dash > 0 && dash + 1 < suffix.size();
  if (coordinate_like) {
    for (const char c : suffix) {
      if (c != '-' && (c < '0' || c > '9')) {
        coordinate_like = false;
        break;
      }
    }
  }
  if (!coordinate_like) {
    const covidx::ContigInfo* contig = index.contig(text);
    if (contig == nullptr) {
      throw covidx::QueryError("malformed --region '" + text +
                               "'; expected chrom or chrom:a-b");
    }
    return ParsedRegion{text, 0, contig->length};
  }
  if (is_known_contig(text) && is_known_contig(chrom)) {
    throw covidx::QueryError("ambiguous --region '" + text +
                             "'; use --contig or a BED file");
  }
  const covidx::ContigInfo* contig = index.contig(chrom);
  if (contig == nullptr)
    throw covidx::QueryError("unknown contig '" + chrom + "'");
  const uint64_t first = parse_coordinate(suffix.substr(0, dash), "start");
  const uint64_t second = parse_coordinate(suffix.substr(dash + 1), "end");
  if (zero_based) {
    if (first >= second || second > contig->length) {
      throw covidx::UsageError("invalid zero-based --region '" + text + "'");
    }
    return ParsedRegion{chrom, first, second};
  }
  if (first < 1 || second < first || second > contig->length) {
    throw covidx::UsageError("invalid 1-based --region '" + text + "'");
  }
  return ParsedRegion{chrom, first - 1, second};
}

} // namespace

int main(int argc, char** argv) {
  std::setlocale(LC_ALL, "C");
  try {
    cli::Args args(argc, argv,
                   {"index", "unit", "bed", "region", "window", "step",
                    "contig", "output", "validation"},
                   {"one-based", "zero-based", "no-header", "mmap", "no-mmap",
                    "validate-crc", "no-validate-crc", "list-units",
                    "list-contigs", "version", "help"});

    if (args.flag("help")) {
      print_help();
      return 0;
    }
    if (args.flag("version")) {
      std::printf("covidx-query %s (format %s)\n", covidx::kBuilderVersion,
                  covidx::kFormatVersion);
      return 0;
    }
    if (!args.has("index"))
      throw covidx::UsageError("--index is required");
    if (args.get("index").empty())
      throw covidx::UsageError("--index must not be empty");
    if (args.flag("one-based") && args.flag("zero-based")) {
      throw covidx::UsageError(
          "--one-based and --zero-based are mutually exclusive");
    }
    if (args.flag("mmap") && args.flag("no-mmap")) {
      throw covidx::UsageError("--mmap and --no-mmap are mutually exclusive");
    }
    if (args.flag("validate-crc") && args.flag("no-validate-crc")) {
      throw covidx::UsageError(
          "--validate-crc and --no-validate-crc are mutually exclusive");
    }
    if (args.flag("list-units") && args.flag("list-contigs")) {
      throw covidx::UsageError(
          "--list-units and --list-contigs are mutually exclusive");
    }
    const bool listing = args.flag("list-units") || args.flag("list-contigs");
    if (listing) {
      if (args.has("bed") || args.has("region") || args.has("window")) {
        throw covidx::UsageError(
            "list modes cannot be combined with a query mode");
      }
    } else {
      if (args.has("bed") && (args.has("region") || args.has("window"))) {
        throw covidx::UsageError(
            "--bed cannot be combined with --region or --window");
      }
      if (!args.has("bed") && !args.has("region") && !args.has("window")) {
        throw covidx::UsageError(
            "one of --bed, --region, or --window is required");
      }
    }
    if (args.has("step") && !args.has("window")) {
      throw covidx::UsageError("--step requires --window");
    }
    if (args.has("contig") && (args.has("bed") || args.has("region"))) {
      throw covidx::UsageError("--contig is only valid with --window");
    }

    uint64_t window_size = 0;
    uint64_t step = 0;
    if (args.has("window")) {
      window_size = cli::parse_u64(args.get("window"), "--window");
      if (window_size == 0)
        throw covidx::UsageError("--window must be positive");
      step = args.has("step") ? cli::parse_u64(args.get("step"), "--step")
                              : window_size;
      if (step == 0)
        throw covidx::UsageError("--step must be positive");
    }

    covidx::OpenOptions open;
    open.mmap = !args.flag("no-mmap");
    open.validate_crc = !args.flag("no-validate-crc");
    if (args.has("validation")) {
      const std::string value = args.get("validation");
      if (value == "none")
        open.validation = covidx::Validation::None;
      else if (value == "structural")
        open.validation = covidx::Validation::Structural;
      else if (value == "full")
        open.validation = covidx::Validation::Full;
      else
        throw covidx::UsageError("invalid --validation '" + value + "'");
    }
    open.log = [](covidx::LogLevel level, std::string_view message) {
      std::fprintf(stderr, "covidx-query: %s: %.*s\n",
                   level == covidx::LogLevel::Warning ? "warning" : "info",
                   static_cast<int>(message.size()), message.data());
    };

    covidx::CovIndex index;
    index.load(args.get("index"), open);

    const std::string output_path = args.get("output", "-");
    const bool header = !args.flag("no-header");

    if (args.flag("list-units")) {
      Output output(output_path);
      std::ostream& out = output.stream();
      if (header)
        out << "# unit\tkind\n";
      for (const covidx::UnitInfo& unit : index.units()) {
        out << unit.name << '\t' << kind_name(unit.kind) << '\n';
      }
      output.finish();
      return 0;
    }
    if (args.flag("list-contigs")) {
      Output output(output_path);
      std::ostream& out = output.stream();
      if (header)
        out << "# contig\tlength\n";
      for (const covidx::ContigInfo& contig : index.contigs()) {
        out << contig.name << '\t' << contig.length << '\n';
      }
      output.finish();
      return 0;
    }

    std::string unit;
    if (args.has("unit")) {
      unit = args.get("unit");
      if (index.unit(unit) == nullptr) {
        std::string available;
        for (const covidx::UnitInfo& item : index.units()) {
          if (!available.empty())
            available += ", ";
          available += item.name;
        }
        throw covidx::QueryError("unknown unit '" + unit +
                                 "' (available: " + available + ")");
      }
    } else {
      if (index.units().empty()) {
        throw covidx::QueryError("index contains no units");
      }
      if (index.units().size() != 1) {
        std::string available;
        for (const covidx::UnitInfo& item : index.units()) {
          if (!available.empty())
            available += ", ";
          available += item.name;
        }
        throw covidx::UsageError("--unit is required; available: " + available);
      }
      unit = index.units().front().name;
    }

    std::vector<Window> windows;
    if (args.has("bed")) {
      const auto records = covidx::detail::read_bed(args.get("bed"));
      for (const auto& record : records) {
        try {
          (void)index.covered_count(unit, record.contig, record.start,
                                    record.end);
        } catch (const covidx::QueryError& error) {
          throw covidx::QueryError("BED line " + std::to_string(record.line) +
                                   ": " + error.what());
        }
        windows.push_back(Window{record.contig, record.start, record.end});
      }
    } else if (args.has("region")) {
      const ParsedRegion region =
          parse_region(index, args.get("region"), args.flag("zero-based"));
      if (index.contig(region.contig) == nullptr) {
        throw covidx::QueryError("unknown contig '" + region.contig + "'");
      }
      if (!window_inside_scope(index.scope(region.contig), region.start,
                               region.end)) {
        throw covidx::QueryError("--region is outside the indexed scope");
      }
      windows =
          tile(region.start, region.end,
               args.has("window") ? window_size : region.end - region.start,
               args.has("window") ? step : region.end - region.start);
      for (Window& window : windows)
        window.contig = region.contig;
    } else {
      if (args.has("contig")) {
        const std::string contig = args.get("contig");
        if (index.contig(contig) == nullptr)
          throw covidx::QueryError("unknown contig '" + contig + "'");
        for (const covidx::Region& interval : index.scope(contig)) {
          for (Window& window :
               tile(interval.start, interval.end, window_size, step)) {
            window.contig = contig;
            windows.push_back(window);
          }
        }
      } else {
        for (const covidx::ContigInfo& contig : index.contigs()) {
          if (contig.length == 0)
            continue;
          for (const covidx::Region& interval : index.scope(contig.name)) {
            for (Window& window :
                 tile(interval.start, interval.end, window_size, step)) {
              window.contig = contig.name;
              windows.push_back(window);
            }
          }
        }
      }
    }

    // Create/truncate the output only after all query input has been read and
    // validated. In particular, --bed FILE --output FILE must not truncate the
    // BED before it is parsed, and a bad unit/window must not leave a partial
    // file behind.
    Output output(output_path);
    std::ostream& out = output.stream();
    if (header)
      out << "#chrom\tstart\tend\tunit\tcovered\tuncovered\tcovered_fraction\n";
    for (const Window& window : windows)
      write_row(out, index, unit, window);
    output.finish();
    return 0;
  } catch (const std::exception& error) {
    std::fprintf(stderr, "covidx-query: %s: %s\n",
                 cli::exception_category(error).c_str(), error.what());
    return cli::exception_exit_code(error);
  }
}
