// Stress-test BAM generator with controllable alignment parameters.
//
// Writes a coordinate-sorted BAM whose reads vary by CIGAR profile, MAPQ,
// flags, read-group assignment, unmapped tail, and spatial clustering. The
// generated file is indexed with BAI or CSI so it can drive covidx-build.
//
// Example:
//   covidx_stress_bam --out /tmp/stress.bam --contigs 4 --contig-length
//   10000000 \
//     --reads-per-contig 200000 --read-length 150 --cigar mixed --mapq varied \
//     --flags mixed --samples 8 --cluster-percent 40
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <htslib/sam.h>

namespace {

struct Options {
  std::string out;
  uint64_t contigs = 4;
  uint64_t contig_length = 10000000;
  uint64_t reads_per_contig = 200000;
  uint64_t read_length = 150;
  uint64_t samples = 1;
  uint64_t unmapped_reads = 0;
  std::string cigar = "m";
  std::string mapq = "constant";
  std::string flags = "clean";
  std::string index = "auto";
  double cluster_percent = 0.0;
  double unknown_rg_percent = 0.0;
  double no_rg_percent = 0.0;
  double n_base_percent = 0.0;
  uint64_t seed = 20261024;
};

[[noreturn]] void usage(const char* message = nullptr) {
  if (message != nullptr)
    std::fprintf(stderr, "stress_bam: %s\n", message);
  std::fprintf(
      stderr,
      "usage: covidx_stress_bam --out <file.bam> [options]\n"
      "  --contigs N                 number of contigs (default 4)\n"
      "  --contig-length N           length of every contig (default "
      "10000000)\n"
      "  --reads-per-contig N        mapped reads per contig (default 200000)\n"
      "  --read-length N             target query length (default 150)\n"
      "  --cigar m|deletions|skips|indels|softclip|hardclip|eqx|manyops|mixed\n"
      "  --mapq constant|varied|unavailable\n"
      "  --flags clean|mixed\n"
      "  --samples N                 number of @RG samples (default 1, 0 = no "
      "RG)\n"
      "  --unmapped-reads N          unmapped records appended at the end\n"
      "  --cluster-percent P         percent of reads in local hotspots\n"
      "  --unknown-rg-percent P      percent tagged with an undeclared RG\n"
      "  --no-rg-percent P           percent with no RG tag\n"
      "  --n-base-percent P          percent of read bases replaced with N\n"
      "  --index auto|bai|csi        coordinate index type (default auto)\n"
      "  --seed N                    deterministic generator seed\n");
  std::exit(2);
}

uint64_t parse_u64(const char* text) {
  char* end = nullptr;
  const unsigned long long value = std::strtoull(text, &end, 10);
  if (text[0] == '\0' || *end != '\0')
    usage("invalid integer");
  return static_cast<uint64_t>(value);
}

double parse_percent(const char* text) {
  char* end = nullptr;
  const double value = std::strtod(text, &end);
  if (text[0] == '\0' || *end != '\0' || value < 0.0 || value > 100.0)
    usage("invalid percentage");
  return value;
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    std::string name = argv[i];
    std::string value;
    if (name == "--help" || name == "-h")
      usage();
    std::string* target = nullptr;
    double* percent_target = nullptr;
    uint64_t* number_target = nullptr;
    const std::size_t equals = name.find('=');
    if (equals != std::string::npos) {
      value = name.substr(equals + 1);
      name = name.substr(0, equals);
    } else {
      if (i + 1 >= argc)
        usage("missing option value");
      value = argv[++i];
    }
    if (name == "--out")
      target = &options.out;
    else if (name == "--cigar")
      target = &options.cigar;
    else if (name == "--mapq")
      target = &options.mapq;
    else if (name == "--flags")
      target = &options.flags;
    else if (name == "--index")
      target = &options.index;
    else if (name == "--contigs")
      number_target = &options.contigs;
    else if (name == "--contig-length")
      number_target = &options.contig_length;
    else if (name == "--reads-per-contig")
      number_target = &options.reads_per_contig;
    else if (name == "--read-length")
      number_target = &options.read_length;
    else if (name == "--samples")
      number_target = &options.samples;
    else if (name == "--unmapped-reads")
      number_target = &options.unmapped_reads;
    else if (name == "--seed")
      number_target = &options.seed;
    else if (name == "--cluster-percent")
      percent_target = &options.cluster_percent;
    else if (name == "--unknown-rg-percent")
      percent_target = &options.unknown_rg_percent;
    else if (name == "--no-rg-percent")
      percent_target = &options.no_rg_percent;
    else if (name == "--n-base-percent")
      percent_target = &options.n_base_percent;
    else
      usage("unknown option");
    if (target != nullptr)
      *target = value;
    else if (percent_target != nullptr)
      *percent_target = parse_percent(value.c_str());
    else
      *number_target = parse_u64(value.c_str());
  }
  if (options.out.empty() || options.contigs == 0 ||
      options.contig_length == 0 || options.reads_per_contig == 0 ||
      options.read_length < 2) {
    usage("invalid required parameter");
  }
  if (options.contig_length > 2147483647ULL)
    usage("contig length exceeds INT32_MAX");
  const bool known_cigar =
      options.cigar == "m" || options.cigar == "deletions" ||
      options.cigar == "skips" || options.cigar == "indels" ||
      options.cigar == "softclip" || options.cigar == "hardclip" ||
      options.cigar == "eqx" || options.cigar == "manyops" ||
      options.cigar == "mixed";
  const bool known_mapq = options.mapq == "constant" ||
                          options.mapq == "varied" ||
                          options.mapq == "unavailable";
  const bool known_flags = options.flags == "clean" || options.flags == "mixed";
  const bool known_index = options.index == "auto" || options.index == "bai" ||
                           options.index == "csi";
  if (!known_cigar || !known_mapq || !known_flags || !known_index)
    usage("unknown enum value");
  return options;
}

struct CigarPlan {
  std::vector<std::pair<char, uint32_t>> ops;
  uint64_t ref_span = 0;
  uint64_t query_len = 0;
};

void add_op(CigarPlan& plan, char op, uint32_t length) {
  if (length == 0)
    return;
  if (!plan.ops.empty() && plan.ops.back().first == op) {
    plan.ops.back().second += length;
  } else {
    plan.ops.push_back({op, length});
  }
  switch (op) {
  case 'M':
  case '=':
  case 'X':
    plan.ref_span += length;
    plan.query_len += length;
    break;
  case 'D':
  case 'N':
    plan.ref_span += length;
    break;
  case 'I':
  case 'S':
    plan.query_len += length;
    break;
  default:
    break;
  }
}

CigarPlan base_plan(const std::string& profile, uint64_t read_length) {
  CigarPlan plan;
  const uint32_t length = static_cast<uint32_t>(read_length);
  if (profile == "m") {
    add_op(plan, 'M', length);
  } else if (profile == "deletions") {
    const uint32_t chunk = std::max<uint32_t>(1, length / 4);
    add_op(plan, 'M', chunk);
    add_op(plan, 'D', 2);
    add_op(plan, 'M', chunk);
    add_op(plan, 'D', 3);
    add_op(plan, 'M', chunk);
    add_op(plan, 'M', length - 3 * chunk);
  } else if (profile == "skips") {
    const uint32_t half = std::max<uint32_t>(1, length / 2);
    add_op(plan, 'M', half);
    add_op(plan, 'N', 1000);
    add_op(plan, 'M', length - half);
  } else if (profile == "indels") {
    const uint32_t chunk = std::max<uint32_t>(1, length / 3);
    add_op(plan, 'M', chunk);
    add_op(plan, 'I', 2);
    add_op(plan, 'M', chunk);
    add_op(plan, 'D', 3);
    add_op(plan, 'M', length - 2 * chunk);
  } else if (profile == "softclip") {
    const uint32_t clip = std::max<uint32_t>(1, length / 8);
    add_op(plan, 'S', clip);
    add_op(plan, 'M', length - 2 * clip);
    add_op(plan, 'S', clip);
  } else if (profile == "hardclip") {
    const uint32_t clip = std::max<uint32_t>(1, length / 8);
    add_op(plan, 'H', clip);
    add_op(plan, 'M', length);
    add_op(plan, 'H', clip);
  } else if (profile == "eqx") {
    const uint32_t half = std::max<uint32_t>(1, length / 2);
    add_op(plan, '=', half);
    add_op(plan, 'X', length - half);
  } else if (profile == "manyops") {
    // Hundreds of operations per read: stresses CIGAR decoding and the
    // per-read block buffer without creating an enormous reference span.
    const uint32_t repeats = std::max<uint32_t>(1, length / 5);
    for (uint32_t i = 0; i < repeats; ++i) {
      add_op(plan, 'M', 2);
      add_op(plan, 'I', 1);
      add_op(plan, 'D', 1);
    }
  } else {
    // Mixed: alternate a few CIGAR elements while keeping a bounded reference
    // span. The profile is deterministic so position selection and record
    // emission agree on the reference span. Every element is legal and query
    // length is derived, not guessed.
    const uint32_t chunk = std::max<uint32_t>(2, length / 6);
    add_op(plan, 'M', chunk);
    add_op(plan, 'I', 1);
    add_op(plan, 'M', chunk);
    add_op(plan, 'D', 2);
    add_op(plan, 'M', chunk);
    add_op(plan, 'N', 50);
    add_op(plan, 'M', chunk);
    add_op(plan, 'S', 3);
    add_op(plan, 'M', length - 3 * chunk);
  }
  if (plan.query_len == 0 || plan.ref_span == 0)
    usage("generated CIGAR has an empty span");
  return plan;
}

uint8_t mapq_for(const std::string& profile, uint64_t index,
                 std::mt19937_64& rng) {
  if (profile == "constant")
    return 60;
  if (profile == "varied") {
    static constexpr uint8_t kMapq[] = {0, 1, 10, 20, 30, 42, 60, 60, 60, 255};
    return kMapq[index % (sizeof(kMapq) / sizeof(kMapq[0]))];
  }
  // "unavailable": most reads map well, some carry MAPQ 255, some are low.
  const uint64_t draw = rng() % 100;
  if (draw < 10)
    return 255;
  if (draw < 20)
    return 10;
  if (draw < 25)
    return 0;
  return 60;
}

uint16_t flag_for(const std::string& profile, uint64_t index) {
  if (profile == "clean")
    return 0;
  const uint64_t bucket = index % 100;
  if (bucket < 10)
    return BAM_FDUP;
  if (bucket < 18)
    return BAM_FSECONDARY;
  if (bucket < 23)
    return BAM_FSUPPLEMENTARY;
  if (bucket < 26)
    return BAM_FQCFAIL;
  return 0;
}

std::string sequence(std::size_t length, std::mt19937_64& rng,
                     double n_base_percent) {
  static constexpr char kBase[] = "ACGT";
  std::string out(length, 'A');
  for (std::size_t i = 0; i < length; ++i) {
    const double draw = static_cast<double>(rng() % 1000000) / 10000.0;
    out[i] = draw < n_base_percent ? 'N' : kBase[rng() % 4];
  }
  return out;
}

std::string read_group(uint64_t index, const Options& options,
                       std::mt19937_64& rng) {
  if (options.samples == 0)
    return {};
  const double draw = static_cast<double>(rng() % 1000000) / 10000.0;
  if (draw < options.unknown_rg_percent)
    return "rgUnknown";
  if (draw < options.unknown_rg_percent + options.no_rg_percent)
    return {};
  return "rg" + std::to_string(index % options.samples);
}

void check(int status, const char* what) {
  if (status < 0) {
    std::fprintf(stderr, "stress_bam: htslib failure: %s\n", what);
    std::exit(1);
  }
}

} // namespace

int main(int argc, char** argv) {
  const Options options = parse_options(argc, argv);
  std::mt19937_64 rng(options.seed);

  std::string header = "@HD\tVN:1.6\tSO:coordinate\n";
  for (uint64_t contig = 0; contig < options.contigs; ++contig) {
    header += "@SQ\tSN:chr" + std::to_string(contig + 1) +
              "\tLN:" + std::to_string(options.contig_length) + "\n";
  }
  for (uint64_t sample = 0; sample < options.samples; ++sample) {
    header += "@RG\tID:rg" + std::to_string(sample) + "\tSM:S" +
              std::to_string(sample) + "\n";
  }

  samFile* out = sam_open(options.out.c_str(), "wb");
  if (out == nullptr) {
    std::fprintf(stderr, "stress_bam: cannot create %s\n", options.out.c_str());
    return 1;
  }
  sam_hdr_t* hdr = sam_hdr_parse(header.size(), header.data());
  if (hdr == nullptr) {
    std::fprintf(stderr, "stress_bam: cannot parse header\n");
    sam_close(out);
    return 1;
  }
  check(sam_hdr_write(out, hdr), "sam_hdr_write");
  bam1_t* record = bam_init1();

  const uint64_t hotspot_count = 256;
  std::vector<uint64_t> hotspot_centers(hotspot_count);
  for (uint64_t& center : hotspot_centers)
    center = 1 + rng() % (options.contig_length - 2);

  for (uint64_t contig = 0; contig < options.contigs; ++contig) {
    std::vector<std::pair<uint64_t, uint64_t>> positions; // (pos, read index)
    positions.reserve(options.reads_per_contig);
    for (uint64_t i = 0; i < options.reads_per_contig; ++i) {
      const CigarPlan plan = base_plan(options.cigar, options.read_length);
      const uint64_t max_pos = options.contig_length - plan.ref_span;
      uint64_t pos = 0;
      const double draw = static_cast<double>(rng() % 1000000) / 10000.0;
      if (draw < options.cluster_percent) {
        const uint64_t center = hotspot_centers[(i + contig) % hotspot_count];
        const uint64_t jitter = rng() % 2001;
        pos = center > jitter ? center - jitter : 0;
        pos = std::min<uint64_t>(pos, max_pos);
      } else {
        pos = rng() % (max_pos + 1);
      }
      positions.emplace_back(pos, i);
    }
    std::sort(positions.begin(), positions.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    for (const auto& [pos, i] : positions) {
      const CigarPlan plan = base_plan(options.cigar, options.read_length);
      if (pos + plan.ref_span > options.contig_length)
        continue;
      std::vector<uint32_t> cigar;
      cigar.reserve(plan.ops.size());
      for (const auto& [op, length] : plan.ops) {
        uint32_t code = 0;
        switch (op) {
        case 'M':
          code = 0;
          break;
        case 'I':
          code = 1;
          break;
        case 'D':
          code = 2;
          break;
        case 'N':
          code = 3;
          break;
        case 'S':
          code = 4;
          break;
        case 'H':
          code = 5;
          break;
        case 'P':
          code = 6;
          break;
        case '=':
          code = 7;
          break;
        case 'X':
          code = 8;
          break;
        default:
          code = 0;
          break;
        }
        cigar.push_back((length << 4) | code);
      }
      const std::string seq = sequence(static_cast<std::size_t>(plan.query_len),
                                       rng, options.n_base_percent);
      const std::string qual(seq.size(), 'I');
      const std::string name =
          "r" + std::to_string(contig) + "_" + std::to_string(i);
      const uint16_t flag = flag_for(options.flags, i);
      const uint8_t mapq = mapq_for(options.mapq, i, rng);
      check(bam_set1(record, name.size(), name.c_str(), flag,
                     static_cast<int32_t>(contig), static_cast<int64_t>(pos),
                     mapq, cigar.size(), cigar.data(), -1, -1, 0, seq.size(),
                     seq.c_str(), qual.c_str(), 0),
            "bam_set1");
      const std::string rg = read_group(i, options, rng);
      if (!rg.empty()) {
        check(bam_aux_append(record, "RG", 'Z', static_cast<int>(rg.size() + 1),
                             reinterpret_cast<const uint8_t*>(rg.c_str())),
              "bam_aux_append");
      }
      check(sam_write1(out, hdr, record), "sam_write1");
    }
  }

  // Unmapped reads belong at the end of a coordinate-sorted BAM.
  for (uint64_t i = 0; i < options.unmapped_reads; ++i) {
    const std::string seq =
        sequence(static_cast<std::size_t>(options.read_length), rng,
                 options.n_base_percent);
    const std::string qual(seq.size(), 'I');
    const std::string name = "unmapped_" + std::to_string(i);
    check(bam_set1(record, name.size(), name.c_str(), BAM_FUNMAP, -1, -1, 0, 0,
                   nullptr, -1, -1, 0, seq.size(), seq.c_str(), qual.c_str(),
                   0),
          "bam_set1 unmapped");
    check(sam_write1(out, hdr, record), "sam_write1 unmapped");
  }

  bam_destroy1(record);
  sam_hdr_destroy(hdr);
  check(sam_close(out), "sam_close");

  int min_shift = 0;
  if (options.index == "csi")
    min_shift = 14;
  else if (options.index == "auto" && options.contig_length > (1ULL << 29))
    min_shift = 14;
  if (sam_index_build3(options.out.c_str(), nullptr, min_shift, 1) < 0) {
    std::fprintf(stderr, "stress_bam: index build failed (index=%s)\n",
                 min_shift == 0 ? "bai" : "csi");
    return 1;
  }
  std::printf(
      "wrote %s: contigs=%llu length=%llu reads_per_contig=%llu cigar=%s "
      "mapq=%s flags=%s "
      "samples=%llu index=%s\n",
      options.out.c_str(), static_cast<unsigned long long>(options.contigs),
      static_cast<unsigned long long>(options.contig_length),
      static_cast<unsigned long long>(options.reads_per_contig),
      options.cigar.c_str(), options.mapq.c_str(), options.flags.c_str(),
      static_cast<unsigned long long>(options.samples),
      min_shift == 0 ? "bai" : "csi");
  return 0;
}
