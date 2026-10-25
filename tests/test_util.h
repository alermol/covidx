#pragma once
//
// Test-only helpers: deterministic BAM construction and process helpers.
//
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

#include <htslib/sam.h>

namespace tests {

namespace fs = std::filesystem;

struct CigarOp {
  char op = 'M';
  uint32_t len = 0;
};

struct ReadSpec {
  std::string name;
  uint16_t flag = 0;
  int32_t tid = 0;
  int64_t pos = 0;
  uint8_t mapq = 60;
  std::vector<CigarOp> cigar;
  std::string seq;
  std::string rg; // empty = no RG aux tag
};

inline uint32_t cigar_code(char op) {
  switch (op) {
  case 'M':
    return 0;
  case 'I':
    return 1;
  case 'D':
    return 2;
  case 'N':
    return 3;
  case 'S':
    return 4;
  case 'H':
    return 5;
  case 'P':
    return 6;
  case '=':
    return 7;
  case 'X':
    return 8;
  case 'B':
    return 9;
  default:
    throw std::runtime_error(std::string("bad test CIGAR op: ") + op);
  }
}

inline std::vector<uint32_t> encode_cigar(const std::vector<CigarOp>& ops) {
  std::vector<uint32_t> out;
  out.reserve(ops.size());
  for (const CigarOp& op : ops)
    out.push_back((op.len << 4) | cigar_code(op.op));
  return out;
}

inline std::string repeat_sequence(std::size_t length,
                                   std::string_view unit = "ACGT") {
  std::string out;
  out.reserve(length);
  for (std::size_t i = 0; i < length; ++i)
    out.push_back(unit[i % unit.size()]);
  return out;
}

inline void check_hts(int status, const char* what) {
  if (status < 0)
    throw std::runtime_error(std::string("htslib failure: ") + what);
}

// Writes a coordinate-sorted BAM and builds its BAI index. Reads are sorted by
// (tid, pos); equal positions keep their input order.
inline void make_bam(const fs::path& path, std::string_view header_text,
                     std::vector<ReadSpec> reads) {
  std::stable_sort(reads.begin(), reads.end(),
                   [](const ReadSpec& a, const ReadSpec& b) {
                     if (a.tid != b.tid)
                       return a.tid < b.tid;
                     return a.pos < b.pos;
                   });

  samFile* out = sam_open(path.string().c_str(), "wb");
  if (out == nullptr)
    throw std::runtime_error("cannot create " + path.string());
  sam_hdr_t* hdr = sam_hdr_parse(header_text.size(), header_text.data());
  if (hdr == nullptr) {
    sam_close(out);
    throw std::runtime_error("cannot parse test header");
  }
  check_hts(sam_hdr_write(out, hdr), "sam_hdr_write");

  bam1_t* record = bam_init1();
  try {
    for (const ReadSpec& spec : reads) {
      const std::vector<uint32_t> cigar = encode_cigar(spec.cigar);
      const int32_t l_seq = static_cast<int32_t>(spec.seq.size());
      const char* seq = l_seq == 0 ? nullptr : spec.seq.c_str();
      std::string qual(static_cast<std::size_t>(l_seq), 'I');
      const char* qual_ptr = l_seq == 0 ? nullptr : qual.c_str();
      const int rc = bam_set1(
          record, spec.name.size(), spec.name.c_str(), spec.flag, spec.tid,
          spec.pos, spec.mapq, cigar.size(), cigar.data(), -1, -1, 0,
          static_cast<std::size_t>(l_seq), seq, qual_ptr, 0);
      if (rc < 0)
        throw std::runtime_error("bam_set1 failed for " + spec.name);
      if (!spec.rg.empty()) {
        const std::size_t len = spec.rg.size() + 1;
        check_hts(
            bam_aux_append(record, "RG", 'Z', static_cast<int>(len),
                           reinterpret_cast<const uint8_t*>(spec.rg.c_str())),
            "bam_aux_append");
      }
      check_hts(sam_write1(out, hdr, record), "sam_write1");
    }
  } catch (...) {
    bam_destroy1(record);
    sam_hdr_destroy(hdr);
    sam_close(out);
    throw;
  }
  bam_destroy1(record);
  sam_hdr_destroy(hdr);
  check_hts(sam_close(out), "sam_close");
  check_hts(sam_index_build(path.string().c_str(), 0), "sam_index_build");
}

inline std::string read_file(const fs::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot read " + path.string());
  return std::string(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
}

inline void write_file(const fs::path& path, std::string_view contents) {
  std::ofstream out(path, std::ios::binary);
  if (!out)
    throw std::runtime_error("cannot write " + path.string());
  out.write(contents.data(), static_cast<std::streamsize>(contents.size()));
}

// Runs a shell command and returns its exit status. Captures stdout/stderr
// into output when requested.
inline int run_command(const std::string& command,
                       std::string* output = nullptr) {
  FILE* pipe = ::popen(command.c_str(), "r");
  if (pipe == nullptr)
    return -1;
  std::string captured;
  char buffer[4096];
  while (std::fgets(buffer, sizeof(buffer), pipe) != nullptr)
    captured += buffer;
  const int status = ::pclose(pipe);
  if (output != nullptr)
    *output = std::move(captured);
  return status;
}

inline int run_command_status(const std::string& command,
                              std::string* output = nullptr) {
  const int status = run_command(command, output);
  if (status == -1)
    return -1;
  if (WIFEXITED(status))
    return WEXITSTATUS(status);
  if (WIFSIGNALED(status))
    return 128 + WTERMSIG(status);
  return -1;
}

inline bool command_exists(const char* name) {
  return run_command(std::string("command -v ") + name + " >/dev/null 2>&1") ==
         0;
}

inline std::vector<std::string> split_tabs(const std::string& line) {
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
  return fields;
}

// A unique scratch directory removed on destruction.
class TempDir {
public:
  explicit TempDir(std::string_view prefix = "covidx-test-") {
    static unsigned counter = 0;
    path_ = fs::temp_directory_path() /
            (std::string(prefix) + std::to_string(::getpid()) + "-" +
             std::to_string(counter++));
    fs::remove_all(path_);
    fs::create_directories(path_);
  }
  ~TempDir() {
    std::error_code ec;
    fs::remove_all(path_, ec);
  }
  TempDir(const TempDir&) = delete;
  TempDir& operator=(const TempDir&) = delete;

  [[nodiscard]] const fs::path& path() const {
    return path_;
  }

private:
  fs::path path_;
};

// The canonical tiny fixture header used by tests and tests/data/tiny.bam.
inline constexpr std::string_view kTinyHeader = "@HD\tVN:1.6\tSO:coordinate\n"
                                                "@SQ\tSN:chr1\tLN:1000\n"
                                                "@SQ\tSN:chr2\tLN:500\n"
                                                "@RG\tID:rg1\tSM:S1\n"
                                                "@RG\tID:rg2\tSM:S2\n"
                                                "@RG\tID:rg3\tSM:S2\n"
                                                "@RG\tID:rg4\tSM:S3\n"
                                                "@RG\tID:rg5\tSM:S4\n";

// Deterministic reads covering every CIGAR/base-policy case, four samples,
// empty sample S3, and one SEQ=* read.
inline std::vector<ReadSpec> tiny_reads() {
  std::vector<ReadSpec> reads;
  auto add = [&reads](std::string name, uint16_t flag, int64_t pos,
                      uint8_t mapq, std::vector<CigarOp> cigar, std::string seq,
                      std::string rg) {
    reads.push_back(ReadSpec{std::move(name), flag, 0, pos, mapq,
                             std::move(cigar), std::move(seq), std::move(rg)});
  };

  add("r_start", 0, 0, 60, {{'M', 10}}, repeat_sequence(10), "rg1");
  add("r_i", 0, 100, 60, {{'M', 5}, {'I', 2}, {'M', 5}}, repeat_sequence(12),
      "rg2");
  add("r_d", 0, 200, 60, {{'M', 5}, {'D', 2}, {'M', 5}}, repeat_sequence(10),
      "rg2");
  add("r_n", 0, 300, 60, {{'M', 5}, {'N', 10}, {'M', 5}}, repeat_sequence(10),
      "rg3");
  add("r_s", 0, 400, 60, {{'S', 5}, {'M', 10}, {'S', 5}}, repeat_sequence(20),
      "rg1");
  add("r_h", 0, 500, 60, {{'H', 5}, {'M', 10}, {'H', 5}}, repeat_sequence(10),
      "rg1");
  add("r_p", 0, 600, 60, {{'M', 5}, {'P', 1}, {'M', 5}}, repeat_sequence(10),
      "rg1");
  add("r_eqx", 0, 700, 60, {{'M', 5}, {'=', 1}, {'X', 5}}, repeat_sequence(11),
      "rg2");
  add("r_ambig", 0, 800, 60, {{'M', 10}}, "ACGTNNACGT", "rg1");
  add("r_dup", BAM_FDUP, 900, 60, {{'M', 10}}, repeat_sequence(10), "rg1");
  add("r_secondary", BAM_FSECONDARY, 910, 60, {{'M', 10}}, repeat_sequence(10),
      "rg1");
  add("r_qcfail", BAM_FQCFAIL, 920, 60, {{'M', 10}}, repeat_sequence(10),
      "rg1");
  add("r_supp", BAM_FSUPPLEMENTARY, 930, 60, {{'M', 10}}, repeat_sequence(10),
      "rg1");
  add("r_lowq", 0, 950, 10, {{'M', 10}}, repeat_sequence(10), "rg1");
  add("r_star", 0, 960, 60, {{'M', 10}}, "", "rg5");
  add("r_end", 0, 990, 60, {{'M', 10}}, repeat_sequence(10), "rg5");
  return reads;
}

} // namespace tests
