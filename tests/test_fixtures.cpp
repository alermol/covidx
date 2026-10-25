#include "doctest.h"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "test_util.h"

namespace fs = std::filesystem;

namespace {

struct BamSummary {
  int contigs = 0;
  int read_groups = 0;
  std::vector<std::string> names;
  std::vector<uint64_t> lengths;
  std::vector<uint64_t> reads_per_contig;
  bool saw_star_seq = false;
  bool saw_ambig_seq = false;
};

BamSummary summarize(const fs::path& path) {
  BamSummary out;
  samFile* fp = sam_open(path.string().c_str(), "r");
  if (fp == nullptr)
    throw std::runtime_error("cannot open test BAM");
  sam_hdr_t* hdr = sam_hdr_read(fp);
  if (hdr == nullptr) {
    sam_close(fp);
    throw std::runtime_error("cannot read test BAM header");
  }
  out.contigs = sam_hdr_nref(hdr);
  out.read_groups = sam_hdr_count_lines(hdr, "RG");
  for (int i = 0; i < out.contigs; ++i) {
    out.names.emplace_back(sam_hdr_tid2name(hdr, i));
    out.lengths.push_back(static_cast<uint64_t>(sam_hdr_tid2len(hdr, i)));
  }
  out.reads_per_contig.assign(static_cast<std::size_t>(out.contigs), 0);

  hts_idx_t* idx = sam_index_load(fp, path.string().c_str());
  if (idx == nullptr) {
    sam_hdr_destroy(hdr);
    sam_close(fp);
    throw std::runtime_error("cannot load test BAM index");
  }
  bam1_t* record = bam_init1();
  for (int tid = 0; tid < out.contigs; ++tid) {
    hts_itr_t* itr = sam_itr_queryi(idx, tid, 0, sam_hdr_tid2len(hdr, tid));
    if (itr == nullptr)
      throw std::runtime_error("cannot create test iterator");
    int rc = 0;
    while ((rc = sam_itr_next(fp, itr, record)) >= 0) {
      ++out.reads_per_contig[static_cast<std::size_t>(tid)];
      if (record->core.l_qseq == 0)
        out.saw_star_seq = true;
      const uint8_t* seq = bam_get_seq(record);
      for (int i = 0; i < record->core.l_qseq; ++i) {
        if (bam_seqi(seq, i) == 15)
          out.saw_ambig_seq = true;
      }
    }
    if (rc != -1)
      throw std::runtime_error("test BAM iterator failed");
    hts_itr_destroy(itr);
  }
  bam_destroy1(record);
  hts_idx_destroy(idx);
  sam_hdr_destroy(hdr);
  sam_close(fp);
  return out;
}

} // namespace

TEST_CASE("make_bam writes a readable indexed BAM") {
  tests::TempDir temp;
  const fs::path bam = temp.path() / "tiny.bam";
  tests::make_bam(bam, tests::kTinyHeader, tests::tiny_reads());
  CHECK(fs::exists(bam));
  CHECK(fs::exists(fs::path(bam.string() + ".bai")));

  const BamSummary summary = summarize(bam);
  CHECK(summary.contigs == 2);
  CHECK(summary.read_groups == 5);
  CHECK(summary.names.at(0) == "chr1");
  CHECK(summary.names.at(1) == "chr2");
  CHECK(summary.lengths.at(0) == 1000);
  CHECK(summary.lengths.at(1) == 500);
  CHECK(summary.reads_per_contig.at(0) == tests::tiny_reads().size());
  CHECK(summary.reads_per_contig.at(1) == 0);
  CHECK(summary.saw_star_seq);
  CHECK(summary.saw_ambig_seq);
}

TEST_CASE("committed tiny fixture is present and stable") {
  const fs::path bam = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  REQUIRE(fs::exists(bam));
  REQUIRE(fs::exists(fs::path(bam.string() + ".bai")));
  const BamSummary summary = summarize(bam);
  CHECK(summary.contigs == 2);
  CHECK(summary.read_groups == 5);
  CHECK(summary.reads_per_contig.at(0) == tests::tiny_reads().size());
  CHECK(summary.reads_per_contig.at(1) == 0);
}
