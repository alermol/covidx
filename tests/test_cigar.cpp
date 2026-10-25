#include "doctest.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "covidx/detail/cigar.h"
#include "covidx/error.h"
#include "test_util.h"

using covidx::AmbiguousBasePolicy;
using covidx::detail::Block;
using covidx::detail::ExtractPolicy;

namespace {

struct BamDeleter {
  void operator()(bam1_t* b) const {
    bam_destroy1(b);
  }
};
using BamPtr = std::unique_ptr<bam1_t, BamDeleter>;

BamPtr make_read(int64_t pos, const std::vector<tests::CigarOp>& cigar_ops,
                 const std::string& seq, uint16_t flag = 0, uint8_t mapq = 60) {
  const std::vector<uint32_t> cigar = tests::encode_cigar(cigar_ops);
  BamPtr read(bam_init1());
  const int32_t l_seq = static_cast<int32_t>(seq.size());
  std::string qual(seq.size(), 'I');
  const int rc =
      bam_set1(read.get(), 4, "read", flag, 0, pos, mapq, cigar.size(),
               cigar.data(), -1, -1, 0, static_cast<std::size_t>(l_seq),
               l_seq == 0 ? nullptr : seq.c_str(),
               l_seq == 0 ? nullptr : qual.c_str(), 0);
  if (rc < 0)
    throw std::runtime_error("bam_set1 failed in test");
  return read;
}

std::vector<Block>
extract(int64_t pos, const std::vector<tests::CigarOp>& cigar,
        const std::string& seq, bool count_deletions = true,
        AmbiguousBasePolicy ambiguous = AmbiguousBasePolicy::Cover,
        uint64_t length = 1000) {
  BamPtr read = make_read(pos, cigar, seq);
  return covidx::detail::extract_blocks(
      *read, ExtractPolicy{count_deletions, ambiguous}, length);
}

void check_blocks(const std::vector<Block>& got,
                  const std::vector<std::pair<uint64_t, uint64_t>>& want) {
  REQUIRE(got.size() == want.size());
  for (std::size_t i = 0; i < want.size(); ++i) {
    CHECK(got[i].start == want[i].first);
    CHECK(got[i].end == want[i].second);
  }
}

} // namespace

TEST_CASE("CIGAR policy table") {
  check_blocks(extract(0, {{'M', 10}}, tests::repeat_sequence(10)), {{0, 10}});
  check_blocks(
      extract(0, {{'M', 5}, {'I', 2}, {'M', 5}}, tests::repeat_sequence(12)),
      {{0, 10}});
  check_blocks(extract(0, {{'M', 5}, {'D', 2}, {'M', 5}},
                       tests::repeat_sequence(10), true),
               {{0, 12}});
  check_blocks(extract(0, {{'M', 5}, {'D', 2}, {'M', 5}},
                       tests::repeat_sequence(10), false),
               {{0, 5}, {7, 12}});
  check_blocks(
      extract(0, {{'M', 5}, {'N', 10}, {'M', 5}}, tests::repeat_sequence(10)),
      {{0, 5}, {15, 20}});
  check_blocks(
      extract(0, {{'D', 2}, {'M', 10}}, tests::repeat_sequence(10), true),
      {{0, 12}});
  check_blocks(
      extract(0, {{'D', 2}, {'M', 10}}, tests::repeat_sequence(10), false),
      {{2, 12}});
  check_blocks(
      extract(0, {{'M', 10}, {'D', 2}}, tests::repeat_sequence(10), true),
      {{0, 12}});
  check_blocks(
      extract(0, {{'M', 10}, {'D', 2}}, tests::repeat_sequence(10), false),
      {{0, 10}});
  check_blocks(
      extract(0, {{'M', 5}, {'P', 1}, {'M', 5}}, tests::repeat_sequence(10)),
      {{0, 10}});
  check_blocks(
      extract(0, {{'S', 5}, {'M', 10}, {'S', 5}}, tests::repeat_sequence(20)),
      {{0, 10}});
  check_blocks(
      extract(0, {{'H', 5}, {'M', 10}, {'H', 5}}, tests::repeat_sequence(10)),
      {{0, 10}});
  check_blocks(
      extract(0, {{'M', 5}, {'=', 1}, {'X', 5}}, tests::repeat_sequence(11)),
      {{0, 11}});
  check_blocks(extract(0, {{'N', 5}, {'M', 10}}, tests::repeat_sequence(10)),
               {{5, 15}});
  check_blocks(extract(0, {{'M', 10}, {'N', 5}}, tests::repeat_sequence(10)),
               {{0, 10}});
}

TEST_CASE("CIGAR rejection cases") {
  CHECK_THROWS_AS(
      extract(0, {{'M', 5}, {'B', 1}, {'M', 5}}, tests::repeat_sequence(10)),
      covidx::InputError);
  CHECK_THROWS_AS(extract(0, {{'M', 0}}, ""), covidx::InputError);
  CHECK_THROWS_AS(extract(-1, {{'M', 10}}, tests::repeat_sequence(10)),
                  covidx::InputError);
  CHECK_THROWS_AS(extract(95, {{'M', 10}}, tests::repeat_sequence(10), true,
                          AmbiguousBasePolicy::Cover, 100),
                  covidx::InputError);
  CHECK_THROWS_AS(extract(0, {}, ""), covidx::InputError);
}

TEST_CASE("Ambiguous read bases cover by default") {
  check_blocks(extract(0, {{'M', 10}}, "ACGTNNACGT"), {{0, 10}});
}

TEST_CASE("Ambiguous read bases exclude splits blocks") {
  check_blocks(
      extract(0, {{'M', 10}}, "ACGTNNACGT", true, AmbiguousBasePolicy::Exclude),
      {{0, 4}, {6, 10}});
  check_blocks(
      extract(0, {{'M', 5}}, "NACGT", true, AmbiguousBasePolicy::Exclude),
      {{1, 5}});
  check_blocks(
      extract(0, {{'M', 10}}, "ACGTACGTRY", true, AmbiguousBasePolicy::Exclude),
      {{0, 8}});
  check_blocks(
      extract(0, {{'M', 8}}, "ACGT====", true, AmbiguousBasePolicy::Exclude),
      {{0, 4}});
  check_blocks(extract(0, {{'M', 4}, {'D', 2}, {'M', 4}}, "ACGNACGT", true,
                       AmbiguousBasePolicy::Exclude),
               {{0, 3}, {4, 10}});
}

TEST_CASE("SEQ=* is accepted under cover and rejected under exclude") {
  check_blocks(extract(0, {{'M', 10}}, ""), {{0, 10}});
  CHECK_THROWS_AS(
      extract(0, {{'M', 10}}, "", true, AmbiguousBasePolicy::Exclude),
      covidx::InputError);
}

TEST_CASE("ReadFilter applies flags, MAPQ, and extraction") {
  covidx::detail::ReadFilter filter;
  filter.contig_length = 1000;
  filter.effective_exclude = 0x0F04;
  std::vector<Block> blocks;

  BamPtr normal = make_read(10, {{'M', 10}}, tests::repeat_sequence(10));
  CHECK(filter.accept(*normal, blocks));
  REQUIRE(blocks.size() == 1);
  CHECK(blocks[0].start == 10);

  BamPtr unmapped =
      make_read(10, {{'M', 10}}, tests::repeat_sequence(10), BAM_FUNMAP);
  CHECK_FALSE(filter.accept(*unmapped, blocks));
  BamPtr secondary =
      make_read(10, {{'M', 10}}, tests::repeat_sequence(10), BAM_FSECONDARY);
  CHECK_FALSE(filter.accept(*secondary, blocks));
  BamPtr supplementary = make_read(10, {{'M', 10}}, tests::repeat_sequence(10),
                                   BAM_FSUPPLEMENTARY);
  CHECK_FALSE(filter.accept(*supplementary, blocks));

  filter.effective_exclude = 0x0704;
  CHECK(filter.accept(*supplementary, blocks));

  filter.min_mapq = 20;
  BamPtr low = make_read(10, {{'M', 10}}, tests::repeat_sequence(10), 0, 10);
  CHECK_FALSE(filter.accept(*low, blocks));
  BamPtr unavailable =
      make_read(10, {{'M', 10}}, tests::repeat_sequence(10), 0, 255);
  CHECK(filter.accept(*unavailable, blocks));
}
