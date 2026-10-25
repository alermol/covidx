#include "doctest.h"

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>

#include "covidx/error.h"
#include "covidx/manifest.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::Manifest;

namespace {

Manifest valid_manifest() {
  Manifest manifest;
  manifest.created = "2026-10-24T00:00:00Z";
  manifest.layout = "both";
  manifest.contigs = {
      covidx::ContigInfo{0, "chr1", 1000},
      covidx::ContigInfo{1, "chr2", 500},
  };
  manifest.scope = {
      covidx::ScopeInfo{"chr1", 0, 1000},
      covidx::ScopeInfo{"chr2", 0, 500},
  };
  manifest.samples = {
      covidx::SampleInfo{"S1", {"rg1"}},
      covidx::SampleInfo{"S2", {"rg2"}},
  };
  covidx::UnitInfo s1;
  s1.id = 0;
  s1.name = "S1";
  s1.kind = covidx::UnitInfo::Kind::Sample;
  s1.members = {"S1"};
  s1.file = "data/u000000.bin";
  s1.file_size = 64;
  s1.total_runs = 2;
  s1.total_covered = 10;
  s1.contigs_with_runs = 1;
  covidx::UnitInfo s2 = s1;
  s2.id = 1;
  s2.name = "S2";
  s2.members = {"S2"};
  s2.file = "data/u000001.bin";
  covidx::UnitInfo pooled = s1;
  pooled.id = 2;
  pooled.name = "pooled";
  pooled.kind = covidx::UnitInfo::Kind::Pooled;
  pooled.members = {"S1", "S2"};
  pooled.file = "data/u000002.bin";
  manifest.units = {s1, s2, pooled};
  manifest.files = {
      {"data/u000000.bin", covidx::FileInfo{64, 0}},
      {"data/u000001.bin", covidx::FileInfo{64, 0}},
      {"data/u000002.bin", covidx::FileInfo{64, 0}},
      {"contigs.tsv", covidx::FileInfo{20, 0}},
  };
  manifest.index.contig_count = 2;
  manifest.index.unit_count = 3;
  manifest.index.total_runs = 6;
  manifest.index.partial_scope = false;
  manifest.source.index_path = "in.bam.bai";
  manifest.source.path = "in.bam";
  manifest.source.records_accepted = 12;
  manifest.source.records_examined = 16;
  manifest.source.size = 528;
  manifest.warnings = {"b warning", "a warning"};
  return manifest;
}

std::string replace_once(std::string text, std::string_view from,
                         std::string_view to) {
  const std::size_t pos = text.find(from);
  if (pos == std::string::npos)
    throw std::runtime_error("test replacement did not match");
  text.replace(pos, from.size(), to);
  return text;
}

} // namespace

TEST_CASE("manifest round-trips with full and partial scope") {
  Manifest manifest = valid_manifest();
  const std::string text = manifest.to_json();
  const Manifest parsed = Manifest::from_json(text, "test");
  CHECK(parsed.contigs.size() == 2);
  CHECK(parsed.scope.size() == 2);
  CHECK(parsed.samples.size() == 2);
  CHECK(parsed.units.size() == 3);
  CHECK(parsed.layout == "both");
  CHECK(parsed.policy.effective_exclude_flags == 0x0F04);
  CHECK(parsed.warnings.size() == 2);
  CHECK(parsed.to_json() == text);

  Manifest partial = valid_manifest();
  partial.scope = {covidx::ScopeInfo{"chr1", 100, 500},
                   covidx::ScopeInfo{"chr2", 0, 500}};
  partial.index.partial_scope = true;
  const Manifest parsed_partial =
      Manifest::from_json(partial.to_json(), "test");
  CHECK(parsed_partial.index.partial_scope);
  CHECK(parsed_partial.scope[0].start == 100);
}

TEST_CASE("manifest is canonical") {
  const std::string text = valid_manifest().to_json();
  CHECK(text.starts_with(
      "{\n  \"builder_version\": \"0.2.0\",\n  \"contigs\": [\n"));
  CHECK(text.ends_with("}\n"));
  CHECK(text.find("\"generated\"") == std::string::npos);
  CHECK(text.find("\n\n") == std::string::npos);
  const std::size_t warning_a = text.find("a warning");
  const std::size_t warning_b = text.find("b warning");
  CHECK(warning_a < warning_b);
}

TEST_CASE("manifest ignores unknown keys") {
  std::string text = valid_manifest().to_json();
  text =
      replace_once(text, "{\n", "{\n  \"future_key\": {\"nested\": [1, 2]},\n");
  CHECK_NOTHROW((void)Manifest::from_json(text, "test"));
}

TEST_CASE(
    "manifest rejects missing keys, duplicate keys, and version mismatches") {
  CHECK_THROWS_AS((void)Manifest::from_json("{}", "test"), covidx::FormatError);
  std::string duplicate = valid_manifest().to_json();
  duplicate = replace_once(duplicate, "{\n", "{\n  \"format\": \"covidx\",\n");
  CHECK_THROWS_AS((void)Manifest::from_json(duplicate, "test"),
                  covidx::FormatError);

  std::string wrong_major = valid_manifest().to_json();
  wrong_major = replace_once(wrong_major, "\"format_version\": \"1.1.0\"",
                             "\"format_version\": \"2.0.0\"");
  CHECK_THROWS_AS((void)Manifest::from_json(wrong_major, "test"),
                  covidx::FormatError);

  std::string wrong_minor = valid_manifest().to_json();
  wrong_minor = replace_once(wrong_minor, "\"format_version\": \"1.1.0\"",
                             "\"format_version\": \"1.0.0\"");
  CHECK_THROWS_AS((void)Manifest::from_json(wrong_minor, "test"),
                  covidx::FormatError);

  std::string no_units = valid_manifest().to_json();
  no_units = replace_once(no_units, "\"units\": [", "\"not_units\": [");
  CHECK_THROWS_AS((void)Manifest::from_json(no_units, "test"),
                  covidx::FormatError);

  std::string bad_type = valid_manifest().to_json();
  bad_type = replace_once(bad_type, "\"format\": \"covidx\"", "\"format\": 5");
  CHECK_THROWS_AS((void)Manifest::from_json(bad_type, "test"),
                  covidx::FormatError);
}

TEST_CASE("manifest validates scope ordering, merging, and ranges") {
  Manifest manifest = valid_manifest();
  manifest.scope = {covidx::ScopeInfo{"chr2", 0, 500},
                    covidx::ScopeInfo{"chr1", 0, 1000}};
  CHECK_THROWS_AS((void)Manifest::from_json(manifest.to_json(), "test"),
                  covidx::FormatError);

  manifest = valid_manifest();
  manifest.scope = {covidx::ScopeInfo{"chr1", 0, 600},
                    covidx::ScopeInfo{"chr1", 500, 1000},
                    covidx::ScopeInfo{"chr2", 0, 500}};
  CHECK_THROWS_AS((void)Manifest::from_json(manifest.to_json(), "test"),
                  covidx::FormatError);

  manifest = valid_manifest();
  manifest.scope = {covidx::ScopeInfo{"chr1", 0, 500},
                    covidx::ScopeInfo{"chr1", 500, 1000},
                    covidx::ScopeInfo{"chr2", 0, 500}};
  CHECK_THROWS_AS((void)Manifest::from_json(manifest.to_json(), "test"),
                  covidx::FormatError);

  manifest = valid_manifest();
  manifest.scope = {covidx::ScopeInfo{"chr1", 0, 1001},
                    covidx::ScopeInfo{"chr2", 0, 500}};
  CHECK_THROWS_AS((void)Manifest::from_json(manifest.to_json(), "test"),
                  covidx::FormatError);

  manifest = valid_manifest();
  manifest.scope = {covidx::ScopeInfo{"chr1", 0, 1000}};
  CHECK_THROWS_AS((void)Manifest::from_json(manifest.to_json(), "test"),
                  covidx::FormatError);

  manifest = valid_manifest();
  manifest.scope = {covidx::ScopeInfo{"chr1", 100, 500},
                    covidx::ScopeInfo{"chr2", 0, 500}};
  manifest.index.partial_scope = false;
  CHECK_THROWS_AS((void)Manifest::from_json(manifest.to_json(), "test"),
                  covidx::FormatError);
}

TEST_CASE("manifest feature-count metadata is mandatory") {
  Manifest manifest = valid_manifest();
  manifest.units[0].total_features = 2;
  manifest.units[1].total_features = 3;
  manifest.units[2].total_features = 5;
  const std::string text = manifest.to_json();
  CHECK(text.find("\"total_features\": 5") != std::string::npos);
  CHECK(text.find("\"features\":") != std::string::npos);
  const Manifest parsed = Manifest::from_json(text, "test");
  CHECK(parsed.units[0].total_features == 2);
  CHECK(parsed.units[2].total_features == 5);
  CHECK(parsed.to_json() == text);

  Manifest old_format = manifest;
  old_format.format_version = "1.0.0";
  CHECK_THROWS_AS((void)Manifest::from_json(old_format.to_json(), "test"),
                  covidx::FormatError);

  Manifest mismatch = manifest;
  mismatch.units[2].total_features = 4;
  CHECK_THROWS_AS((void)Manifest::from_json(mismatch.to_json(), "test"),
                  covidx::FormatError);
}

TEST_CASE("manifest validates names, CRCs, and source hashes") {
  Manifest manifest = valid_manifest();
  manifest.samples[0].name = "pooled";
  CHECK_THROWS_AS((void)Manifest::from_json(manifest.to_json(), "test"),
                  covidx::FormatError);

  manifest = valid_manifest();
  manifest.samples[0].name = "bad\tname";
  CHECK_THROWS_AS((void)Manifest::from_json(manifest.to_json(), "test"),
                  covidx::FormatError);

  manifest = valid_manifest();
  manifest.samples[0].name = std::string("bad\xff", 4);
  CHECK_THROWS_AS((void)Manifest::from_json(manifest.to_json(), "test"),
                  covidx::FormatError);

  manifest = valid_manifest();
  manifest.files["data/u000000.bin"].crc32 = 0xABCDEF01;
  std::string text = manifest.to_json();
  text =
      replace_once(text, "\"crc32\": \"abcdef01\"", "\"crc32\": \"ABCDEF01\"");
  CHECK_THROWS_AS((void)Manifest::from_json(text, "test"), covidx::FormatError);

  manifest = valid_manifest();
  text = manifest.to_json();
  text = replace_once(text, "\"sha256\": null", "\"sha256\": \"abc\"");
  CHECK_THROWS_AS((void)Manifest::from_json(text, "test"), covidx::FormatError);
}

TEST_CASE("manifest rejects contigs beyond INT32_MAX") {
  Manifest manifest = valid_manifest();
  manifest.contigs[0].length = covidx::kMaxContigLength + 1;
  manifest.scope[0].end = covidx::kMaxContigLength + 1;
  CHECK_THROWS_AS((void)Manifest::from_json(manifest.to_json(), "test"),
                  covidx::FormatError);
}

TEST_CASE("manifest save and load") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "manifest.json";
  const Manifest manifest = valid_manifest();
  manifest.save(path.string());
  const Manifest loaded = Manifest::load(path.string());
  CHECK(loaded.to_json() == manifest.to_json());
}

TEST_CASE("SOURCE_DATE_EPOCH handling") {
  ::setenv("SOURCE_DATE_EPOCH", "0", 1);
  CHECK(covidx::created_from_env() == "1970-01-01T00:00:00Z");
  ::setenv("SOURCE_DATE_EPOCH", "1", 1);
  CHECK(covidx::created_from_env() == "1970-01-01T00:00:01Z");
  ::setenv("SOURCE_DATE_EPOCH", "1700000000", 1);
  CHECK(covidx::created_from_env() == "2023-11-14T22:13:20Z");
  ::setenv("SOURCE_DATE_EPOCH", "-1", 1);
  CHECK_THROWS_AS((void)covidx::created_from_env(), covidx::UsageError);
  ::setenv("SOURCE_DATE_EPOCH", "abc", 1);
  CHECK_THROWS_AS((void)covidx::created_from_env(), covidx::UsageError);
  ::setenv("SOURCE_DATE_EPOCH", "", 1);
  CHECK_THROWS_AS((void)covidx::created_from_env(), covidx::UsageError);
  ::unsetenv("SOURCE_DATE_EPOCH");
  CHECK(covidx::created_from_env().ends_with("Z"));
}
