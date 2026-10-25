#include "doctest.h"

#include "covidx/error.h"
#include "covidx/format.h"

TEST_CASE("Format constants match the contract") {
  CHECK(std::string(covidx::kFormatVersion) == "1.0.0");
  CHECK(std::string(covidx::kBuilderVersion) == "0.1.1");
  CHECK(covidx::kMetaRecordSize == 48);
  CHECK(covidx::kPackAlignment == 8);
  CHECK(covidx::kMaxContigLength == 2147483647ULL);
  CHECK(std::string(covidx::kUnitMagic, 8) == "COVIDXU1");
}

TEST_CASE("Layout names round-trip") {
  for (covidx::Layout layout : {covidx::Layout::PerSample,
                                covidx::Layout::Pooled, covidx::Layout::Both}) {
    CHECK(covidx::parse_layout(covidx::layout_name(layout)) == layout);
  }
  CHECK_THROWS_AS((void)covidx::parse_layout("nope"), covidx::UsageError);
}

TEST_CASE("Ambiguous base policy names round-trip") {
  CHECK(covidx::parse_ambiguous_base_policy("cover") ==
        covidx::AmbiguousBasePolicy::Cover);
  CHECK(covidx::parse_ambiguous_base_policy("exclude") ==
        covidx::AmbiguousBasePolicy::Exclude);
  CHECK_THROWS_AS((void)covidx::parse_ambiguous_base_policy("nope"),
                  covidx::UsageError);
}

TEST_CASE("Error hierarchy relationships") {
  CHECK_THROWS_AS(throw covidx::InputError("x"), covidx::IoError);
  CHECK_THROWS_AS(throw covidx::OutputError("x"), covidx::IoError);
  CHECK_THROWS_AS(throw covidx::FormatError("x"), covidx::Error);
  CHECK_THROWS_AS(throw covidx::QueryError("x"), covidx::Error);
  CHECK_THROWS_AS(throw covidx::UsageError("x"), covidx::Error);
}
