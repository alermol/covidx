#include "covidx/format.h"

namespace covidx {

Layout parse_layout(const std::string& value) {
  if (value == "per-sample")
    return Layout::PerSample;
  if (value == "pooled")
    return Layout::Pooled;
  if (value == "both")
    return Layout::Both;
  throw UsageError("invalid layout '" + value +
                   "'; expected per-sample, pooled, or both");
}

const char* layout_name(Layout layout) {
  switch (layout) {
  case Layout::PerSample:
    return "per-sample";
  case Layout::Pooled:
    return "pooled";
  case Layout::Both:
    return "both";
  }
  return "per-sample";
}

AmbiguousBasePolicy parse_ambiguous_base_policy(const std::string& value) {
  if (value == "cover")
    return AmbiguousBasePolicy::Cover;
  if (value == "exclude")
    return AmbiguousBasePolicy::Exclude;
  throw UsageError("invalid --ambiguous-read-bases value '" + value +
                   "'; expected cover or exclude");
}

const char* ambiguous_base_policy_name(AmbiguousBasePolicy policy) {
  return policy == AmbiguousBasePolicy::Cover ? "cover" : "exclude";
}

} // namespace covidx
