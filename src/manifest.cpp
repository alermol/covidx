#include "covidx/manifest.h"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "covidx/detail/crc32.h"
#include "covidx/detail/json.h"
#include "covidx/detail/text.h"
#include "covidx/error.h"

namespace covidx {

namespace {

using detail::JsonValue;

[[nodiscard]] const JsonValue& require_object(const JsonValue& value,
                                              std::string_view what) {
  if (!value.is_object())
    throw FormatError(std::string(what) + " must be a JSON object");
  return value;
}

[[nodiscard]] const JsonValue& require_field(const JsonValue& object,
                                             std::string_view key,
                                             const std::string& context) {
  const JsonValue* value = object.find(key);
  if (value == nullptr) {
    throw FormatError(context + ": missing required key '" + std::string(key) +
                      "'");
  }
  return *value;
}

[[nodiscard]] std::string require_string(const JsonValue& object,
                                         std::string_view key,
                                         const std::string& context) {
  const JsonValue& value = require_field(object, key, context);
  if (!value.is_string()) {
    throw FormatError(context + ": key '" + std::string(key) +
                      "' must be a string");
  }
  return value.text;
}

[[nodiscard]] uint64_t require_uint(const JsonValue& object,
                                    std::string_view key,
                                    const std::string& context) {
  const JsonValue& value = require_field(object, key, context);
  if (!value.is_number() || !value.is_integer) {
    throw FormatError(context + ": key '" + std::string(key) +
                      "' must be an unsigned integer");
  }
  return value.integer;
}

[[nodiscard]] bool require_bool(const JsonValue& object, std::string_view key,
                                const std::string& context) {
  const JsonValue& value = require_field(object, key, context);
  if (!value.is_bool()) {
    throw FormatError(context + ": key '" + std::string(key) +
                      "' must be a boolean");
  }
  return value.boolean;
}

[[nodiscard]] const JsonValue& require_array(const JsonValue& object,
                                             std::string_view key,
                                             const std::string& context) {
  const JsonValue& value = require_field(object, key, context);
  if (!value.is_array()) {
    throw FormatError(context + ": key '" + std::string(key) +
                      "' must be an array");
  }
  return value;
}

[[nodiscard]] std::vector<std::string>
string_array(const JsonValue& object, std::string_view key,
             const std::string& context) {
  const JsonValue& array = require_array(object, key, context);
  std::vector<std::string> out;
  out.reserve(array.array.size());
  for (const JsonValue& item : array.array) {
    if (!item.is_string()) {
      throw FormatError(context + ": '" + std::string(key) +
                        "' entries must be strings");
    }
    out.push_back(item.text);
  }
  return out;
}

[[nodiscard]] std::string validate_name_field(const std::string& name,
                                              const char* kind) {
  try {
    detail::validate_unit_name(name, kind, true);
  } catch (const UsageError& e) {
    throw FormatError(std::string("invalid ") + kind +
                      " name in manifest: " + e.what());
  }
  return name;
}

[[nodiscard]] UnitInfo::Kind parse_kind(const std::string& value,
                                        const std::string& context) {
  if (value == "sample")
    return UnitInfo::Kind::Sample;
  if (value == "group")
    return UnitInfo::Kind::Group;
  if (value == "pooled")
    return UnitInfo::Kind::Pooled;
  throw FormatError(context + ": invalid unit kind '" + value + "'");
}

[[nodiscard]] std::string kind_name(UnitInfo::Kind kind) {
  switch (kind) {
  case UnitInfo::Kind::Sample:
    return "sample";
  case UnitInfo::Kind::Group:
    return "group";
  case UnitInfo::Kind::Pooled:
    return "pooled";
  }
  return "sample";
}

[[nodiscard]] bool is_hex8(std::string_view text) {
  if (text.size() != 8)
    return false;
  for (const char c : text) {
    if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')))
      return false;
  }
  return true;
}

[[nodiscard]] bool safe_relative_path(std::string_view path) {
  if (path.empty() || path.front() == '/' ||
      path.find('\0') != std::string_view::npos)
    return false;
  std::size_t start = 0;
  while (start <= path.size()) {
    const std::size_t slash = path.find('/', start);
    const std::string_view part = path.substr(
        start, slash == std::string_view::npos ? std::string_view::npos
                                               : slash - start);
    if (part == "..")
      return false;
    if (slash == std::string_view::npos)
      break;
    start = slash + 1;
  }
  return true;
}

[[nodiscard]] bool is_hex64(std::string_view text) {
  if (text.size() != 64)
    return false;
  for (const char c : text) {
    if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')))
      return false;
  }
  return true;
}

[[nodiscard]] uint32_t parse_hex8(const std::string& text,
                                  const std::string& context) {
  if (!is_hex8(text))
    throw FormatError(context + ": CRC32 must be 8 lowercase hex digits");
  uint32_t value = 0;
  for (const char c : text) {
    value <<= 4;
    value |= static_cast<uint32_t>(c <= '9' ? c - '0' : c - 'a' + 10);
  }
  return value;
}

[[nodiscard]] JsonValue make_object() {
  return JsonValue::make_object();
}

void set_string(JsonValue& object, const char* key, const std::string& value) {
  object.set(key, JsonValue::make_string(value));
}

void set_uint(JsonValue& object, const char* key, uint64_t value) {
  object.set(key, JsonValue::make_uint(value));
}

void set_bool(JsonValue& object, const char* key, bool value) {
  object.set(key, JsonValue::make_bool(value));
}

[[nodiscard]] JsonValue
string_array_value(const std::vector<std::string>& values) {
  JsonValue array = JsonValue::make_array();
  for (const std::string& value : values)
    array.array.push_back(JsonValue::make_string(value));
  return array;
}

} // namespace

std::string Manifest::to_json() const {
  JsonValue root = make_object();
  set_string(root, "format", format);
  set_string(root, "format_version", format_version);
  set_string(root, "builder_version", builder_version);
  set_string(root, "created", created);
  set_string(root, "endianness", endianness);
  set_string(root, "layout", layout);

  JsonValue source_value = make_object();
  set_string(source_value, "format", source.format);
  set_string(source_value, "header_so", source.header_so);
  set_string(source_value, "index_path", source.index_path);
  set_string(source_value, "path", source.path);
  set_uint(source_value, "records_accepted", source.records_accepted);
  set_uint(source_value, "records_examined", source.records_examined);
  if (source.sha256.has_value()) {
    source_value.set("sha256", JsonValue::make_string(*source.sha256));
  } else {
    source_value.set("sha256", JsonValue::make_null());
  }
  set_uint(source_value, "size", source.size);
  root.set("source", std::move(source_value));

  JsonValue policy_value = make_object();
  set_string(policy_value, "ambiguous_read_bases", policy.ambiguous_read_bases);
  policy_value.set("cigar_covered_ops", string_array_value({"M", "=", "X"}));
  policy_value.set("cigar_ignored_ops",
                   string_array_value({"I", "S", "H", "P"}));
  policy_value.set("cigar_split_ops", string_array_value({"N"}));
  set_bool(policy_value, "count_deletions", policy.count_deletions);
  set_string(policy_value, "deletion_policy",
             policy.count_deletions ? "covered" : "gap");
  set_uint(policy_value, "effective_exclude_flags",
           policy.effective_exclude_flags);
  set_uint(policy_value, "exclude_flags", policy.exclude_flags);
  set_bool(policy_value, "include_supplementary", policy.include_supplementary);
  set_uint(policy_value, "min_mapq", policy.min_mapq);
  set_string(policy_value, "reference_base_identity",
             policy.reference_base_identity);
  set_string(policy_value, "unmapped_policy", policy.unmapped_policy);
  root.set("policy", std::move(policy_value));

  JsonValue contigs_value = JsonValue::make_array();
  for (const ContigInfo& contig : contigs) {
    JsonValue item = make_object();
    set_uint(item, "id", contig.id);
    set_uint(item, "length", contig.length);
    set_string(item, "name", contig.name);
    contigs_value.array.push_back(std::move(item));
  }
  root.set("contigs", std::move(contigs_value));

  JsonValue scope_value = JsonValue::make_array();
  for (const ScopeInfo& interval : scope) {
    JsonValue item = make_object();
    set_string(item, "contig", interval.contig);
    set_uint(item, "start", interval.start);
    set_uint(item, "end", interval.end);
    scope_value.array.push_back(std::move(item));
  }
  root.set("scope", std::move(scope_value));

  JsonValue samples_value = JsonValue::make_array();
  for (const SampleInfo& sample : samples) {
    JsonValue item = make_object();
    set_string(item, "name", sample.name);
    item.set("read_groups", string_array_value(sample.read_groups));
    samples_value.array.push_back(std::move(item));
  }
  root.set("samples", std::move(samples_value));

  JsonValue units_value = JsonValue::make_array();
  for (const UnitInfo& unit : units) {
    JsonValue item = make_object();
    set_uint(item, "contigs_with_runs", unit.contigs_with_runs);
    set_string(item, "file", unit.file);
    set_string(item, "file_crc32", detail::crc32_hex(unit.file_crc32));
    set_uint(item, "file_size", unit.file_size);
    set_uint(item, "id", unit.id);
    set_string(item, "kind", kind_name(unit.kind));
    item.set("members", string_array_value(unit.members));
    set_string(item, "name", unit.name);
    set_uint(item, "total_covered", unit.total_covered);
    set_uint(item, "total_features", unit.total_features);
    set_uint(item, "total_runs", unit.total_runs);
    units_value.array.push_back(std::move(item));
  }
  root.set("units", std::move(units_value));

  JsonValue files_value = make_object();
  for (const auto& [path, info] : files) {
    JsonValue item = make_object();
    set_string(item, "crc32", detail::crc32_hex(info.crc32));
    set_uint(item, "size", info.size);
    files_value.set(path, std::move(item));
  }
  root.set("files", std::move(files_value));

  JsonValue index_value = make_object();
  set_uint(index_value, "contig_count", index.contig_count);
  set_uint(index_value, "meta_record_size", index.meta_record_size);
  set_uint(index_value, "pack_alignment", index.pack_alignment);
  set_bool(index_value, "partial_scope", index.partial_scope);
  set_uint(index_value, "total_runs", index.total_runs);
  set_uint(index_value, "unit_count", index.unit_count);
  root.set("index", std::move(index_value));

  JsonValue semantics_value = make_object();
  set_string(semantics_value, "covered", semantics_covered);
  set_string(semantics_value, "features", semantics_features);
  root.set("semantics", std::move(semantics_value));

  std::vector<std::string> sorted_warnings = warnings;
  std::sort(sorted_warnings.begin(), sorted_warnings.end());
  root.set("warnings", string_array_value(sorted_warnings));

  return detail::to_canonical_json(root);
}

Manifest Manifest::from_json(std::string_view text, std::string_view origin) {
  const JsonValue root = detail::parse_json(text, origin);
  if (!root.is_object())
    throw FormatError("manifest must be a JSON object");
  const std::string context = "manifest " + std::string(origin);

  Manifest manifest;
  manifest.format = require_string(root, "format", context);
  if (manifest.format != "covidx") {
    throw FormatError(context + ": unsupported format '" + manifest.format +
                      "'");
  }
  manifest.format_version = require_string(root, "format_version", context);
  if (manifest.format_version != kFormatVersion) {
    throw FormatError(context + ": unsupported format_version '" +
                      manifest.format_version + "'");
  }
  manifest.builder_version = require_string(root, "builder_version", context);
  manifest.created = require_string(root, "created", context);
  manifest.endianness = require_string(root, "endianness", context);
  if (manifest.endianness != "little") {
    throw FormatError(context + ": only little-endian indexes are supported");
  }
  manifest.layout = require_string(root, "layout", context);
  try {
    (void)parse_layout(manifest.layout);
  } catch (const UsageError& e) {
    throw FormatError(context + ": " + e.what());
  }

  const JsonValue& source =
      require_object(require_field(root, "source", context), "source");
  manifest.source.format =
      require_string(source, "format", context + ".source");
  if (manifest.source.format != "bam") {
    throw FormatError(context + ": unsupported source format '" +
                      manifest.source.format + "'");
  }
  manifest.source.header_so =
      require_string(source, "header_so", context + ".source");
  if (manifest.source.header_so != "coordinate") {
    throw FormatError(context + ": source header_so must be coordinate");
  }
  manifest.source.index_path =
      require_string(source, "index_path", context + ".source");
  manifest.source.path = require_string(source, "path", context + ".source");
  manifest.source.records_accepted =
      require_uint(source, "records_accepted", context + ".source");
  manifest.source.records_examined =
      require_uint(source, "records_examined", context + ".source");
  const JsonValue* sha = source.find("sha256");
  if (sha == nullptr)
    throw FormatError(context + ": missing required key 'source.sha256'");
  if (sha->is_string()) {
    if (!is_hex64(sha->text)) {
      throw FormatError(
          context + ": source.sha256 must be 64 lowercase hex digits or null");
    }
    manifest.source.sha256 = sha->text;
  } else if (!sha->is_null()) {
    throw FormatError(context + ": source.sha256 must be a string or null");
  }
  manifest.source.size = require_uint(source, "size", context + ".source");

  const JsonValue& policy =
      require_object(require_field(root, "policy", context), "policy");
  manifest.policy.ambiguous_read_bases =
      require_string(policy, "ambiguous_read_bases", context + ".policy");
  if (manifest.policy.ambiguous_read_bases != "cover" &&
      manifest.policy.ambiguous_read_bases != "exclude") {
    throw FormatError(context + ": invalid ambiguous_read_bases policy");
  }
  manifest.policy.count_deletions =
      require_bool(policy, "count_deletions", context + ".policy");
  const std::string deletion_policy =
      require_string(policy, "deletion_policy", context + ".policy");
  if (deletion_policy !=
      (manifest.policy.count_deletions ? "covered" : "gap")) {
    throw FormatError(context +
                      ": deletion_policy disagrees with count_deletions");
  }
  const uint64_t effective =
      require_uint(policy, "effective_exclude_flags", context + ".policy");
  const uint64_t exclude =
      require_uint(policy, "exclude_flags", context + ".policy");
  if (effective > 0xFFFF || exclude > 0xFFFF) {
    throw FormatError(context + ": flag masks must fit in 16 bits");
  }
  manifest.policy.effective_exclude_flags = static_cast<uint16_t>(effective);
  manifest.policy.exclude_flags = static_cast<uint16_t>(exclude);
  manifest.policy.include_supplementary =
      require_bool(policy, "include_supplementary", context + ".policy");
  const uint64_t min_mapq =
      require_uint(policy, "min_mapq", context + ".policy");
  if (min_mapq > 255)
    throw FormatError(context + ": min_mapq must be <= 255");
  manifest.policy.min_mapq = static_cast<uint8_t>(min_mapq);
  manifest.policy.reference_base_identity =
      require_string(policy, "reference_base_identity", context + ".policy");
  manifest.policy.unmapped_policy =
      require_string(policy, "unmapped_policy", context + ".policy");
  const auto check_ops = [&](const char* key,
                             std::initializer_list<const char*> expected) {
    const std::vector<std::string> actual =
        string_array(policy, key, context + ".policy");
    if (actual.size() != expected.size() ||
        !std::equal(actual.begin(), actual.end(), expected.begin())) {
      throw FormatError(context + ": policy." + key +
                        " is inconsistent with the format");
    }
  };
  check_ops("cigar_covered_ops", {"M", "=", "X"});
  check_ops("cigar_ignored_ops", {"I", "S", "H", "P"});
  check_ops("cigar_split_ops", {"N"});
  if (manifest.policy.reference_base_identity != "ignored" ||
      manifest.policy.unmapped_policy != "skip") {
    throw FormatError(context +
                      ": unsupported reference-base or unmapped policy");
  }
  const uint16_t expected_effective = static_cast<uint16_t>(
      (manifest.policy.exclude_flags | 0x0800) &
      ~(manifest.policy.include_supplementary ? 0x0800 : 0));
  if (manifest.policy.effective_exclude_flags != expected_effective) {
    throw FormatError(
        context + ": effective_exclude_flags disagrees with the flag policy");
  }

  const JsonValue& contigs = require_array(root, "contigs", context);
  std::unordered_set<std::string> contig_names;
  for (std::size_t i = 0; i < contigs.array.size(); ++i) {
    const JsonValue& item = require_object(contigs.array[i], "contigs[]");
    const uint64_t id = require_uint(item, "id", context + ".contigs");
    if (id != i)
      throw FormatError(context + ": contig IDs must be contiguous from zero");
    const uint64_t length = require_uint(item, "length", context + ".contigs");
    if (length > kMaxContigLength) {
      throw FormatError(context + ": contig length exceeds INT32_MAX");
    }
    const std::string name = require_string(item, "name", context + ".contigs");
    if (name.empty() || !detail::is_valid_utf8(name) ||
        name.find('\0') != std::string::npos ||
        name.find('\t') != std::string::npos ||
        name.find('\r') != std::string::npos ||
        name.find('\n') != std::string::npos) {
      throw FormatError(context + ": invalid contig name");
    }
    if (!contig_names.insert(name).second) {
      throw FormatError(context + ": duplicate contig name '" + name + "'");
    }
    manifest.contigs.push_back(
        ContigInfo{static_cast<uint32_t>(id), name, length});
  }

  // Scope has one or more intervals per contig. Build the name index once;
  // scanning all contigs per interval made load O(contigs^2) for references
  // with thousands of contigs.
  std::unordered_map<std::string, uint32_t> contig_ids;
  contig_ids.reserve(manifest.contigs.size());
  for (const ContigInfo& contig : manifest.contigs)
    contig_ids.emplace(contig.name, contig.id);

  const JsonValue& scope = require_array(root, "scope", context);
  std::vector<std::vector<ScopeInfo>> scope_by_contig(manifest.contigs.size());
  std::size_t previous_contig_id = 0;
  bool first_scope = true;
  for (const JsonValue& item_value : scope.array) {
    const JsonValue& item = require_object(item_value, "scope[]");
    ScopeInfo interval;
    interval.contig = require_string(item, "contig", context + ".scope");
    interval.start = require_uint(item, "start", context + ".scope");
    interval.end = require_uint(item, "end", context + ".scope");
    if (interval.start >= interval.end) {
      throw FormatError(context + ": scope interval is empty or reversed");
    }
    const auto found = contig_ids.find(interval.contig);
    if (found == contig_ids.end()) {
      throw FormatError(context + ": scope refers to unknown contig '" +
                        interval.contig + "'");
    }
    const uint32_t contig_id = found->second;
    if (interval.end > manifest.contigs.at(contig_id).length) {
      throw FormatError(context + ": scope interval ends beyond contig length");
    }
    if (!first_scope && contig_id < previous_contig_id) {
      throw FormatError(context + ": scope is not sorted by contig order");
    }
    if (!first_scope && contig_id == previous_contig_id) {
      const auto& list = scope_by_contig[contig_id];
      if (!list.empty() && list.back().end >= interval.start) {
        throw FormatError(
            context + ": scope intervals overlap or touch and are not merged");
      }
    }
    previous_contig_id = contig_id;
    first_scope = false;
    scope_by_contig[contig_id].push_back(interval);
    manifest.scope.push_back(interval);
  }
  for (std::size_t i = 0; i < manifest.contigs.size(); ++i) {
    if (scope_by_contig[i].empty()) {
      throw FormatError(context + ": indexed contig '" +
                        manifest.contigs[i].name + "' has no scope interval");
    }
  }

  const JsonValue& samples = require_array(root, "samples", context);
  std::unordered_set<std::string> sample_names;
  std::unordered_set<std::string> rg_ids;
  for (const JsonValue& item_value : samples.array) {
    const JsonValue& item = require_object(item_value, "samples[]");
    SampleInfo sample;
    sample.name = validate_name_field(
        require_string(item, "name", context + ".samples"), "sample");
    if (sample.name == "pooled") {
      throw FormatError(context + ": sample name 'pooled' is reserved");
    }
    if (!sample_names.insert(sample.name).second) {
      throw FormatError(context + ": duplicate sample name '" + sample.name +
                        "'");
    }
    sample.read_groups =
        string_array(item, "read_groups", context + ".samples");
    for (const std::string& rg : sample.read_groups) {
      if (rg.empty())
        throw FormatError(context + ": empty read-group ID");
      if (!rg_ids.insert(rg).second) {
        throw FormatError(context + ": duplicate read-group ID '" + rg + "'");
      }
    }
    manifest.samples.push_back(std::move(sample));
  }

  const JsonValue& units = require_array(root, "units", context);
  std::unordered_set<std::string> unit_names;
  for (std::size_t i = 0; i < units.array.size(); ++i) {
    const JsonValue& item = require_object(units.array[i], "units[]");
    UnitInfo unit;
    const uint64_t id = require_uint(item, "id", context + ".units");
    if (id != i)
      throw FormatError(context + ": unit IDs must be contiguous from zero");
    unit.id = static_cast<uint32_t>(id);
    unit.name = validate_name_field(
        require_string(item, "name", context + ".units"), "unit");
    if (!unit_names.insert(unit.name).second) {
      throw FormatError(context + ": duplicate unit name '" + unit.name + "'");
    }
    unit.kind =
        parse_kind(require_string(item, "kind", context + ".units"), context);
    unit.members = string_array(item, "members", context + ".units");
    unit.file = require_string(item, "file", context + ".units");
    if (!safe_relative_path(unit.file)) {
      throw FormatError(context + ": unit file path must be relative and safe");
    }
    unit.file_crc32 =
        parse_hex8(require_string(item, "file_crc32", context + ".units"),
                   context + ".units");
    unit.file_size = require_uint(item, "file_size", context + ".units");
    unit.total_runs = require_uint(item, "total_runs", context + ".units");
    unit.total_covered =
        require_uint(item, "total_covered", context + ".units");
    unit.contigs_with_runs =
        require_uint(item, "contigs_with_runs", context + ".units");
    const JsonValue* total_features = item.find("total_features");
    if (total_features == nullptr) {
      throw FormatError(context +
                        ": missing required key 'units[].total_features'");
    }
    if (!total_features->is_number() || !total_features->is_integer) {
      throw FormatError(context +
                        ": unit total_features must be an unsigned integer");
    }
    unit.total_features = total_features->integer;
    if (unit.kind == UnitInfo::Kind::Pooled) {
      if (unit.name != "pooled")
        throw FormatError(context + ": pooled unit must be named pooled");
    } else if (unit.name == "pooled") {
      throw FormatError(context + ": 'pooled' is reserved for the pooled unit");
    }
    if (unit.kind == UnitInfo::Kind::Sample) {
      if (unit.members.size() != 1 || unit.members.front() != unit.name) {
        throw FormatError(context + ": sample unit members must be [name]");
      }
      if (!sample_names.contains(unit.name)) {
        throw FormatError(context + ": sample unit '" + unit.name +
                          "' has no sample entry");
      }
    }
    if (unit.kind == UnitInfo::Kind::Group) {
      if (unit.members.empty()) {
        throw FormatError(context + ": group unit '" + unit.name +
                          "' has no members");
      }
      for (const std::string& member : unit.members) {
        if (!sample_names.contains(member)) {
          throw FormatError(context + ": group '" + unit.name +
                            "' has unknown member '" + member + "'");
        }
      }
    }
    manifest.units.push_back(std::move(unit));
  }

  // Layout consistency.
  bool has_pooled = false;
  for (const UnitInfo& unit : manifest.units) {
    if (unit.kind == UnitInfo::Kind::Pooled)
      has_pooled = true;
  }
  if (manifest.layout == "pooled" &&
      (!has_pooled || manifest.units.size() != 1)) {
    throw FormatError(context +
                      ": pooled layout must contain exactly one pooled unit");
  }
  if (manifest.layout == "per-sample" && has_pooled) {
    throw FormatError(context +
                      ": per-sample layout must not contain a pooled unit");
  }
  if (manifest.layout == "both" && !has_pooled) {
    throw FormatError(context + ": both layout must contain a pooled unit");
  }

  const JsonValue& files =
      require_object(require_field(root, "files", context), "files");
  for (std::size_t i = 0; i < files.keys.size(); ++i) {
    const std::string& path = files.keys[i];
    if (!safe_relative_path(path)) {
      throw FormatError(context +
                        ": files keys must be relative and safe paths");
    }
    const JsonValue& item = require_object(files.values[i], "files[]");
    FileInfo info;
    info.crc32 = parse_hex8(require_string(item, "crc32", context + ".files"),
                            context + ".files");
    info.size = require_uint(item, "size", context + ".files");
    manifest.files.emplace(path, info);
  }

  const JsonValue& index_value =
      require_object(require_field(root, "index", context), "index");
  manifest.index.contig_count =
      require_uint(index_value, "contig_count", context + ".index");
  const uint64_t meta_size =
      require_uint(index_value, "meta_record_size", context + ".index");
  if (meta_size != kMetaRecordSize) {
    throw FormatError(context + ": unsupported meta_record_size");
  }
  const uint64_t pack =
      require_uint(index_value, "pack_alignment", context + ".index");
  if (pack != kPackAlignment)
    throw FormatError(context + ": unsupported pack_alignment");
  manifest.index.meta_record_size = static_cast<uint32_t>(meta_size);
  manifest.index.pack_alignment = static_cast<uint32_t>(pack);
  manifest.index.partial_scope =
      require_bool(index_value, "partial_scope", context + ".index");
  manifest.index.total_runs =
      require_uint(index_value, "total_runs", context + ".index");
  manifest.index.unit_count =
      require_uint(index_value, "unit_count", context + ".index");
  if (manifest.index.contig_count != manifest.contigs.size()) {
    throw FormatError(context + ": index.contig_count disagrees with contigs");
  }
  if (manifest.index.unit_count != manifest.units.size()) {
    throw FormatError(context + ": index.unit_count disagrees with units");
  }
  std::unordered_map<std::string, uint64_t> sample_features;
  for (const UnitInfo& unit : manifest.units) {
    if (unit.kind == UnitInfo::Kind::Sample)
      sample_features.emplace(unit.name, unit.total_features);
  }
  for (const UnitInfo& unit : manifest.units) {
    if (unit.kind == UnitInfo::Kind::Sample || unit.members.empty())
      continue;
    uint64_t expected = 0;
    for (const std::string& member : unit.members) {
      const auto found = sample_features.find(member);
      if (found == sample_features.end()) {
        throw FormatError(context + ": unit '" + unit.name +
                          "' has no feature-count sample member '" + member +
                          "'");
      }
      expected += found->second;
    }
    if (unit.total_features != expected) {
      throw FormatError(context + ": unit '" + unit.name +
                        "' total_features disagrees with its members");
    }
  }
  uint64_t total_runs = 0;
  for (const UnitInfo& unit : manifest.units) {
    const auto file = manifest.files.find(unit.file);
    if (file == manifest.files.end()) {
      throw FormatError(context + ": unit file '" + unit.file +
                        "' is not listed in files");
    }
    if (file->second.size != unit.file_size ||
        file->second.crc32 != unit.file_crc32) {
      throw FormatError(context +
                        ": unit file metadata disagrees with files for '" +
                        unit.file + "'");
    }
    if (unit.contigs_with_runs > manifest.contigs.size()) {
      throw FormatError(context +
                        ": unit contigs_with_runs exceeds the contig count");
    }
    total_runs += unit.total_runs;
  }
  if (total_runs != manifest.index.total_runs) {
    throw FormatError(context + ": index.total_runs disagrees with the units");
  }

  bool partial = false;
  for (std::size_t i = 0; i < manifest.contigs.size(); ++i) {
    const auto& intervals = scope_by_contig[i];
    if (intervals.size() != 1 || intervals.front().start != 0 ||
        intervals.front().end != manifest.contigs[i].length) {
      partial = true;
      break;
    }
  }
  if (manifest.index.partial_scope != partial) {
    throw FormatError(
        context + ": index.partial_scope disagrees with the scope intervals");
  }

  const JsonValue& semantics =
      require_object(require_field(root, "semantics", context), "semantics");
  manifest.semantics_covered =
      require_string(semantics, "covered", context + ".semantics");
  if (manifest.semantics_covered.empty()) {
    throw FormatError(context + ": semantics.covered must not be empty");
  }
  const JsonValue* feature_semantics = semantics.find("features");
  if (feature_semantics == nullptr) {
    throw FormatError(context + ": missing required key 'semantics.features'");
  }
  if (!feature_semantics->is_string() || feature_semantics->text.empty()) {
    throw FormatError(context +
                      ": semantics.features must be a non-empty string");
  }
  manifest.semantics_features = feature_semantics->text;

  const JsonValue& warnings = require_array(root, "warnings", context);
  for (const JsonValue& item : warnings.array) {
    if (!item.is_string())
      throw FormatError(context + ": warnings must be strings");
    manifest.warnings.push_back(item.text);
  }
  return manifest;
}

void Manifest::save(const std::string& path) const {
  std::ofstream out(path, std::ios::binary);
  if (!out)
    throw OutputError("cannot open manifest for writing: '" + path + "'");
  const std::string text = to_json();
  out.write(text.data(), static_cast<std::streamsize>(text.size()));
  if (!out)
    throw OutputError("write error on manifest '" + path + "'");
  out.close();
  if (!out)
    throw OutputError("close error on manifest '" + path + "'");
}

Manifest Manifest::load(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw IoError("cannot open manifest '" + path + "'");
  const std::string text((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
  if (in.bad())
    throw IoError("read error on manifest '" + path + "'");
  return from_json(text, path);
}

std::string created_from_env() {
  uint64_t epoch = 0;
  const char* value = std::getenv("SOURCE_DATE_EPOCH");
  if (value == nullptr) {
    const std::time_t now = std::time(nullptr);
    if (now < 0)
      throw UsageError("system clock is before the Unix epoch");
    epoch = static_cast<uint64_t>(now);
  } else {
    const std::string text(value);
    if (text.empty())
      throw UsageError("SOURCE_DATE_EPOCH is empty");
    const auto parsed =
        std::from_chars(text.data(), text.data() + text.size(), epoch);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()) {
      throw UsageError("SOURCE_DATE_EPOCH must be a non-negative integer");
    }
  }
  if (epoch > 253402300799ULL) {
    throw UsageError(
        "SOURCE_DATE_EPOCH is outside the representable timestamp range");
  }
  std::time_t seconds = static_cast<std::time_t>(epoch);
  std::tm tm{};
  if (::gmtime_r(&seconds, &tm) == nullptr) {
    throw UsageError("SOURCE_DATE_EPOCH cannot be converted to UTC");
  }
  char buffer[32];
  const int written =
      std::snprintf(buffer, sizeof(buffer), "%04d-%02d-%02dT%02d:%02d:%02dZ",
                    tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,
                    tm.tm_min, tm.tm_sec);
  if (written <= 0 || static_cast<std::size_t>(written) >= sizeof(buffer)) {
    throw UsageError("failed to format the timestamp");
  }
  return std::string(buffer, static_cast<std::size_t>(written));
}

} // namespace covidx
