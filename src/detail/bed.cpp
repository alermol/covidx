#include "covidx/detail/bed.h"

#include <charconv>
#include <fstream>
#include <iostream>
#include <string_view>

#include "covidx/error.h"

namespace covidx::detail {

namespace {

[[nodiscard]] uint64_t parse_uint(const std::string& text,
                                  const std::string& context) {
  if (text.empty())
    throw InputError(context + ": empty integer field");
  uint64_t value = 0;
  const auto result =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    throw InputError(context + ": invalid unsigned integer '" + text + "'");
  }
  return value;
}

[[nodiscard]] std::vector<std::string> split_tabs(const std::string& line) {
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

} // namespace

std::vector<BedRecord> read_bed_stream(std::istream& in,
                                       const std::string& origin) {
  std::vector<BedRecord> records;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(in, line)) {
    ++line_number;
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    if (line.empty())
      continue;
    if (line.starts_with("#") || line.starts_with("track") ||
        line.starts_with("browser")) {
      continue;
    }
    const std::vector<std::string> fields = split_tabs(line);
    if (fields.size() < 3) {
      throw InputError(origin + ":" + std::to_string(line_number) +
                       ": BED row needs at least three tab-separated fields");
    }
    if (fields[0].empty()) {
      throw InputError(origin + ":" + std::to_string(line_number) +
                       ": empty contig name");
    }
    BedRecord record;
    record.contig = fields[0];
    record.start =
        parse_uint(fields[1], origin + ":" + std::to_string(line_number));
    record.end =
        parse_uint(fields[2], origin + ":" + std::to_string(line_number));
    if (record.start >= record.end) {
      throw InputError(origin + ":" + std::to_string(line_number) +
                       ": zero-length or reversed BED interval");
    }
    record.line = line_number;
    records.push_back(std::move(record));
  }
  if (in.bad())
    throw IoError("read error on BED input '" + origin + "'");
  return records;
}

std::vector<BedRecord> read_bed(const std::string& path) {
  if (path == "-")
    return read_bed_stream(std::cin, "stdin");
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw IoError("cannot open BED file '" + path + "'");
  return read_bed_stream(in, path);
}

} // namespace covidx::detail
