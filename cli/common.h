#pragma once
//
// Small hand-rolled CLI helpers shared by covidx-build and covidx-query.
//
#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include "covidx/cov_index.h"
#include "covidx/error.h"
#include "covidx/format.h"

namespace cli {

class Args {
public:
  Args(int argc, char** argv, std::set<std::string> value_options,
       std::set<std::string> flag_options,
       std::set<std::string> repeatable = {})
      : value_options_(std::move(value_options)),
        flag_options_(std::move(flag_options)),
        repeatable_(std::move(repeatable)) {
    for (int i = 1; i < argc; ++i) {
      std::string token = argv[i];
      if (!token.starts_with("--"))
        throw covidx::UsageError("unexpected argument '" + token + "'");
      const std::size_t equals = token.find('=');
      std::string name = token.substr(
          2, equals == std::string::npos ? std::string::npos : equals - 2);
      std::string inline_value;
      const bool has_inline = equals != std::string::npos;
      if (has_inline)
        inline_value = token.substr(equals + 1);
      if (value_options_.contains(name)) {
        if (!repeatable_.contains(name) && values_.contains(name)) {
          throw covidx::UsageError("option --" + name +
                                   " was given more than once");
        }
        std::string value;
        if (has_inline) {
          value = std::move(inline_value);
        } else {
          if (i + 1 >= argc)
            throw covidx::UsageError("option --" + name + " needs a value");
          value = argv[++i];
        }
        values_[name].push_back(std::move(value));
      } else if (flag_options_.contains(name)) {
        if (has_inline)
          throw covidx::UsageError("option --" + name +
                                   " does not take a value");
        flags_.insert(name);
      } else {
        throw covidx::UsageError("unknown option --" + name);
      }
    }
  }

  [[nodiscard]] bool has(std::string_view name) const {
    return flags_.contains(std::string(name)) ||
           values_.contains(std::string(name));
  }
  [[nodiscard]] bool flag(std::string_view name) const {
    return flags_.contains(std::string(name));
  }
  [[nodiscard]] std::string get(std::string_view name,
                                std::string fallback = {}) const {
    const auto it = values_.find(std::string(name));
    if (it == values_.end() || it->second.empty())
      return fallback;
    return it->second.back();
  }
  [[nodiscard]] std::vector<std::string> all(std::string_view name) const {
    const auto it = values_.find(std::string(name));
    return it == values_.end() ? std::vector<std::string>{} : it->second;
  }

private:
  std::set<std::string> value_options_;
  std::set<std::string> flag_options_;
  std::set<std::string> repeatable_;
  std::map<std::string, std::vector<std::string>> values_;
  std::set<std::string> flags_;
};

[[nodiscard]] inline uint64_t parse_u64(const std::string& text,
                                        const std::string& option) {
  if (text.empty())
    throw covidx::UsageError("option " + option + " needs an integer value");
  uint64_t value = 0;
  const auto result =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    throw covidx::UsageError("option " + option +
                             " has an invalid integer value '" + text + "'");
  }
  return value;
}

[[nodiscard]] inline uint64_t parse_flags(const std::string& text,
                                          const std::string& option) {
  if (text.empty())
    throw covidx::UsageError("option " + option + " needs a value");
  int base = 10;
  std::string_view digits = text;
  if (digits.size() > 2 && digits[0] == '0' &&
      (digits[1] == 'x' || digits[1] == 'X')) {
    base = 16;
    digits.remove_prefix(2);
  }
  if (digits.empty() || digits.front() == '-' || digits.front() == '+') {
    throw covidx::UsageError("option " + option + " has an invalid value '" +
                             text + "'");
  }
  uint64_t value = 0;
  const auto result = std::from_chars(
      digits.data(), digits.data() + digits.size(), value, base);
  if (result.ec != std::errc{} || result.ptr != digits.data() + digits.size()) {
    throw covidx::UsageError("option " + option + " has an invalid value '" +
                             text + "'");
  }
  if (value > 0xFFFF) {
    throw covidx::UsageError("option " + option +
                             " must be in the range 0..0xFFFF");
  }
  return value;
}

[[nodiscard]] inline std::vector<std::string>
split_csv(const std::string& text) {
  std::vector<std::string> parts;
  std::size_t start = 0;
  while (start <= text.size()) {
    const std::size_t comma = text.find(',', start);
    parts.push_back(text.substr(
        start, comma == std::string::npos ? std::string::npos : comma - start));
    if (comma == std::string::npos)
      break;
    start = comma + 1;
  }
  return parts;
}

[[nodiscard]] inline std::vector<std::string>
read_lines(const std::string& path, const std::string& what) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw covidx::IoError("cannot open " + what + " '" + path + "'");
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    if (line.empty())
      continue;
    lines.push_back(std::move(line));
  }
  if (in.bad())
    throw covidx::IoError("read error on " + what + " '" + path + "'");
  return lines;
}

[[nodiscard]] inline unsigned default_threads() {
  const unsigned hardware = std::thread::hardware_concurrency();
  return std::min<unsigned>(8, hardware == 0 ? 1 : hardware);
}

[[nodiscard]] inline int exception_exit_code(const std::exception& error) {
  if (dynamic_cast<const covidx::UsageError*>(&error) != nullptr)
    return 1;
  if (dynamic_cast<const covidx::OutputError*>(&error) != nullptr)
    return 3;
  if (dynamic_cast<const covidx::InputError*>(&error) != nullptr)
    return 2;
  if (dynamic_cast<const covidx::FormatError*>(&error) != nullptr)
    return 2;
  if (dynamic_cast<const covidx::QueryError*>(&error) != nullptr)
    return 2;
  if (dynamic_cast<const covidx::IoError*>(&error) != nullptr)
    return 2;
  return 4;
}

[[nodiscard]] inline std::string
exception_category(const std::exception& error) {
  if (dynamic_cast<const covidx::UsageError*>(&error) != nullptr)
    return "usage error";
  if (dynamic_cast<const covidx::OutputError*>(&error) != nullptr)
    return "output error";
  if (dynamic_cast<const covidx::InputError*>(&error) != nullptr)
    return "input error";
  if (dynamic_cast<const covidx::FormatError*>(&error) != nullptr)
    return "format error";
  if (dynamic_cast<const covidx::QueryError*>(&error) != nullptr)
    return "query error";
  if (dynamic_cast<const covidx::IoError*>(&error) != nullptr)
    return "I/O error";
  return "internal error";
}

} // namespace cli
