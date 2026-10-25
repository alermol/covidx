#pragma once
//
// Canonical JSON value, parser, and writer used by manifest.json.
//
// The parser rejects duplicate object keys and trailing content. The writer
// sorts object keys by unsigned byte order, indents with two spaces, uses LF,
// and terminates the document with exactly one newline.
//
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace covidx::detail {

struct JsonValue {
  enum class Type { Null, Bool, Number, String, Array, Object };

  Type type = Type::Null;
  bool boolean = false;
  bool is_integer = false;
  uint64_t integer = 0;
  double number = 0.0;
  std::string text;
  std::vector<JsonValue> array;
  // Object members are parallel key/value arrays kept in insertion order.
  std::vector<std::string> keys;
  std::vector<JsonValue> values;

  [[nodiscard]] const JsonValue* find(std::string_view key) const;
  [[nodiscard]] bool is_object() const {
    return type == Type::Object;
  }
  [[nodiscard]] bool is_array() const {
    return type == Type::Array;
  }
  [[nodiscard]] bool is_string() const {
    return type == Type::String;
  }
  [[nodiscard]] bool is_number() const {
    return type == Type::Number;
  }
  [[nodiscard]] bool is_bool() const {
    return type == Type::Bool;
  }
  [[nodiscard]] bool is_null() const {
    return type == Type::Null;
  }

  static JsonValue make_string(std::string value);
  static JsonValue make_uint(uint64_t value);
  static JsonValue make_bool(bool value);
  static JsonValue make_null();
  static JsonValue make_object();
  static JsonValue make_array();

  // Adds or replaces an object member. Duplicate keys are rejected when
  // parsing; programmatic construction replaces the existing value.
  void set(const std::string& key, JsonValue value);
};

// Parses one JSON document. Throws FormatError on malformed input.
[[nodiscard]] JsonValue parse_json(std::string_view text,
                                   std::string_view origin);

// Serializes a value canonically, including exactly one trailing newline.
[[nodiscard]] std::string to_canonical_json(const JsonValue& value);

} // namespace covidx::detail
