#include "covidx/detail/json.h"

#include <algorithm>
#include <charconv>
#include <format>
#include <limits>
#include <utility>

#include "covidx/error.h"

namespace covidx::detail {

const JsonValue* JsonValue::find(std::string_view key) const {
  for (std::size_t i = 0; i < values.size(); ++i)
    if (keys[i] == key)
      return &values[i];
  return nullptr;
}

JsonValue JsonValue::make_string(std::string value) {
  JsonValue v;
  v.type = Type::String;
  v.text = std::move(value);
  return v;
}

JsonValue JsonValue::make_uint(uint64_t value) {
  JsonValue v;
  v.type = Type::Number;
  v.is_integer = true;
  v.integer = value;
  return v;
}

JsonValue JsonValue::make_bool(bool value) {
  JsonValue v;
  v.type = Type::Bool;
  v.boolean = value;
  return v;
}

JsonValue JsonValue::make_null() {
  return JsonValue{};
}

JsonValue JsonValue::make_object() {
  JsonValue v;
  v.type = Type::Object;
  return v;
}

JsonValue JsonValue::make_array() {
  JsonValue v;
  v.type = Type::Array;
  return v;
}

void JsonValue::set(const std::string& key, JsonValue value) {
  for (std::size_t i = 0; i < keys.size(); ++i) {
    if (keys[i] == key) {
      values[i] = std::move(value);
      return;
    }
  }
  keys.push_back(key);
  values.push_back(std::move(value));
}

namespace {

class JsonParser {
public:
  JsonParser(std::string_view text, std::string origin)
      : text_(text), origin_(std::move(origin)) {}

  JsonValue parse_document() {
    skip_whitespace();
    JsonValue value = parse_value();
    skip_whitespace();
    if (pos_ != text_.size())
      fail("trailing content after the top-level value");
    return value;
  }

private:
  [[noreturn]] void fail(const std::string& message) const {
    throw FormatError("malformed JSON in " + origin_ + " at byte " +
                      std::to_string(pos_) + ": " + message);
  }

  struct DepthGuard {
    std::size_t& depth;
    ~DepthGuard() {
      --depth;
    }
  };

  void skip_whitespace() {
    while (pos_ < text_.size()) {
      const char c = text_[pos_];
      if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
        ++pos_;
      else
        break;
    }
  }

  char peek() {
    if (pos_ >= text_.size())
      fail("unexpected end of input");
    return text_[pos_];
  }

  void expect(char c) {
    if (peek() != c)
      fail(std::string("expected '") + c + "'");
    ++pos_;
  }

  void expect_literal(std::string_view literal) {
    if (text_.substr(pos_, literal.size()) != literal)
      fail("expected '" + std::string(literal) + "'");
    pos_ += literal.size();
  }

  JsonValue parse_value() {
    switch (peek()) {
    case '{':
      return parse_object();
    case '[':
      return parse_array();
    case '"': {
      JsonValue v;
      v.type = JsonValue::Type::String;
      v.text = parse_string();
      return v;
    }
    case 't':
      expect_literal("true");
      return JsonValue::make_bool(true);
    case 'f':
      expect_literal("false");
      return JsonValue::make_bool(false);
    case 'n':
      expect_literal("null");
      return JsonValue();
    default:
      return parse_number();
    }
  }

  JsonValue parse_object() {
    if (depth_ >= 64)
      fail("nesting is too deep");
    DepthGuard guard{depth_};
    ++depth_;
    JsonValue v = JsonValue::make_object();
    expect('{');
    skip_whitespace();
    if (peek() == '}') {
      ++pos_;
      return v;
    }
    for (;;) {
      skip_whitespace();
      std::string key = parse_string();
      if (v.find(key) != nullptr)
        fail("duplicate object key '" + key + "'");
      skip_whitespace();
      expect(':');
      skip_whitespace();
      JsonValue member = parse_value();
      v.keys.push_back(std::move(key));
      v.values.push_back(std::move(member));
      skip_whitespace();
      const char c = peek();
      if (c == ',') {
        ++pos_;
        continue;
      }
      if (c == '}') {
        ++pos_;
        return v;
      }
      fail("expected ',' or '}' in object");
    }
  }

  JsonValue parse_array() {
    if (depth_ >= 64)
      fail("nesting is too deep");
    DepthGuard guard{depth_};
    ++depth_;
    JsonValue v = JsonValue::make_array();
    expect('[');
    skip_whitespace();
    if (peek() == ']') {
      ++pos_;
      return v;
    }
    for (;;) {
      v.array.push_back(parse_value());
      skip_whitespace();
      const char c = peek();
      if (c == ',') {
        ++pos_;
        skip_whitespace();
        continue;
      }
      if (c == ']') {
        ++pos_;
        return v;
      }
      fail("expected ',' or ']' in array");
    }
  }

  unsigned parse_hex4() {
    if (pos_ + 4 > text_.size())
      fail("truncated \\u escape");
    unsigned value = 0;
    for (int i = 0; i < 4; ++i) {
      const char c = text_[pos_++];
      unsigned digit = 0;
      if (c >= '0' && c <= '9')
        digit = static_cast<unsigned>(c - '0');
      else if (c >= 'a' && c <= 'f')
        digit = static_cast<unsigned>(c - 'a' + 10);
      else if (c >= 'A' && c <= 'F')
        digit = static_cast<unsigned>(c - 'A' + 10);
      else
        fail("invalid hex digit in \\u escape");
      value = (value << 4) | digit;
    }
    return value;
  }

  static void append_utf8(std::string& out, unsigned cp) {
    if (cp < 0x80) {
      out.push_back(static_cast<char>(cp));
    } else if (cp < 0x800) {
      out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else if (cp < 0x10000) {
      out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else {
      out.push_back(static_cast<char>(0xF0 | (cp >> 18)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    }
  }

  std::string parse_string() {
    expect('"');
    std::string out;
    for (;;) {
      if (pos_ >= text_.size())
        fail("unterminated string");
      const auto c = static_cast<unsigned char>(text_[pos_++]);
      if (c == '"')
        return out;
      if (c < 0x20)
        fail("unescaped control character in string");
      if (c != '\\') {
        out.push_back(static_cast<char>(c));
        continue;
      }
      if (pos_ >= text_.size())
        fail("unterminated escape");
      const char e = text_[pos_++];
      switch (e) {
      case '"':
        out.push_back('"');
        break;
      case '\\':
        out.push_back('\\');
        break;
      case '/':
        out.push_back('/');
        break;
      case 'b':
        out.push_back('\b');
        break;
      case 'f':
        out.push_back('\f');
        break;
      case 'n':
        out.push_back('\n');
        break;
      case 'r':
        out.push_back('\r');
        break;
      case 't':
        out.push_back('\t');
        break;
      case 'u': {
        unsigned cp = parse_hex4();
        if (cp >= 0xDC00 && cp <= 0xDFFF)
          fail("lone low surrogate in \\u escape");
        if (cp >= 0xD800 && cp <= 0xDBFF) {
          if (pos_ + 1 < text_.size() && text_[pos_] == '\\' &&
              text_[pos_ + 1] == 'u') {
            pos_ += 2;
            const unsigned low = parse_hex4();
            if (low >= 0xDC00 && low <= 0xDFFF)
              cp = 0x10000 + ((cp - 0xD800) << 10) + (low - 0xDC00);
            else
              fail("invalid low surrogate in \\u escape");
          } else {
            fail("missing low surrogate in \\u escape");
          }
        }
        append_utf8(out, cp);
        break;
      }
      default:
        fail("unknown escape sequence");
      }
    }
  }

  JsonValue parse_number() {
    const std::size_t start = pos_;
    if (pos_ < text_.size() && text_[pos_] == '-')
      ++pos_;
    if (pos_ >= text_.size() || text_[pos_] < '0' || text_[pos_] > '9')
      fail("invalid value");
    const std::size_t first_digit = pos_;
    while (pos_ < text_.size() && text_[pos_] >= '0' && text_[pos_] <= '9')
      ++pos_;
    if (pos_ - first_digit > 1 && text_[first_digit] == '0')
      fail("leading zeros are not allowed in numbers");
    bool is_integer = true;
    if (pos_ < text_.size() && text_[pos_] == '.') {
      is_integer = false;
      ++pos_;
      if (pos_ >= text_.size() || text_[pos_] < '0' || text_[pos_] > '9')
        fail("expected a digit after '.'");
      while (pos_ < text_.size() && text_[pos_] >= '0' && text_[pos_] <= '9')
        ++pos_;
    }
    if (pos_ < text_.size() && (text_[pos_] == 'e' || text_[pos_] == 'E')) {
      is_integer = false;
      ++pos_;
      if (pos_ < text_.size() && (text_[pos_] == '+' || text_[pos_] == '-'))
        ++pos_;
      if (pos_ >= text_.size() || text_[pos_] < '0' || text_[pos_] > '9')
        fail("expected a digit in the exponent");
      while (pos_ < text_.size() && text_[pos_] >= '0' && text_[pos_] <= '9')
        ++pos_;
    }
    const std::string_view token = text_.substr(start, pos_ - start);
    JsonValue v;
    v.type = JsonValue::Type::Number;
    const char* first = token.data();
    const char* last = first + token.size();
    if (is_integer && token[0] != '-') {
      uint64_t value = 0;
      const std::from_chars_result parsed = std::from_chars(first, last, value);
      if (parsed.ec == std::errc{} && parsed.ptr == last) {
        v.is_integer = true;
        v.integer = value;
        return v;
      }
    }
    double number = 0.0;
    const std::from_chars_result parsed =
        std::from_chars(first, last, number, std::chars_format::general);
    if (parsed.ec == std::errc::result_out_of_range)
      number = token.starts_with('-') ? -std::numeric_limits<double>::infinity()
                                      : std::numeric_limits<double>::infinity();
    v.is_integer = false;
    v.number = number;
    return v;
  }

  const std::string_view text_;
  std::string origin_;
  std::size_t pos_ = 0;
  std::size_t depth_ = 0;
};

void write_json_string(std::string_view s, std::string& out) {
  out.push_back('"');
  for (const char c : s) {
    const auto u = static_cast<unsigned char>(c);
    switch (u) {
    case '"':
      out += "\\\"";
      break;
    case '\\':
      out += "\\\\";
      break;
    case '\b':
      out += "\\b";
      break;
    case '\f':
      out += "\\f";
      break;
    case '\n':
      out += "\\n";
      break;
    case '\r':
      out += "\\r";
      break;
    case '\t':
      out += "\\t";
      break;
    default:
      if (u < 0x20) {
        std::format_to(std::back_inserter(out), "\\u{:04x}", u);
      } else {
        out.push_back(c);
      }
    }
  }
  out.push_back('"');
}

void write_json(const JsonValue& value, std::string& out, int depth) {
  const std::string pad(static_cast<std::size_t>(depth) * 2, ' ');
  const std::string pad_inner(static_cast<std::size_t>(depth + 1) * 2, ' ');
  switch (value.type) {
  case JsonValue::Type::Null:
    out += "null";
    break;
  case JsonValue::Type::Bool:
    out += value.boolean ? "true" : "false";
    break;
  case JsonValue::Type::Number:
    if (value.is_integer) {
      out += std::to_string(value.integer);
    } else {
      std::format_to(std::back_inserter(out), "{:.17g}", value.number);
    }
    break;
  case JsonValue::Type::String:
    write_json_string(value.text, out);
    break;
  case JsonValue::Type::Array: {
    if (value.array.empty()) {
      out += "[]";
      break;
    }
    out += "[\n";
    for (std::size_t i = 0; i < value.array.size(); ++i) {
      out += pad_inner;
      write_json(value.array[i], out, depth + 1);
      out += (i + 1 < value.array.size()) ? ",\n" : "\n";
    }
    out += pad;
    out += "]";
    break;
  }
  case JsonValue::Type::Object: {
    if (value.values.empty()) {
      out += "{}";
      break;
    }
    std::vector<std::size_t> order(value.values.size());
    for (std::size_t i = 0; i < order.size(); ++i)
      order[i] = i;
    std::ranges::sort(order, std::less<>{},
                      [&value](std::size_t i) -> const std::string& {
                        return value.keys[i];
                      });
    out += "{\n";
    for (std::size_t i = 0; i < order.size(); ++i) {
      const std::size_t at = order[i];
      out += pad_inner;
      write_json_string(value.keys[at], out);
      out += ": ";
      write_json(value.values[at], out, depth + 1);
      out += (i + 1 < order.size()) ? ",\n" : "\n";
    }
    out += pad;
    out += "}";
    break;
  }
  }
}

} // namespace

JsonValue parse_json(std::string_view text, std::string_view origin) {
  return JsonParser(text, std::string(origin)).parse_document();
}

std::string to_canonical_json(const JsonValue& value) {
  std::string out;
  write_json(value, out, 0);
  out.push_back('\n');
  return out;
}

} // namespace covidx::detail
