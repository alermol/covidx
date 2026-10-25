#include "doctest.h"

#include <string>

#include "covidx/detail/json.h"
#include "covidx/error.h"

using covidx::detail::JsonValue;
using covidx::detail::parse_json;
using covidx::detail::to_canonical_json;

TEST_CASE("JSON parses scalar and nested values") {
  const JsonValue v =
      parse_json(R"({"a":1,"b":[true,null,"x"],"c":{"d":false}})", "test");
  REQUIRE(v.is_object());
  REQUIRE(v.find("a") != nullptr);
  CHECK(v.find("a")->integer == 1);
  CHECK(v.find("b")->array.size() == 3);
  CHECK(v.find("c")->find("d")->boolean == false);
}

TEST_CASE("JSON rejects trailing content and duplicate keys") {
  CHECK_THROWS_AS((void)parse_json("{} garbage", "test"), covidx::FormatError);
  CHECK_THROWS_AS((void)parse_json(R"({"a":1,"a":2})", "test"),
                  covidx::FormatError);
  CHECK_THROWS_AS((void)parse_json(R"({"a":})", "test"), covidx::FormatError);
  CHECK_THROWS_AS((void)parse_json(R"({"a":"\q"})", "test"),
                  covidx::FormatError);
  CHECK_THROWS_AS((void)parse_json("\"unterminated", "test"),
                  covidx::FormatError);
}

TEST_CASE(
    "Canonical JSON sorts keys, indents with two spaces, and ends with LF") {
  JsonValue root = JsonValue::make_object();
  root.set("z", JsonValue::make_uint(1));
  root.set("a", JsonValue::make_string("v"));
  root.set("m", JsonValue::make_bool(true));
  const std::string text = to_canonical_json(root);
  CHECK(text == "{\n"
                "  \"a\": \"v\",\n"
                "  \"m\": true,\n"
                "  \"z\": 1\n"
                "}\n");
}

TEST_CASE("Canonical JSON round-trips") {
  JsonValue root = JsonValue::make_object();
  JsonValue array = JsonValue::make_array();
  array.array.push_back(JsonValue::make_uint(0));
  array.array.push_back(JsonValue::make_string("a\nb"));
  root.set("array", std::move(array));
  root.set("empty", JsonValue::make_object());
  const std::string text = to_canonical_json(root);
  const JsonValue parsed = parse_json(text, "roundtrip");
  CHECK(parsed.find("array")->array.at(1).text == "a\nb");
  CHECK(parsed.find("empty")->values.empty());
  CHECK(to_canonical_json(parsed) == text);
}

TEST_CASE("JSON string escapes are canonical") {
  JsonValue root = JsonValue::make_object();
  root.set("s", JsonValue::make_string(std::string("a\tb\001c")));
  CHECK(to_canonical_json(root) == "{\n  \"s\": \"a\\tb\\u0001c\"\n}\n");
}
