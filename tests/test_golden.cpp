#include "doctest.h"

#include <cstdint>
#include <filesystem>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "covidx/detail/sha256.h"
#include "test_util.h"

namespace fs = std::filesystem;

namespace {

fs::path source_dir() {
  return fs::path(COVIDX_SOURCE_DIR);
}

int build_golden(const fs::path& out, unsigned threads) {
  const std::string command =
      "cd " + source_dir().string() + " && SOURCE_DATE_EPOCH=0 " +
      COVIDX_BUILD_BIN + " --bam tests/data/tiny.bam --out '" + out.string() +
      "' --threads " + std::to_string(threads) + " >/dev/null 2>&1";
  return tests::run_command_status(command);
}

std::map<std::string, std::string> snapshot(const fs::path& root) {
  std::map<std::string, std::string> files;
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file())
      continue;
    files[fs::relative(entry.path(), root).generic_string()] =
        tests::read_file(entry.path());
  }
  return files;
}

std::map<std::string, std::string> load_golden_hashes() {
  const fs::path path = fs::path(COVIDX_TEST_DATA_DIR) / "golden_index.sha256";
  const std::string text = tests::read_file(path);
  std::map<std::string, std::string> hashes;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    REQUIRE(line.size() > 66);
    const std::string hash = line.substr(0, 64);
    const std::string file = line.substr(66);
    hashes.emplace(file, hash);
  }
  return hashes;
}

} // namespace

TEST_CASE("golden index is byte-identical for 1/2/8 threads") {
  const auto hashes = load_golden_hashes();
  const auto expected_files = [&] {
    std::vector<std::string> names;
    for (const auto& [name, hash] : hashes) {
      (void)hash;
      names.push_back(name);
    }
    return names;
  }();

  tests::TempDir temp;
  std::map<std::string, std::string> first;
  for (unsigned threads : {1u, 2u, 8u}) {
    const fs::path out = temp.path() / ("out-" + std::to_string(threads));
    REQUIRE(build_golden(out, threads) == 0);
    const fs::path index = out / "tiny.covidx";
    const auto files = snapshot(index);
    REQUIRE(files.size() == hashes.size());
    for (const std::string& name : expected_files) {
      REQUIRE(files.contains(name));
      CHECK(covidx::detail::sha256_hex(files.at(name).data(),
                                       files.at(name).size()) ==
            hashes.at(name));
    }
    if (first.empty()) {
      first = files;
    } else {
      CHECK(files == first);
    }
  }
}
