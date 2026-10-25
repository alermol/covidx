#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <sys/stat.h>

#include "covidx/cov_index.h"
#include "covidx/error.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::BuildOptions;
using covidx::CovIndex;

namespace {

fs::path tiny_bam() {
  return fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
}

fs::path find_work_dir(const fs::path& parent) {
  for (const auto& entry : fs::directory_iterator(parent)) {
    if (entry.is_directory() &&
        entry.path().filename().string().starts_with("covidx-work-")) {
      return entry.path();
    }
  }
  return {};
}

ino_t inode_of(const fs::path& path) {
  struct stat st = {};
  REQUIRE(::stat(path.c_str(), &st) == 0);
  return st.st_ino;
}

std::vector<std::string> relative_files(const fs::path& root) {
  std::vector<std::string> files;
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (entry.is_regular_file())
      files.push_back(fs::relative(entry.path(), root).generic_string());
  }
  return files;
}

} // namespace

TEST_CASE("save uses hardlinks on the same filesystem") {
  tests::TempDir temp;
  const fs::path work_parent = temp.path() / "work";
  fs::create_directories(work_parent);
  BuildOptions options;
  options.work_dir = work_parent.string();
  options.keep_work = true;
  CovIndex index;
  index.build(tiny_bam().string(), options);
  index.save(temp.path().string(), "tiny");

  const fs::path work = find_work_dir(work_parent);
  REQUIRE(!work.empty());
  const fs::path source = work / "data" / "u000000.bin";
  const fs::path target = temp.path() / "tiny.covidx" / "data" / "u000000.bin";
  REQUIRE(fs::exists(source));
  REQUIRE(fs::exists(target));
  CHECK(inode_of(source) == inode_of(target));

  const std::vector<std::string> files =
      relative_files(temp.path() / "tiny.covidx");
  CHECK(std::find(files.begin(), files.end(), "manifest.json") != files.end());
  CHECK(std::find(files.begin(), files.end(), "contigs.tsv") != files.end());
  CHECK(std::find(files.begin(), files.end(), "contigs.tsv.crc32") !=
        files.end());
  CHECK(std::find(files.begin(), files.end(), "data/u000000.bin") !=
        files.end());
  CHECK(std::find(files.begin(), files.end(), "data/u000000.bin.crc32") !=
        files.end());
  CHECK(std::find(files.begin(), files.end(), "spools") == files.end());
}

TEST_CASE("save falls back to copying across filesystems") {
  if (!fs::exists("/dev/shm")) {
    MESSAGE("/dev/shm is not available; skipping cross-filesystem save test");
    return;
  }
  tests::TempDir temp;
  const fs::path work_parent = temp.path() / "work";
  fs::create_directories(work_parent);
  BuildOptions options;
  options.work_dir = work_parent.string();
  options.keep_work = true;
  CovIndex index;
  index.build(tiny_bam().string(), options);

  const fs::path other =
      fs::path("/dev/shm") / ("covidx-save-" + std::to_string(::getpid()));
  fs::remove_all(other);
  fs::create_directories(other);
  try {
    index.save(other.string(), "tiny");
    const fs::path work = find_work_dir(work_parent);
    REQUIRE(!work.empty());
    CHECK(inode_of(work / "data" / "u000000.bin") !=
          inode_of(other / "tiny.covidx" / "data" / "u000000.bin"));
    CovIndex loaded;
    loaded.load((other / "tiny.covidx").string());
    CHECK(loaded.covered_count("S1", "chr1", 0, 1000) == 60);
  } catch (...) {
    fs::remove_all(other);
    throw;
  }
  fs::remove_all(other);
}

TEST_CASE("failed save leaves no temp directory") {
  tests::TempDir temp;
  BuildOptions options;
  options.work_dir = temp.path().string();
  options.keep_work = true;
  CovIndex index;
  index.build(tiny_bam().string(), options);
  const fs::path work = find_work_dir(temp.path());
  REQUIRE(!work.empty());
  {
    std::ofstream out(work / "units.tsv", std::ios::binary | std::ios::app);
    out.put('x');
  }
  const fs::path out_dir = temp.path() / "out";
  CHECK_THROWS_AS(index.save(out_dir.string(), "tiny"), covidx::OutputError);
  CHECK_FALSE(fs::exists(out_dir / "tiny.covidx"));
  if (fs::exists(out_dir)) {
    for (const auto& entry : fs::directory_iterator(out_dir)) {
      CHECK(entry.path().filename().string().find(".tmp.") ==
            std::string::npos);
    }
  }
}
