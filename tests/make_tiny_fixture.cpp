// Regenerates the committed tiny fixture. Run from the repository root:
//   ./build/covidx_make_fixture
#include <filesystem>
#include <iostream>
#include <string_view>

#include "test_util.h"

int main(int argc, char** argv) {
  namespace fs = std::filesystem;
  fs::path data_dir =
      argc > 1 ? fs::path(argv[1]) : fs::path(COVIDX_TEST_DATA_DIR);
  fs::create_directories(data_dir);

  const fs::path bam = data_dir / "tiny.bam";
  tests::make_bam(bam, tests::kTinyHeader, tests::tiny_reads());

  tests::write_file(data_dir / "tiny_scope.bed", "chr1\t100\t350\n"
                                                 "chr1\t340\t650\n"
                                                 "chr1\t900\t1000\n"
                                                 "chr2\t0\t500\n");

  std::cout << "wrote " << bam << " and scope BED\n";
  return 0;
}
