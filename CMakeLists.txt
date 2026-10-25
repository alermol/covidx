cmake_minimum_required(VERSION 3.20)

project(covidx VERSION 0.2.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(COVIDX_TOP_LEVEL ON)
else()
  set(COVIDX_TOP_LEVEL OFF)
endif()
option(COVIDX_BUILD_TESTS "Build covidx tests" ${COVIDX_TOP_LEVEL})
option(COVIDX_BUILD_CLI "Build covidx command-line tools" ON)
option(COVIDX_BUILD_BENCH "Build benchmarks" OFF)
option(COVIDX_BUILD_FUZZERS "Build fuzzing targets" OFF)
option(COVIDX_ASAN "Enable AddressSanitizer" OFF)
option(COVIDX_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(COVIDX_TSAN "Enable ThreadSanitizer" OFF)
option(COVIDX_WERROR "Treat warnings as errors" OFF)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(HTSLIB QUIET IMPORTED_TARGET htslib>=1.16)
endif()
# Keep one target name in the exported interface regardless of how htslib
# was discovered, so consumers do not depend on the build-time method.
if(NOT TARGET HTSLIB::HTSLIB)
  if(TARGET PkgConfig::HTSLIB)
    add_library(HTSLIB::HTSLIB INTERFACE IMPORTED)
    set_target_properties(HTSLIB::HTSLIB PROPERTIES
      INTERFACE_LINK_LIBRARIES PkgConfig::HTSLIB)
  else()
    find_path(HTSLIB_INCLUDE_DIR htslib/sam.h REQUIRED)
    find_library(HTSLIB_LIBRARY hts REQUIRED)
    add_library(HTSLIB::HTSLIB UNKNOWN IMPORTED)
    set_target_properties(HTSLIB::HTSLIB PROPERTIES
      IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}")
  endif()
endif()

# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------

set(COVIDX_SOURCES
  src/format.cpp
  src/bucket_table.cpp
  src/build.cpp
  src/cigar.cpp
  src/cov_index.cpp
  src/manifest.cpp
  src/detail/bam_reader.cpp
  src/detail/bed.cpp
  src/detail/merge.cpp
  src/detail/scope.cpp
  src/detail/spool.cpp
  src/detail/text.cpp
  src/run_index.cpp
  src/run_sweep.cpp
  src/unit_file.cpp
  src/detail/crc32.cpp
  src/detail/json.cpp
  src/detail/sha256.cpp
  src/detail/thread_pool.cpp
)

add_library(covidx ${COVIDX_SOURCES})
add_library(covidx::covidx ALIAS covidx)

target_include_directories(covidx
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_definitions(covidx PUBLIC _FILE_OFFSET_BITS=64)
target_link_libraries(covidx PRIVATE ZLIB::ZLIB Threads::Threads)
target_link_libraries(covidx PRIVATE HTSLIB::HTSLIB)

set(COVIDX_WARNINGS -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
target_compile_options(covidx PRIVATE ${COVIDX_WARNINGS})
if(COVIDX_WERROR)
  target_compile_options(covidx PRIVATE -Werror)
endif()

# ---------------------------------------------------------------------------
# Sanitizers
# ---------------------------------------------------------------------------

if(COVIDX_ASAN AND COVIDX_TSAN)
  message(FATAL_ERROR "COVIDX_ASAN and COVIDX_TSAN are mutually exclusive")
endif()
if(COVIDX_ASAN)
  add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address)
endif()
if(COVIDX_UBSAN)
  add_compile_options(-fsanitize=undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=undefined)
endif()
if(COVIDX_TSAN)
  add_compile_options(-fsanitize=thread -fno-omit-frame-pointer)
  add_link_options(-fsanitize=thread)
endif()

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if(COVIDX_BUILD_CLI)
  add_executable(covidx-build cli/covidx_build.cpp)
  add_executable(covidx-query cli/covidx_query.cpp)
  foreach(target covidx-build covidx-query)
    target_link_libraries(${target} PRIVATE covidx)
    target_compile_options(${target} PRIVATE ${COVIDX_WARNINGS})
    if(COVIDX_WERROR)
      target_compile_options(${target} PRIVATE -Werror)
    endif()
  endforeach()
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

if(COVIDX_BUILD_TESTS AND NOT COVIDX_BUILD_CLI)
  message(FATAL_ERROR "COVIDX_BUILD_TESTS=ON requires COVIDX_BUILD_CLI=ON")
endif()

if(COVIDX_BUILD_TESTS)
  enable_testing()
  add_library(covidx_doctest INTERFACE)
  target_include_directories(covidx_doctest INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/tests)

  set(COVIDX_TEST_SOURCES
    tests/test_main.cpp
    tests/test_crc32.cpp
    tests/test_json.cpp
    tests/test_thread_pool.cpp
    tests/test_format.cpp
    tests/test_fixtures.cpp
    tests/test_cigar.cpp
    tests/test_run_sweep.cpp
    tests/test_bucket_table.cpp
    tests/test_bucket_edges.cpp
    tests/test_run_index.cpp
    tests/test_bam_reader.cpp
    tests/test_build_single.cpp
    tests/test_parity_samtools.cpp
    tests/test_parity_bedtools.cpp
    tests/test_build_scope.cpp
    tests/test_group_merge.cpp
    tests/test_build_multi.cpp
    tests/test_spool.cpp
    tests/test_determinism.cpp
    tests/test_manifest.cpp
    tests/test_unit_file.cpp
    tests/test_query.cpp
    tests/test_feature_counts.cpp
    tests/test_atomic_save.cpp
    tests/test_corruption.cpp
    tests/test_concurrency.cpp
    tests/test_mmap.cpp
    tests/test_cli_build.cpp
    tests/test_cli_query.cpp
    tests/test_golden.cpp
    tests/test_end_to_end.cpp
  )

  add_executable(covidx_tests ${COVIDX_TEST_SOURCES})
  target_link_libraries(covidx_tests PRIVATE
    covidx covidx_doctest ZLIB::ZLIB Threads::Threads)
  target_include_directories(covidx_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_compile_definitions(covidx_tests PRIVATE
    COVIDX_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data"
    COVIDX_BUILD_BIN="$<TARGET_FILE:covidx-build>"
    COVIDX_QUERY_BIN="$<TARGET_FILE:covidx-query>"
    COVIDX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(covidx_tests PRIVATE HTSLIB::HTSLIB)
  target_compile_options(covidx_tests PRIVATE ${COVIDX_WARNINGS})
  if(COVIDX_WERROR)
    target_compile_options(covidx_tests PRIVATE -Werror)
  endif()

  if(COVIDX_BUILD_CLI)
    add_dependencies(covidx_tests covidx-build covidx-query)
  endif()
  add_test(NAME covidx_unit COMMAND covidx_tests)
  set_tests_properties(covidx_unit PROPERTIES LABELS "unit")

  configure_file(tests/install_consumer_test.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/install_consumer_test.cmake @ONLY)
  add_test(NAME covidx_install_consumer
    COMMAND ${CMAKE_COMMAND} -P
      ${CMAKE_CURRENT_BINARY_DIR}/install_consumer_test.cmake)
  set_tests_properties(covidx_install_consumer PROPERTIES LABELS "install")

  add_executable(covidx_make_fixture tests/make_tiny_fixture.cpp)
  target_include_directories(covidx_make_fixture PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  target_compile_definitions(covidx_make_fixture PRIVATE
    COVIDX_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
  target_link_libraries(covidx_make_fixture PRIVATE HTSLIB::HTSLIB)
endif()

# ---------------------------------------------------------------------------
# Benchmarks and fuzzers (optional)
# ---------------------------------------------------------------------------

if(COVIDX_BUILD_BENCH)
  add_executable(covidx_bench_query bench/bench_query.cpp)
  add_executable(covidx_bench_build bench/bench_build.cpp)
  add_executable(covidx_stress_bam bench/stress_bam.cpp)
  target_link_libraries(covidx_bench_query PRIVATE covidx)
  target_link_libraries(covidx_bench_build PRIVATE covidx)
  target_link_libraries(covidx_stress_bam PRIVATE covidx)
  target_link_libraries(covidx_stress_bam PRIVATE HTSLIB::HTSLIB)
endif()

if(COVIDX_BUILD_FUZZERS)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "COVIDX_BUILD_FUZZERS requires Clang and libFuzzer")
  endif()
  foreach(fuzz_target fuzz_json fuzz_bed fuzz_manifest fuzz_cigar fuzz_query)
    add_executable(${fuzz_target} tests/fuzz/${fuzz_target}.cpp)
    target_link_libraries(${fuzz_target} PRIVATE covidx)
    target_compile_options(${fuzz_target} PRIVATE
      -fsanitize=fuzzer,address,undefined -fno-omit-frame-pointer)
    target_link_options(${fuzz_target} PRIVATE
      -fsanitize=fuzzer,address,undefined)
    target_link_libraries(${fuzz_target} PRIVATE HTSLIB::HTSLIB)
  endforeach()
endif()

# ---------------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------------

install(TARGETS covidx
  EXPORT covidxTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(FILES
  include/covidx/cov_index.h
  include/covidx/error.h
  include/covidx/format.h
  include/covidx/manifest.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/covidx
)
if(COVIDX_BUILD_CLI)
  install(TARGETS covidx-build covidx-query
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

install(EXPORT covidxTargets
  FILE covidxTargets.cmake
  NAMESPACE covidx::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/covidx
)
configure_package_config_file(
  cmake/covidxConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/covidxConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/covidx
)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/covidxConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/covidxConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/covidxConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/covidx
)
