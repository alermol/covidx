#pragma once
//
// Typed exception hierarchy.
//
// Every exception thrown by the covidx library derives from covidx::Error.
// The CLI maps these types onto the documented exit codes.
//
#include <stdexcept>
#include <string>

namespace covidx {

// Base class of all covidx errors.
class Error : public std::runtime_error {
public:
  explicit Error(const std::string& what) : std::runtime_error(what) {}
};

// Low-level I/O failure (open/read/write/mmap).
class IoError : public Error {
public:
  using Error::Error;
};

// Malformed or unreadable input (BAM, BED, malformed values).
class InputError : public IoError {
public:
  using IoError::IoError;
};

// Failure while writing or replacing an index on disk.
class OutputError : public IoError {
public:
  using IoError::IoError;
};

// Malformed, corrupt, or unsupported index data.
class FormatError : public Error {
public:
  using Error::Error;
};

// Invalid query arguments (bad range, unknown unit/contig, out of scope).
class QueryError : public Error {
public:
  using Error::Error;
};

// Invalid API/CLI usage (bad option, invalid argument value).
class UsageError : public Error {
public:
  using Error::Error;
};

} // namespace covidx
