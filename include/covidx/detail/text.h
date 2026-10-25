#pragma once
//
// Small UTF-8 and name-validation helpers shared by the manifest and BAM
// header parsing.
//
#include <string>
#include <string_view>

namespace covidx::detail {

[[nodiscard]] bool is_valid_utf8(std::string_view text);

// Validates a sample, group, or unit name: non-empty, <= 255 bytes, valid
// UTF-8, and free of NUL, tab, CR, and LF. Throws UsageError.
void validate_unit_name(std::string_view name, std::string_view kind,
                        bool allow_pooled);

// Validates an index name: [A-Za-z0-9._-]+, <= 128 bytes, not "." or "..".
[[nodiscard]] bool valid_index_name(std::string_view name);

} // namespace covidx::detail
