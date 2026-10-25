# Third-party notices

`covidx` is distributed under the MIT license (see `LICENSE`). It contains code adapted from the following MIT-licensed projects.

## gcrank

- Source: https://github.com/alermol/gcrank
- Pinned commit: `704eee87135ff221ec3c4d60068c39c48bd0bfaa`
- License: MIT
- Adapted files:
  - `include/covidx/detail/crc32.h`, `src/detail/crc32.cpp`
  - `include/covidx/detail/thread_pool.h`, `src/detail/thread_pool.cpp`
  - `include/covidx/detail/sha256.h`, `src/detail/sha256.cpp`
  - `include/covidx/detail/json.h`, `src/detail/json.cpp`
  - The atomic-save design pattern in `src/cov_index.cpp`
  - The mmap loading pattern in `src/unit_file.cpp`
  - The exception hierarchy shape in `include/covidx/error.h`

Copyright (c) Aleksey Ermolaev. Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## htslib

- Source: https://github.com/samtools/htslib
- License: MIT/BSD-style (see the htslib distribution)
- Use: BAM reading, indexing, and header parsing. htslib is linked as an external dependency and is not vendored in this repository.

## doctest

- Source: https://github.com/doctest/doctest
- Version: 2.5.0 single header (`tests/doctest.h`)
- License: MIT
- Copyright (c) 2016-2023 Viktor Kirilov

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
