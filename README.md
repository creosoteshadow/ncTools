# ncTools

A small collection of high-quality, non-cryptographic utilities in a single header.

- **wyrand** — ultra-fast 64-bit PRNG (~13 GB/s), stream splitting
- **RNG256** — 256-bit state PRNG, 2^256 period, 2^128 streams, ~4.2 GB/s
- **CompactHash** — ~10 GB/s 128-bit hasher, SMHasher clean
- **CompactHash_streaming** — incremental version (identical output)

## Features

- 100% C++, Header-only, no dependencies beyond C++17/20 intrinsics
- Passes PractRand (multi-TB) and SMHasher (full suite)
- Deterministic seeding resistant to poor inputs
- Minimal code size (~few hundred LOC core)

## Usage

```cpp
#include "ncTools.h"
using namespace ncTools;

wyrand rng(12345);
uint64_t r = rng();

RNG256 big_rng;
double d = big_rng.uni();

auto h = CompactHash(&data[0], data.size(), 42);
auto [lo, hi] = h;
```

## Quick Demo

```cpp
#define ENABLE_NCTOOLS_DEMO
#include "ncTools.h"
#include <iostream>

int main() {
	ncTools_demo();
}
```

See comments in ncTools.h for full API and test results.

## License:
MIT.
