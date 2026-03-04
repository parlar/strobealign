# Performance Optimizations

Benchmarked on real E. coli K-12 paired-end reads (DRR217225, ~6.8M read pairs, 150bp, 4 threads).

## Summary

| Stage | Baseline | Optimized | Improvement |
|---|---|---|---|
| Creating strobemers | 8.27s | 5.31s | **-35.8%** |
| Finding hits | 7.79s | 7.55s | -3.1% |
| Chaining | 2.70s | 2.64s | -2.2% |
| Extending & pairing | 16.22s | 15.70s | -3.2% |
| **Total mapping** | **37.23s** | **33.45s** | **-10.2%** |
| **Wall clock** | **37.56s** | **33.75s** | **-10.1%** |

## Changes

### 1. AVX2 Smith-Waterman alignment (extension -5.9%)

Added 256-bit AVX2 implementation of the striped Smith-Waterman kernel alongside the existing SSE2 code. Enabled at compile time with `-DENABLE_AVX=ON`. Doubles the SIMD width from 16 to 32 byte-lanes per vector.

**Files:** `cpp/ext/ssw/ssw_avx2.c`, `cpp/ext/ssw/ssw_avx2.h`, `cpp/ext/ssw/ssw.c`, `CMakeLists.txt`

### 2. Software prefetching in hit finding (hit finding -2.4%)

Prefetch the next query strobe's hash bucket while processing the current one, hiding memory access latency. Expect larger gains on human genome where the index doesn't fit in cache.

**Files:** `cpp/index.hpp`, `cpp/hits.cpp`

### 3. Hardware POPCNT instruction (strobemer creation -8.7%)

Replaced `std::bitset<64>::count()` with `__builtin_popcountll()` and added `-mpopcnt` compile flag. The compiler was emitting software `__popcountdi2` calls despite the CPU having hardware POPCNT support.

**Files:** `cpp/randstrobes.cpp`, `CMakeLists.txt`

### 4. Eliminate allocations in has_shared_substring()

Replaced `std::string::substr()` (heap allocation per iteration) with `std::string_view::substr()` (zero-copy). Bigger impact expected on human genome where rescue alignment is triggered more frequently.

**Files:** `cpp/aln.cpp`

### 5. Circular buffer for SyncmerIterator (strobemer creation -25.6%)

Replaced `std::deque<uint64_t>` with an inline fixed-capacity circular buffer (`SmerBuffer`). The deque was heap-allocating for a 5-element, 40-byte sliding window — massive overhead relative to data size. The inline buffer keeps everything on the stack.

**Files:** `cpp/randstrobes.hpp`

### 6. Pre-allocate QueryRandstrobe vectors

Added `reserve(syncmers.size())` to avoid ~7 reallocations per read per strand during randstrobe generation.

**Files:** `cpp/randstrobes.cpp`

## Build

```bash
# With AVX2 (recommended on x86-64 with AVX2 support):
cmake -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_AVX=ON

# Without AVX2 (default, SSE2 only):
cmake -B build -DCMAKE_BUILD_TYPE=Release
```

## Correctness

- All 50 unit tests pass in both SSE2-only and AVX2 configurations
- Mapping results are bit-identical between SSE2 and AVX2 builds
- All changes are safe for human genome mapping
