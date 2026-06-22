#pragma once

// =============================================================================
// java_compat.h
//
// Exact C++ replicas of the Java primitives that drive the active-subnetwork
// search's deterministic ordering and random-number generation.  Pulled into
// a single header for re-usability.
//
// see: https://tinyurl.com/java-implementation
//
// Every routine here has been validated bit-for-bit against
// ActiveSubnetworkSearch.jar (OpenJDK 21):
//
//   * java_string_hashcode  — Java String.hashCode()
//   * java_spread           — HashMap internal hash spreading (h ^ h>>>16)
//   * java_cap_for          — HashMap/HashSet table capacity for a given size
//   * JavaRandom            — java.util.Random (48-bit LCG)
//   * java_shuffle          — Collections.shuffle(list, rnd)
//
// =============================================================================

#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// JAVA COLLECTION ORDERING
//
// Java's HashMap/HashSet iterates keys in bucket order:
//   bucket = spread(hashCode) & (capacity - 1)
// where capacity grows as the smallest power of two above size/0.75.
// Within a bucket, elements appear in insertion order (a stable sort on the
// bucket index therefore reproduces the full iteration order exactly).
// -----------------------------------------------------------------------------

// Java String.hashCode(): h = 31*h + c  (signed 32-bit int, wraps on overflow)
static inline int java_string_hashcode(const std::string& s) {
  int h = 0;
  for (unsigned char c : s) h = 31 * h + static_cast<int>(c);
  return h;
}

// HashMap internal spreading: h ^ (h >>> 16)  (unsigned right-shift)
static inline int java_spread(int h) {
  unsigned uh = static_cast<unsigned>(h);
  return static_cast<int>(uh ^ (uh >> 16));
}

// HashMap/HashSet table capacity for `size` elements
//   cap = 16;  while (size > cap * 0.75) cap *= 2;
// NOTE: parity holds only for the linked-list bucket model.
// Java treeifies a bucket when it reaches >= 8 entries *and* the table
// capacity is >= 64, at which point iteration order becomes tree order and
// would diverge from this reconstruction.  At load factor 0.75 this
// requires a bucket collision run of 8+, which is astronomically unlikely
// for String.hashCode values on realistic gene/protein names, so in
// practice this is never triggered.
static inline int java_cap_for(int size) {
  int cap = 16;
  while (size > static_cast<int>(cap * 0.75)) cap <<= 1;
  return cap;
}

// -----------------------------------------------------------------------------
// JAVA RANDOM NUMBER GENERATION
//
// Exact replica of java.util.Random (48-bit linear congruential generator).
// The constants and algorithm are specified in the Java API docs and are
// stable across JDK versions.
// -----------------------------------------------------------------------------

struct JavaRandom {
  uint64_t seed;

  static const uint64_t MULT = 0x5DEECE66DULL;
  static const uint64_t ADD  = 0xBULL;
  static const uint64_t MASK = (1ULL << 48) - 1;

  explicit JavaRandom(int64_t s) {
    seed = (static_cast<uint64_t>(s) ^ MULT) & MASK;
  }

  // Mirrors java.util.Random.next(bits)
  int next(int bits) {
    seed = (seed * MULT + ADD) & MASK;
    return static_cast<int>(static_cast<int64_t>(seed) >> (48 - bits));
  }

  // Mirrors java.util.Random.nextInt(bound)
  int nextInt(int bound) {
    if ((bound & -bound) == bound)   // power of two: fast path
      return static_cast<int>(
        (static_cast<int64_t>(bound) * static_cast<int64_t>(next(31))) >> 31);
    int bits, val;
    do {
      bits = next(31);
      val  = bits % bound;
    } while (bits - val + (bound - 1) < 0);
    return val;
  }

  // Mirrors java.util.Random.nextDouble()
  double nextDouble() {
    return (  (static_cast<int64_t>(next(26)) << 27)
                + static_cast<int64_t>(next(27)) )
    * (1.0 / (1LL << 53));
  }
};

// Mirrors java.util.Collections.shuffle(list, rnd):
//   for (int i = size; i > 1; i--) swap(list[i-1], list[rnd.nextInt(i)]);
// Works for any element type T.
template <typename T>
static inline void java_shuffle(std::vector<T>& v, JavaRandom& rng) {
  for (int i = static_cast<int>(v.size()); i > 1; --i) {
    int j = rng.nextInt(i);
    std::swap(v[i - 1], v[j]);
  }
}
