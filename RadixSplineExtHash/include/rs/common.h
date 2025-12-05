#pragma once

#include <cstddef>
#include <cstdint>

namespace rs {

// A CDF coordinate.
template <class KeyType>
struct Coord {
  KeyType x;
  double y;

  // Overload the equality operator for std::remove
  bool operator==(const Coord& other) const {
      return x == other.x;
  }
  
  bool operator<(const Coord& other) const {
    return x < other.x;
    }
};

struct SearchBound {
  size_t begin;
  size_t end;  // Exclusive.
};

}  // namespace rs
