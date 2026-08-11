#pragma once

#include <array>
#include <complex>
#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "Body.h"
#include "Vec2.h"

namespace fmm {

// Number of multipole/local expansion coefficients retained per box
constexpr int kExpansionTerms = 25;

class Box {
public:
  Vec2 root;                    // Bottom-left corner of this box.
  double extent;                // Side length of this box.
  std::vector<Box *> ancestors; // Root-most-first chain of parent boxes.

  bool has_child_boxes;
  // Child quadrant order:
  // 0 = (root),                   1 = (root + (0, split_size)),
  // 2 = (root + (split_size, 0)), 3 = (root + (split_size, split_size)).
  // A null entry means that quadrant held no bodies and was pruned
  std::array<std::unique_ptr<Box>, 4> child_boxes;

  std::vector<Body>
      bodies_in_box; // Only meaningful when has_child_boxes is false.

  // Same-level neighbour boxes.
  std::array<Box *, 8> colleagues;

  // FMM interaction lists
  std::set<Box *> U;
  std::set<Box *> V;
  std::set<Box *> W;
  std::set<Box *> X;

  Vec2 centre;
  std::array<std::complex<double>, kExpansionTerms> multipole_expansion;
  std::array<std::complex<double>, kExpansionTerms> local_expansion;
  std::vector<Vec2> forces;

  Box();
  Box(Vec2 root_, double extent_, std::vector<Box *> ancestors_ = {},
      Vec2 centre_ = Vec2());
  ~Box();

  // Boxes own their children, so they are move-only
  Box(const Box &) = delete;
  Box &operator=(const Box &) = delete;
  Box(Box &&) noexcept;
  Box &operator=(Box &&) noexcept;

  bool operator==(const Box &other) const;
  bool operator!=(const Box &other) const;

  int num_child_boxes() const;

  std::string repr() const;
};

std::ostream &operator<<(std::ostream &os, const Box &box);

} // namespace fmm
