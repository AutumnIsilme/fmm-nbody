#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "Body.h"
#include "Box.h"

namespace fmm {

// Splits `box` into four children and distributes `bodies` among them,
// recursing into any child that ends up with more than `bodies_per_box`
// bodies.
void split(Box &box, const std::vector<Body> &bodies,
           std::size_t bodies_per_box);

// Builds a quadtree over `bodies` covering the fixed domain [-1, 1] x
// [-1, 1], splitting any box that contains
// more than `bodies_per_box` bodies and then linking colleague lists.
// Returns the owning root box.
std::unique_ptr<Box> create_quadtree(const std::vector<Body> &bodies,
                                     std::size_t bodies_per_box);

} // namespace fmm
