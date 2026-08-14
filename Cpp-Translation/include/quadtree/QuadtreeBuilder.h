#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "Body.h"
#include "Box.h"

namespace fmm {

void split(Box &box, const std::vector<Body> &bodies,
           std::size_t bodies_per_box);

std::unique_ptr<Box> create_quadtree(const std::vector<Body> &bodies,
                                     std::size_t bodies_per_box);

} // namespace fmm
