#pragma once

#include <cstddef>
#include <vector>

#include "Body.h"
#include "Box.h"

namespace fmm {

void split(Box &box, const std::vector<Body> &bodies,
           std::size_t bodies_per_box, BoxAllocator &box_alloc);

Box *create_quadtree(std::vector<Body> &bodies,
                     std::size_t bodies_per_box, BoxAllocator &box_alloc);

} // namespace fmm
