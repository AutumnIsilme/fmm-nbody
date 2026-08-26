#pragma once

#include <vector>

#include "quadtree/Body.h"

namespace fmm {

std::vector<Body> generate_2d_bodies_uniform_random(int n_points,
                                                    double max_mass);

} // namespace fmm
