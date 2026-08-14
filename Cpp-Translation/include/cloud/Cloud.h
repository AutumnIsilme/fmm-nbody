#pragma once

#include <string>
#include <vector>

namespace fmm {

std::vector<std::vector<double>>
generate_2d_bodies_uniform_random(int n_points, double max_mass);

} // namespace fmm
