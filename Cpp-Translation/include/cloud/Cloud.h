#pragma once

#include <string>
#include <vector>

namespace fmm {

std::vector<std::string> lines(const std::vector<std::vector<double>> &bodies);

void write_points_and_masses(const std::vector<std::vector<double>> &bodies,
                             const std::string &filename);

// std::vector<std::vector<double>>
// generate_points_and_masses_uniform_random_distribution(int n_points,
//                                                        double max_mass);

std::vector<std::vector<double>>
generate_2d_bodies_uniform_random(int n_points, double max_mass);

} // namespace fmm
