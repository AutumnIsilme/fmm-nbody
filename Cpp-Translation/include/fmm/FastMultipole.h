#pragma once

#include <vector>

#include "../quadtree/Body.h"
#include "../quadtree/Box.h"
#include "../quadtree/Vec2.h"

namespace fmm {

std::vector<Body>
bodies_from_rows(const std::vector<std::vector<double>> &rows);

// This is the force ON `right` FROM `left`
Vec2 pairwise_force(const Body &left, const Body &right);

void basic(int N);

void run_fmm(int N = 500, int bodies_per_box = 5, double epsilon = 1e-7);

} // namespace fmm
