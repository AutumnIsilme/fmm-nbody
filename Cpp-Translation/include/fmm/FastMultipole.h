#pragma once
#include "../live_view/LiveSettings.h"
#include "../quadtree/Body.h"
#include "../quadtree/Box.h"
#include "../quadtree/Strata.h"
#include "../quadtree/Vec2.h"
#include <string>
#include <vector>
namespace fmm {

std::vector<Body>
bodies_from_rows(const std::vector<std::vector<double>> &rows);

Vec2 pairwise_force(const Body &left, const Body &right);

void solve_fmm_forces(Strata &strata, std::vector<Box *> &leaf_boxes,
                      double epsilon = 1e-7);

void run_fmm_simulation(int N = 500, int bodies_per_box = 25,
                        double epsilon = 1e-7, int num_steps = 500,
                        double dt = 1e-4, int rebuild_every = 1,
                        bool show_boxes = false,
                        const LiveViewSettings *live_view = nullptr);

} // namespace fmm
