#include "../../include/fmm/FastMultipole.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "../../include/cloud/Cloud.h"
#include "../../include/live_view/LiveSettings.h"
#include "../../include/quadtree/Colleagues.h"
#include "../../include/quadtree/InteractionLists.h"
#include "../../include/quadtree/QuadtreeBuilder.h"
#include "../../include/quadtree/Strata.h"
#include "../../include/quadtree/Visualisation.h"

namespace fmm {

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(const Clock::time_point &start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

double binomial(int n, int k) {
  if (k < 0 || n < 0 || k > n)
    return 0.0;
  if (k > n - k)
    k = n - k;
  double result = 1.0;
  for (int i = 0; i < k; ++i) {
    result *= static_cast<double>(n - i);
    result /= static_cast<double>(i + 1);
  }
  return result;
}

// I saw a post saying this was sensible for more than one time step to reduce
// things blowing up: FIX this might need tuning
constexpr double kSofteningSquared = 1e-4;

void accumulate_cross_pairwise_forces(const std::vector<Body> &left_bodies,
                                      std::vector<Vec2> &left_forces,
                                      const std::vector<Body> &right_bodies,
                                      std::vector<Vec2> &right_forces) {
  for (std::size_t i = 0; i < left_bodies.size(); ++i) {
    for (std::size_t j = 0; j < right_bodies.size(); ++j) {
      Vec2 force = pairwise_force(left_bodies[i], right_bodies[j]);
      left_forces[i] = left_forces[i] + force;
      right_forces[j] = right_forces[j] - force;
    }
  }
}

void accumulate_self_pairwise_forces(const std::vector<Body> &bodies,
                                     std::vector<Vec2> &forces) {
  if (bodies.empty())
    return;
  for (std::size_t i = 0; i + 1 < bodies.size(); ++i) {
    for (std::size_t j = i + 1; j < bodies.size(); ++j) {
      const Body &left = bodies[i];
      const Body &right = bodies[j];
      Vec2 force = pairwise_force(left, right);
      forces[i] = forces[i] + force;
      forces[j] = forces[j] - force;
    }
  }
}

// Allow for multiple steps by clearing tree
void reset_expansions(Strata &strata) {
  for (auto &level : strata) {
    for (Box *box : level) {
      std::fill(std::begin(box->multipole_expansion),
                std::end(box->multipole_expansion),
                std::complex<double>(0.0, 0.0));
      std::fill(std::begin(box->local_expansion),
                std::end(box->local_expansion), std::complex<double>(0.0, 0.0));
    }
  }
}

class LiveWriter {
public:
  LiveWriter(std::string path, std::string box_colour,
             std::string particle_colour)
      : path_(std::move(path)), tmp_path_(path_ + ".tmp"),
        box_colour_(std::move(box_colour)),
        particle_colour_(std::move(particle_colour)) {
    auto parent = std::filesystem::path(path_).parent_path();
    if (!parent.empty()) {
      std::filesystem::create_directories(parent);
    }
  }

  void write(const Box &root) {
    graph_quadtree(root, tmp_path_, box_colour_, particle_colour_);
    std::filesystem::rename(tmp_path_, path_);
  }

private:
  std::string path_;
  std::string tmp_path_;
  std::string box_colour_;
  std::string particle_colour_;
};
} // namespace

std::vector<Body>
bodies_from_rows(const std::vector<std::vector<double>> &rows) {
  std::vector<Body> result;
  result.reserve(rows.size());
  for (const auto &row : rows) {
    result.emplace_back(row);
  }
  return result;
}

// I added a softening factor in here
Vec2 pairwise_force(const Body &left, const Body &right) {
  double dx = right.x() - left.x();
  double dy = right.y() - left.y();
  double norm_sq = dx * dx + dy * dy + kSofteningSquared;
  double denom = norm_sq * std::sqrt(norm_sq); // (r^2 + eps^2)^{3/2}
  double mass_product = left.mass() * right.mass();
  double scale_factor = mass_product / denom;
  return Vec2(dx * scale_factor, dy * scale_factor);
}

void basic(int N) {
  auto start = Clock::now();

  std::vector<Body> bodies =
      bodies_from_rows(fmm::generate_2d_bodies_uniform_random(N, 10.0));

  Box box(Vec2(-1.0, -1.0), 2.0);
  box.bodies_in_box = bodies;
  box.forces.assign(box.bodies_in_box.size(), Vec2(0.0, 0.0));

  auto constructed = Clock::now();

  accumulate_self_pairwise_forces(box.bodies_in_box, box.forces);

  auto done = Clock::now();

  std::cout << "Total time: " << seconds_since(start) << "\n";
  std::cout << "Construction time: "
            << std::chrono::duration<double>(constructed - start).count()
            << "\n";
  std::cout << "Processing time: "
            << std::chrono::duration<double>(done - constructed).count()
            << "\n";
}

void solve_fmm_forces(Strata &strata, std::vector<Box *> &leaf_boxes,
                      double epsilon) {
  int p = static_cast<int>(std::ceil(-std::log2(epsilon)));
  if (p + 1 > kExpansionTerms) {
    throw std::runtime_error(
        "solve_fmm_forces: epsilon requires p=" + std::to_string(p) +
        " expansion terms, exceeding the fixed Box array size "
        "(kExpansionTerms=" +
        std::to_string(kExpansionTerms) + ").");
  }

  reset_expansions(strata);

  // --- Step 2.1 ---
  for (Box *leaf : leaf_boxes) {
    for (const Body &body : leaf->bodies_in_box) {
      std::complex<double> z_i(body.x() - leaf->centre.x,
                               body.y() - leaf->centre.y);
      double mass = body.mass();
      leaf->multipole_expansion[0] += mass;
      std::complex<double> z_pow(1.0, 0.0);
      for (int k = 1; k <= p; ++k) {
        z_pow *= z_i;
        leaf->multipole_expansion[static_cast<std::size_t>(k)] +=
            (-mass / static_cast<double>(k)) * z_pow;
      }
    }
  }

  // --- Step 2.2 ---
  for (int i = static_cast<int>(strata.size()) - 2; i >= 2; --i) {
    for (Box *box : strata[static_cast<std::size_t>(i)]) {
      if (!box->has_child_boxes)
        continue;
      for (auto &child_ptr : box->child_boxes) {
        Box *child = child_ptr.get();
        if (child == nullptr)
          continue;
        std::complex<double> z_0(child->centre.x - box->centre.x,
                                 child->centre.y - box->centre.y);
        box->multipole_expansion[0] += child->multipole_expansion[0];
        for (int l = 1; l <= p; ++l) {
          for (int k = 1; k <= l - 1; ++k) {
            box->multipole_expansion[static_cast<std::size_t>(l)] +=
                child->multipole_expansion[static_cast<std::size_t>(k)] *
                std::pow(z_0, l - k) * binomial(l - 1, k - 1);
          }
          box->multipole_expansion[static_cast<std::size_t>(l)] -=
              (child->multipole_expansion[0] / static_cast<double>(l)) *
              std::pow(z_0, l);
        }
      }
    }
  }

  // --- Step 3 ---

  // Zero-init forces on every leaf
  for (Box *leaf : leaf_boxes) {
    leaf->forces.assign(leaf->bodies_in_box.size(), Vec2(0.0, 0.0));
  }

  std::less<Box *> box_order;
  for (Box *leaf : leaf_boxes) {
    for (Box *adjacent : leaf->U) {
      if (!box_order(leaf, adjacent))
        continue; // wait for the other side of the pair to process it
      accumulate_cross_pairwise_forces(leaf->bodies_in_box, leaf->forces,
                                       adjacent->bodies_in_box,
                                       adjacent->forces);
    }
    accumulate_self_pairwise_forces(leaf->bodies_in_box, leaf->forces);
  }

  // --- Step 4 ---
  for (std::size_t i = 2; i < strata.size(); ++i) {
    for (Box *box : strata[i]) {
      for (Box *b_j : box->V) {
        std::complex<double> z_0(b_j->centre.x - box->centre.x,
                                 b_j->centre.y - box->centre.y);
        std::complex<double> recip_z_0 = 1.0 / z_0;
        std::complex<double> neg_recip_z_0 = -1.0 / z_0;

        std::complex<double> sum0(0.0, 0.0);
        std::complex<double> neg_pow(1.0, 0.0);
        for (int k = 1; k <= p; ++k) {
          neg_pow *= neg_recip_z_0;
          sum0 +=
              b_j->multipole_expansion[static_cast<std::size_t>(k)] * neg_pow;
        }
        box->local_expansion[0] +=
            b_j->multipole_expansion[0] * std::log(-z_0) + sum0;

        for (int l = 1; l <= p; ++l) {
          std::complex<double> inner =
              -b_j->multipole_expansion[0] / static_cast<double>(l);
          std::complex<double> neg_pow_k(1.0, 0.0);
          for (int k = 1; k <= p; ++k) {
            neg_pow_k *= neg_recip_z_0;
            inner += b_j->multipole_expansion[static_cast<std::size_t>(k)] *
                     neg_pow_k * binomial(l + k - 1, k - 1);
          }
          box->local_expansion[static_cast<std::size_t>(l)] +=
              std::pow(recip_z_0, l) * inner;
        }
      }
    }
  }

  // --- Step 5 ---
  for (Box *leaf : leaf_boxes) {
    for (std::size_t i = 0; i < leaf->bodies_in_box.size(); ++i) {
      const Body &body = leaf->bodies_in_box[i];
      std::complex<double> z(body.x(), body.y());
      for (Box *box : leaf->W) {
        std::complex<double> relative =
            z - std::complex<double>(box->centre.x, box->centre.y);
        std::complex<double> force = box->multipole_expansion[0] / relative;
        std::complex<double> inv_rel_pow = 1.0 / relative;
        for (int k = 1; k <= p; ++k) {
          inv_rel_pow /= relative; // relative^-(k+1)
          force += (-static_cast<double>(k)) *
                   box->multipole_expansion[static_cast<std::size_t>(k)] *
                   inv_rel_pow;
        }
        leaf->forces[i] = leaf->forces[i] + Vec2(force.real(), force.imag());
      }
    }
  }

  // --- Step 6 ---
  for (std::size_t i = 2; i < strata.size(); ++i) {
    for (Box *box : strata[i]) {
      if (box->X.empty())
        continue;
      std::complex<double> centre(box->centre.x, box->centre.y);
      for (Box *big_leaf : box->X) {
        for (const Body &body : big_leaf->bodies_in_box) {
          std::complex<double> z =
              centre - std::complex<double>(body.x(), body.y());
          double mass = body.mass();
          box->local_expansion[0] += mass * std::log(-z);
          for (int l = 1; l <= p; ++l) {
            box->local_expansion[static_cast<std::size_t>(l)] -=
                mass / (static_cast<double>(l) * std::pow(z, l));
          }
        }
      }
    }
  }

  // --- Step 7 ---
  for (std::size_t i = 2; i < strata.size(); ++i) {
    for (Box *box : strata[i]) {
      if (!box->has_child_boxes)
        continue;
      std::complex<double> centre(box->centre.x, box->centre.y);
      for (auto &child_ptr : box->child_boxes) {
        Box *child = child_ptr.get();
        if (child == nullptr)
          continue;
        std::complex<double> z =
            centre - std::complex<double>(child->centre.x, child->centre.y);
        for (int l = 0; l <= p; ++l) {
          for (int k = l; k <= p; ++k) {
            child->local_expansion[static_cast<std::size_t>(l)] +=
                box->local_expansion[static_cast<std::size_t>(k)] *
                std::pow(-z, k - l) * binomial(k, l);
          }
        }
      }
    }
  }

  // --- Step 8 ---
  for (Box *leaf : leaf_boxes) {
    std::complex<double> centre(leaf->centre.x, leaf->centre.y);
    for (std::size_t i = 0; i < leaf->bodies_in_box.size(); ++i) {
      const Body &body = leaf->bodies_in_box[i];
      std::complex<double> z =
          std::complex<double>(body.x(), body.y()) - centre;
      std::complex<double> force(0.0, 0.0);
      for (int l = 0; l <= p; ++l) {
        force += static_cast<double>(l) *
                 leaf->local_expansion[static_cast<std::size_t>(l)] *
                 std::pow(z, l - 1);
      }
      leaf->forces[i] = leaf->forces[i] + Vec2(force.real(), force.imag());
    }
  }
}

void run_fmm(int N, int bodies_per_box, double epsilon) {
  auto start = Clock::now();

  std::vector<Body> bodies =
      bodies_from_rows(generate_2d_bodies_uniform_random(N, 10.0));
  std::unique_ptr<Box> bodies_tree =
      create_quadtree(bodies, static_cast<std::size_t>(bodies_per_box));
  Strata strata = stratify_quadtree(*bodies_tree);
  std::vector<Box *> leaf_boxes = leaves(strata);
  populate_list_1(leaf_boxes);
  populate_list_2(strata);
  populate_list_3_and_4(leaf_boxes);

  auto populated = Clock::now();

  solve_fmm_forces(strata, leaf_boxes, epsilon);

  auto solved = Clock::now();

  std::cout << "Total time: " << seconds_since(start) << "\n";
  std::cout << "Tree build + interaction lists: "
            << std::chrono::duration<double>(populated - start).count() << "\n";
  std::cout << "Force solve (steps 2.1-8): "
            << std::chrono::duration<double>(solved - populated).count()
            << "\n";
}

// Each step:
//   1. solve_fmm_forces on the current tree
//   2. kick-drift: v += (F/m)*dt, then x += v*dt
//   3. draw a frame
//   4. rebuild the tree once `rebuild_every` steps have passed since the
//      last rebuild
void run_fmm_simulation(int N, int bodies_per_box, double epsilon,
                        int num_steps, double dt, int rebuild_every,
                        const LiveViewSettings *live_view) {
  if (rebuild_every <= 0) {
    throw std::invalid_argument(
        "run_fmm_simulation: rebuild_every must be positive");
  }

  auto start = Clock::now();

  std::vector<Body> bodies =
      bodies_from_rows(generate_2d_bodies_uniform_random(N, 10.0));

  for (std::size_t i = 0; i < bodies.size(); ++i) {
    bodies[i].id = i;
  }

  std::unique_ptr<Box> bodies_tree =
      create_quadtree(bodies, static_cast<std::size_t>(bodies_per_box));
  Strata strata = stratify_quadtree(*bodies_tree);
  std::vector<Box *> leaf_boxes = leaves(strata);
  populate_list_1(leaf_boxes);
  populate_list_2(strata);
  populate_list_3_and_4(leaf_boxes);

  std::vector<Vec2> velocities(static_cast<std::size_t>(N), Vec2(0.0, 0.0));

  std::vector<double> masses(static_cast<std::size_t>(N), 0.0);
  for (Box *leaf : leaf_boxes) {
    for (const Body &body : leaf->bodies_in_box) {
      masses[body.id] = body.mass();
    }
  }

  // rebuild tree based on new particle locations
  auto rebuild_tree = [&]() {
    std::vector<Body> flat_bodies;
    flat_bodies.reserve(static_cast<std::size_t>(N));
    for (Box *leaf : leaf_boxes) {
      for (const Body &body : leaf->bodies_in_box) {
        flat_bodies.push_back(body);
      }
    }
    bodies_tree =
        create_quadtree(flat_bodies, static_cast<std::size_t>(bodies_per_box));
    strata = stratify_quadtree(*bodies_tree);
    leaf_boxes = leaves(strata);
    populate_list_1(leaf_boxes);
    populate_list_2(strata);
    populate_list_3_and_4(leaf_boxes);
  };

  // animations
  std::unique_ptr<LiveWriter> live_writer;
  if (live_view != nullptr) {
    live_writer = std::make_unique<LiveWriter>(live_view->output_path,
                                               live_view->box_colour,
                                               live_view->particle_colour);
    live_writer->write(*bodies_tree); // initial configuration, t = 0
  }

  int steps_since_rebuild = 0;
  int steps_since_live = 0;

  for (int step = 0; step < num_steps; ++step) {
    solve_fmm_forces(strata, leaf_boxes, epsilon);

    for (Box *leaf : leaf_boxes) {
      for (std::size_t i = 0; i < leaf->bodies_in_box.size(); ++i) {
        Body &body = leaf->bodies_in_box[i];
        std::size_t id = body.id;

        Vec2 accel = leaf->forces[i] * (1.0 / masses[id]);
        velocities[id] = velocities[id] + (accel * dt); // kick
        Vec2 displacement = velocities[id] * dt;        // drift
        body.set_position(body.x() + displacement.x, body.y() + displacement.y);
      }
    }

    ++steps_since_live;
    if (live_writer != nullptr && steps_since_live >= live_view->frame_stride) {
      live_writer->write(*bodies_tree);
      steps_since_live = 0;
    }

    ++steps_since_rebuild;
    if (steps_since_rebuild >= rebuild_every) {
      rebuild_tree();
      steps_since_rebuild = 0;
    }
  }

  std::cout << "Total simulation time (" << num_steps
            << " steps): " << seconds_since(start) << "\n";
}

} // namespace fmm
