#include "../../include/fmm/FastMultipole.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include "../../include/SimConstants.h"
#include "../../include/cloud/Cloud.h"
#include "../../include/live_view/LiveSettings.h"
#include "../../include/quadtree/InteractionLists.h"
#include "../../include/quadtree/QuadtreeBuilder.h"
#include "../../include/quadtree/Strata.h"
#include "../../include/quadtree/Visualisation.h"

// --- Helper function to log iteration state --- TODO remove
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

void log_step_diagnostics(int step, double current_time,
                          const std::vector<fmm::Body> &bodies,
                          const std::vector<fmm::Vec2> &velocities,
                          const std::vector<fmm::Vec2> &forces, double eps_sq,
                          std::ostream &out = std::cout) {
  return; // NOTE
  const std::size_t num_bodies = bodies.size();
  if (num_bodies == 0) {
    out << "==================== STEP " << step << " (t = " << std::fixed
        << std::setprecision(4) << current_time << ") ====================\n";
    out << "No active bodies remaining.\n";
    out << "==================================================================="
           "\n\n";
    return;
  }

  const double G = 1.0;

  double total_mass = 0.0;
  double cm_x = 0.0, cm_y = 0.0;
  double px = 0.0, py = 0.0;
  double angular_momentum_z = 0.0;
  double kinetic_energy = 0.0;
  double potential_energy = 0.0;

  double max_force_mag = 0.0;
  double max_speed = 0.0;

  // --- 1. Single-Body Accumulations ---
  for (std::size_t i = 0; i < num_bodies; ++i) {
    const auto &b = bodies[i];
    double m = b.mass();
    double x = b.x();
    double y = b.y();

    // Safely pull velocity and force with fallback if vectors are mismatched
    double vx = (i < velocities.size()) ? velocities[i].x : 0.0;
    double vy = (i < velocities.size()) ? velocities[i].y : 0.0;
    double fx = (i < forces.size()) ? forces[i].x : 0.0;
    double fy = (i < forces.size()) ? forces[i].y : 0.0;

    total_mass += m;
    cm_x += m * x;
    cm_y += m * y;

    px += m * vx;
    py += m * vy;

    angular_momentum_z += m * (x * vy - y * vx);

    double speed_sq = vx * vx + vy * vy;
    kinetic_energy += 0.5 * m * speed_sq;

    double force_mag = std::sqrt(fx * fx + fy * fy);
    max_force_mag = std::max(max_force_mag, force_mag);
    max_speed = std::max(max_speed, std::sqrt(speed_sq));
  }

  if (total_mass > 0.0) {
    cm_x /= total_mass;
    cm_y /= total_mass;
  }

  // --- 2. Pairwise Potential Energy ---
  for (std::size_t i = 0; i < num_bodies; ++i) {
    for (std::size_t j = i + 1; j < num_bodies; ++j) {
      double dx = bodies[j].x() - bodies[i].x();
      double dy = bodies[j].y() - bodies[i].y();
      double r2_soft = dx * dx + dy * dy + eps_sq;

      potential_energy +=
          0.5 * G * bodies[i].mass() * bodies[j].mass() * std::log(r2_soft);
    }
  }

  double total_energy = kinetic_energy + potential_energy;

  // --- 3. Structured Printout ---
  out << "==================== STEP " << step << " (t = " << std::fixed
      << std::setprecision(4) << current_time << ") ====================\n";
  out << std::scientific << std::setprecision(6);
  out << "Total Mass       : " << total_mass << "\n";
  out << "Center of Mass   : (" << cm_x << ", " << cm_y << ")\n";
  out << "Linear Momentum  : (" << px << ", " << py
      << ")  [|P| = " << std::sqrt(px * px + py * py) << "]\n";
  out << "Angular Mom (L)  : " << angular_momentum_z << "\n";
  out << "Kinetic Energy   : " << kinetic_energy << "\n";
  out << "Potential Energy : " << potential_energy << "\n";
  out << "Total Energy     : " << total_energy << "\n";
  out << "Max Speed        : " << max_speed << "\n";
  out << "Max Force Mag    : " << max_force_mag << "\n";

  // --- 4. Sample Tracking (Safely bounded) ---
  if (num_bodies > 0) {
    double v0_x = (!velocities.empty()) ? velocities[0].x : 0.0;
    double v0_y = (!velocities.empty()) ? velocities[0].y : 0.0;
    double f0_x = (!forces.empty()) ? forces[0].x : 0.0;
    double f0_y = (!forces.empty()) ? forces[0].y : 0.0;

    out << "--- Sample Body 0 State ---\n";
    out << "  pos: (" << bodies[0].x() << ", " << bodies[0].y() << ")\n";
    out << "  vel: (" << v0_x << ", " << v0_y << ")\n";
    out << "  f_ext: (" << f0_x << ", " << f0_y << ")\n";
  }

  if (num_bodies > 1) {
    double v1_x = (velocities.size() > 1) ? velocities[1].x : 0.0;
    double v1_y = (velocities.size() > 1) ? velocities[1].y : 0.0;
    double f1_x = (forces.size() > 1) ? forces[1].x : 0.0;
    double f1_y = (forces.size() > 1) ? forces[1].y : 0.0;

    out << "--- Sample Body 1 State ---\n";
    out << "  pos: (" << bodies[1].x() << ", " << bodies[1].y() << ")\n";
    out << "  vel: (" << v1_x << ", " << v1_y << ")\n";
    out << "  f_ext: (" << f1_x << ", " << f1_y << ")\n";
  }
  out << "==================================================================="
         "\n\n";
}

namespace fmm {

namespace {

std::vector<std::size_t> remove_escaped_bodies(std::vector<Box *> &leaf_boxes,
                                               const Vec2 &domain_centre,
                                               double escape_radius_sq) {
  std::vector<std::size_t> removed_ids;
  for (Box *leaf : leaf_boxes) {
    std::vector<Body> kept;
    kept.reserve(leaf->bodies_in_box.size());
    for (const Body &body : leaf->bodies_in_box) {
      double dx = body.x() - domain_centre.x;
      double dy = body.y() - domain_centre.y;
      bool finite = std::isfinite(body.x()) && std::isfinite(body.y());
      bool escaped = finite && (dx * dx + dy * dy > escape_radius_sq);
      if (!finite || escaped) {
        removed_ids.push_back(body.id);
      } else {
        kept.push_back(body);
      }
    }
    if (kept.size() != leaf->bodies_in_box.size()) {
      leaf->bodies_in_box = std::move(kept);
    }
  }
  return removed_ids;
}

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

double get_cubic_spline_force_factor(std::complex<double> dz, double h) {
  double r_sq = std::norm(dz); // |dz|^2
  double h_sq = h * h;

  // 1. Outside softening radius: Exact 2D Newtonian Gravity
  if (r_sq >= h_sq) {
    return 1.0 / r_sq;
  }

  // 2. Inside softening radius: Smooth C^1 Cubic Spline
  double r = std::sqrt(r_sq);
  if (r < 1e-12) {
    return 0.0; // Force vanishes at r = 0
  }

  return (3.0 / h_sq) - (2.0 * r / (h_sq * h));
}

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
  std::complex<double> dz(right.x() - left.x(), right.y() - left.y());

  // Derive the softening radius 'h' from your existing squared constant
  double h = std::sqrt(SimConstants::kSofteningSquared);
  double k = get_cubic_spline_force_factor(dz, h);

  double mass_product = left.mass() * right.mass();

  return Vec2(dz.real() * k * mass_product, dz.imag() * k * mass_product);
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

  for (Box *leaf : leaf_boxes) {
    std::vector<Body> finite_only;
    finite_only.reserve(leaf->bodies_in_box.size());
    for (const Body &body : leaf->bodies_in_box) {
      if (std::isfinite(body.x()) && std::isfinite(body.y())) {
        finite_only.push_back(body);
      } else {
        std::cerr << "solve_fmm_forces: body id=" << body.id
                  << " has a non-finite position and was excluded from "
                     "this solve (it should already have been removed by "
                     "remove_escaped_bodies)\n";
      }
    }
    if (finite_only.size() != leaf->bodies_in_box.size()) {
      leaf->bodies_in_box = std::move(finite_only);
    }
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
          for (int k = 1; k <= l; ++k) {
            box->multipole_expansion[static_cast<std::size_t>(l)] +=
                child->multipole_expansion[static_cast<std::size_t>(k)] *
                std::pow(z_0, l - k) * binomial(l - 1, k - 1);
          }
          double sign = (l % 2 == 0) ? 1.0 : -1.0;
          box->multipole_expansion[static_cast<std::size_t>(l)] +=
              sign * (child->multipole_expansion[0] / static_cast<double>(l)) *
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
        // Inside Step 5 (evaluating List 3 multipoles directly on leaf
        // particles)
        std::complex<double> relative =
            z - std::complex<double>(box->centre.x, box->centre.y);

        std::complex<double> force = box->multipole_expansion[0] / relative;
        std::complex<double> inv_rel_pow = 1.0 / relative;
        for (int k = 1; k <= p; ++k) {
          inv_rel_pow /= relative;
          force += (-static_cast<double>(k)) *
                   box->multipole_expansion[static_cast<std::size_t>(k)] *
                   inv_rel_pow;
        }

        leaf->forces[i] =
            leaf->forces[i] + Vec2(-force.real(), force.imag()) * body.mass();
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
              std::complex<double>(body.x(), body.y()) - centre;
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
      for (int l = 1; l <= p; ++l) {
        force += static_cast<double>(l) *
                 leaf->local_expansion[static_cast<std::size_t>(l)] *
                 std::pow(z, l - 1);
      }
      leaf->forces[i] =
          leaf->forces[i] + Vec2(-force.real(), force.imag()) * body.mass();
      ;
    }
  }
}

// Each step:
//   1. solve_fmm_forces on the current tree
//   2. kick-drift: v += (F/m)*dt, then x += v*dt
//   3. draw a frame
//   4. rebuild the tree once `rebuild_every` substeps have passed since
//      the last rebuild
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
  symmetrize_v_list(strata);
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
    symmetrize_v_list(strata);
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
    double time_remaining = dt;
    int substeps_taken = 0;

    while (time_remaining > 0.0) {
      solve_fmm_forces(strata, leaf_boxes, epsilon);

      double max_accel_sq = 0.0;
      bool saw_nonfinite_accel = false;
      for (Box *leaf : leaf_boxes) {
        for (std::size_t i = 0; i < leaf->bodies_in_box.size(); ++i) {
          const Body &body = leaf->bodies_in_box[i];
          Vec2 accel = leaf->forces[i] * (1.0 / masses[body.id]);
          double accel_sq = accel.x * accel.x + accel.y * accel.y;
          if (std::isfinite(accel_sq)) {
            max_accel_sq = std::max(max_accel_sq, accel_sq);
          } else {
            std::cerr << "step " << step << ": body id=" << body.id
                      << " has a non-finite acceleration this substep; "
                         "forcing the minimum substep\n";
            saw_nonfinite_accel = true;
          }
        }
      }
      if (saw_nonfinite_accel) {
        max_accel_sq = std::numeric_limits<double>::max();
      }

      ++substeps_taken;

      double dt_sub = time_remaining;
      if (max_accel_sq > 0.0) {
        double max_accel = std::sqrt(max_accel_sq);
        double softening_length = std::sqrt(SimConstants::kSofteningSquared);
        double dt_stable = SimConstants::kTimestepSafetyFactor *
                           std::sqrt(softening_length / max_accel);
        dt_stable = std::max(dt_stable, SimConstants::kMinSubstepDt);
        dt_sub = std::min(dt_sub, dt_stable);
      }

      if (substeps_taken == SimConstants::kMaxSubstepsPerFrame) {
        std::cerr << "step " << step << ": exceeded " << substeps_taken
                  << " substeps (max_accel_sq=" << max_accel_sq
                  << ") -- continuing with the floored step size\n";
      }

      for (Box *leaf : leaf_boxes) {
        for (std::size_t i = 0; i < leaf->bodies_in_box.size(); ++i) {
          Body &body = leaf->bodies_in_box[i];
          std::size_t id = body.id;
          Vec2 accel = leaf->forces[i] * (1.0 / masses[id]);
          double accel_sq = accel.x * accel.x + accel.y * accel.y;

          if (!std::isfinite(accel_sq)) {
            body.set_position(std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN());
            continue;
          }

          // TODO
          // if (accel_sq > SimConstants::kMaxAccelSquared) {
          //   double scale = std::sqrt(SimConstants::kMaxAccelSquared /
          //   accel_sq); accel = accel * scale;
          // }

          velocities[id] = velocities[id] + (accel * dt_sub);
          Vec2 displacement = velocities[id] * dt_sub;
          body.set_position(body.x() + displacement.x,
                            body.y() + displacement.y);
        }
      }

      Vec2 domain_centre = bodies_tree->root +
                           Vec2(bodies_tree->extent, bodies_tree->extent) / 2.0;
      double escape_radius =
          SimConstants::kEscapeRadiusMultiplier * bodies_tree->extent;
      auto removed = remove_escaped_bodies(leaf_boxes, domain_centre,
                                           escape_radius * escape_radius);
      for (std::size_t id : removed) {
        std::cerr << "body id=" << id << " escaped the domain at step " << step
                  << " - removing it from the simulation\n";
      }

      time_remaining -= dt_sub;

      ++steps_since_rebuild;
      if (steps_since_rebuild >= rebuild_every) {
        rebuild_tree();
        steps_since_rebuild = 0;
      }
    }

    // TODO remove
    std::size_t alive_count = 0;
    for (Box *leaf : leaf_boxes) {
      alive_count += leaf->bodies_in_box.size();
    }

    std::vector<Body> current_bodies;
    std::vector<Vec2> current_velocities;
    std::vector<Vec2> current_forces;
    current_bodies.reserve(alive_count);
    current_velocities.reserve(alive_count);
    current_forces.reserve(alive_count);

    for (Box *leaf : leaf_boxes) {
      for (std::size_t i = 0; i < leaf->bodies_in_box.size(); ++i) {
        const Body &body = leaf->bodies_in_box[i];
        current_bodies.push_back(body);

        // Guard global velocity vector lookup by ID
        if (body.id < velocities.size()) {
          current_velocities.push_back(velocities[body.id]);
        } else {
          current_velocities.push_back(Vec2(0.0, 0.0));
        }

        // Guard leaf force vector indexing
        if (i < leaf->forces.size()) {
          current_forces.push_back(leaf->forces[i]);
        } else {
          current_forces.push_back(Vec2(0.0, 0.0));
        }
      }
    }

    log_step_diagnostics(step, step * dt, current_bodies, current_velocities,
                         current_forces, SimConstants::kSofteningSquared);

    ++steps_since_live;
    if (live_writer != nullptr && steps_since_live >= live_view->frame_stride) {
      live_writer->write(*bodies_tree);
      steps_since_live = 0;
    }
  }

  std::cout << "Total simulation time (" << num_steps
            << " steps): " << seconds_since(start) << "\n";
}

} // namespace fmm
