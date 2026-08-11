#include "../../include/fmm/FastMultipole.h"

#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

#include "../../include/cloud/Cloud.h"
#include "../../include/quadtree/Colleagues.h"
#include "../../include/quadtree/InteractionLists.h"
#include "../../include/quadtree/QuadtreeBuilder.h"
#include "../../include/quadtree/Strata.h"

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

void accumulate_self_pairwise_forces_as_written(const std::vector<Body> &bodies,
                                                std::vector<Vec2> &forces) {
  if (bodies.empty())
    return;
  for (std::size_t i = 0; i + 1 < bodies.size(); ++i) {
    for (std::size_t j = 0; j < bodies.size() - (i + 1); ++j) {
      const Body &left = bodies[i];
      const Body &right = bodies[i + 1 + j];
      Vec2 force = pairwise_force(left, right);
      forces[i] = forces[i] + force;
      forces[i + 1 + j] = forces[i + 1 + j] - force;
    }
  }
}

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

Vec2 pairwise_force(const Body &left, const Body &right) {
  double dx = right.x() - left.x();
  double dy = right.y() - left.y();
  double norm = std::sqrt(dx * dx + dy * dy);
  double ux = dx / norm;
  double uy = dy / norm;
  double mass_product = left.mass() * right.mass();
  double scale = mass_product / (norm * norm);
  return Vec2(ux * scale, uy * scale);
}

void basic(int N) {
  auto start = Clock::now();

  std::vector<Body> bodies =
      bodies_from_rows(generate_2d_bodies_uniform_random(N, 10.0));

  Box box(Vec2(-1.0, -1.0), 2.0);
  box.bodies_in_box = bodies;
  box.forces.assign(box.bodies_in_box.size(), Vec2(0.0, 0.0));

  auto constructed = Clock::now();

  accumulate_self_pairwise_forces_as_written(box.bodies_in_box, box.forces);

  auto done = Clock::now();

  std::cout << "Total time: " << seconds_since(start) << "\n";
  std::cout << "Construction time: "
            << std::chrono::duration<double>(constructed - start).count()
            << "\n";
  std::cout << "Processing time: "
            << std::chrono::duration<double>(done - constructed).count()
            << "\n";
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

  int p = static_cast<int>(std::ceil(-std::log2(epsilon)));
  std::cout << "p: " << p << "\n";
  if (p + 1 > kExpansionTerms) {
    throw std::runtime_error(
        "run_fmm: epsilon requires p=" + std::to_string(p) +
        " expansion terms, exceeding the fixed Box array size "
        "(kExpansionTerms=" +
        std::to_string(kExpansionTerms) +
        "). Python's equivalent (a fixed-size-25 numpy array) "
        "would raise an IndexError here; use an epsilon >= 1e-7 or resize "
        "kExpansionTerms in Box.h.");
  }

  // --- Step 2.1: P2M (multipole expansion of each leaf's own bodies) ---
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

  auto step_2_1 = Clock::now();

  // --- Step 2.2: M2M (upward pass, strata[-2:1:-1] i.e. size-2 .. 2 inclusive)
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

  auto step_2_2 = Clock::now();

  // --- Step 3: near field, direct pairwise via U list, plus leaf
  // self-interaction --- See the KNOWN BUG note at the top of this file for
  // both the adjacent->forces reset-on-every-visit behaviour and the
  // self-interaction indexing bug, both reproduced faithfully.

  // Zero-init forces on every leaf
  for (Box *leaf : leaf_boxes) {
    leaf->forces.assign(leaf->bodies_in_box.size(), Vec2(0.0, 0.0));
  }

  for (Box *leaf : leaf_boxes) {
    std::vector<Box *> adjacents(leaf->U.begin(), leaf->U.end());
    for (Box *adjacent : adjacents) {
      adjacent->forces.assign(adjacent->bodies_in_box.size(),
                              Vec2(0.0, 0.0)); // faithful reset, see note above
      accumulate_cross_pairwise_forces(leaf->bodies_in_box, leaf->forces,
                                       adjacent->bodies_in_box,
                                       adjacent->forces);
      adjacent->U.erase(leaf);
    }
    accumulate_self_pairwise_forces_as_written(leaf->bodies_in_box,
                                               leaf->forces);
  }

  auto step_3 = Clock::now();

  // --- Step 4: M2L (V list -> local expansion) ---
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

  auto step_4 = Clock::now();

  // --- Step 5: W list, direct multipole evaluation added straight into leaf
  // forces ---
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

  auto step_5 = Clock::now();

  // --- Step 6: P2L (X list -> local expansion, direct body evaluation) ---
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

  auto step_6 = Clock::now();

  // --- Step 7: L2L (downward pass, local expansion shift) ---
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

  auto step_7 = Clock::now();

  // --- Step 8: evaluate local expansion at each leaf body, add to forces ---
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

  auto step_8 = Clock::now();
  auto end = Clock::now();

  std::cout << "Total time: " << seconds_since(start) << "\n";
  std::cout << "Step 1: "
            << std::chrono::duration<double>(populated - start).count() << "\n";
  std::cout << "Step 2.1: "
            << std::chrono::duration<double>(step_2_1 - populated).count()
            << "\n";
  std::cout << "Step 2.2: "
            << std::chrono::duration<double>(step_2_2 - step_2_1).count()
            << "\n";
  std::cout << "Step 3: "
            << std::chrono::duration<double>(step_3 - step_2_2).count() << "\n";
  std::cout << "Step 4: "
            << std::chrono::duration<double>(step_4 - step_3).count() << "\n";
  std::cout << "Step 5: "
            << std::chrono::duration<double>(step_5 - step_4).count() << "\n";
  std::cout << "Step 6: "
            << std::chrono::duration<double>(step_6 - step_5).count() << "\n";
  std::cout << "Step 7: "
            << std::chrono::duration<double>(step_7 - step_6).count() << "\n";
  std::cout << "Step 8: "
            << std::chrono::duration<double>(step_8 - step_7).count() << "\n";
  (void)end;
}

} // namespace fmm
