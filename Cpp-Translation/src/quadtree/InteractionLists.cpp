#include "../../include/quadtree/InteractionLists.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <set>
namespace fmm {
void populate_list_1(const std::vector<Box *> &leaf_boxes) {
  for (Box *box : leaf_boxes) {
    for (Box *colleague : box->colleagues) {
      if (colleague != nullptr && !colleague->has_child_boxes) {
        box->U.insert(colleague);
      }
    }
    for (Box *cobox : box->U) {
      cobox->U.insert(box);
    }
  }
}
void populate_list_2(const Strata &strata) {
  for (std::size_t i = 0; i + 1 < strata.size(); ++i) {
    for (Box *box : strata[i]) {
      if (!box->has_child_boxes)
        continue;
      std::vector<Box *> colleague_children;
      for (Box *colleague : box->colleagues) {
        if (colleague != nullptr && colleague->has_child_boxes) {
          for (auto &cc : colleague->child_boxes) {
            colleague_children.push_back(cc.get());
          }
        } else {
          colleague_children.push_back(nullptr);
        }
      }
      for (auto &child_ptr : box->child_boxes) {
        Box *child = child_ptr.get();
        if (child == nullptr)
          continue;
        for (Box *candidate : colleague_children) {
          if (candidate == nullptr)
            continue;
          if (child->U.find(candidate) != child->U.end())
            continue;
          if (std::find(child->colleagues.begin(), child->colleagues.end(),
                        candidate) != child->colleagues.end()) {
            continue;
          }
          child->V.insert(candidate);
        }
      }
    }
  }
}

namespace {
// V-membership is supposed to be a symmetric relation: if box A is
// well-separated from box B, B is equally well-separated from A. Unlike U
// (which populate_list_1 explicitly patches for symmetry a few lines up),
// nothing forces this for V -- populate_list_2 only ever inserts into
// child->V from the perspective of a single ancestor's colleague set. If
// any tree-level irregularity ever leaves that one-sided (A gets B in its
// V, but B never gets A), B silently never applies the reciprocal M2L
// force back onto A. That's a real force on one side of the pair with no
// equal-and-opposite reaction anywhere -- a small, systematic, per-step
// momentum leak that predates and does not require any close encounter,
// which is exactly the drift-before-fireworks pattern being chased here.
//
// This is a diagnostic/defensive pass, not a substitute for verifying the
// tree itself never produces an orphaned (un-mirrored) V relationship in
// the first place -- if this measurably fixes the drift, that confirms
// the mechanism and the next step is enforcing 2:1 level-balancing in the
// quadtree so the underlying asymmetry can't occur.
} // namespace

void symmetrize_v_list(const Strata &strata) {
  for (const auto &level : strata) {
    for (Box *box : level) {
      for (Box *other : box->V) {
        other->V.insert(box);
      }
    }
  }
}

namespace {
std::vector<Box *> collect_w_descendants(Box *origin_leaf, Box *box,
                                         const std::set<int> &qualifying) {
  std::vector<Box *> result;
  for (int i = 0; i < 4; ++i) {
    Box *child = box->child_boxes[static_cast<std::size_t>(i)].get();
    if (child == nullptr)
      continue;
    if (qualifying.count(i) || !child->has_child_boxes) {
      result.push_back(child);
      child->X.insert(origin_leaf);
    } else {
      std::vector<Box *> nested =
          collect_w_descendants(origin_leaf, child, qualifying);
      result.insert(result.end(), nested.begin(), nested.end());
    }
  }
  return result;
}
} // namespace
void populate_list_3_and_4(const std::vector<Box *> &leaf_boxes) {
  static const std::array<std::set<int>, 8> kQualifyingSets = {
      std::set<int>{1, 3},    std::set<int>{1, 2, 3}, std::set<int>{2, 3},
      std::set<int>{0, 2, 3}, std::set<int>{0, 2},    std::set<int>{0, 1, 2},
      std::set<int>{0, 1},    std::set<int>{0, 1, 3}};
  for (Box *box : leaf_boxes) {
    for (int c = 0; c < 8; ++c) {
      Box *colleague = box->colleagues[static_cast<std::size_t>(c)];
      if (colleague != nullptr && colleague->has_child_boxes) {
        std::vector<Box *> found = collect_w_descendants(
            box, colleague, kQualifyingSets[static_cast<std::size_t>(c)]);
        box->W.insert(found.begin(), found.end());
      }
    }
  }
}
} // namespace fmm
