#pragma once
//
//   U -- for a leaf box, its neighbouring leaf boxes (near interactions,
//        handled directly rather than via expansions).
//   V -- for a childless-parent's children, same-level boxes that are
//        well separated from it but not from its parent (the core
//        multipole-to-local translation list).
//   W / X -- the auxiliary lists used when adjacent boxes are leaves of
//        different sizes: W holds a leaf's well-separated descendants of
//        larger colleagues, and X (populated as a side effect of
//        building W) holds, for each such descendant, the leaf that is
//        well separated from it.

#include <vector>

#include "Box.h"
#include "Strata.h"

namespace fmm {

void symmetrize_v_list(const Strata &strata);

void populate_list_1(const std::vector<Box *> &leaf_boxes);

void populate_list_2(const Strata &strata);

void populate_list_3_and_4(const std::vector<Box *> &leaf_boxes);

} // namespace fmm
