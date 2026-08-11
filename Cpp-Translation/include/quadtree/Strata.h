#pragma once

#include <vector>

#include "Box.h"

namespace fmm {

using Strata = std::vector<std::vector<Box *>>;

// Groups all boxes in the tree rooted at `root` into levels, with the
// root itself as stratum 0.
Strata stratify_quadtree(Box &root);

// Collects every leaf (childless) box across all strata.
std::vector<Box *> leaves(const Strata &strata);

} // namespace fmm
