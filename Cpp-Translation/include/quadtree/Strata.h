#pragma once

#include <vector>

#include "Box.h"

namespace fmm {

using Strata = std::vector<std::vector<Box *>>;

Strata stratify_quadtree(Box &root);

// Collects every leaf (childless) box across all strata.
std::vector<Box *> leaves(const Strata &strata);

} // namespace fmm
