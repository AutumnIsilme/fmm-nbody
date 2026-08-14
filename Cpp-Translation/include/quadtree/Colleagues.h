#pragma once

#include "Box.h"

namespace fmm {

// Walks the tree rooted at `box` top-down and, for each box that has
// children, works out each child's 8 same-level neighbours
// ("colleagues") from the parent's own colleague list.
void colleagify_quadtree(Box &box);

} // namespace fmm
