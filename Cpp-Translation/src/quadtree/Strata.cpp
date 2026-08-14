#include "../../include/quadtree/Strata.h"

#include <cstddef>
#include <utility>

namespace fmm {

Strata stratify_quadtree(Box &root) {
  Strata strata;
  strata.push_back({&root});

  std::size_t i = 0;
  while (i < strata.size()) {
    for (Box *box : strata[i]) {
      if (box != nullptr && box->has_child_boxes) {
        if (strata.size() == i + 1) {
          strata.push_back({});
        }
        for (auto &child : box->child_boxes) {
          strata[i + 1].push_back(child.get());
        }
      }
    }

    std::vector<Box *> filtered;
    filtered.reserve(strata[i].size());
    for (Box *box : strata[i]) {
      if (box != nullptr)
        filtered.push_back(box);
    }
    strata[i] = std::move(filtered);

    ++i;
  }

  return strata;
}

std::vector<Box *> leaves(const Strata &strata) {
  std::vector<Box *> result;
  for (const auto &stratum : strata) {
    for (Box *box : stratum) {
      if (!box->has_child_boxes) {
        result.push_back(box);
      }
    }
  }
  return result;
}

} // namespace fmm
