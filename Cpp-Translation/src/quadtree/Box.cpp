#include "../../include/quadtree/Box.h"

#include <sstream>
#include <utility>

namespace fmm {

Box::Box()
    : root(), extent(0.0), ancestors(), has_child_boxes(false), child_boxes(),
      bodies_in_box(), colleagues{}, U(), V(), W(), X(), centre(),
      multipole_expansion{}, local_expansion{}, forces() {}

Box::Box(Vec2 root_, double extent_, std::vector<Box *> ancestors_,
         Vec2 centre_)
    : root(root_), extent(extent_), ancestors(std::move(ancestors_)),
      has_child_boxes(false), child_boxes(), bodies_in_box(), colleagues{}, U(),
      V(), W(), X(), centre(centre_), multipole_expansion{}, local_expansion{},
      forces() {}

Box::~Box() = default;

Box::Box(Box &&) noexcept = default;
Box &Box::operator=(Box &&) noexcept = default;

bool Box::operator==(const Box &other) const {
  return root == other.root && extent == other.extent;
}

bool Box::operator!=(const Box &other) const { return !(*this == other); }

int Box::num_child_boxes() const {
  if (!has_child_boxes)
    return 0;
  int count = 0;
  for (const auto &child : child_boxes) {
    if (child)
      ++count;
  }
  return count;
}

std::string Box::repr() const {
  std::ostringstream oss;
  oss << "Box[root=(" << root.x << ", " << root.y << ")"
      << ", extent=" << extent << ", parent=";
  if (ancestors.empty()) {
    oss << "None";
  } else {
    const Box *parent = ancestors.back();
    oss << "(" << parent->root.x << ", " << parent->root.y << ")x"
        << parent->extent;
  }
  oss << ", has_child_boxes=" << (has_child_boxes ? "true" : "false")
      << ", num_child_boxes=" << num_child_boxes()
      << ", num_bodies_in_box=" << bodies_in_box.size() << "]";
  return oss.str();
}

std::ostream &operator<<(std::ostream &os, const Box &box) {
  return os << box.repr();
}

} // namespace fmm
