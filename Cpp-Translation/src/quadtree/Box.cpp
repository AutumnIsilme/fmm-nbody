#include "quadtree/Box.h"

#include <sstream>

namespace fmm {

Box::Box()
    : root(), extent(0.0), has_child_boxes(false), child_boxes(),
      bodies_in_box(), colleagues{}, U(), V(), W(), X(), centre(),
      multipole_expansion{}, local_expansion{}, forces() {}

Box::Box(Vec2 root_, double extent_, Vec2 centre_)
    : root(root_), extent(extent_),
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
    oss << ", has_child_boxes=" << (has_child_boxes ? "true" : "false")
        << ", num_child_boxes=" << num_child_boxes()
        << ", num_bodies_in_box=" << bodies_in_box.size() << "]";
    return oss.str();
}

std::ostream &operator<<(std::ostream &os, const Box &box) {
    return os << box.repr();
}

Box *BoxAllocator::make_box() {
    Box *box = &boxes[alloc_point];
    *box = {};
    alloc_point++;
    return box;
}

Box *BoxAllocator::make_box(Box data) {
    Box *box = &boxes[alloc_point];
    *box = std::move(data);
    alloc_point++;
    return box;
}

void BoxAllocator::clear_allocator() {
    alloc_point = 0;
}

} // namespace fmm
