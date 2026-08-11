#include "../../include/quadtree/QuadtreeBuilder.h"

#include <array>
#include <utility>

#include "../../include/quadtree/Colleagues.h"

namespace fmm {

void split(Box &box, const std::vector<Body> &bodies,
           std::size_t bodies_per_box) {
  const double split_size = box.extent / 2.0;

  std::vector<Box *> ancestors = box.ancestors;
  ancestors.push_back(&box);

  // Quadrant roots
  // 0 = bottom-left, 1 = top-left, 2 = bottom-right, 3 = top-right.
  const std::array<Vec2, 4> roots = {box.root, box.root + Vec2(0.0, split_size),
                                     box.root + Vec2(split_size, 0.0),
                                     box.root + Vec2(split_size, split_size)};

  std::array<std::unique_ptr<Box>, 4> owned_children;
  std::array<Box *, 4> raw_children{};
  for (int i = 0; i < 4; ++i) {
    Vec2 child_centre =
        roots[static_cast<std::size_t>(i)] + Vec2(split_size, split_size) / 2.0;
    owned_children[static_cast<std::size_t>(i)] =
        std::unique_ptr<Box>(new Box(roots[static_cast<std::size_t>(i)],
                                     split_size, ancestors, child_centre));
    raw_children[static_cast<std::size_t>(i)] =
        owned_children[static_cast<std::size_t>(i)].get();
  }

  // Partition bodies by quadrant.
  std::array<std::vector<Body>, 4> sub_box_bodies;
  for (const Body &body : bodies) {
    const double nx = body.x() - box.root.x;
    const double ny = body.y() - box.root.y;
    if (nx < split_size && ny < split_size) {
      sub_box_bodies[0].push_back(body);
    } else if (nx < split_size && ny >= split_size) {
      sub_box_bodies[1].push_back(body);
    } else if (nx >= split_size && ny < split_size) {
      sub_box_bodies[2].push_back(body);
    } else if (nx >= split_size && ny >= split_size) {
      sub_box_bodies[3].push_back(body);
    }
  }

  std::array<bool, 4> keep_child = {true, true, true, true};
  for (int i = 0; i < 4; ++i) {
    Box *sub_box = raw_children[static_cast<std::size_t>(i)];

    // Temporary placeholder colleagues (siblings only).
    sub_box->colleagues = {raw_children[0], raw_children[1], raw_children[2],
                           raw_children[3], nullptr,         nullptr,
                           nullptr,         nullptr};

    auto &bucket = sub_box_bodies[static_cast<std::size_t>(i)];
    if (bucket.size() > bodies_per_box) {
      split(*sub_box, bucket, bodies_per_box);
    } else if (bucket.empty()) {
      keep_child[static_cast<std::size_t>(i)] = false;
    } else {
      sub_box->bodies_in_box = std::move(bucket);
    }
  }

  box.has_child_boxes = true;
  for (int i = 0; i < 4; ++i) {
    box.child_boxes[static_cast<std::size_t>(i)] =
        keep_child[static_cast<std::size_t>(i)]
            ? std::move(owned_children[static_cast<std::size_t>(i)])
            : nullptr;
  }
}

std::unique_ptr<Box> create_quadtree(const std::vector<Body> &bodies,
                                     std::size_t bodies_per_box) {
  const Vec2 box_0_root(-1.0, -1.0);
  const double box_0_extent = 2.0;
  auto box_0 = std::unique_ptr<Box>(
      new Box(box_0_root, box_0_extent, std::vector<Box *>{},
              box_0_root + Vec2(box_0_extent, box_0_extent) / 2.0));

  if (bodies.size() <= bodies_per_box) {
    box_0->bodies_in_box = bodies;
  } else {
    split(*box_0, bodies, bodies_per_box);
    colleagify_quadtree(*box_0);
  }

  return box_0;
}

} // namespace fmm
