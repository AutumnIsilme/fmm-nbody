#include "../../include/quadtree/QuadtreeBuilder.h"
#include "../../include/quadtree/Colleagues.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
namespace fmm {

namespace {
constexpr double kMinBoxExtent = 1e-10;
} // namespace

void split(Box &box, const std::vector<Body> &bodies,
           std::size_t bodies_per_box) {
    const double split_size = box.extent / 2.0;
    std::vector<Box *> ancestors = box.ancestors;
    ancestors.push_back(&box);
    // Quadrant roots
    // 0 = bottom-left, 1 = top-left, 2 = bottom-right, 3 = top-right.
    const std::array<Vec2, 4> roots = {box.root,
                                       box.root + Vec2(0.0, split_size),
                                       box.root + Vec2(split_size, 0.0),
                                       box.root + Vec2(split_size, split_size)};
    std::array<std::unique_ptr<Box>, 4> owned_children;
    std::array<Box *, 4> raw_children{};
    for (int i = 0; i < 4; ++i) {
        Vec2 child_centre = roots[static_cast<std::size_t>(i)] +
                            Vec2(split_size, split_size) / 2.0;
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
        if (!std::isfinite(nx) || !std::isfinite(ny)) {
            throw std::runtime_error(
                "split: body id=" + std::to_string(body.id) +
                " has a non-finite position (x=" + std::to_string(body.x()) +
                ", y=" + std::to_string(body.y()) +
                ") -- likely a force blowup from an unresolved close "
                "encounter. "
                "Consider increasing softening (kSofteningSquared in "
                "SimConstants.h) or reducing dt.");
        }
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
        auto &bucket = sub_box_bodies[static_cast<std::size_t>(i)];
        if (bucket.size() > bodies_per_box) {
            if (split_size > kMinBoxExtent) {
                split(*sub_box, bucket, bodies_per_box);
            } else {
                std::cerr << "split: " << bucket.size()
                          << " bodies could not be separated below extent "
                          << split_size
                          << " (likely near-coincident positions from an "
                             "unresolved close encounter) -- keeping them in a "
                             "single leaf instead of recursing further\n";
                sub_box->bodies_in_box = std::move(bucket);
            }
        } else if (bucket.empty()) {
            keep_child[static_cast<std::size_t>(i)] = false;
        } else {
            sub_box->bodies_in_box = std::move(bucket);
        }
    }

    for (int i = 0; i < 4; ++i) {
        if (!keep_child[static_cast<std::size_t>(i)])
            continue;
        Box *sub_box = raw_children[static_cast<std::size_t>(i)];
        sub_box->colleagues = {keep_child[0] ? raw_children[0] : nullptr,
                               keep_child[1] ? raw_children[1] : nullptr,
                               keep_child[2] ? raw_children[2] : nullptr,
                               keep_child[3] ? raw_children[3] : nullptr,
                               nullptr,
                               nullptr,
                               nullptr,
                               nullptr};
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
    if (bodies.empty()) {
        const Vec2 box_0_root(-1.0, -1.0);
        const double box_0_extent = 2.0;
        return std::unique_ptr<Box>(
            new Box(box_0_root, box_0_extent, std::vector<Box *>{},
                    box_0_root + Vec2(box_0_extent, box_0_extent) / 2.0));
    }

    double min_x = bodies[0].x();
    double max_x = bodies[0].x();
    double min_y = bodies[0].y();
    double max_y = bodies[0].y();

    for (const auto &body : bodies) {
        if (std::isfinite(body.x()) && std::isfinite(body.y())) {
            min_x = std::min(min_x, body.x());
            max_x = std::max(max_x, body.x());
            min_y = std::min(min_y, body.y());
            max_y = std::max(max_y, body.y());
        }
    }

    double width = max_x - min_x;
    double height = max_y - min_y;
    double extent = std::max(width, height);

    if (extent < 1e-6) {
        extent = 2.0; // Default size if bodies are at a single point
    } else {
        extent *= 1.10;
    }

    static double smooth_extent = 2.0;

    extent = std::max(extent, 2.0);

    if (extent > smooth_extent) {
        smooth_extent = extent; // Instantly grow to fit escaping particles
    } else {
        smooth_extent = 0.98 * smooth_extent + 0.02 * extent; // Smooth shrink
    }

    extent = smooth_extent; // Use the smoothed extent for this frame

    double center_x = 0.5 * (min_x + max_x);
    double center_y = 0.5 * (min_y + max_y);

    Vec2 box_0_root(center_x - 0.5 * extent, center_y - 0.5 * extent);
    double box_0_extent = extent;

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
