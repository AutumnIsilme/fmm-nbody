#include "../../include/quadtree/Colleagues.h"

#include <cstddef>

namespace fmm {

namespace {

Box *child_or_self(Box *colleague, int child_index) {
  if (colleague != nullptr && colleague->has_child_boxes) {
    return colleague->child_boxes[static_cast<std::size_t>(child_index)].get();
  }
  return colleague;
}

} // namespace

void colleagify_quadtree(Box &box) {
  if (box.child_boxes[0]) {
    box.child_boxes[0]->colleagues = {box.child_boxes[1].get(),
                                      box.child_boxes[3].get(),
                                      box.child_boxes[2].get(),
                                      child_or_self(box.colleagues[4], 3),
                                      child_or_self(box.colleagues[4], 1),
                                      child_or_self(box.colleagues[5], 3),
                                      child_or_self(box.colleagues[6], 2),
                                      child_or_self(box.colleagues[6], 3)};
  }
  if (box.child_boxes[1]) {
    box.child_boxes[1]->colleagues = {child_or_self(box.colleagues[0], 0),
                                      child_or_self(box.colleagues[0], 2),
                                      box.child_boxes[3].get(),
                                      box.child_boxes[2].get(),
                                      box.child_boxes[0].get(),
                                      child_or_self(box.colleagues[6], 2),
                                      child_or_self(box.colleagues[6], 3),
                                      child_or_self(box.colleagues[7], 2)};
  }
  if (box.child_boxes[2]) {
    box.child_boxes[2]->colleagues = {box.child_boxes[3].get(),
                                      child_or_self(box.colleagues[2], 1),
                                      child_or_self(box.colleagues[2], 0),
                                      child_or_self(box.colleagues[3], 1),
                                      child_or_self(box.colleagues[4], 3),
                                      child_or_self(box.colleagues[4], 1),
                                      box.child_boxes[0].get(),
                                      box.child_boxes[1].get()};
  }
  if (box.child_boxes[3]) {
    box.child_boxes[3]->colleagues = {child_or_self(box.colleagues[0], 2),
                                      child_or_self(box.colleagues[1], 0),
                                      child_or_self(box.colleagues[2], 1),
                                      child_or_self(box.colleagues[2], 0),
                                      box.child_boxes[2].get(),
                                      box.child_boxes[0].get(),
                                      box.child_boxes[1].get(),
                                      child_or_self(box.colleagues[0], 0)};
  }

  for (auto &child : box.child_boxes) {
    if (child && child->has_child_boxes) {
      colleagify_quadtree(*child);
    }
  }
}

} // namespace fmm
