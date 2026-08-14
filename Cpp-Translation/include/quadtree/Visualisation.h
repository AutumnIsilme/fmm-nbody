#pragma once

#include <sstream>
#include <string>

#include "Box.h"

namespace fmm {

// Recursively draws the outline of `box` and every descendant box.
void draw_box_square(const Box &box, std::ostringstream &svg,
                     const std::string &colour = "blue",
                     bool show_boxes = false);

// Recursively draws leaf outlines and the bodies contained in each leaf
// (stem boxes contribute only their descendants' outlines, never their
// own).
void graph_box(const Box &box, std::ostringstream &svg,
               const std::string &box_colour = "blue",
               const std::string &particle_colour = "dimgray",
               bool show_boxes = false);

// Renders the full quadtree (leaf outlines + bodies) to an SVG file at
// `filepath`.
// The view box assumes the fixed [-1, 1] x [-1, 1] domain
// used by `create_quadtree.
void graph_quadtree(const Box &root, const std::string &filepath,
                    const std::string &box_colour = "blue",
                    const std::string &particle_colour = "dimgray",
                    bool show_boxes = false);

} // namespace fmm
