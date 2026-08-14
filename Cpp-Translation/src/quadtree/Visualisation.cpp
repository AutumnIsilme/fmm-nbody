#include "../../include/quadtree/Visualisation.h"

#include <fstream>
#include <iostream>

namespace fmm {

namespace {

void write_box_outline(std::ostringstream &svg, const Box &box,
                       const std::string &colour, bool show_boxes) {
  if (!show_boxes)
    return;
  svg << "  <rect x=\"" << box.root.x << "\" y=\"" << box.root.y
      << "\" width=\"" << box.extent << "\" height=\"" << box.extent
      << "\" fill=\"none\" stroke=\"" << colour
      << "\" stroke-width=\"0.0015\" />\n";
}

} // namespace

void draw_box_square(const Box &box, std::ostringstream &svg,
                     const std::string &colour, bool show_boxes) {
  if (box.has_child_boxes) {
    for (const auto &child : box.child_boxes) {
      if (child)
        draw_box_square(*child, svg, colour, show_boxes);
    }
  }
  write_box_outline(svg, box, colour, show_boxes);
}

void graph_box(const Box &box, std::ostringstream &svg,
               const std::string &box_colour,
               const std::string &particle_colour, bool show_boxes) {
  if (box.has_child_boxes) {
    for (const auto &child : box.child_boxes) {
      if (child)
        graph_box(*child, svg, box_colour, particle_colour, show_boxes);
    }
    return;
  }

  write_box_outline(svg, box, box_colour, show_boxes);
  for (const Body &body : box.bodies_in_box) {
    svg << "  <circle cx=\"" << body.x() << "\" cy=\"" << body.y()
        << "\" r=\"0.003\" fill=\"" << particle_colour << "\" />\n";
  }
}

void graph_quadtree(const Box &root, const std::string &filepath,
                    const std::string &box_colour,
                    const std::string &particle_colour, bool show_boxes) {
  std::ostringstream svg;
  svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"-1.05 -1.05 2.1 "
         "2.1\">\n";
  svg << "<g transform=\"scale(1,-1)\">\n";
  graph_box(root, svg, box_colour, particle_colour, show_boxes);
  svg << "</g>\n</svg>\n";

  std::ofstream out(filepath);
  out << svg.str();
}

} // namespace fmm
