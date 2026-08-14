#pragma once
#include <string>

struct LiveViewSettings {
  std::string output_path = "live.svg";
  int frame_stride = 1;
  std::string box_colour = "#333333";
  std::string particle_colour = "#e63946";
};
