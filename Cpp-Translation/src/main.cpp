#include <cstdlib>
#include <iostream>
#include <string>

#include "../include/fmm/FastMultipole.h"
#include "../include/live_view/LiveSettings.h"

void print_usage(const char *prog_name) {
  std::cout
      << "Usage: " << prog_name << " [options]\n\n"
      << "Options:\n"
      << "  -n, --bodies <int>          Number of particles N (default: 500)\n"
      << "  -b, --bodies-per-box <int>  Max particles per leaf box (default: "
         "25)\n"
      << "  -e, --epsilon <double>      Precision tolerance (default: 1e-7)\n"
      << "  -s, --steps <int>           Number of integration steps (default: "
         "5000)\n"
      << "      --dt <double>           Time step delta t (default: 1e-4)\n"
      << "  -r, --rebuild-every <int>   Quadtree rebuild frequency in steps "
         "(default: 1)\n"
      << "      --no-boxes              Disable background quadtree box "
         "rendering (default)\n"
      << "      --show-boxes            Enable background quadtree box "
         "rendering \n"
      << "  -h, --help                  Show this help message and exit\n";
}

int main(int argc, char **argv) {
  // Default values
  int N = 500;
  int bodies_per_box = 25;
  double epsilon = 1e-7;
  int num_steps = 5000;
  double dt = 1e-4;
  int rebuild_every = 1;
  bool show_boxes = false;

  // Parse command line arguments
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg == "-h" || arg == "--help") {
      print_usage(argv[0]);
      return 0;
    } else if ((arg == "-n" || arg == "--bodies") && i + 1 < argc) {
      N = std::stoi(argv[++i]);
    } else if ((arg == "-b" || arg == "--bodies-per-box") && i + 1 < argc) {
      bodies_per_box = std::stoi(argv[++i]);
    } else if ((arg == "-e" || arg == "--epsilon") && i + 1 < argc) {
      epsilon = std::stod(argv[++i]);
    } else if ((arg == "-s" || arg == "--steps") && i + 1 < argc) {
      num_steps = std::stoi(argv[++i]);
    } else if (arg == "--dt" && i + 1 < argc) {
      dt = std::stod(argv[++i]);
    } else if ((arg == "-r" || arg == "--rebuild-every") && i + 1 < argc) {
      rebuild_every = std::stoi(argv[++i]);
    } else if (arg == "--no-boxes") {
      show_boxes = false;
    } else if (arg == "--show-boxes") {
      show_boxes = true;
    } else {
      std::cerr << "Unrecognized or incomplete argument: " << arg << "\n";
      print_usage(argv[0]);
      return 1;
    }
  }

  // Print configuration summary
  std::cout << "===========================================\n"
            << "       2D Fast Multipole Simulation        \n"
            << "===========================================\n"
            << "  Particles (N)      : " << N << "\n"
            << "  Bodies per Box     : " << bodies_per_box << "\n"
            << "  Epsilon            : " << epsilon << "\n"
            << "  Simulation Steps   : " << num_steps << "\n"
            << "  Time Step (dt)     : " << dt << "\n"
            << "  Rebuild Frequency  : Every " << rebuild_every << " step(s)\n"
            << "  Render Quadtree    : "
            << (show_boxes ? "Enabled" : "Disabled") << "\n"
            << "===========================================\n\n";

  LiveViewSettings live;
  live.output_path = "live.svg";
  live.frame_stride = 1;
  live.show_boxes = show_boxes;

  fmm::run_fmm_simulation(N, bodies_per_box, epsilon, num_steps, dt,
                          rebuild_every, show_boxes, &live);

  return 0;
}
