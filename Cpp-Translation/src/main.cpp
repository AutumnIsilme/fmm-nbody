#include <cstdlib>
#include <iostream>

#include "../include/fmm/FastMultipole.h"
#include "../include/live_view/LiveSettings.h"

int main(int argc, char **argv) {
  int N = (argc > 1) ? std::atoi(argv[1]) : 500;
  int bodies_per_box = (argc > 2) ? std::atoi(argv[2]) : 25;
  double epsilon = (argc > 3) ? std::atof(argv[3]) : 1e-7;

  std::cout << "Running FMM: N=" << N << " bodies_per_box=" << bodies_per_box
            << " epsilon=" << epsilon << "\n";

  fmm::run_fmm(N, bodies_per_box, epsilon);
  LiveViewSettings live;
  live.output_path = "live/live.svg";
  live.frame_stride = 1;
  fmm::run_fmm_simulation(N, bodies_per_box, epsilon, 100000, 1e-3, 10, &live);

  return 0;
}
