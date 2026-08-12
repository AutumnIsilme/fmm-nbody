#include <cstdlib>
#include <iostream>

#include "../include/fmm/FastMultipole.h"
#include "../include/live_view/LiveSettings.h"

int main(int argc, char **argv) {
  int N = (argc > 1) ? std::atoi(argv[1]) : 500;
  int bodies_per_box = (argc > 2) ? std::atoi(argv[2]) : 25;
  double epsilon = (argc > 3) ? std::atof(argv[3]) : 1e-7;
  int num_steps = (argc > 4) ? std::atoi(argv[4]) : 5000;
  double dt = (argc > 5) ? std::atof(argv[5]) : 1e-3;
  int rebuild_every = (argc > 6) ? std::atoi(argv[6]) : 10;

  std::cout << "N = " << N << ", bodies_per_box = " << bodies_per_box << "\n";
  std::cout << "Running FMM: N=" << N << " bodies_per_box=" << bodies_per_box
            << " epsilon=" << epsilon << "\n";

  // TODO debugging settings
  N = 50;
  bodies_per_box = 5;
  num_steps = 100000;
  dt = 1e-5;

  LiveViewSettings live;
  live.output_path = "live.svg";
  live.frame_stride = 1;
  fmm::run_fmm_simulation(N, bodies_per_box, epsilon, num_steps, dt,
                          rebuild_every, &live);

  return 0;
}
