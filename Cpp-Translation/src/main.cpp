#include <cstdlib>
#include <iostream>

#include "../include/fmm/FastMultipole.h"

int main(int argc, char **argv) {
  int N = (argc > 1) ? std::atoi(argv[1]) : 500;
  int bodies_per_box = (argc > 2) ? std::atoi(argv[2]) : 25;
  double epsilon = (argc > 3) ? std::atof(argv[3]) : 1e-7;

  std::cout << "Running FMM: N=" << N << " bodies_per_box=" << bodies_per_box
            << " epsilon=" << epsilon << "\n";

  fmm::run_fmm(N, bodies_per_box, epsilon);

  return 0;
}
