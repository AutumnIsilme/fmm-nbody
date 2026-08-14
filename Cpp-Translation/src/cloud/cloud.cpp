#include "../../include/cloud/Cloud.h"

#include <fstream>
#include <random>
#include <string>
#include <vector>

using namespace std;

namespace fmm {

vector<string> lines(const vector<vector<double>> &bodies) {
  vector<string> result;

  for (const auto &body : bodies) {
    result.push_back(to_string(body[0]) + ", " + to_string(body[1]) + ", " +
                     to_string(body[2]) + ", " + to_string(body[3]) + ", " +
                     to_string(body[4]) + "\n");
  }

  return result;
}

void write_points_and_masses(const vector<vector<double>> &bodies,
                             const string &filename) {
  ofstream file(filename);

  for (const auto &line : lines(bodies)) {
    file << line;
  }
}

vector<vector<double>> generate_2d_bodies_uniform_random(int n_points,
                                                         double max_mass) {
  vector<vector<double>> data;

  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> random(0.0, 1.0);

  for (int i = 0; i < n_points; i++) {
    data.push_back({random(gen) * 2 - 1, random(gen) * 2 - 1, 0.0, 0.0, 0.0,
                    random(gen) * max_mass});
  }

  return data;
}
} // namespace fmm
