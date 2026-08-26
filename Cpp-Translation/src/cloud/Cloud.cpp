#include "cloud/Cloud.h"

#include <random>
#include <vector>

using namespace std;

namespace fmm {
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
