#pragma once
#include <cstddef>
#include <vector>
namespace fmm {
struct Body {
  public:
    double x;
    double y;
    double vx;
    double vy;
    double mass;

    // stable per-body id. Allows particles to be matched between steps
    std::size_t id = 0;

    Body();
    explicit Body(std::vector<double> data_);
    Body(double x, double y);
    Body(double x, double y, double vx, double vy, double mass);

    // double x() const;
    // double y() const;
    // double mass() const;
    // std::array<double, 2> position() const;
    // double operator[](std::size_t index) const;

    void set_position(double new_x, double new_y) {
        x = new_x;
        y = new_y;
    }
};
} // namespace fmm
