#pragma once
#include <array>
#include <cstddef>
#include <vector>
namespace fmm {
class Body {
public:
  std::vector<double> data;

  // stable per-body id. Allows particles to be matched between steps
  std::size_t id = 0;

  Body();
  explicit Body(std::vector<double> data_);
  Body(double x, double y);
  Body(double x, double y, double charge);
  double x() const;
  double y() const;
  // Returns 0.0 if no charge/mass column was supplied.
  double mass() const;
  double charge() const;
  std::array<double, 2> position() const;
  double operator[](std::size_t index) const;

  void set_position(double new_x, double new_y) {
    data[0] = new_x;
    data[1] = new_y;
  }
};
} // namespace fmm
