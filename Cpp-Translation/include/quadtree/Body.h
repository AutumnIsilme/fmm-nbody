#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace fmm {

class Body {
public:
  // Raw row of data: [x, y, charge/mass, ...]. Mirrors a numpy row.
  std::vector<double> data;

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
};

} // namespace fmm
