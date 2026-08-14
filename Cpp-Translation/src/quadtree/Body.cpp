#include "../../include/quadtree/Body.h"

#include <stdexcept>
#include <utility>

namespace fmm {

Body::Body() = default;

Body::Body(std::vector<double> data_) : data(std::move(data_)) {
  if (data.size() < 2) {
    throw std::invalid_argument(
        "Body must have at least an x and y coordinate");
  }
}

Body::Body(double x, double y) : data{x, y} {}

double Body::x() const { return data[0]; }

double Body::y() const { return data[1]; }

double Body::mass() const { return data.size() > 2 ? data.back() : 0.0; }

std::array<double, 2> Body::position() const { return {data[0], data[1]}; }

double Body::operator[](std::size_t index) const { return data[index]; }

} // namespace fmm
