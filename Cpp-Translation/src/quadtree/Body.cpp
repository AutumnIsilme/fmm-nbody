#include "quadtree/Body.h"

namespace fmm {

Body::Body() = default;

// Body::Body(std::vector<double> data_) {
//     if (data_.size() < 2) {
//         throw std::invalid_argument(
//             "Body must have at least an x and y coordinate");
//     }
//
//     x = data_[0];
//     y = data_[1];
//     if (data_.size() >= 4) {
//         vx = data_[2];
//         vy = data_[3];
//     }
//     if (data_.size() >= 5) {
//         mass = data_[4];
//     }
// }

// Body::Body(double x, double y) : data{x, y} {}

Body::Body(double x, double y, double vx, double vy, double mass)
    : x(x), y(y), vx(vx), vy(vy), mass(mass) {}

// double Body::x() const { return data[0]; }

// double Body::y() const { return data[1]; }

// double Body::mass() const { return data.size() > 2 ? data.back() : 0.0; }

// std::array<double, 2> Body::position() const { return {data[0], data[1]}; }

// double Body::operator[](std::size_t index) const { return data[index]; }

} // namespace fmm
