#include "../../include/quadtree/Vec2.h"

namespace fmm {

Vec2::Vec2() : x(0.0), y(0.0) {}

Vec2::Vec2(double x_, double y_) : x(x_), y(y_) {}

Vec2 Vec2::operator+(const Vec2 &other) const {
  return Vec2(x + other.x, y + other.y);
}

Vec2 Vec2::operator-(const Vec2 &other) const {
  return Vec2(x - other.x, y - other.y);
}

Vec2 Vec2::operator*(double scalar) const {
  return Vec2(x * scalar, y * scalar);
}

Vec2 Vec2::operator/(double scalar) const {
  return Vec2(x / scalar, y / scalar);
}

bool Vec2::operator==(const Vec2 &other) const {
  return x == other.x && y == other.y;
}

bool Vec2::operator!=(const Vec2 &other) const { return !(*this == other); }

} // namespace fmm
