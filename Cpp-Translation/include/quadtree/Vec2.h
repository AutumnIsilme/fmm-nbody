#pragma once

// Implements a 2D vector with pairwise operations

namespace fmm {

struct Vec2 {
  double x;
  double y;

  Vec2();
  Vec2(double x_, double y_);

  Vec2 operator+(const Vec2 &other) const;
  Vec2 operator-(const Vec2 &other) const;
  Vec2 operator*(double scalar) const;
  Vec2 operator/(double scalar) const;

  bool operator==(const Vec2 &other) const;
  bool operator!=(const Vec2 &other) const;
};

} // namespace fmm
