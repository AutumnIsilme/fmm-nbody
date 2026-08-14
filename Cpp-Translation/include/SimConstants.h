#pragma once

namespace SimConstants {
// smaller = more accurate / slower
constexpr double kTimestepSafetyFactor = 0.05;
constexpr int kMaxSubstepsPerFrame =
    1000; // bail out of subdividing after this many
constexpr double kEscapeRadiusMultiplier =
    1.5; // after a body has moved this many domains-widths away, treat it as
         // not coming back
constexpr double kSofteningSquared = 0.000225; // 0.015^2
constexpr double kMinSubstepDt = 1e-7;

} // namespace SimConstants
