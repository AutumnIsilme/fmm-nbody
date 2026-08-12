#pragma once

namespace SimConstants {
// smaller = more accurate / slower
constexpr double kTimestepSafetyFactor = 0.2;
constexpr int kMaxSubstepsPerFrame =
    200; // bail out of subdividing after this many
constexpr double kEscapeRadiusMultiplier =
    5.0; // after a body has moved this many domains-widths away, treat it as
         // not coming back
// I saw a post saying this was sensible for more than one time step to reduce
// things blowing up
constexpr double kSofteningSquared = 1e-4;

} // namespace SimConstants
