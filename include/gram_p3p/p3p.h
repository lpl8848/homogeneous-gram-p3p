#pragma once
#include "detail/gram_normal_solver.h"
namespace gram_p3p {
using Pose = ::MatrixPose;
// Exactly three finite world points and three nonzero calibrated bearings.
// Camera-from-world convention: R * X[i] + t = positive_scale[i] * bearing[i].
// Output is overwritten. Returns its size, at most four.
inline int solve(const std::vector<Eigen::Vector3d>& bearings,
                 const std::vector<Eigen::Vector3d>& world,
                 std::vector<Pose>* poses) {
    return gram_normal_solver::solve(bearings, world, poses);
}
}
