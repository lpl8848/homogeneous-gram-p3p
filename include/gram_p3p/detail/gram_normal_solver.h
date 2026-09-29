#pragma once
#include "gram_normal_p3p.h"
// Research candidate frozen after horizontal testing on 2026-09-29.
// Camera-from-world R,t; arbitrary nonzero calibrated ray magnitudes;
// no front-hemisphere or small-angle assumption. Keep author baselines intact.
namespace gram_normal_solver {
inline int solve(const std::vector<Eigen::Vector3d>& rays,const std::vector<Eigen::Vector3d>& world,std::vector<MatrixPose>* poses){
    return gram_normal_p3p::solve<6,3,false,true,true,1>(rays,world,poses);
}
}
