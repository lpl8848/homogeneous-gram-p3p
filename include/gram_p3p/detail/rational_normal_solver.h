#pragma once
#include "rational_normal_p3p.h"

// RN-P3P research candidate, 2026-09-29.
// All calibrated bearing directions are supported, including negative z.
// World points must be distinct and noncollinear. Output is camera-from-world.
// Numerical norm failures use the complementary circle chart. An ordinary
// empty physical solution set is returned directly.
namespace rational_normal_solver {
using Vector=Eigen::Vector3d;
inline int solve_unit(const std::vector<Vector>& bearings,const std::vector<Vector>& points,std::vector<MatrixPose>* poses){
    return rational_normal_p3p::solve<true,2,1,0,false,false>(bearings,points,poses);
}
// Convenience API for arbitrary nonzero bearing vectors. Normalization is
// performed inside this call and must be included when timing this API.
inline int solve(const std::vector<Vector>& bearings,const std::vector<Vector>& points,std::vector<MatrixPose>* poses){
    return rational_normal_p3p::solve<true,2,1,0,false,true>(bearings,points,poses);
}
}
