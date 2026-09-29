#pragma once
#include <Eigen/Dense>
#include <vector>
struct MatrixPose {
    Eigen::Matrix3d rotation;
    Eigen::Vector3d t;
    MatrixPose(const Eigen::Matrix3d& r,const Eigen::Vector3d& v):rotation(r),t(v){}
    const Eigen::Matrix3d& R() const {return rotation;}
};
