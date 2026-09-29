#include <gram_p3p/p3p.h>
#include <iostream>

int main() {
    const std::vector<Eigen::Vector3d> world = {{0,0,0}, {1,0,0}, {.2,.8,.1}};
    const Eigen::Matrix3d truth_R = Eigen::AngleAxisd(.3,
        Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
    const Eigen::Vector3d truth_t(.1,-.2,3);
    std::vector<Eigen::Vector3d> bearings;
    for (const auto& point : world) bearings.push_back(truth_R * point + truth_t);
    std::vector<gram_p3p::Pose> poses;
    gram_p3p::solve(bearings, world, &poses);
    std::cout << "Returned " << poses.size() << " camera-from-world poses\n";
    for (const auto& pose : poses)
        std::cout << "R:\n" << pose.R() << "\nt: " << pose.t.transpose() << "\n";
    return poses.empty() ? 1 : 0;
}
