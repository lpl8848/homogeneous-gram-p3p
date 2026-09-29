#include <gram_p3p/p3p.h>
#include "../simulations/data.h"
#include <iostream>
#include <limits>

bool check(const simulation::Scene& scene, bool scale) {
    auto bearings=scene.bearings;
    if(scale) {bearings[0]*=.2;bearings[1]*=3.;bearings[2]*=1.5;}
    std::vector<gram_p3p::Pose> poses;
    gram_p3p::solve(bearings,scene.world,&poses);
    if(poses.empty() || poses.size()>4) return false;
    double best=std::numeric_limits<double>::infinity();
    for(const auto& p:poses) {
        if(!p.R().allFinite() || !p.t.allFinite()) return false;
        if((p.R().transpose()*p.R()-Eigen::Matrix3d::Identity()).norm()>1e-7 ||
            std::abs(p.R().determinant()-1)>1e-7) return false;
        for(int i=0;i<3;++i) {
            const Eigen::Vector3d y=p.R()*scene.world[i]+p.t;
            const double s=y.dot(bearings[i])/bearings[i].squaredNorm();
            if(s<=0 || (y-s*bearings[i]).norm()>1e-7*std::max(1.,y.norm())) return false;
        }
        best=std::min(best,(p.R()-scene.R).norm()+(p.t-scene.t).norm());
    }
    return best<1e-6;
}
int main() {
    simulation::Generator generator(17);
    for(int i=0;i<1000;++i) {
        const auto scene=generator.draw(0,i%2!=0);
        if(!check(scene,false) || !check(scene,true)) {
            std::cerr<<"Geometry test failed at input "<<i<<'\n';return 1;
        }
    }
    std::cout<<"Passed 1000 perspective/full-sphere scenes with independent ray scales.\n";
    return 0;
}
