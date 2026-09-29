#pragma once
#include <Eigen/Dense>
#include <random>
#include <vector>

namespace simulation {
using V = Eigen::Vector3d;
struct Scene {
    Eigen::Matrix3d R;
    V t;
    std::vector<V> bearings, world;
};
class Generator {
    std::mt19937_64 random;
    std::normal_distribution<double> normal{0,1};
    std::uniform_real_distribution<double> uniform{-1,1};
    V gaussian() {
        // Explicit draw order avoids unspecified function-argument evaluation.
        const double x=normal(random), y=normal(random), z=normal(random);
        return V(x,y,z);
    }
public:
    explicit Generator(std::uint64_t seed) : random(seed) {}
    Scene draw(double epsilon=0, bool sphere=false) {
        for (;;) {
            Scene scene;
            const double w=normal(random), x=normal(random),
                         y=normal(random), z=normal(random);
            Eigen::Quaterniond q(w,x,y,z); q.normalize();
            scene.R=q.toRotationMatrix(); scene.t=gaussian();
            scene.bearings.resize(3); scene.world.resize(3);
            for (int i=0;i<3;++i) {
                V ray;
                if (sphere) ray=gaussian().normalized();
                else {const double u=uniform(random), v=uniform(random);
                    ray=V(u,v,1).normalized();}
                const double range=.1+9.9*(uniform(random)+1)*.5;
                scene.bearings[i]=ray;
                scene.world[i]=scene.R.transpose()*(range*ray-scene.t);
            }
            if (epsilon>0) {
                const V baseline=scene.world[1]-scene.world[0];
                const V perpendicular=baseline.cross(gaussian());
                if (perpendicular.squaredNorm()==0) continue;
                const double alpha=.2+.6*(uniform(random)+1)*.5;
                scene.world[2]=scene.world[0]+alpha*baseline+
                    epsilon*baseline.norm()*perpendicular.normalized();
                const V point=scene.R*scene.world[2]+scene.t;
                if (!sphere && point.z()<=0) continue;
                scene.bearings[2]=point.normalized();
            }
            if ((scene.world[1]-scene.world[0]).cross(
                    scene.world[2]-scene.world[0]).squaredNorm()==0) continue;
            if (!sphere) {
                const Eigen::Vector2d a=scene.bearings[0].hnormalized(),
                    b=scene.bearings[1].hnormalized(),c=scene.bearings[2].hnormalized();
                if ((b.x()-a.x())*(c.y()-a.y())-
                    (b.y()-a.y())*(c.x()-a.x())==0) continue;
            }
            return scene;
        }
    }
};
}
