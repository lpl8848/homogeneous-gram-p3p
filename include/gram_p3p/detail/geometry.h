#pragma once
#include <Eigen/Dense>
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include "matrix_pose.h"
namespace gram_geometry {
using EV=Eigen::Vector3d;
struct V{
    double x,y,z;
    V():x(0),y(0),z(0){} V(double a,double b,double c):x(a),y(b),z(c){}
    V(const EV& v):x(v[0]),y(v[1]),z(v[2]){}
    double operator[](int i)const{return i==0?x:i==1?y:z;}
    V operator+(const V& v)const{return V(x+v.x,y+v.y,z+v.z);}
    V operator-(const V& v)const{return V(x-v.x,y-v.y,z-v.z);}
    V operator-()const{return V(-x,-y,-z);}
    V operator*(double a)const{return V(x*a,y*a,z*a);}
    V operator/(double a)const{double inv=1/a;return V(x*inv,y*inv,z*inv);}
    double dot(const V& v)const{return x*v.x+y*v.y+z*v.z;}
    double squaredNorm()const{return x*x+y*y+z*z;}
    V cross(const V& v)const{return V(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x);}
    V normalized()const{return *this/std::sqrt(squaredNorm());}
    EV eigen()const{return EV(x,y,z);}
};
inline V operator*(double a,const V& v){return v*a;}
}
