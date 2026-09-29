#pragma once
#include <cmath>
#include <cstdint>
#include <cstring>
namespace fast_radical {
// Binary exponent seed followed by Halley and fourth-order Householder.
// Full double refinement, not a reduced-accuracy pose approximation.
// Very small/large inputs retain the platform routine to avoid squaring
// underflow/overflow in the homogeneous fourth-order expression.
inline double cbrt(double value){
    double a=std::abs(value);if(!(a>=1e-100&&a<=1e100))return std::cbrt(value);
    std::uint64_t bits;std::memcpy(&bits,&a,8);bits=bits/3+UINT64_C(0x2a9f76251691a600);
    double y;std::memcpy(&y,&bits,8);double z=y*y*y;
    y*= (z+2*a)/(2*z+a);z=y*y*y;
    double aa=a*a,az=a*z,zz=z*z;
    y*= (4*aa+19*az+4*zz)/(aa+16*az+10*zz);
    return std::copysign(y,value);
}
}
