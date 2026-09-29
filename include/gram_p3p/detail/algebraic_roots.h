#pragma once
#include <algorithm>
#include <cmath>
#include "PoseLib/misc/univariate.h"
namespace chordal {
inline int algebraic_roots(const double *c,double *out) {
    double scale=0;for(int i=0;i<5;++i)scale=std::max(scale,std::abs(c[i]));
    if(std::abs(c[4])<1e-14*scale) {
        if(std::abs(c[3])>1e-14*scale) return poselib::univariate::solve_cubic_real(c[2]/c[3],c[1]/c[3],c[0]/c[3],out);
        return poselib::univariate::solve_quadratic_real(c[2],c[1],c[0],out);
    }
    double inv=1/c[4],b=c[3]*inv,d=c[1]*inv,e=c[0]*inv,a=c[2]*inv;
    // Ferrari's resolvent must not divide by zero for an exact biquadratic.
    double pq=d-.5*b*a+.125*b*b*b;
    if(pq==0) {
        double pp=a-.375*b*b,rr=e-.25*b*d+.0625*b*b*a-3*b*b*b*b/256;
        double z[2];int nz=poselib::univariate::solve_quadratic_real(1,pp,rr,z),n=0;
        for(int i=0;i<nz;++i)if(z[i]>=0){double s=std::sqrt(z[i]);out[n++]=s-.25*b;if(s>0)out[n++]=-s-.25*b;}
        return n;
    }
    return poselib::univariate::solve_quartic_real(b,a,d,e,out);
}
}
