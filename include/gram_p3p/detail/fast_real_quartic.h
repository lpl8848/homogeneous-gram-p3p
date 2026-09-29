#pragma once
#include <cmath>
#include <algorithm>
#include "algebraic_roots.h"
namespace fast_real_quartic {
inline double largest_cubic(double b,double c,double d){
    double shift=b/3,p=c-b*shift,q=(2*shift*shift-c)*shift+d,z;
    if(p<0){
        double r=std::sqrt(-p/3),h=q-2*r*r*r;
        if(h<=0)z=r+std::sqrt(-h/(3*r));
        else z=-r-std::sqrt((q+2*r*r*r)/(3*r));
    }else{
        double bound=std::sqrt(1+std::abs(q));if(p>0)bound=std::min(bound,std::abs(q)/p);
        z=std::copysign(bound,-q);
    }
    for(int i=0;i<10;++i){
        double f=(z*z+p)*z+q,df=3*z*z+p;
        double bound=std::abs(z*z*z)+std::abs(p*z)+std::abs(q);
        if(std::abs(f)<1e-15*bound)return z-shift;
        double den=df*df-3*z*f,step=den>0?f*df/den:f/df;
        z-=step;if(std::abs(step)<2e-16*(1+std::abs(z)))return z-shift;
    }
    double rs[3];int n=poselib::univariate::solve_cubic_real(b,c,d,rs);
    return *std::max_element(rs,rs+n);
}
inline int quadratic(double b,double c,double* out){
    double D=b*b-4*c;if(D<0)return 0;
    double q=-.5*(b+std::copysign(std::sqrt(D),b));
    if(q==0){out[0]=0;return 1;}out[0]=q;out[1]=c/q;return D==0?1:2;
}
inline int solve(const double* c,double* out){
    double scale=0;for(int i=0;i<5;++i)scale=std::max(scale,std::abs(c[i]));
    if(std::abs(c[4])<1e-14*scale)return chordal::algebraic_roots(c,out);
    double inv=1/c[4],b=c[3]*inv,a=c[2]*inv,d=c[1]*inv,e=c[0]*inv,shift=b*.25;
    double p=a-.375*b*b,q=d-.5*b*a+.125*b*b*b,r=e-.25*b*d+.0625*b*b*a-3*b*b*b*b/256;
    if(q==0){
        // Keep the zero root of z^2+p*z+r and its companion explicitly.
        // The legacy quadratic implementation computes 0/0 when r==0.
        double zz[2];int nz=quadratic(p,r,zz),n=0;
        for(int i=0;i<nz;++i)if(zz[i]>=0){double v=std::sqrt(zz[i]);out[n++]=v-shift;if(v>0)out[n++]=-v-shift;}
        return n;
    }
    double m=largest_cubic(2*p,p*p-4*r,-q*q);
    if(!(m>0)||!std::isfinite(m))return chordal::algebraic_roots(c,out);
    double s=std::sqrt(m),qs=q/s,u=.5*(p+m-qs),v=.5*(p+m+qs);
    int n=quadratic(s,u,out);n+=quadratic(-s,v,out+n);
    for(int i=0;i<n;++i)out[i]-=shift;
    return n;
}
}
