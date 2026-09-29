#pragma once
#include "sparse_householder.h"
namespace rotor_householder {
using V2=Eigen::Vector2d;
inline void circle_to_rotor(double x,double y,double& h,double& k){
    if(x>=0){h=std::sqrt(.5*(1+x));k=.5*y/h;}
    else{k=std::copysign(std::sqrt(.5*(1-x)),y);h=k!=0?.5*y/k:0;}
}
// Return (h,k) directly: x=h^2-k^2,y=2hk. A line of unit normal
// (cos beta,sin beta) gives cos(2theta-beta)=offset, so compose the
// two half-angle rotors rather than forming a tangent root and reprojecting it.
inline int intersect(double p,double q,double r,double s,double plus,double minus,std::array<V2,4>& hk){
    auto fallback=[&](){std::array<V2,4> xy;int n=sparse_householder::intersect(p,q,r,s,plus,minus,xy),nn=0;
        for(int i=0;i<n;++i){double h,k;circle_to_rotor(xy[i][0],xy[i][1],h,k);if(h>0)hk[nn++]=V2(h,k);}return nn;};
    if(p==0||q==0||r==0)return fallback();
    double t=sparse_householder::pencil_root(p,q,r,s),a=1-t,b=t;if(!(a>0&&b>0))return fallback();
    double sa=std::sqrt(a),sb=std::sqrt(b),qa=q*sb,rb=r*sa,offsets[2];
    if(qa*rb>0){offsets[0]=-2*sa*sb*(s+p*t)/(qa+rb);offsets[1]=.5*(q/sa+r/sb)/(-p);}
    else{offsets[1]=-2*sa*sb*(s+p*t)/(qa-rb);offsets[0]=.5*(q/sa-r/sb)/(-p);}
    // Preserve the existing cancellation-safe endpoint lift at tangent lines.
    for(double offset:offsets)if(std::abs(1-std::abs(offset))<64*std::numeric_limits<double>::epsilon())return fallback();
    double ct=std::sqrt(.5*(1+sa)),st=.5*sb/ct;int n=0;
    auto add=[&](double h,double k){if(h<0)h=-h,k=-k;if(!(h>0))return;
        for(int i=0;i<n;++i)if((hk[i]-V2(h,k)).squaredNorm()<1e-24)return;if(n<4)hk[n++]=V2(h,k);};
    for(int j=0;j<2;++j){double off=offsets[j];if(!(std::abs(off)<1))continue;
        double ca=std::sqrt(.5*(1+off)),ss=std::sqrt(.5*(1-off)),sn=j==0?st:-st;
        add(ct*ca-sn*ss,sn*ca+ct*ss);add(ct*ca+sn*ss,sn*ca-ct*ss);
    }return n;
}
}
