#pragma once
#include "fast_real_quartic.h"
namespace sparse_householder {
using V=Eigen::Vector3d;using M=Eigen::Matrix3d;
inline int quadratic(double a,double b,double c,double* out){
    if(a==0){if(b==0)return 0;out[0]=-c/b;return 1;}
    double D=b*b-4*a*c;if(D<0)return 0;
    double q=-.5*(b+std::copysign(std::sqrt(D),b));
    if(q==0){out[0]=0;return 1;}
    out[0]=q/a;out[1]=c/q;return D==0?1:2;
}
// Exact geometry: a unit circle intersects p*x*x+q*x+r*y+s=0, p<=0.
// Select an indefinite singular pencil member with a bracketed cubic in [0,1].
inline double pencil_root(double p,double q,double r,double s){
    double qq=q*q,rr=r*r,lo=0,hi=1,t=rr/(qq+rr);
    for(int i=0;i<40;++i){
        double w=s+p*t,tt=t*(1-t);
        double f=t*qq-(1-t)*rr-4*p*tt*w;
        double df=qq+rr-4*p*((1-2*t)*w+p*tt);
        double ddf=8*p*s-8*p*p+24*p*p*t;
        if(f>0)hi=t;else lo=t;
        double error=std::abs(t*qq)+std::abs((1-t)*rr)+std::abs(4*p*tt*w);
        if(std::abs(f)<=4e-16*error)return t;
        // D is cubic, D'''=24*p*p. Householder order four has one
        // division per step, like Halley, but uses this known third derivative.
        double ff=df*df,divisor=df*(ff-f*ddf)+4*p*p*f*f;
        double next=t-(divisor!=0?f*(ff-.5*f*ddf)/divisor:f/df);
        if(!(next>lo&&next<hi))next=(lo+hi)*.5;
        if(next==t)return t;
        t=next;
    }
    return t;
}
inline int intersect(double p,double q,double r,double s,double at_plus_one,double at_minus_one,std::array<Eigen::Vector2d,4>& xy){
    int n=0;
    auto add=[&](double x,double y){
        if(!std::isfinite(x+y))return;
        for(int j=0;j<n;++j)if((xy[j]-Eigen::Vector2d(x,y)).squaredNorm()<1e-24)return;
        if(n<4)xy[n++]=Eigen::Vector2d(x,y);
    };
    auto line=[&](double nx,double ny,double offset){
        double norm2=nx*nx+ny*ny,disc=norm2-offset*offset;
        if(!(norm2>0)||disc<0)return;
        double inv=1/norm2,step=std::sqrt(disc)*inv,c=offset*inv;
        add(c*nx-step*ny,c*ny+step*nx);
        if(disc>0)add(c*nx+step*ny,c*ny-step*nx);
    };
    if(p==0){line(q,r,-s);return n;}
    if(r==0){double zs[2];int nr=quadratic(p,-2*p-q,at_plus_one,zs);
        for(int j=0;j<nr;++j)if(zs[j]>=0&&zs[j]<=2){double x=1-zs[j],y=std::sqrt(zs[j]*(2-zs[j]));add(x,y);if(y>0)add(x,-y);}return n;}
    if(q==0){double ys[2];int nr=quadratic(-p,r,s+p,ys);
        for(int j=0;j<nr;++j)if(std::abs(ys[j])<=1){double x=std::sqrt((1-ys[j])*(1+ys[j]));add(x,ys[j]);if(x>0)add(-x,ys[j]);}return n;}
    double t=pencil_root(p,q,r,s),a=1-t,b=t;
    if(!(a>0&&b>0))return n;
    double sa=std::sqrt(a),sb=std::sqrt(b),qa=q*sb,rb=r*sa;
    // Rationalize the cancellation-prone offset using det(C)=0.
    double offsets[2];
    if(qa*rb>0){offsets[0]=-2*sa*sb*(s+p*t)/(qa+rb);offsets[1]=.5*(q/sa+r/sb)/(-p);}
    else{offsets[1]=-2*sa*sb*(s+p*t)/(qa-rb);offsets[0]=.5*(q/sa-r/sb)/(-p);}
    double leading[2]={offsets[0]+sa,offsets[1]+sa},constant[2]={offsets[0]-sa,offsets[1]-sa};
    // At circle endpoints the product of the factors is known directly from
    // the metric equations. Preserve it instead of subtracting near equals.
    int ic=std::abs(constant[0])>std::abs(constant[1])?1:0;
    int il=std::abs(leading[0])>std::abs(leading[1])?1:0;
    if(constant[1-ic]!=0)constant[ic]=at_plus_one/(p*constant[1-ic]);
    if(leading[1-il]!=0)leading[il]=at_minus_one/(p*leading[1-il]);
    for(int j=0;j<2;++j){double roots[2];int nn=quadratic(leading[j],(j==0?-2:2)*sb,constant[j],roots);
        for(int i=0;i<nn;++i){double v=roots[i];if(std::abs(v)<=1){double vv=v*v,inv=1/(1+vv);add((1-vv)*inv,2*v*inv);}
            else{double u=1/v,uu=u*u,inv=1/(1+uu);add((uu-1)*inv,2*u*inv);}}
    }
    return n;
}
}
