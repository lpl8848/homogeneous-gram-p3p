#pragma once
#include "fast_real_quartic.h"
#include "fast_cuberoot.h"
namespace monic_quartic {
inline long long cubic_counts[2]{};
// Classical real arithmetic only. The two policies select different real
// cubic branches; no camera geometry approximations or angular restrictions.
template<int AlgebraicTrig=0,bool Count=false>
inline double cubic(double b,double c,double d,bool positive){
    double shift=b/3,P=c-b*shift,Q=(2*shift*shift-c)*shift+d,T=P/3,H=Q/2,D=H*H+T*T*T;
    double x;
    if constexpr(Count)++cubic_counts[D>=0?0:1];
    if(D>=0){double v=-H-std::copysign(std::sqrt(D),H),u;
        if constexpr(AlgebraicTrig>=3)u=fast_radical::cbrt(v);else u=std::cbrt(v);
        x=u-(u!=0?T/u:0);}
    else{double radius=std::sqrt(-T),arg=std::clamp(-H/(radius*radius*radius),-1.,1.);
        if constexpr(AlgebraicTrig==2||AlgebraicTrig==4){
            // Shift to the positive stationary point: z^2*(z+3*radius)
            // =2*radius^3-Q. With t=sqrt(2*(1+arg)/27), z=3*radius*t*w,
            // the three-real-root branch becomes w^2*(1+t*w)=1.
            // t is in [0,sqrt(4/27)] and w in [sqrt(3)/2,1], uniformly
            // well conditioned. A degree-five seed plus ONE H4 correction
            // replaces acos/cos; this is not an approximate pose output.
            double t=std::sqrt(std::max(0.,(2./27.)*(1+arg)));
            double w=((((-.5876944488286939*t+1.0018166431130315)*t-.8728298040856343)*t+.6149699613296116)*t-.4996963574003894)*t+.9999984699611322;
            double f=w*w*(1+t*w)-1,df=2*w+3*t*w*w,ddf=2+6*t*w,ff=df*df;
            w-=f*(ff-.5*f*ddf)/(df*(ff-f*ddf)+t*f*f);
            x=radius*(1+3*t*w);
        }else if constexpr(AlgebraicTrig==1){
            // |arg|=4*y^3-3*y has ONE root in [sqrt(3)/2,1].
            // Unlike the direct largest-root chart, its derivative is >=6,
            // including repeated-root limits of the original cubic.
            double target=std::abs(arg),y=.86602540378443864676+.13397459621556135324*target;
            for(int i=0;i<2;++i){double f=(4*y*y-3)*y-target,df=12*y*y-3,ddf=24*y,ff=df*df;
                y-=f*(ff-.5*f*ddf)/(df*(ff-f*ddf)+4*f*f);}
            // 1-y=(1-target)/(2*y+1)^2 preserves a repeated-root distance
            // smaller than one ulp of y. Never subtract the rounded y from 1.
            if(positive)x=arg>=0?2*radius*y:radius*(y+std::sqrt(3*(1-target)*(1+y))/(2*y+1));
            else x=std::copysign(2*radius*y,arg);
        }else{
            if(positive)x=2*radius*std::cos(std::acos(arg)/3);
            else x=std::copysign(2*radius*std::cos(std::acos(std::abs(arg))/3),arg);
        }
    }return x-shift;
}
template<int Mode,bool Count=false,bool ProvidedInverse=false>
inline int solve(const double* cs,double* roots,double supplied_inverse=0){
    if constexpr(Mode==0)return fast_real_quartic::solve(cs,roots);
    double inv;if constexpr(ProvidedInverse)inv=supplied_inverse;else inv=1/cs[4];double a=cs[3]*inv,b=cs[2]*inv,c=cs[1]*inv,d=cs[0]*inv;
    if constexpr(Mode==1){
        double y=cubic(-b,a*c-4*d,4*b*d-a*a*d-c*c,false);
        double D=y*y-4*d;if(!(D>0))return fast_real_quartic::solve(cs,roots);
        double sq=std::sqrt(D),q1=.5*(y+std::copysign(sq,y)),q2=q1!=0?d/q1:0;
        double div=q1-q2;if(div==0)return fast_real_quartic::solve(cs,roots);
        double iv=1/div,p1=(a*q1-c)*iv,p2=(c-a*q2)*iv;
        int n=fast_real_quartic::quadratic(p1,q1,roots);return n+fast_real_quartic::quadratic(p2,q2,roots+n);
    }else{
        double shift=.25*a,ss=shift*shift,p=b-6*ss,q=c+(8*ss-2*b)*shift,r=d-c*shift+(b-3*ss)*ss;
        if(q==0)return fast_real_quartic::solve(cs,roots);
        double m=cubic<Mode==6?4:Mode==5?3:Mode==4?2:Mode==3?1:0,Count>(2*p,p*p-4*r,-q*q,true);
        if(!(m>0&&std::isfinite(m)))return fast_real_quartic::solve(cs,roots);
        double sm=std::sqrt(m),qm=q/sm,u=.5*(p+m-qm),v=.5*(p+m+qm);
        int n=fast_real_quartic::quadratic(sm,u,roots);n+=fast_real_quartic::quadratic(-sm,v,roots+n);
        if(n==0&&r<0){
            // A depressed quartic with r<0 MUST have at least two real roots.
            // Empty quadratic factors therefore signal numerical failure.
            // Vieta: m*((m+p)^2-4*r)=q^2. For r<0 its denominator is a
            // sum of positive terms; recover a tiny m without undo-shift loss.
            double mp=m+p;m=q*q/(mp*mp-4*r);
            if(m>0){sm=std::sqrt(m);qm=q/sm;u=.5*(p+m-qm);v=.5*(p+m+qm);
                n=fast_real_quartic::quadratic(sm,u,roots);n+=fast_real_quartic::quadratic(-sm,v,roots+n);}
            if(n==0)return fast_real_quartic::solve(cs,roots);
        }
        for(int i=0;i<n;++i)roots[i]-=shift;return n;
    }
}
}
