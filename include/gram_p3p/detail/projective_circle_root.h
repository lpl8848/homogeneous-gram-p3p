#pragma once
#include "rotor_householder.h"

// Research candidate. Both homogeneous pencil coordinates are retained;
// a small endpoint distance is never reconstructed by subtracting from 1.
namespace projective_circle_root {
inline long long histogram[5][41]{};
inline long long rejected_steps[5]{};
struct Pair { double a,b; };

template<int Mode, bool Count=false>
inline Pair root(double p,double q,double r,double s) {
    double qq=q*q,rr=r*r;
    if constexpr(Mode==0) {
        double lo=0,hi=1,t=rr/(qq+rr);int steps=0;
        for(;steps<40;++steps) {
            double w=s+p*t,tt=t*(1-t);
            double f=t*qq-(1-t)*rr-4*p*tt*w;
            double df=qq+rr-4*p*((1-2*t)*w+p*tt);
            double ddf=8*p*s-8*p*p+24*p*p*t;
            if(f>0)hi=t;else lo=t;
            double mag=std::abs(t*qq)+std::abs((1-t)*rr)+std::abs(4*p*tt*w);
            if(std::abs(f)<=4e-16*mag)break;
            double ff=df*df,den=df*(ff-f*ddf)+4*p*p*f*f;
            double next=t-(den!=0?f*(ff-.5*f*ddf)/den:f/df);
            if(!(next>lo&&next<hi)){next=(lo+hi)*.5;if constexpr(Count)++rejected_steps[Mode];}
            if(next==t)break;t=next;
        }
        if constexpr(Count)++histogram[Mode][std::min(steps+1,40)];
        return {1-t,t};
    } else {
        // D(0)<0 and D(1)>0. The midpoint chooses a half interval containing
        // a root. Reflection sends that interval to [0,1/2].
        bool flip=.5*(qq-rr)-p*(s+.5*p)<0;
        double ss=flip?-s-p:s,Q=flip?rr:qq,R=flip?qq:rr;
        double a=4*p*p,b=4*p*(ss-p),c=Q+R-4*p*ss;
        double lo=0,hi=.5,t=R/(Q+R);int steps=0;
        if constexpr(Mode>=2) {
            // Endpoint osculating quadratic. Select its smaller positive root
            // by rationalization; the cubic remainder is then polished.
            double disc=c*c+4*b*R;
            if(disc>=0) {double den=c+std::sqrt(disc);if(den>0)t=2*R/den;}
        }
        if(!(t>0&&t<hi))t=.25;
        for(;steps<40;++steps) {
            double f=((a*t+b)*t+c)*t-R;
            double df=(3*a*t+2*b)*t+c,ddf=6*a*t+2*b;
            if(f>0)hi=t;else lo=t;
            double mag=((a*t+std::abs(b))*t+std::abs(c))*t+R;
            if(std::abs(f)<=4e-16*mag)break;
            double ff=df*df,den=df*(ff-f*ddf)+a*f*f;
            double next=t-(den!=0?f*(ff-.5*f*ddf)/den:f/df);
            if(!(next>lo&&next<hi)){next=(lo+hi)*.5;if constexpr(Count)++rejected_steps[Mode];}
            if(next==t)break;t=next;
        }
        if constexpr(Count)++histogram[Mode][std::min(steps+1,40)];
        return flip?Pair{t,1-t}:Pair{1-t,t};
    }
}

template<int Mode,bool Count=false>
inline int intersect(double p,double q,double r,double s,double plus,double minus,std::array<Eigen::Vector2d,4>& hk) {
    if(p==0||q==0||r==0)return rotor_householder::intersect(p,q,r,s,plus,minus,hk);
    auto ab=root<Mode,Count>(p,q,r,s);double a=ab.a,b=ab.b;
    if(!(a>0&&b>0))return rotor_householder::intersect(p,q,r,s,plus,minus,hk);
    double sa=std::sqrt(a),sb=std::sqrt(b),qa=q*sb,rb=r*sa,offsets[2];
    // Express s+p*b in the closer endpoint chart, too.
    double w=b<=.5?s+p*b:(s+p)-p*a;
    unsigned mask=3;
    if constexpr(Mode>=3){
        int big=qa*rb>0?1:0;
        double N=big==1?qa+rb:qa-rb,K=-2*p*sa*sb,C=-2*sa*sb*w;
        mask=0;
        // Homogeneous line offsets. Test circle intersection before division;
        // this also avoids constructing an enormous nonintersecting offset.
        if(std::abs(N)<std::abs(K)){offsets[big]=N/K;mask|=1u<<big;}
        if(std::abs(C)<std::abs(N)){offsets[1-big]=C/N;mask|=1u<<(1-big);}
        double eps=64*std::numeric_limits<double>::epsilon();
        if(std::abs(std::abs(N)-std::abs(K))<=eps*std::max(std::abs(N),std::abs(K))||std::abs(std::abs(C)-std::abs(N))<=eps*std::max(std::abs(C),std::abs(N)))return rotor_householder::intersect(p,q,r,s,plus,minus,hk);
        if(mask==0)return 0;
    }else{
        if(qa*rb>0){offsets[0]=-2*sa*sb*w/(qa+rb);offsets[1]=.5*(q/sa+r/sb)/(-p);}
        else{offsets[1]=-2*sa*sb*w/(qa-rb);offsets[0]=.5*(q/sa-r/sb)/(-p);}
        for(double off:offsets)if(std::abs(1-std::abs(off))<64*std::numeric_limits<double>::epsilon())return rotor_householder::intersect(p,q,r,s,plus,minus,hk);
    }
    double ct=std::sqrt(.5*(1+sa)),st=.5*sb/ct;int n=0;
    auto add=[&](double h,double k){if(h<0)h=-h,k=-k;if(!(h>0))return;
        for(int i=0;i<n;++i)if((hk[i]-Eigen::Vector2d(h,k)).squaredNorm()<1e-24)return;if(n<4)hk[n++]=Eigen::Vector2d(h,k);};
    for(int j=0;j<2;++j){if(!(mask&(1u<<j)))continue;double off=offsets[j];if(!(std::abs(off)<1))continue;
        double ca=std::sqrt(.5*(1+off)),ss=std::sqrt(.5*(1-off)),sn=j==0?st:-st;
        add(ct*ca-sn*ss,sn*ca+ct*ss);add(ct*ca+sn*ss,sn*ca-ct*ss);
    }return n;
}
}
