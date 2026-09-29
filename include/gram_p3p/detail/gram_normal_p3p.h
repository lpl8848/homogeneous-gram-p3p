#pragma once
#include "rational_normal_solver.h"

// Homogeneous Gram/normal-coordinate P3P. Nonzero bearing vectors may have
// arbitrary positive magnitudes. The generic path has no ray normalization.
namespace gram_normal_p3p {
using EV=Eigen::Vector3d;using V=gram_geometry::V;
inline volatile double profile_sink=0;
inline long long polish_counts[6]{},physical_rejects=0,fallback_calls=0;
template<int RootMode=2,int Stage=3,bool Count=false,bool FastGram=false,bool Dedup=true,int Pivot=0>
inline int solve(const std::vector<EV>& rays,const std::vector<EV>& points,std::vector<MatrixPose>* poses){
    poses->clear();poses->reserve(4);std::array<V,3> f={rays[0],rays[1],rays[2]},X={points[0],points[1],points[2]};
    double n0=f[0].squaredNorm(),n1=f[1].squaredNorm(),n2=f[2].squaredNorm();if(!(n0>0&&n1>0&&n2>0))return 0;
    double c=f[0].dot(f[1]),c02=f[0].dot(f[2]),c12=f[1].dot(f[2]);
    V N;double G;
    if constexpr(Pivot==0){
        V normal01,normal02,normal12;double G01,G02,G12;
        if constexpr(FastGram){G01=n0*n1-c*c;G02=n0*n2-c02*c02;G12=n1*n2-c12*c12;}
        else{normal01=f[0].cross(f[1]-f[0]);normal02=f[0].cross(f[2]-f[0]);normal12=f[1].cross(f[2]-f[1]);G01=normal01.squaredNorm();G02=normal02.squaredNorm();G12=normal12.squaredNorm();}
        double score01=G01*n2,score02=G02*n1,score12=G12*n0;G=G01;if constexpr(!FastGram)N=normal01;
        if(score02>score01&&score02>score12){std::swap(f[1],f[2]);std::swap(X[1],X[2]);std::swap(n1,n2);std::swap(c,c02);if constexpr(!FastGram)N=normal02;G=G02;}
        else if(score12>score01){std::swap(f[0],f[2]);std::swap(X[0],X[2]);std::swap(n0,n2);std::swap(c,c12);if constexpr(!FastGram)N=-normal12;G=G12;}
        if constexpr(FastGram){N=f[0].cross(f[1]-f[0]);G=N.squaredNorm();}
    }else{N=f[0].cross(f[1]-f[0]);G=N.squaredNorm();}
    if(!(G>0))return rational_normal_solver::solve(rays,points,poses);
    double tau;if constexpr(FastGram)tau=N.dot(f[2]);else tau=f[0].dot((f[1]-f[0]).cross(f[2]-f[0]));double tau2=tau*tau;
    V wb=X[0]-X[1],ws=X[0]-X[2],wn=wb.cross(ws);double side=wb.squaredNorm(),area=wn.squaredNorm();if(!(side>0&&area>0))return 0;
    double invside=1/side,eta2=area*invside*invside,gamma=wb.dot(ws)*invside;
    double alpha=n1*c02-c*c12,beta=n0*c12-c*c02,l0=gamma*alpha,l1=(1-gamma)*beta,D0=G*gamma*(1-gamma),ee=eta2*G;
    double pole0=l0*c02,pole1=l1*c12;
    if(std::abs(pole0-pole1)<=64*std::numeric_limits<double>::epsilon()*(std::abs(pole0)+std::abs(pole1))){
        double cp=std::abs(c12)>std::abs(c02)?-l0/c12:-l1/c02,dp=(c-cp)*cp+D0;
        if(std::isfinite(cp)&&std::abs(dp)<=128*std::numeric_limits<double>::epsilon()*(std::abs(c*cp)+cp*cp+std::abs(D0)))return rational_normal_solver::solve(rays,points,poses);
    }
    V e0=l0*f[0]-l1*f[1],e1=c12*f[0]-c02*f[1];
    double E0=e0.squaredNorm(),E1=2*alpha*beta,E2=e1.squaredNorm();
    double poly[5]={tau2*D0*D0-ee*E0,-ee*E1+2*tau2*c*D0,E0-ee*E2+tau2*(c*c-2*D0),E1-2*tau2*c,G*n2};
    if constexpr(Stage==0){profile_sink=poly[0]+poly[1]+poly[2]+poly[3]+poly[4];return 0;}
    double invleading=1/poly[4],roots[4];int nr=monic_quartic::solve<RootMode,Count,true>(poly,roots,invleading),n=0;bool failed=false;
    if constexpr(Stage==1){double ss=0;for(int i=0;i<nr;++i)ss+=roots[i];profile_sink=ss;return nr;}
    double local[4][6],depths[4][3];
    for(int i=0;i<nr;++i){double chi=roots[i];if(!std::isfinite(chi))continue;double a0=0,a1=0,D=0,E=0;V ev;bool fit=false,physical=true;
        for(int it=0;it<5;++it){a0=l0+chi*c12;a1=l1+chi*c02;D=(c-chi)*chi+D0;
            if(it==0&&!(D>0?(a0>0&&a1>0):(a0<0&&a1<0))){physical=false;if constexpr(Count)++physical_rejects;break;}
            ev=e0+chi*e1;E=ev.squaredNorm();double a=chi*chi-ee,tD=tau*D,val=a*E+tD*tD,mag=(chi*chi+ee)*E+tD*tD;
            if(std::abs(val)<=64*std::numeric_limits<double>::epsilon()*mag){fit=true;if constexpr(Count)++polish_counts[it];break;}
            double Ep=2*ev.dot(e1),Dp=c-2*chi,der=2*chi*E+a*Ep+2*tau2*D*Dp;
            if(der==0)break;double next=chi-val/der;if(next==chi)break;chi=next;
        }
        if(!physical)continue;
        if(!fit||!(E>0)){failed=true;continue;}
        double inv=std::copysign(1/std::sqrt(E),D),d0=inv*a0,d1=inv*a1,d2=inv*D;
        if(!(d0>0&&d1>0&&d2>0))continue;
        bool duplicate=false;if constexpr(Dedup)for(int j=0;j<n;++j){double x=d0-depths[j][0],y=d1-depths[j][1],z=d2-depths[j][2];if(x*x+y*y+z*z<1e-20*(d0*d0+d1*d1+d2*d2))duplicate=true;}
        if(!duplicate&&n<4){local[n][0]=chi;local[n][1]=D;local[n][2]=inv;local[n][3]=ev.x;local[n][4]=ev.y;local[n][5]=ev.z;depths[n][0]=d0;depths[n][1]=d1;depths[n++][2]=d2;}
    }
    if(failed){if constexpr(Count)++fallback_calls;return rational_normal_solver::solve(rays,points,poses);}
    if constexpr(Stage==2){profile_sink=n?depths[n-1][0]:0;return n;}
    if(n==0)return 0;
    double L=std::sqrt(side),iG=invleading*n2;V w0=wb*(L*invside),w2=wn*(side/area),w1=w2.cross(w0),wx;
    if constexpr(!FastGram)wx=V(w0.dot(X[0]),w1.dot(X[0]),w2.dot(X[0]));
    for(int i=0;i<n;++i){double chi=local[i][0],D=local[i][1],inv=local[i][2];
        V v0=inv*V(local[i][3],local[i][4],local[i][5]),cross=N.cross(v0);double K=chi*iG,T=tau*D*inv*iG;
        V v1=K*cross-T*N,v2=T*cross+K*N;Eigen::Matrix3d R;
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)R(a,b)=v0[a]*w0[b]+v1[a]*w1[b]+v2[a]*w2[b];
        V t;
        if constexpr(FastGram){V rx(R(0,0)*X[0][0]+R(0,1)*X[0][1]+R(0,2)*X[0][2],R(1,0)*X[0][0]+R(1,1)*X[0][1]+R(1,2)*X[0][2],R(2,0)*X[0][0]+R(2,1)*X[0][1]+R(2,2)*X[0][2]);t=(L*depths[i][0])*f[0]-rx;}
        else t=(L*depths[i][0])*f[0]-(wx[0]*v0+wx[1]*v1+wx[2]*v2);poses->emplace_back(R,t.eigen());
    }return n;
}
}
