#pragma once
#include "geometry.h"
#include "experimental_circle_p3p.h"
#include "monic_quartic.h"

// Rational normal-rotation elimination. Baseline angular square roots cancel
// from BOTH coefficient construction and complete matrix pose recovery.
namespace rational_normal_p3p {
using EV=Eigen::Vector3d;using V=gram_geometry::V;
template<bool Fallback=true,int RootMode=0,int LiftMode=0,int CoefficientMode=0,bool EmptyFallback=true,bool NormalizeInput=false>
inline int solve(const std::vector<EV>& rays,const std::vector<EV>& points,std::vector<MatrixPose>* poses){
    poses->clear();poses->reserve(4);std::array<V,3> f={rays[0],rays[1],rays[2]},X={points[0],points[1],points[2]};
    if constexpr(NormalizeInput)for(auto& ray:f)ray=ray.normalized();
    double a01=(f[0]-f[1]).squaredNorm(),a02=(f[0]-f[2]).squaredNorm(),a12=(f[1]-f[2]).squaredNorm();
    double score01=a01*(4-a01),score02=a02*(4-a02),score12=a12*(4-a12),A=a01;
    if(score02>score01&&score02>score12){std::swap(f[1],f[2]);std::swap(X[1],X[2]);A=a02;a02=a01;}
    else if(score12>score01){std::swap(f[0],f[2]);std::swap(X[0],X[2]);A=a12;a12=a01;}
    double B=4-A,AB=A*B;if(!(AB>0))return 0;
    double g=2-.5*(a02+a12),v=.5*(a12-a02),z=-2*f[0].dot((f[1]-f[0]).cross(f[2]-f[0]));
    V wb=X[0]-X[1],ws=X[0]-X[2],wn=wb.cross(ws);double side=wb.squaredNorm(),area=wn.squaredNorm();if(!(side>0&&area>0))return 0;
    double invside=1/side,eta2=area*invside*invside,delta=wb.dot(ws+X[1]-X[2])*invside,ee=eta2*AB;
    double kappa=.5*(B-A),D0=.25*AB*(1-delta*delta),l0=.5*(A*g+delta*B*v),l1=.5*(B*v+delta*A*g);
    // The adjugate chart collapses when both numerators vanish together at
    // a singular linear system. Detect this algebraically before solving a
    // repeated quartic root; root polishing cannot recover lost multiplicity.
    double pole0=l0*v,pole1=l1*g;
    if(std::abs(pole0+pole1)<=64*std::numeric_limits<double>::epsilon()*(std::abs(pole0)+std::abs(pole1))){
        double mp=std::abs(g)>std::abs(v)?-l0/g:l1/v;
        double dp=(kappa-mp)*mp+D0;
        if(std::isfinite(mp)&&std::abs(dp)<=128*std::numeric_limits<double>::epsilon()*(std::abs(kappa*mp)+mp*mp+std::abs(D0))){
            if constexpr(NormalizeInput){std::vector<EV> normalized={rays[0].normalized(),rays[1].normalized(),rays[2].normalized()};return experimental_circle::solve<2>(normalized,points,poses);}
            else return experimental_circle::solve<2>(rays,points,poses);
        }
    }
    double zz=z*z,c[5],invAB=0;
    if constexpr(CoefficientMode==1){
        // A*g^2+B*v^2+z^2=A*B. Apply the unit-bearing identity BEFORE
        // forming the coefficients: eliminate E1,E2 and directly form monic P.
        invAB=1/AB;double G=.25*(1-delta*delta),e0=(A*l0*l0+B*l1*l1)*invAB;
        c[3]=g*g-v*v-2*kappa;c[4]=1;
        c[2]=e0-ee+zz*(eta2-2*G)+zz*invAB*kappa*kappa;
        c[1]=-ee*c[3]+2*zz*kappa*(G-eta2);
        c[0]=zz*AB*G*G-ee*e0;
    }else{
        double E0=A*l0*l0+B*l1*l1,E1=A*A*g*g-B*B*v*v,E2=A*g*g+B*v*v;
        c[0]=zz*D0*D0-ee*E0;c[1]=-ee*E1+2*zz*kappa*D0;c[2]=E0-ee*E2+zz*(kappa*kappa-2*D0);c[3]=E1-2*zz*kappa;c[4]=AB;
    }
    double roots[4];int nr=monic_quartic::solve<RootMode>(c,roots),n=0;bool failed=false;
    double local[4][5],depths[4][3];
    for(int i=0;i<nr;++i){double m=roots[i];if(!std::isfinite(m))continue;bool fit=false,physical=true;
        double N0=0,N1=0,D=0,E=0;
        for(int j=0;j<5;++j){N0=l0+m*g;N1=l1-m*v;D=(kappa-m)*m+D0;
            if constexpr(LiftMode==1)if(j==0&&!((D>0?N0:-N0)>std::abs(N1))){physical=false;break;}
            E=A*N0*N0+B*N1*N1;
            double a=m*m-ee,zD=z*D,val=a*E+zD*zD,mag=(m*m+ee)*E+zD*zD;
            if(std::abs(val)<=64*std::numeric_limits<double>::epsilon()*mag){fit=true;break;}
            double Ep=2*(A*N0*g-B*N1*v),Dp=kappa-2*m,der=2*m*E+a*Ep+2*zz*D*Dp;
            if(der==0)break;double next=m-val/der;if(next==m){if constexpr(LiftMode==0)fit=true;break;}m=next;
        }
        if(!physical)continue;
        if constexpr(LiftMode==0){N0=l0+m*g;N1=l1-m*v;D=(kappa-m)*m+D0;E=A*N0*N0+B*N1*N1;}
        if(!((D>0?N0:-N0)>std::abs(N1)))continue;
        if(!(E>0)||!fit){failed=true;continue;}
        double inv=std::copysign(1/std::sqrt(E),D);
        if constexpr(LiftMode==0){double norm=m*m+zz*D*D/E;
            if(std::abs(norm-ee)>128*std::numeric_limits<double>::epsilon()*(norm+ee)){failed=true;continue;}}
        double d0=(N0+N1)*inv,d1=(N0-N1)*inv,d2=D*inv;
        bool duplicate=false;for(int j=0;j<n;++j){double x=d0-depths[j][0],y=d1-depths[j][1],z=d2-depths[j][2];if(x*x+y*y+z*z<1e-20*(d0*d0+d1*d1+d2*d2))duplicate=true;}
        if(!duplicate&&n<4){local[n][0]=m;local[n][1]=N0;local[n][2]=N1;local[n][3]=D;local[n][4]=inv;depths[n][0]=d0;depths[n][1]=d1;depths[n++][2]=d2;}
    }
    if constexpr(Fallback)if(failed||(EmptyFallback&&n==0)){
        if constexpr(NormalizeInput){std::vector<EV> normalized={rays[0].normalized(),rays[1].normalized(),rays[2].normalized()};return experimental_circle::solve<2>(normalized,points,poses);}
        else return experimental_circle::solve<2>(rays,points,poses);
    }
    if(n==0)return 0;
    if constexpr(CoefficientMode==0)invAB=1/AB;
    double L=std::sqrt(side),invA=B*invAB,invB=A*invAB;
    V sum=f[0]+f[1],dif=f[0]-f[1],ns=sum.cross(dif)*invAB,es=sum*invB,ds=dif*invA;
    V w0=wb*(L*invside),w2=wn*(side/area),w1=w2.cross(w0);
    V wx(w0.dot(X[0]),w1.dot(X[0]),w2.dot(X[0]));
    for(int i=0;i<n;++i){double m=local[i][0],N0=local[i][1],N1=local[i][2],D=local[i][3],inv=local[i][4];
        double h=N0*inv,k=N1*inv,S=-z*D*inv;
        V perp=h*es-k*ds,v0=h*dif+k*sum,v1=m*perp+S*ns,v2=S*perp-m*ns;
        Eigen::Matrix3d R;for(int a=0;a<3;++a)for(int b=0;b<3;++b)R(a,b)=v0[a]*w0[b]+v1[a]*w1[b]+v2[a]*w2[b];
        V t=(L*depths[i][0])*f[0]-(wx[0]*v0+wx[1]*v1+wx[2]*v2);poses->emplace_back(R,t.eigen());
    }return n;
}
}
