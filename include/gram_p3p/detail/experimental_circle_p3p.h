// Generated benchmark derivative; production solver unchanged.
#pragma once
#include "fast_real_quartic.h"
#include "projective_circle_root.h"
#include <limits>
#include "matrix_pose.h"
namespace experimental_circle {
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
inline void quaternion(const double R[3][3],double* q){
    double trace=R[0][0]+R[1][1]+R[2][2];
    if(trace>0){q[0]=.5*std::sqrt(1+trace);double inv=.25/q[0];q[1]=(R[2][1]-R[1][2])*inv;q[2]=(R[0][2]-R[2][0])*inv;q[3]=(R[1][0]-R[0][1])*inv;}
    else{int i=R[1][1]>R[0][0]?1:0;if(R[2][2]>R[i][i])i=2;int j=(i+1)%3,k=(i+2)%3;
        q[i+1]=.5*std::sqrt(1+R[i][i]-R[j][j]-R[k][k]);double inv=.25/q[i+1];
        q[0]=(R[k][j]-R[j][k])*inv;q[j+1]=(R[j][i]+R[i][j])*inv;q[k+1]=(R[k][i]+R[i][k])*inv;
    }
}
inline bool bounded_pair(const double* c,double epsilon,double* roots){
    if(!(epsilon<1&&c[4]>0&&c[0]<0))return false;
    double bound=1/std::sqrt((1-epsilon)*(1+epsilon));
    double vertex=std::clamp(-c[3]/(4*c[4]),-bound,bound),curvature=(12*c[4]*vertex+6*c[3])*vertex+2*c[2];
    // Convexity on the PROVED physical interval, rather than subtraction of
    // two nearly equal quartic-discriminant invariants, isolates both roots.
    if(!(curvature>0))return false;
    auto eval=[&](double t){return (((c[4]*t+c[3])*t+c[2])*t+c[1])*t+c[0];};
    if(!(eval(bound)>0&&eval(-bound)>0))return false;
    double disc=std::sqrt(c[2]*c[2]-4*c[4]*c[0]),zz=c[2]>=0?-2*c[0]/(c[2]+disc):(-c[2]+disc)/(2*c[4]);
    for(int j=0;j<2;++j){double sign=j==0?1:-1,lo=0,hi=bound,t=std::sqrt(zz);if(!(t<bound))t=bound*.5;
        for(int it=0;it<40;++it){double x=sign*t,f=eval(x),df=sign*(((4*c[4]*x+3*c[3])*x+2*c[2])*x+c[1]),ddf=(12*c[4]*x+6*c[3])*x+2*c[2];
            if(f>0)hi=t;else lo=t;double ax=t,mag=(((std::abs(c[4])*ax+std::abs(c[3]))*ax+std::abs(c[2]))*ax+std::abs(c[1]))*ax+std::abs(c[0]);
            if(std::abs(f)<4e-16*mag)break;double den=2*df*df-f*ddf,next=t-2*f*df/den;
            if(!(next>lo&&next<hi))next=(lo+hi)*.5;if(next==t)break;t=next;
        }roots[j]=sign*t;
    }return true;
}
template<int Mode=1, bool Count=false, int Stage=3> inline int solve(const std::vector<EV>& rays,const std::vector<EV>& points,std::vector<MatrixPose>* poses){
    poses->clear();poses->reserve(4);std::array<V,3> f={rays[0],rays[1],rays[2]},X={points[0],points[1],points[2]};
    double a01=(f[0]-f[1]).squaredNorm(),a02pair=(f[0]-f[2]).squaredNorm(),a12pair=(f[1]-f[2]).squaredNorm();
    double score01=a01*(4-a01),score02=a02pair*(4-a02pair),score12=a12pair*(4-a12pair);
    // Reorder the already computed three chord metrics with the input pair.
    double A=a01,a02=a02pair,a12=a12pair;
    if(score02>score01&&score02>score12){std::swap(f[1],f[2]);std::swap(X[1],X[2]);A=a02pair;a02=a01;}
    else if(score12>score01){std::swap(f[0],f[2]);std::swap(X[0],X[2]);A=a12pair;a12=a01;}
    V dif=f[0]-f[1],sum=f[0]+f[1];double B=4-A;if(!(A>0&&B>0))return 0;
    double sa=std::sqrt(A),sb=std::sqrt(B),sab=sa*sb,invab=1/sab,ia=sb*invab,ib=sa*invab;
    // Unit-ray identities are applied BEFORE projecting onto a difference axis.
    // The volume uses two short differences and preserves its small value.
    double U=(2-.5*(a02+a12))*ib,Vv=.5*(a12-a02)*ia;
    double Z=-2*f[0].dot((f[1]-f[0]).cross(f[2]-f[0]))*invab;V wb=X[0]-X[1],ws=X[0]-X[2],wn=wb.cross(ws);
    double side=wb.squaredNorm(),area=wn.squaredNorm();if(!(side>0&&area>0))return 0;
    double invside=1/side,height2=area*invside*invside,height=0,delta=wb.dot(ws+X[1]-X[2])*invside;
    double alpha=.5*(B*Vv*invab+delta*U),beta=-.5*(A*U*invab+delta*Vv),m2=alpha*alpha+beta*beta;
    constexpr bool scaled=false;
    double mm=0,H=1,K=0,H1=0,K1=1,eps=1,rho=scaled?Z/height:Z,mu=2*invab,nu=-.5*delta,c[5];
    if(scaled){
        mm=std::sqrt(m2);H=-beta/mm;K=alpha/mm;H1=alpha/mm;K1=beta/mm;eps=height/mm;
        double g=Vv*H+U*K,z=Vv*H1+U*K1;
        double T0=mu*H*K+nu,T1=mu*(H*K1+H1*K),T2=mu*H1*K1+nu;
        double e2=eps*eps,rr=rho*rho;
        c[0]=-g*g+rr*T0*T0;c[1]=2*eps*(-g*z+rr*T0*T1);
        c[2]=1-e2*(g*g+z*z)+rr*e2*(T1*T1+2*T0*T2);
        c[3]=2*eps*e2*(-g*z+rr*T1*T2);c[4]=e2*(1-e2*z*z)+rr*e2*e2*T2*T2;
    }else{
        // h>0 for positive depths, so t=k/h has no physical point at infinity.
        // The equal odd coefficients are an exact consequence of h^2+k^2=1.
        double al=alpha,be=beta,rr=rho*rho;
        c[0]=al*al+rr*nu*nu-height2*Vv*Vv;
        c[4]=be*be+rr*nu*nu-height2*U*U;
    }
    double roots[4];int nr=0;std::array<Eigen::Vector2d,4> rotor;bool direct_rotor=false;
    if(U==0&&Vv==0){
        // The squared eliminant has a repeated factor. Use T=0 unsquared.
        double y=delta*sab*.5;if(std::abs(y)>1)return 0;
        double x=std::sqrt((1-y)*(1+y));
        if(x>-1)roots[nr++]=y/(1+x);
        if(x<1)roots[nr++]=y/(1-x);
    }
    else if(scaled)nr=bounded_pair(c,eps,roots)?2:fast_real_quartic::solve(c,roots);
    else{
        double al=alpha,be=beta,rr=rho*rho;
        double p=-.25*rr*mu*mu,q=.5*(al*al-be*be+height2*(U*U-Vv*Vv));
        double r=al*be+rr*nu*mu-height2*U*Vv,s=.5*(al*al+be*be-height2*(U*U+Vv*Vv))+rr*nu*nu-p;
        if constexpr(Stage==0){poses->emplace_back(Eigen::Matrix3d::Identity(),EV(p+q+r+s,c[0],c[4]));return 1;}
        nr=projective_circle_root::intersect<Mode,Count>(p,q,r,s,c[0],c[4],rotor);direct_rotor=true;

    }
    if constexpr(Stage==1){if(nr)poses->emplace_back(Eigen::Matrix3d::Identity(),EV(nr,direct_rotor?rotor[0][0]:roots[0],0));return nr;}
    // A local rotation determines both baseline depths, hence translation too.
    // Keep only these compact hypotheses while solving and lifting the roots.
    double local_q[4][4],local_depth[4];int nsol=0;
    for(int j=0;j<nr;++j){
        double t=direct_rotor?0:roots[j],v=eps*t,I=1,h,k;
        if(direct_rotor){h=rotor[j][0];k=rotor[j][1];}
        else{if(!std::isfinite(t))continue;I=std::copysign(1/std::sqrt(1+v*v),H+v*H1);h=I*(H+v*H1);k=I*(K+v*K1);}
        if(!(h*ia>std::abs(k)*ib))continue;
        double F=k*U+h*Vv,T=mu*h*k+nu;
        if(F*T<0&&std::abs(F)>32*std::numeric_limits<double>::epsilon()*(std::abs(k*U)+std::abs(h*Vv))&&std::abs(T)>32*std::numeric_limits<double>::epsilon()*(std::abs(mu*h*k)+std::abs(nu)))continue;
        double raw_co=scaled?I*t:alpha*h+beta*k,raw_si=-rho*T;
        double norm=raw_co*raw_co+raw_si*raw_si,target=(scaled?1:height2)*F*F;
        bool fit=std::abs(norm-target)<=64*std::numeric_limits<double>::epsilon()*(norm+target);
        if(!fit){
            // P(theta)=C^2+Z^2*T^2-eta^2*F^2. Use its tangent derivative
            // without entering a monomial quartic or tangent parameter chart.
            double cd=-alpha*k+beta*h,td=mu*(h*h-k*k),fd=h*U-k*Vv;
            double half_der=raw_co*cd+Z*Z*T*td-height2*F*fd;
            double u=-(norm-target)/(4*half_der);
            if(std::isfinite(u)&&std::abs(u)<1){
                // Cayley retraction preserves h^2+k^2=1 in exact arithmetic.
                double uu=u*u,inv=1/(1+uu),hh=((1-uu)*h-2*u*k)*inv;
                k=(2*u*h+(1-uu)*k)*inv;h=hh;
                F=k*U+h*Vv;T=mu*h*k+nu;
            }else{
            if(direct_rotor)t=k/h;
        // Preserve the small incidence factors instead of a dense Horner sum.
        double al=scaled?alpha/height:alpha,be=scaled?beta/height:beta;
        double gg=Vv*H+U*K,zz=Vv*H1+U*K1;
        double T0=mu*H*K+nu,T1=mu*(H*K1+H1*K),T2=mu*H1*K1+nu;
        for(int it=0;it<4&&!(U==0&&Vv==0);++it){
            double l=scaled?t:al+be*t,ld=scaled?1:be;
            double fn=scaled?gg+eps*zz*t:Vv+U*t,fd=scaled?eps*zz:U;
            double tn=scaled?T0+eps*T1*t+eps*eps*T2*t*t:nu*(1+t*t)+mu*t;
            double td=scaled?eps*T1+2*eps*eps*T2*t:2*nu*t+mu;
            double e=scaled?eps*eps:1,D=1+e*t*t,ll=l*l,weight=scaled?1:height2,ff=weight*fn*fn,rt=rho*tn;
            double val=(ll-ff)*D+rt*rt,mag=(ll+ff)*D+rt*rt;
            if(std::abs(val)<=4*std::numeric_limits<double>::epsilon()*mag)break;
            double der=2*(l*ld-weight*fn*fd)*D+2*e*t*(ll-ff)+2*rho*rho*tn*td;
            if(der==0)break;double next=t-val/der;if(next==t)break;t=next;
        }

            v=eps*t;I=std::copysign(1/std::sqrt(1+v*v),H+v*H1);h=I*(H+v*H1);k=I*(K+v*K1);
            F=k*U+h*Vv;T=mu*h*k+nu;
            }
        }
        double d0=h*ia+k*ib,d1=h*ia-k*ib;if(!(d0>0&&d1>0))continue;
        auto emit=[&](double co,double si,bool check=true)->bool{
            double r2=co*co+si*si;if(!(r2>0))return false;
            // Store eta*cos(phi), eta*sin(phi), not the unit rotation pair.
            // On the exact incidence equation the root gives this length
            // already. The altitude cancels with the dual world frame below.
            if(check){if(height==0)height=std::sqrt(height2);
            double inv=height/std::sqrt(r2);co*=inv;si*=inv;
            double ye=.5*h*ia/ib-.5*delta*k-co*h;
            double yd=.5*k*ib/ia-.5*delta*h+co*k;
            double yn=-si,qd=ye*U+yd*Vv+yn*Z;
            double re=ye-qd*U,rd=yd-qd*Vv,rn=yn-qd*Z;
            if(!(qd>0)||re*re+rd*rd+rn*rn>1e-14*(ye*ye+yd*yd+yn*yn))return false;}
            double local[4]={h,k,co,si};bool duplicate=false;
            for(int n=0;n<nsol;++n){double dist=0;for(int i=0;i<4;++i){double diff=local_q[n][i]-local[i];dist+=diff*diff*(i<2?height2:1);}if(dist<height2*1e-20)duplicate=true;}
            if(!duplicate&&nsol<4){for(int i=0;i<4;++i)local_q[nsol][i]=local[i];local_depth[nsol++]=d0;}
            return true;
        };
        auto emit_normal=[&](){double co=scaled?I*t:(alpha*h+beta*k),si=-rho*T;
            if(!(F*T>0))return false;
            double norm=co*co+si*si,target=(scaled?1:height2)*F*F;
            // The unexpanded original incidence norm is the acceptance equation.
            // It replaces a redundant third-point reprojection calculation.
            bool root_fit=std::abs(norm-target)<=64*std::numeric_limits<double>::epsilon()*(norm+target);
            double invF=1/F;return emit(co*invF,si*invF,!root_fit);
        };
        if(F*F>64*std::numeric_limits<double>::epsilon()*(U*U+Vv*Vv)&&emit_normal())continue;
        // Choose a single lift chart per root. Two sphere branches are retained
        // at a projection pole, while one pose is not counted twice across charts.
        double dl=d0,side02=ws.squaredNorm()/side;
        double discr=side02-dl*dl*a02*(1-.25*a02);
        if(discr<0){emit_normal();continue;}
        double center=dl*(1-.5*a02),offset=std::sqrt(discr),qs[2];
        qs[0]=center+std::copysign(offset,center);
        qs[1]=qs[0]!=0?(dl*dl-side02)/qs[0]:center-offset;
        bool sphere_fit=false;
        for(int i=0;i<2;++i)if(qs[i]>0){
            // A backward-error check of the baseline-plane equation.
            // Its scale is arithmetic in the dimensionless geometric equation,
            // rather than an image residual multiplied by a very large depth.
            double error_budget=32*std::numeric_limits<double>::epsilon()*(std::abs(qs[i]*F)+std::abs(mu*h*k)+std::abs(nu)+1+height);
            if(std::abs(qs[i]*F-T)>error_budget)continue;
            double co=.5*(h*h*ia/ib-k*k*ib/ia)-qs[i]*(h*U-k*Vv),si=-qs[i]*Z;sphere_fit|=emit(co,si);
        }
        if(!sphere_fit)emit_normal();
    }
    if constexpr(Stage==2){if(nsol)poses->emplace_back(Eigen::Matrix3d::Identity(),EV(nsol,local_depth[0],local_q[0][2]));return nsol;}
    if(nsol==0)return 0;
    double L=std::sqrt(side);V ce=sum*ib,cn=sum.cross(dif)*invab,cd=cn.cross(ce);
    V w0=wb*(L*invside),w2=wn*(side/area),w1=w2.cross(w0);
    V wx(w0.dot(X[0]),w1.dot(X[0]),w2.dot(X[0]));
    for(int n=0;n<nsol;++n){
        double h=local_q[n][0],k=local_q[n][1],co=local_q[n][2],si=local_q[n][3];
        V v0=h*cd+k*ce,perp=h*ce-k*cd,v1=co*perp+si*cn,v2=si*perp-co*cn;
        double R[3][3];for(int i=0;i<3;++i)for(int j=0;j<3;++j)R[i][j]=v0[i]*w0[j]+v1[i]*w1[j]+v2[i]*w2[j];
        V tr=(L*local_depth[n])*f[0]-(wx[0]*v0+wx[1]*v1+wx[2]*v2);
        Eigen::Matrix3d rotation;for(int i=0;i<3;++i)for(int j=0;j<3;++j)rotation(i,j)=R[i][j];
        poses->emplace_back(rotation,tr.eigen());
    }
    return nsol;
}
}
