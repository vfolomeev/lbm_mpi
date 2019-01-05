#pragma once
#ifndef PHYSICS_H
#define PHYSICS_H
#include "math_api.h"

class Physics
{
public:
     double rho0;
     double G;
     double alpha;
     double omega;
    Physics(){
        std::cout<<"created Physics\n";
        rho0=1.;
        G=2.25e-4/120;
        alpha=0.01;
        omega=1.0/(3.*alpha+0.5);

    }
    Physics(Physics* p){
        std::cout<<"copied Physics\n";
        rho0=p->rho0;
        G=p->G;
        alpha=p->alpha;
        omega=p->omega;

    }

    Vector3d updatef(double* f,Vector3d* c,double* w){
        //compute rho and u
            double sum=0.;
            Vector3d sumu=Vector3d();

            for(int k=0;k<27;k++){
                sum=sum+f[k];

                sumu=sumu+(f[k])*c[k];


            }
            double rhon=sum;

            Vector3d u=(sumu+Vector3d(G,0,0)*(1./omega))*(1/rhon); //Application of forcing;

        //collision
            double cu,feq;
            for(int k=0;k<27;k++){
                cu=c[k]*u;
                feq=w[k]*rhon*(1+3*cu+4.5*cu*cu-1.5*(u*u));
                f[k]=omega*feq+(1.-omega)*f[k];
            }
            return u;
    }

};
std::ostream& operator<<(std::ostream& cout_,const Physics& p){
    cout_<<"rho0= "<<p.rho0<<"\n";
    cout_<<"G= "<<p.G<<"\n";
    cout_<<"alpha= "<<p.alpha<<"\n";
    cout_<<"omega= "<<p.omega<<"\n";

    return cout_;
}

#endif // PHYSICS_H
