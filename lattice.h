#pragma once
#ifndef LATTICE_H
#define LATTICE_H
#include "math_api.h"
#include "physics.h"
#include <iostream>
#include <string>
#include "fstream"
#include "mpi.h"
using namespace std;
template<int nx,int ny,int nz>
class Lattice
{
public:
    Vector3d c[27];
    double w[27];
    double f[nx][ny][nz][27];
    Vector3d u[nx][ny][nz];
    Physics phys;
    MPI_Datatype MpiVectorType, MpiVectorType2;
    MPI_Status status;
    int mpi_grid_rank,mpi_up_z,mpi_down_z;

    Lattice(Physics physics);
    ~Lattice();

    void initSolution();
    void collide();
    void stream();
    //mpi exchange
    void mpi_Z(MPI_Comm grid_comm);
    //boundary conditions
    void bouncebackBC(Vector3d dir);
    void periodicBC(Vector3d dir);
    //results processing
    void collectStat(double* umean,double* usqrt);
    void saveResult(string fname);
    void saveState(string ufname);
    void restoreState(string ufname);
};
template<int nx,int ny,int nz>
Lattice<nx,ny,nz>::Lattice(Physics physics)
{
    phys=physics;

    //Create MPI data types
    MPI_Type_vector(9, 1, 3, MPI_DOUBLE, &MpiVectorType);
    MPI_Type_commit(&MpiVectorType);
    //The first vector structure selects data with z positive or negative c[k] directions
    //There are 9 such vectors equaly spaced with the gap=3
    MPI_Aint sizeofVector;
    MPI_Type_extent(MpiVectorType,&sizeofVector);
    //MPI_Type_vector(nx*ny, 1, nz, MpiVectorType, &MpiVectorType2);
    MPI_Type_hvector(nx*ny, 1, nz*27*sizeof(double), MpiVectorType, &MpiVectorType2);
    MPI_Type_commit(&MpiVectorType2);
    //The second vector structure uses previous one and selects data with prescribed z coordinate
    //There data are equaly spaced with the gap =nz*27, nx*ny is the area

    //Init Ck tensor
    int k;
    for(int lx=-1;lx<2;lx++){
        for(int ly=-1;ly<2;ly++){
            for(int lz=-1;lz<2;lz++){
                k=ind(lx,ly,lz);
                c[k]=Vector3d(lx,ly,lz);
                if((c[k]*c[k])==0) w[k]=8./27;
                if(c[k]*c[k]==1) w[k]=2./27;
                if(c[k]*c[k]==2) w[k]=1./54;
                if(c[k]*c[k]==3) w[k]=1./216;
                //cout<<c[k]<<" ";
            }
        }
    }
    cout<<"Created Lattice\n";
    cout<<phys;
}
template<int nx,int ny,int nz> Lattice<nx,ny,nz>::~Lattice()
{
    MPI_Type_free(&MpiVectorType);
    MPI_Type_free(&MpiVectorType2);
}


template<int nx,int ny,int nz>
void Lattice<nx,ny,nz>::initSolution(){
    int k;
#pragma omp parallel for private(k)
    for(int lx=0;lx<nx;lx++){
        for(int ly=0;ly<ny;ly++){
            for(int lz=0;lz<nz;lz++){
                double u0x=(4*2.7/(ny*ny))*ly*(ny-ly);
                double u1z=(2.7*0.01)*sin(2*3.1415*4*lz/nz)*sin(2*3.1415*4*ly/ny)*sin(2*3.1415*10*lx/nx);
                double u1y=(2.7*0.01)*cos(2*3.1415*4*lz/nz)*cos(2*3.1415*4*ly/ny)*cos(2*3.1415*10*lx/nx);
                //Vector3d u=Vector3d(u0x,u1z,u1y);
                Vector3d u=Vector3d(u0x,0.01,0);

                for(k=0;k<27;k++){
                    double cu=c[k]*u;
                    f[lx][ly][lz][k]=w[k]*1.*(1+3*cu+4.5*cu*cu-1.5*(u*u));
                    //f[lx][ly][lz][k]=w[k]*1.;// zero initial velocity
                    //cout<<f[lx][ly][lz][k]<<" ";
                   }
             }
        }
    }
}
template<int nx,int ny,int nz>
void Lattice<nx,ny,nz>::collide(){
#pragma omp parallel for
        for(int lx=0;lx<nx;lx++){
            for(int ly=0;ly<ny;ly++){
                for(int lz=0;lz<nz;lz++){
                u[lx][ly][lz]=phys.updatef(&f[lx][ly][lz][0],&c[0],&w[0]);
                }
           }

       }
}
template<int nx,int ny,int nz>
void Lattice<nx,ny,nz>::stream(){
#pragma omp parallel for
        for(int k=0;k<27;k++){
            if(k!=ind(0,0,0)){
                auto run_contraction=Vector3d(c[k].x*c[k].x,c[k].y*c[k].y,c[k].z*c[k].z);
                auto run_dir=Vector3d(c[k].x>0?-1:1,c[k].y>0?-1:1,c[k].z>0?-1:1);
                auto run_start=Vector3d(c[k].x>0?nx-1:0,c[k].y>0?ny-1:0,c[k].z>0?nz-1:0);
                auto run_shift=(-1)*c[k];

                for(int lx=0;lx<(nx-run_contraction.x);lx++){
                    for(int ly=0;ly<(ny-run_contraction.y);ly++){
                        for(int lz=0;lz<(nz-run_contraction.z);lz++){

 f[(int)(lx*run_dir.x+run_start.x)][(int)(ly*run_dir.y+run_start.y)][(int)(lz*run_dir.z+run_start.z)][k]=
 f[(int)(run_start.x+run_dir.x*lx+run_shift.x)][(int)(run_start.y+run_dir.y*ly+run_shift.y)][(int)(run_start.z+run_dir.z*lz+run_shift.z)][k];
                         }
                     }

                }
            }


        }

}
template<int nx,int ny,int nz>
void Lattice<nx,ny,nz>::bouncebackBC(Vector3d dir){
    Vector3d l0=Vector3d((dir.x==1)?nx-1:0,(dir.y==1)?ny-1:0,(dir.z==1)?nz-1:0);
    Vector3d ln=Vector3d((dir.x==-1)?1:nx,(dir.y==-1)?1:ny,(dir.z==-1)?1:nz);

#pragma omp parallel for
     for(int lx=l0.x;lx<(int)ln.x;lx++){
        for(int ly=l0.y;ly<(int)ln.y;ly++){
            for(int lz=l0.z;lz<(int)ln.z;lz++){
//
                for(int k=0;k<27;k++){
                    if((c[k]*dir)<0){
                        f[lx][ly][lz][k]=f[lx][ly][lz][ind(-c[k].x,-c[k].y,-c[k].z)];
                    }
                }
            }
        }
     }
}
template<int nx,int ny,int nz>
void Lattice<nx,ny,nz>::periodicBC(Vector3d dir){
    Vector3d l0=Vector3d((dir.x==1)?nx-1:0,(dir.y==1)?ny-1:0,(dir.z==1)?nz-1:0);
    Vector3d ln=Vector3d((dir.x==-1)?1:nx,(dir.y==-1)?1:ny,(dir.z==-1)?1:nz);

#pragma omp parallel for
    for(int lx=l0.x;lx<(int)ln.x;lx++){
        for(int ly=l0.y;ly<(int)ln.y;ly++){
            for(int lz=l0.z;lz<(int)ln.z;lz++){

                for(int k=0;k<27;k++){
                    if((c[k]*dir)<0){
                        f[lx][ly][lz][k]=f[(int)(lx-dir.x*(nx-1))][(int)(ly-dir.y*(ny-1))][(int)(lz-dir.z*(nz-1))][k];
                    }
                    if((c[k]*dir)>0){
                        f[(int)(lx-dir.x*(nx-1))][(int)(ly-dir.y*(ny-1))][(int)(lz-dir.z*(nz-1))][k]=f[lx][ly][lz][k];
                    }
                }
            }
        }
     }
}
template<int nx,int ny,int nz>
void Lattice<nx,ny,nz>::collectStat(double* umean,double* usqrt){
        for(int ly=0;ly<ny;ly++){
            umean[ly]+=u[nx/2][ly][nz/2].x;
            usqrt[ly]+=u[nx/2][ly][nz/2].x*u[nx/2][ly][nz/2].x;
        }
}
template<int nx,int ny,int nz>
void Lattice<nx,ny,nz>::mpi_Z(MPI_Comm grid_comm){
    // mpi exchange in z direction
    MPI_Cart_shift(grid_comm,0,1,&mpi_grid_rank,&mpi_up_z);
    MPI_Sendrecv(&f[0][0][nz-1][ind(-1,-1,1)],1,MpiVectorType2,mpi_up_z,1,&f[0][0][0][ind(-1,-1,1)],1,MpiVectorType2,mpi_grid_rank,1,grid_comm,&status);

    MPI_Cart_shift(grid_comm,0,-1,&mpi_grid_rank,&mpi_down_z);
    MPI_Sendrecv(&f[0][0][0][ind(-1,-1,-1)],1,MpiVectorType2,mpi_down_z,2,&f[0][0][nz-1][ind(-1,-1,-1)],1,MpiVectorType2,mpi_grid_rank,2,grid_comm,&status);

}

template<int nx,int ny,int nz>
void Lattice<nx,ny,nz>::saveResult(string fname){
    ofstream file;
    file.open(fname,ios::out);
    file<<"X  Y   Z   u   v   w \n";
    //cout<<"X  Y   Z   u   v   w \n";

    for(int lx=0;lx<nx;lx++){
        for(int ly=0;ly<ny;ly++){

            file<<lx<<"    "<<ly<<" 0  "<<u[lx][ly][nz/2]<<"\n";
            //cout<<lx<<"    "<<ly<<" 0  "<<u[lx][ly][nz/2]<<"\n";

        }

    }
    file.flush();
    file.close();
}
template<int nx,int ny,int nz>
void Lattice<nx,ny,nz>::saveState(string fbname){
    //Save state
    ofstream bfileOut;
    bfileOut.open(fbname, ios::out | ios::binary);
    bfileOut.write((char *) &f, sizeof f);
    bfileOut.close();
}
template<int nx,int ny,int nz>
void Lattice<nx,ny,nz>::restoreState(string fbname){
    ifstream bfileIn;
    bfileIn.open(fbname, ios::out | ios::binary);
    bfileIn.read((char *) &f, sizeof f);
    bfileIn.close();

}


#endif // LATTICE_H
