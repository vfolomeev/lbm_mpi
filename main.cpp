#include <iostream>

#include "fstream"
#include "math_api.h"
#include "omp.h"
#include <time.h>

using namespace std;

int main()
{
    int const nx=1000; //number of lattice nodes in two dimentions
    int const ny=40;
    int const nz=10;
    //int const n=nx*ny*nz;

    //double f[nx][ny][nz][27],rho[nx][ny][nz],feq,x[nx],y[ny],z[nz],w[27],cu,rhon;
    double w[27];
    double feq,cu,rhon;
    auto f = new double[nx][ny][nz][27];
    Vector3d c[27];
    auto u = new Vector3d[nx][ny][nz];
    Vector3d sumu;
    //double dt=1.0,dy=1.0,dx=dy,dz=dy;
    double u0=0.2,rho0=5.;
    int i,j,k;
        int mstep=1000;//total number of steps

    double alpha=0.02;
    //double Re=u0*nx/alpha;
    double sum;
    //double csq=dx*dx/(dt*dt);

    double omega=1.0/(3.*alpha+0.5);
// start wall clock
   double wall_timer = omp_get_wtime();
  // clock_t clock_timer = clock();

    //init c vectors
    for(int lx=-1;lx<2;lx++){
        for(int ly=-1;ly<2;ly++){
            for(int lz=-1;lz<2;lz++){
                k=ind(lx,ly,lz);
                c[k]=Vector3d(lx,ly,lz);
                if((c[k]*c[k])==0) w[k]=8./27;
                if(c[k]*c[k]==1) w[k]=2./27;
                if(c[k]*c[k]==2) w[k]=1./54;
                if(c[k]*c[k]==3) w[k]=1./216;


            }
        }
    }


omp_set_num_threads(4);

    //initiate solution
#pragma omp parallel for private(k)
    for(int lx=0;lx<nx;lx++){
        for(int ly=0;ly<ny;ly++){
            for(int lz=0;lz<nz;lz++){
                //rho[lx][ly][lz]=rho0;// Can be deleted to reduce memory consumption
                u[lx][ly][lz]=Vector3d();//TODO delete

                for(k=0;k<27;k++){
                    f[lx][ly][lz][k]=w[k]*rho0;// zero initial velocity


                   }
             }
        }
    }

    //main loop

    for(int t=0;t<mstep;t++){
       // qDebug()<<t;
#pragma omp parallel for private(i,j,k,sum,sumu,feq,cu,rhon)

        for(int lx=0;lx<nx;lx++){
            for(int ly=0;ly<ny;ly++){
                for(int lz=0;lz<nz;lz++){
                //compute rho and u
                    sum=0.;
                    sumu=Vector3d();

                    for(k=0;k<27;k++){
                        sum=sum+f[lx][ly][lz][k];

                        sumu=sumu+(f[lx][ly][lz][k])*c[k];


                    }
                    rhon=sum;

                    u[lx][ly][lz]=sumu*(1/rhon);

                //collision
                    for(k=0;k<27;k++){
                        cu=c[k]*u[lx][ly][lz];
                        feq=w[k]*rhon*(1+3*cu+4.5*cu*cu-1.5*(u[lx][ly][lz]*u[lx][ly][lz]));//feq=w[k]*rhon*(1+3*cu+4.5*cu*cu-1.5*(u[lx][ly][lz]*u[lx][ly][lz]));
                        f[lx][ly][lz][k]=omega*feq+(1.-omega)*f[lx][ly][lz][k];
                    }
                }
           }

       }

        //streaming

#pragma omp parallel for private(i,j)
        for(k=0;k<27;k++){
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


        //boundary conditions
// z normal bc
#pragma omp parallel for private(i,j,rhon)
        for(int l1=0;l1<nx;l1++){
            for(int l2=0;l2<ny;l2++){
                for(int k1=-1;k1<2;k1++){
                    for(int k2=-1;k2<2;k2++){
            //z=0 sym
                    f[l1][l2][0][ind(k1,k2,1)]=f[l1][l2][0][ind(k1,k2,-1)];
            //z=l sym
                    f[l1][l2][nz-1][ind(k1,k2,-1)]=f[l1][l2][nz-1][ind(k1,k2,1)];

                    }
               }
           }

        }
// x normal bc
#pragma omp parallel for private(i,j,rhon,k)
        for(int l1=0;l1<ny;l1++){
            for(int l2=0;l2<nz;l2++){
                for(int k1=-1;k1<2;k1++){
                    for(int k2=-1;k2<2;k2++){
                        double rho=0;
                        for(k=0;k<27;k++){
                            if(c[k]*Vector3d(-1,0,0)>=0){
                                rho=rho+f[0][l1][l2][k];
                            }
                            else{
                                rho=rho+f[0][l1][l2][ind(-c[k].x,-c[k].x,-c[k].x)];
                            }
                        }
                //x=0 inlet
                        f[0][l1][l2][ind(1,k1,k2)]=f[0][l1][l2][ind(-1,-k1,-k2)]+3*2*w[ind(1,k1,k2)]*rho*1.25*c[ind(1,k1,k2)]*Vector3d(u0,0,0);
                //x=l outlet
                        f[nx-1][l1][l2][ind(-1,k1,k2)]=2*f[nx-2][l1][l2][ind(-1,k1,k2)]-f[nx-3][l1][l2][ind(-1,k1,k2)];
                    }
                }
            }
        }
// y normal bc
#pragma omp parallel for private(i,j,rhon)
        for(int l1=0;l1<nx;l1++){
            for(int l2=0;l2<nz;l2++){
                for(int k1=-1;k1<2;k1++){
                    for(int k2=-1;k2<2;k2++){
                //y=0 bounce back
                        f[l1][0][l2][ind(k1,1,k2)]=f[l1][0][l2][ind(-k1,-1,-k2)];
                //y=l bounce back
                        //rhon=f[i][m-1][ind(0,0)]+f[i][m-1][ind(1,0)]+f[i][m-1][ind(-1,0)]+2.*(f[i][m-1][ind(0,1)]+f[i][m-1][ind(-1,1)]+f[i][m-1][ind(1,1)]);

                        f[l1][ny-1][l2][ind(k1,-1,k2)]=f[l1][ny-1][l2][ind(-k1,1,-k2)];
                    }
                }
            }
        }

}
    std::cout<<   " time on wall: " <<  omp_get_wtime() - wall_timer << "\n";
    //writing results file
    ofstream file;
    file.open("D:\\Development\\lbm_mpi\\res.csv",ios::out);

    file<<"X  Y   Z   u   v   w \n";
    for(int lx=0;lx<nx;lx++){
        for(int ly=0;ly<ny;ly++){

            file<<lx<<"    "<<ly<<" 0  "<<u[lx][ly][nz/2]<<"\n";

        }

    }
    file.flush();
    file.close();
}
