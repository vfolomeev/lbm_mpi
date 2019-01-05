#include <iostream>
#include <string>

#include "fstream"
#include "math_api.h"
#include "omp.h"
#include <time.h>
#include "mpi.h"
#include <cmath>
#include "lattice.h"
using namespace std;

int main(int argc,char *argv[])
{
//Initialize MPI
    MPI_Comm grid_comm;
    MPI_Init(&argc, &argv);    
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    int dims[1];        //mpi dimensions
    int periodic[1];    //mpi dimention periodicity
    int reorder = 0,  ndims = 1, maxdims = mpi_size;
    int coordinates[1];
    int coords[2];
    MPI_Comm row_comm;
    dims[0] = mpi_size;
    cout<<"mpi size: "<<mpi_size<<"\n";
    periodic[0] = 1;
    coords[0] = 0;
    double mpi_w_time=MPI_Wtime();
//create cartesian topology
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, reorder, &grid_comm);
    MPI_Comm_rank(grid_comm, &mpi_rank);
    MPI_Cart_coords(grid_comm,mpi_rank,ndims,coordinates);

    int const nx=400; //number of lattice nodes in two dimentions
    int const ny=240;
    int const nz=12;

       int mstep=10000;//total number of steps

//results and autosave file
    ofstream file;
    ofstream ufile;
    ofstream bfileOut;
    ifstream bfileIn;
    string fname="D:/Development/lbm_mpi/res";
    fname+=std::to_string(coordinates[0]);
    string fbname=fname;
    string ufname=fname;
    ufname+="vel.csv";
    fbname+=".bin";
    fname+=".csv";
    int stat_start=8.5e5;
    double umean[ny],usqrt[ny];
    for(int ly=0;ly<ny;ly++){
        umean[ly]=usqrt[ly]=0;
    }

// start wall clock
   double wall_timer = omp_get_wtime();
  // clock_t clock_timer = clock();

   //Init lattice
   Physics phys=Physics();

   auto mesh=new Lattice<nx,ny,nz>(phys);



omp_set_num_threads(4);

    //initiate solution

    mesh->initSolution();
    //mesh->restoreState(fbname);



    //main loop

    for(int t=0;t<mstep;t++){

       // Collide

        mesh->collide();

        //streaming

        mesh->stream();

        //boundary conditions


        //mpi exchange in z direction

        mesh->mpi_Z(grid_comm);

        //Boundary conditions

        //Periodic in X direction

        mesh->periodicBC(Vector3d(1,0,0));

        //Wall normal to Y

        mesh->bouncebackBC(Vector3d(0,1,0));
        mesh->bouncebackBC(Vector3d(0,-1,0));

        //Legacy

                //x=0 inlet
                        //f[0][l1][l2][ind(1,k1,k2)]=f[0][l1][l2][ind(-1,-k1,-k2)]+3*2*w[ind(1,k1,k2)]*rho*1.25*c[ind(1,k1,k2)]*Vector3d(u0,0,0);

                //x=l outlet
                        //f[nx-1][l1][l2][ind(-1,k1,k2)]=2*f[nx-2][l1][l2][ind(-1,k1,k2)]-f[nx-3][l1][l2][ind(-1,k1,k2)];

        //Save results
        cout<<t<<"\n";
        if ((t % 10)==0) {
            mesh->saveResult(fname);
            mesh->saveState(fbname);
            cout<<fname;
            }
            //ufile.flush();
            //ufile.close();
            cout.flush();
        }

    std::cout<<   " time on wall: " <<  MPI_Wtime() - mpi_w_time << "\n";

    MPI_Finalize();
}
