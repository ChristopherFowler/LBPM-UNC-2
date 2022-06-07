#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

//#include "common/pmmc.h"
#include "common/Domain.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"

#include "RayTrace3.h"




using namespace std;

int RT_DOUBLE = 8;
int RT_INT = 1;



inline int populate(int ijk ,
             int n  ,
             int m   ,
             int dir ,
             double * LIBBqA, double * LIBBqBC, double * LIBBqD, double * qDistances,
             int N  ,
             int Np
) {
    int count = 0;
    double A,BC,D;
    A = BC = D = 0; // Initialize all weights to zero
    double q = qDistances[ijk+dir*N];
    if (m == 0) { /* if neighbor is a solid site */
        if (q >= 0.5) {
            A = 1.-1./(2*q);
            BC = 1./(2*q);
            D = 0;
        }
        if (q < 0.5) {
            A = 0;
            BC = 2*q;
            D = (1.-2*q);
            
        }
        LIBBqA[n+dir*Np] = A;
        LIBBqBC[n+dir*Np] = BC;
        LIBBqD[n+dir*Np] = D;
        count++;
    }
    else { /* if neighbor is a fluid site */
        LIBBqA[n+dir*Np] = 0;
        LIBBqBC[n+dir*Np] = 1;
        LIBBqD[n+dir*Np] = 0;
        qDistances[ijk+dir*N] = 1;
    }
    return count;
}


inline int LIBB_Init(int * Map,
char *id,
int strideY,
int strideZ,
int start,
int finish,
int N, int Np, double * qDistances, double * LIBBqA, double * LIBBqBC, double * LIBBqD) {
 
 int nn,ijk,m;
 size_t count = 0;
 for (int n=start; n<finish; n++){
     // Get the 1D index based on regular data layout
     ijk = Map[n]; // standard coordinates
     //........................................................................
     nn = ijk-1;                            // neighbor index (get convention)
     m = id[nn];                        // get neighbor for id - 2
     count += populate(ijk,n,m,1,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np); 
//     //........................................................................
     nn = ijk+1;                            // neighbor index (get convention)
     m = id[nn];                        // get neighbor for id - 1
     count += populate(ijk,n,m,0,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk-strideY;                            // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 4
     count += populate(ijk,n,m,3,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk+strideY;                            // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 3
     count += populate(ijk,n,m,2,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk-strideZ;                        // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 6
     count += populate(ijk,n,m,5,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk+strideZ;                        // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 5
     count += populate(ijk,n,m,4,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk-strideY-1;                        // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 8
     count += populate(ijk,n,m,7,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk+strideY+1;                        // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 7
     count += populate(ijk,n,m,6,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk+strideY-1;                        // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 10
     count += populate(ijk,n,m,9,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk-strideY+1;                        // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 9
     count += populate(ijk,n,m,8,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk-strideZ-1;                        // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 12
     count += populate(ijk,n,m,11,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk+strideZ+1;                        // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 11
     count += populate(ijk,n,m,10,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk+strideZ-1;                        // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 14
     count += populate(ijk,n,m,13,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk-strideZ+1;                        // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 13
     count += populate(ijk,n,m,12,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk-strideZ-strideY;                    // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 16
     count += populate(ijk,n,m,15,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk+strideZ+strideY;                    // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 15
     count += populate(ijk,n,m,14,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk+strideZ-strideY;                    // neighbor index (get convention)
     m = id[nn];                   // get neighbor for id - 18
     count += populate(ijk,n,m,17,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     //........................................................................
     nn = ijk-strideZ+strideY;                    // neighbor index (get convention)
     m = id[nn];                    // get neighbor for id - 17
     count += populate(ijk,n,m,16,LIBBqA,LIBBqBC,LIBBqD,qDistances,N,Np);
     
 }
 return 0;
}

inline double calculate_porosity(char * id,
                              int Nx, int Ny, int Nz, Domain &Dm, int rank) {
 int nprocx,nprocy,nprocz;
 
 nprocx = Dm.nprocx();
 nprocy = Dm.nprocy();
 nprocz = Dm.nprocz();
    size_t count = 0;
    size_t totalGlobal = 0;
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                size_t n = k*Nx*Ny+j*Nx+i;
                if (id[n] == 2)  count++;
            }
        }
    }
    
    // total Global is the number of nodes in the pore-space
    MPI_Allreduce(&count,&totalGlobal,1,MPI_INT,MPI_SUM,Dm.Comm);
    double volume=double(nprocx*nprocy*nprocz)*double(Nx-2)*double(Ny-2)*double(Nz-2);
    double porosity=double(totalGlobal)/volume;
    
    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    double dcount = 0;
    double dcountGlobal = 0;
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                for (int p=0;p<8;p++){
                    int n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
                    if (id[n] == 2)  {  dcount += 0.125; }
                }
            }
        }
    }
    MPI_Allreduce(&dcount,&dcountGlobal,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    porosity=double(dcountGlobal)/volume;
 
 return porosity;
}


inline void fillBoolfield(bool* boolfield,
                          int Nx, int Ny, int Nz) {
    
    
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                size_t n = k*Nx*Ny+j*Nx+i;
                boolfield[n] = false;
            }
        }
    }
    
    
    
}

inline void PackID(int *list,
                   int count, char *sendbuf, char *ID) {
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        sendbuf[idx] = ID[n];
    }
}

inline void UnpackID(int *list,
                     int count, char *recvbuf, char *ID) {
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        ID[n] = recvbuf[idx];
    }
}

inline void PackID_Double(int *list,
                          int count, double *sendbuf, double *ID) {
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        sendbuf[idx] = ID[n];
    }
}

inline void UnpackID_Double(int *list,
                            int count, double *recvbuf, double *ID) {
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        ID[n] = recvbuf[idx];
    }
}

inline void PopulateHalo_Double(double *id,
                                Domain &Dm, int nx, int ny, int nz, int rank) {


    int Nx = nx;
    int Ny = ny;
    int Nz = nz;

    double count;
    count = 0.f;

    // Communication buffers
    double *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
    double *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
    double *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
    double *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
    double *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
    double *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
    // send buffers
    sendID_x = new double [Dm.sendCount_x];
    sendID_y = new double [Dm.sendCount_y];
    sendID_z = new double [Dm.sendCount_z];
    sendID_X = new double [Dm.sendCount_X];
    sendID_Y = new double [Dm.sendCount_Y];
    sendID_Z = new double [Dm.sendCount_Z];
    sendID_xy = new double [Dm.sendCount_xy];
    sendID_yz = new double [Dm.sendCount_yz];
    sendID_xz = new double [Dm.sendCount_xz];
    sendID_Xy = new double [Dm.sendCount_Xy];
    sendID_Yz = new double [Dm.sendCount_Yz];
    sendID_xZ = new double [Dm.sendCount_xZ];
    sendID_xY = new double [Dm.sendCount_xY];
    sendID_yZ = new double [Dm.sendCount_yZ];
    sendID_Xz = new double [Dm.sendCount_Xz];
    sendID_XY = new double [Dm.sendCount_XY];
    sendID_YZ = new double [Dm.sendCount_YZ];
    sendID_XZ = new double [Dm.sendCount_XZ];
    //......................................................................................
    // recv buffers
    recvID_x = new double [Dm.recvCount_x];
    recvID_y = new double [Dm.recvCount_y];
    recvID_z = new double [Dm.recvCount_z];
    recvID_X = new double [Dm.recvCount_X];
    recvID_Y = new double [Dm.recvCount_Y];
    recvID_Z = new double [Dm.recvCount_Z];
    recvID_xy = new double [Dm.recvCount_xy];
    recvID_yz = new double [Dm.recvCount_yz];
    recvID_xz = new double [Dm.recvCount_xz];
    recvID_Xy = new double [Dm.recvCount_Xy];
    recvID_xZ = new double [Dm.recvCount_xZ];
    recvID_xY = new double [Dm.recvCount_xY];
    recvID_yZ = new double [Dm.recvCount_yZ];
    recvID_Yz = new double [Dm.recvCount_Yz];
    recvID_Xz = new double [Dm.recvCount_Xz];
    recvID_XY = new double [Dm.recvCount_XY];
    recvID_YZ = new double [Dm.recvCount_YZ];
    recvID_XZ = new double [Dm.recvCount_XZ];
    //......................................................................................
    int sendtag,recvtag;
    sendtag = recvtag = 7;
    
    // Pack and send the updated ID values
    PackID_Double(Dm.sendList_x, Dm.sendCount_x ,sendID_x, id);
    PackID_Double(Dm.sendList_X, Dm.sendCount_X ,sendID_X, id);
    PackID_Double(Dm.sendList_y, Dm.sendCount_y ,sendID_y, id);
    PackID_Double(Dm.sendList_Y, Dm.sendCount_Y ,sendID_Y, id);
    PackID_Double(Dm.sendList_z, Dm.sendCount_z ,sendID_z, id);
    PackID_Double(Dm.sendList_Z, Dm.sendCount_Z ,sendID_Z, id);
    PackID_Double(Dm.sendList_xy, Dm.sendCount_xy ,sendID_xy, id);
    PackID_Double(Dm.sendList_Xy, Dm.sendCount_Xy ,sendID_Xy, id);
    PackID_Double(Dm.sendList_xY, Dm.sendCount_xY ,sendID_xY, id);
    PackID_Double(Dm.sendList_XY, Dm.sendCount_XY ,sendID_XY, id);
    PackID_Double(Dm.sendList_xz, Dm.sendCount_xz ,sendID_xz, id);
    PackID_Double(Dm.sendList_Xz, Dm.sendCount_Xz ,sendID_Xz, id);
    PackID_Double(Dm.sendList_xZ, Dm.sendCount_xZ ,sendID_xZ, id);
    PackID_Double(Dm.sendList_XZ, Dm.sendCount_XZ ,sendID_XZ, id);
    PackID_Double(Dm.sendList_yz, Dm.sendCount_yz ,sendID_yz, id);
    PackID_Double(Dm.sendList_Yz, Dm.sendCount_Yz ,sendID_Yz, id);
    PackID_Double(Dm.sendList_yZ, Dm.sendCount_yZ ,sendID_yZ, id);
    PackID_Double(Dm.sendList_YZ, Dm.sendCount_YZ ,sendID_YZ, id);
    //......................................................................................
    MPI_Sendrecv(sendID_x,Dm.sendCount_x,MPI_DOUBLE,Dm.rank_x(),sendtag,
                 recvID_X,Dm.recvCount_X,MPI_DOUBLE,Dm.rank_X(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_X,Dm.sendCount_X,MPI_DOUBLE,Dm.rank_X(),sendtag,
                 recvID_x,Dm.recvCount_x,MPI_DOUBLE,Dm.rank_x(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_y,Dm.sendCount_y,MPI_DOUBLE,Dm.rank_y(),sendtag,
                 recvID_Y,Dm.recvCount_Y,MPI_DOUBLE,Dm.rank_Y(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Y,Dm.sendCount_Y,MPI_DOUBLE,Dm.rank_Y(),sendtag,
                 recvID_y,Dm.recvCount_y,MPI_DOUBLE,Dm.rank_y(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_z,Dm.sendCount_z,MPI_DOUBLE,Dm.rank_z(),sendtag,
                 recvID_Z,Dm.recvCount_Z,MPI_DOUBLE,Dm.rank_Z(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Z,Dm.sendCount_Z,MPI_DOUBLE,Dm.rank_Z(),sendtag,
                 recvID_z,Dm.recvCount_z,MPI_DOUBLE,Dm.rank_z(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xy,Dm.sendCount_xy,MPI_DOUBLE,Dm.rank_xy(),sendtag,
                 recvID_XY,Dm.recvCount_XY,MPI_DOUBLE,Dm.rank_XY(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_XY,Dm.sendCount_XY,MPI_DOUBLE,Dm.rank_XY(),sendtag,
                 recvID_xy,Dm.recvCount_xy,MPI_DOUBLE,Dm.rank_xy(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Xy,Dm.sendCount_Xy,MPI_DOUBLE,Dm.rank_Xy(),sendtag,
                 recvID_xY,Dm.recvCount_xY,MPI_DOUBLE,Dm.rank_xY(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xY,Dm.sendCount_xY,MPI_DOUBLE,Dm.rank_xY(),sendtag,
                 recvID_Xy,Dm.recvCount_Xy,MPI_DOUBLE,Dm.rank_Xy(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xz,Dm.sendCount_xz,MPI_DOUBLE,Dm.rank_xz(),sendtag,
                 recvID_XZ,Dm.recvCount_XZ,MPI_DOUBLE,Dm.rank_XZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_XZ,Dm.sendCount_XZ,MPI_DOUBLE,Dm.rank_XZ(),sendtag,
                 recvID_xz,Dm.recvCount_xz,MPI_DOUBLE,Dm.rank_xz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Xz,Dm.sendCount_Xz,MPI_DOUBLE,Dm.rank_Xz(),sendtag,
                 recvID_xZ,Dm.recvCount_xZ,MPI_DOUBLE,Dm.rank_xZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xZ,Dm.sendCount_xZ,MPI_DOUBLE,Dm.rank_xZ(),sendtag,
                 recvID_Xz,Dm.recvCount_Xz,MPI_DOUBLE,Dm.rank_Xz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_yz,Dm.sendCount_yz,MPI_DOUBLE,Dm.rank_yz(),sendtag,
                 recvID_YZ,Dm.recvCount_YZ,MPI_DOUBLE,Dm.rank_YZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_YZ,Dm.sendCount_YZ,MPI_DOUBLE,Dm.rank_YZ(),sendtag,
                 recvID_yz,Dm.recvCount_yz,MPI_DOUBLE,Dm.rank_yz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Yz,Dm.sendCount_Yz,MPI_DOUBLE,Dm.rank_Yz(),sendtag,
                 recvID_yZ,Dm.recvCount_yZ,MPI_DOUBLE,Dm.rank_yZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_yZ,Dm.sendCount_yZ,MPI_DOUBLE,Dm.rank_yZ(),sendtag,
                 recvID_Yz,Dm.recvCount_Yz,MPI_DOUBLE,Dm.rank_Yz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    //......................................................................................
    UnpackID_Double(Dm.recvList_x, Dm.recvCount_x ,recvID_x, id);
    UnpackID_Double(Dm.recvList_X, Dm.recvCount_X ,recvID_X, id);
    UnpackID_Double(Dm.recvList_y, Dm.recvCount_y ,recvID_y, id);
    UnpackID_Double(Dm.recvList_Y, Dm.recvCount_Y ,recvID_Y, id);
    UnpackID_Double(Dm.recvList_z, Dm.recvCount_z ,recvID_z, id);
    UnpackID_Double(Dm.recvList_Z, Dm.recvCount_Z ,recvID_Z, id);
    UnpackID_Double(Dm.recvList_xy, Dm.recvCount_xy ,recvID_xy, id);
    UnpackID_Double(Dm.recvList_Xy, Dm.recvCount_Xy ,recvID_Xy, id);
    UnpackID_Double(Dm.recvList_xY, Dm.recvCount_xY ,recvID_xY, id);
    UnpackID_Double(Dm.recvList_XY, Dm.recvCount_XY ,recvID_XY, id);
    UnpackID_Double(Dm.recvList_xz, Dm.recvCount_xz ,recvID_xz, id);
    UnpackID_Double(Dm.recvList_Xz, Dm.recvCount_Xz ,recvID_Xz, id);
    UnpackID_Double(Dm.recvList_xZ, Dm.recvCount_xZ ,recvID_xZ, id);
    UnpackID_Double(Dm.recvList_XZ, Dm.recvCount_XZ ,recvID_XZ, id);
    UnpackID_Double(Dm.recvList_yz, Dm.recvCount_yz ,recvID_yz, id);
    UnpackID_Double(Dm.recvList_Yz, Dm.recvCount_Yz ,recvID_Yz, id);
    UnpackID_Double(Dm.recvList_yZ, Dm.recvCount_yZ ,recvID_yZ, id);
    UnpackID_Double(Dm.recvList_YZ, Dm.recvCount_YZ ,recvID_YZ, id);
    //......................................................................................
    
    /* Correct corner distance map */
    for (int k=0;k<2;k++){
       for (int j=0;j<2;j++){
           for (int i=0;i<2;i++){
               int n = k*Nx*Ny+j*Nx+i;
               if (i == 1 && j == 1 && k == 1) {int m = (Nz-1)*Nx*Ny + (Ny-1)*Nx+(Nz-1);    id[m] = id[n]; }
           }
       }
    }
}

inline void PopulateHalo_Char(char *id,
                              Domain &Dm, int nx, int ny, int nz, int rank) {


    int Nx = nx;
    int Ny = ny;
    int Nz = nz;

    double count;
    count = 0.f;

    // Communication buffers
    char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
    char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
    char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
    char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
    char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
    char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
    // send buffers
    sendID_x = new char [Dm.sendCount_x];
    sendID_y = new char [Dm.sendCount_y];
    sendID_z = new char [Dm.sendCount_z];
    sendID_X = new char [Dm.sendCount_X];
    sendID_Y = new char [Dm.sendCount_Y];
    sendID_Z = new char [Dm.sendCount_Z];
    sendID_xy = new char [Dm.sendCount_xy];
    sendID_yz = new char [Dm.sendCount_yz];
    sendID_xz = new char [Dm.sendCount_xz];
    sendID_Xy = new char [Dm.sendCount_Xy];
    sendID_Yz = new char [Dm.sendCount_Yz];
    sendID_xZ = new char [Dm.sendCount_xZ];
    sendID_xY = new char [Dm.sendCount_xY];
    sendID_yZ = new char [Dm.sendCount_yZ];
    sendID_Xz = new char [Dm.sendCount_Xz];
    sendID_XY = new char [Dm.sendCount_XY];
    sendID_YZ = new char [Dm.sendCount_YZ];
    sendID_XZ = new char [Dm.sendCount_XZ];
    //......................................................................................
    // recv buffers
    recvID_x = new char [Dm.recvCount_x];
    recvID_y = new char [Dm.recvCount_y];
    recvID_z = new char [Dm.recvCount_z];
    recvID_X = new char [Dm.recvCount_X];
    recvID_Y = new char [Dm.recvCount_Y];
    recvID_Z = new char [Dm.recvCount_Z];
    recvID_xy = new char [Dm.recvCount_xy];
    recvID_yz = new char [Dm.recvCount_yz];
    recvID_xz = new char [Dm.recvCount_xz];
    recvID_Xy = new char [Dm.recvCount_Xy];
    recvID_xZ = new char [Dm.recvCount_xZ];
    recvID_xY = new char [Dm.recvCount_xY];
    recvID_yZ = new char [Dm.recvCount_yZ];
    recvID_Yz = new char [Dm.recvCount_Yz];
    recvID_Xz = new char [Dm.recvCount_Xz];
    recvID_XY = new char [Dm.recvCount_XY];
    recvID_YZ = new char [Dm.recvCount_YZ];
    recvID_XZ = new char [Dm.recvCount_XZ];
    //......................................................................................
    int sendtag,recvtag;
    sendtag = recvtag = 7;
    
    // Pack and send the updated ID values
    PackID(Dm.sendList_x, Dm.sendCount_x ,sendID_x, id);
    PackID(Dm.sendList_X, Dm.sendCount_X ,sendID_X, id);
    PackID(Dm.sendList_y, Dm.sendCount_y ,sendID_y, id);
    PackID(Dm.sendList_Y, Dm.sendCount_Y ,sendID_Y, id);
    PackID(Dm.sendList_z, Dm.sendCount_z ,sendID_z, id);
    PackID(Dm.sendList_Z, Dm.sendCount_Z ,sendID_Z, id);
    PackID(Dm.sendList_xy, Dm.sendCount_xy ,sendID_xy, id);
    PackID(Dm.sendList_Xy, Dm.sendCount_Xy ,sendID_Xy, id);
    PackID(Dm.sendList_xY, Dm.sendCount_xY ,sendID_xY, id);
    PackID(Dm.sendList_XY, Dm.sendCount_XY ,sendID_XY, id);
    PackID(Dm.sendList_xz, Dm.sendCount_xz ,sendID_xz, id);
    PackID(Dm.sendList_Xz, Dm.sendCount_Xz ,sendID_Xz, id);
    PackID(Dm.sendList_xZ, Dm.sendCount_xZ ,sendID_xZ, id);
    PackID(Dm.sendList_XZ, Dm.sendCount_XZ ,sendID_XZ, id);
    PackID(Dm.sendList_yz, Dm.sendCount_yz ,sendID_yz, id);
    PackID(Dm.sendList_Yz, Dm.sendCount_Yz ,sendID_Yz, id);
    PackID(Dm.sendList_yZ, Dm.sendCount_yZ ,sendID_yZ, id);
    PackID(Dm.sendList_YZ, Dm.sendCount_YZ ,sendID_YZ, id);
    //......................................................................................
    MPI_Sendrecv(sendID_x,Dm.sendCount_x,MPI_CHAR,Dm.rank_x(),sendtag,
                 recvID_X,Dm.recvCount_X,MPI_CHAR,Dm.rank_X(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_X,Dm.sendCount_X,MPI_CHAR,Dm.rank_X(),sendtag,
                 recvID_x,Dm.recvCount_x,MPI_CHAR,Dm.rank_x(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_y,Dm.sendCount_y,MPI_CHAR,Dm.rank_y(),sendtag,
                 recvID_Y,Dm.recvCount_Y,MPI_CHAR,Dm.rank_Y(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Y,Dm.sendCount_Y,MPI_CHAR,Dm.rank_Y(),sendtag,
                 recvID_y,Dm.recvCount_y,MPI_CHAR,Dm.rank_y(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_z,Dm.sendCount_z,MPI_CHAR,Dm.rank_z(),sendtag,
                 recvID_Z,Dm.recvCount_Z,MPI_CHAR,Dm.rank_Z(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Z,Dm.sendCount_Z,MPI_CHAR,Dm.rank_Z(),sendtag,
                 recvID_z,Dm.recvCount_z,MPI_CHAR,Dm.rank_z(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xy,Dm.sendCount_xy,MPI_CHAR,Dm.rank_xy(),sendtag,
                 recvID_XY,Dm.recvCount_XY,MPI_CHAR,Dm.rank_XY(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_XY,Dm.sendCount_XY,MPI_CHAR,Dm.rank_XY(),sendtag,
                 recvID_xy,Dm.recvCount_xy,MPI_CHAR,Dm.rank_xy(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Xy,Dm.sendCount_Xy,MPI_CHAR,Dm.rank_Xy(),sendtag,
                 recvID_xY,Dm.recvCount_xY,MPI_CHAR,Dm.rank_xY(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xY,Dm.sendCount_xY,MPI_CHAR,Dm.rank_xY(),sendtag,
                 recvID_Xy,Dm.recvCount_Xy,MPI_CHAR,Dm.rank_Xy(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xz,Dm.sendCount_xz,MPI_CHAR,Dm.rank_xz(),sendtag,
                 recvID_XZ,Dm.recvCount_XZ,MPI_CHAR,Dm.rank_XZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_XZ,Dm.sendCount_XZ,MPI_CHAR,Dm.rank_XZ(),sendtag,
                 recvID_xz,Dm.recvCount_xz,MPI_CHAR,Dm.rank_xz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Xz,Dm.sendCount_Xz,MPI_CHAR,Dm.rank_Xz(),sendtag,
                 recvID_xZ,Dm.recvCount_xZ,MPI_CHAR,Dm.rank_xZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_xZ,Dm.sendCount_xZ,MPI_CHAR,Dm.rank_xZ(),sendtag,
                 recvID_Xz,Dm.recvCount_Xz,MPI_CHAR,Dm.rank_Xz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_yz,Dm.sendCount_yz,MPI_CHAR,Dm.rank_yz(),sendtag,
                 recvID_YZ,Dm.recvCount_YZ,MPI_CHAR,Dm.rank_YZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_YZ,Dm.sendCount_YZ,MPI_CHAR,Dm.rank_YZ(),sendtag,
                 recvID_yz,Dm.recvCount_yz,MPI_CHAR,Dm.rank_yz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_Yz,Dm.sendCount_Yz,MPI_CHAR,Dm.rank_Yz(),sendtag,
                 recvID_yZ,Dm.recvCount_yZ,MPI_CHAR,Dm.rank_yZ(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendID_yZ,Dm.sendCount_yZ,MPI_CHAR,Dm.rank_yZ(),sendtag,
                 recvID_Yz,Dm.recvCount_Yz,MPI_CHAR,Dm.rank_Yz(),recvtag,Dm.Comm,MPI_STATUS_IGNORE);
    //......................................................................................
    UnpackID(Dm.recvList_x, Dm.recvCount_x ,recvID_x, id);
    UnpackID(Dm.recvList_X, Dm.recvCount_X ,recvID_X, id);
    UnpackID(Dm.recvList_y, Dm.recvCount_y ,recvID_y, id);
    UnpackID(Dm.recvList_Y, Dm.recvCount_Y ,recvID_Y, id);
    UnpackID(Dm.recvList_z, Dm.recvCount_z ,recvID_z, id);
    UnpackID(Dm.recvList_Z, Dm.recvCount_Z ,recvID_Z, id);
    UnpackID(Dm.recvList_xy, Dm.recvCount_xy ,recvID_xy, id);
    UnpackID(Dm.recvList_Xy, Dm.recvCount_Xy ,recvID_Xy, id);
    UnpackID(Dm.recvList_xY, Dm.recvCount_xY ,recvID_xY, id);
    UnpackID(Dm.recvList_XY, Dm.recvCount_XY ,recvID_XY, id);
    UnpackID(Dm.recvList_xz, Dm.recvCount_xz ,recvID_xz, id);
    UnpackID(Dm.recvList_Xz, Dm.recvCount_Xz ,recvID_Xz, id);
    UnpackID(Dm.recvList_xZ, Dm.recvCount_xZ ,recvID_xZ, id);
    UnpackID(Dm.recvList_XZ, Dm.recvCount_XZ ,recvID_XZ, id);
    UnpackID(Dm.recvList_yz, Dm.recvCount_yz ,recvID_yz, id);
    UnpackID(Dm.recvList_Yz, Dm.recvCount_Yz ,recvID_Yz, id);
    UnpackID(Dm.recvList_yZ, Dm.recvCount_yZ ,recvID_yZ, id);
    UnpackID(Dm.recvList_YZ, Dm.recvCount_YZ ,recvID_YZ, id);
    //......................................................................................
    for (int k=0;k<2;k++){
          for (int j=0;j<2;j++){
              for (int i=0;i<2;i++){
                  int n = k*Nx*Ny+j*Nx+i;
                  if (i == 1 && j == 1 && k == 1) {int m = (Nz-1)*Nx*Ny + (Ny-1)*Nx+(Nz-1);    id[m] = id[n]; }
              }
          }
       }
}

int main(int argc, char **argv)
{
    // Initialize MPI
    int rank,nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);
    {
        if (rank == 0) cout << "========================================\nRAY TRACING PRE-PROCESSOR\n========================================\n";
        size_t n;
        size_t count = 0;
        
        int Nx,Ny,Nz;
        double Lx,Ly,Lz;
        int iproc,jproc,kproc;
        int nprocx,nprocy,nprocz;

        // Filenames used
        char LocalRankString[8];
        char LocalRankFilename[40];
        
        /* Read in input.db information */
        string filename;
        if (argc > 1) filename=argv[1];
        else ERROR("No input database provided\n");

        // read the input database
        auto db = std::make_shared<Database>( filename );
        auto domain_db = db->getDatabase( "Domain" );
        auto Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));
        auto Mask  = std::shared_ptr<Domain>(new Domain(domain_db,comm));
        auto analysis_db = db->getDatabase( "Analysis" );
        
        Nx = Dm->Nx;  Ny = Dm->Ny;  Nz = Dm->Nz;
        Lx = Dm->Lx;  Ly = Dm->Ly;  Lz = Dm->Lz;
        iproc = Dm->iproc();  jproc = Dm->jproc(); kproc = Dm->kproc();
        nprocx = Dm->nprocx();  nprocy = Dm->nprocy();  nprocz = Dm->nprocz();

        int * TmpMap;

        int nspheres = domain_db->getScalar<int>( "nspheres");
        auto ReadValues = domain_db->getVector<char>( "ReadValues" );  // Should this be int??
        auto WriteValues = domain_db->getVector<char>( "WriteValues" );
        
        auto rt_db = db->getDatabase( "RayTrace" );
        int geometry = rt_db->getScalar<int>( "geometry" );
        double offset = rt_db->getScalar<double>( "offset" );
        int BoundaryCondition = domain_db->getScalar<int>( "BC" );

        int rmin = 1;
        int rmax = Nz-1;
        if (BoundaryCondition > 0){
            if ( domain_db->keyExists( "ReservoirMin" ) && kproc == 0){
                rmin = domain_db->getScalar<int>( "ReservoirMin" );
            } 
            if ( domain_db->keyExists( "ReservoirMax" ) && kproc == (nprocz - 1)){
                rmax = domain_db->getScalar<int>( "ReservoirMax" );
                rmax = Nz-1-((Nz-2)*nprocz - rmax);
            }
        }

        int plate_min = 1;
        int plate_max = Nz-1;
        if ( rt_db->keyExists( "PlateMin" ) ){
            plate_min = rt_db->getScalar<int>( "PlateMin" );
        }
        if ( rt_db->keyExists( "PlateMax" ) ){
            plate_max = rt_db->getScalar<int>( "PlateMax" );
        }


        if (rank == 0) { cout << "Geometry:";
            if (geometry == 1 ) cout << "periodic sphere pack" << endl;
            if (geometry == 2 ) cout << "slab" << endl;
            if (geometry == 3 ) cout << "parallel plates" << endl;
        }
        
        /* Set up domain sizes with halo */
        int Nq = 18;

        int N = Nx*Ny*Nz;
       
        /* Initialize domain class communication */
        for (n=0; n<N; n++) Dm->id[n] = 1;  Dm->CommInit();  for (n=0; n<N; n++) Dm->id[n] = 2; // Set wetting phase (id=2)
        
        /* Declare ray tracing data arrays */
        auto boolfield = new bool[N];  for (n=0; n<N; n++) boolfield[n]=true; fillBoolfield(boolfield, Nx, Ny, Nz);
        
        double * TemporaryField = new double[N];
        for (size_t n = 0; n < N; n++) TemporaryField[n] = 0.0;
        
        /* Initialize ray-trace data */
        auto LIBB_id = new char[N];  for (n=0; n<N; n++) LIBB_id[n] = 2; // Set domain to wetting fluid
        auto qList = new double[Nq*N]; for (size_t n = 0; n < Nq*N; n++) qList[n] = -1;
        auto rtDistance = new double[N];  for (size_t n = 0; n < N; n++) rtDistance[n] = 100.0;
        auto offsetDistance = new double[N];  for (size_t n = 0; n < N; n++) offsetDistance[n] = 100.0;
        
        
        std::shared_ptr<TwoPhase> Averages;
        for (size_t i=0; i<N; i++) Dm->id[i] = LIBB_id[i];
        
        // Set up Analysis routines
        std::array<size_t,3> subdivide = { 1, 1, 1 };
        
        if ( analysis_db->keyExists( "subdivide" ) ) { auto tmp = analysis_db->getVector<size_t>( "subdivide" );  subdivide = { tmp[0], tmp[1], tmp[2] }; }
        
        Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm,subdivide) ); // TwoPhase analysis object
        MPI_Barrier(comm);
        
        
        
        
        // Sphere pack
        if ( geometry == 1 ) {
            ComputeLIBBqDistances_sp(LIBB_id,
                                  qList, rtDistance, offsetDistance, Nx, Ny, Nz, Lx, nspheres, iproc,  jproc,
                                  kproc,  nprocx,  nprocy,  nprocz, rank, *Dm,BoundaryCondition,rmin,rmax);
            /* Clear halo and then populate halo */
            for (int n = 0; n < N; n++) if (boolfield[n] == true) LIBB_id[n] = 0;
            for (int n = 0; n < N; n++) if (boolfield[n] == true) rtDistance[n] = 1;
            for (int n = 0; n < N; n++) if (boolfield[n] == true) offsetDistance[n] = 1;
            
            
            PopulateHalo_Char(LIBB_id, *Dm, Nx, Ny, Nz, rank);
            PopulateHalo_Double(rtDistance, *Dm, Nx, Ny, Nz, rank);
            PopulateHalo_Double(offsetDistance, *Dm, Nx, Ny, Nz, rank);
        } else if (geometry == 2 ) {
            ComputeLIBBqDistances_s(LIBB_id,
                                  qList, rtDistance, Nx, Ny, Nz, Lx, nspheres, iproc,  jproc,
                                  kproc,  nprocx,  nprocy,  nprocz, rank, *Dm, plate_min, plate_max);
            

            // for (int k=0;k<Nz;k++){
            //     for (int j=0;j<Ny;j++){
            //         for (int i=0;i<Nx;i++){
            //             int n=k*Nx*Ny+j*Nx+i;
            //             if (k==0){rtDistance[n] = 1.0;}
            //             if (k==1){rtDistance[n] = -1.0;}
            //             if (k==2){rtDistance[n] = -1.0;}
            //             if (k==3){rtDistance[n] = 1.0;}
            //             if (k==4){rtDistance[n] = 2.0;}
            //             if (k==5){rtDistance[n] = 1.0;}
            //             if (k==6){rtDistance[n] = 1.0;}
            //         }
            //     }
            // }

            for (size_t n = 0; n < N; n++) offsetDistance[n] = rtDistance[n];


        } else if (geometry == 3 ) {
            ComputeLIBBqDistances_pp(LIBB_id,
                                  qList, rtDistance, Nx, Ny, Nz, Lx, nspheres, iproc,  jproc,
                                  kproc,  nprocx,  nprocy,  nprocz, rank, *Dm, plate_min, plate_max);

            for (size_t n = 0; n < N; n++) offsetDistance[n] = rtDistance[n];


        }
        

        // if (rank==0){
        //     printf("rtDistance \n");
        // }
        // for (int k=0;k<Nz;k++){
        //     for (int j=0;j<Ny;j++){
        //         for (int i=0;i<Nx;i++){
        //             int n=k*Nx*Ny+j*Nx+i;
        //             printf("%.2f ",rtDistance[n]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n\n");
        // }

        
        for (int n=0; n<N; n++) Averages->SDs(n) = rtDistance[n];
        //if (geometry == 1){
            Averages->ComputeVolumeFraction(Averages->SDs, Averages->SDs_x, Averages->SDs_y, Averages->SDs_z, Averages->VFmask, true, rmin, rmax, geometry);
        //} else if (geometry == 2){
        //     for (int i = 0; i < N; i++) {
        //         if (LIBB_id[i] == 0){
        //             Averages->VFmask(i) = 1.0;
        //         } else {
        //             Averages->VFmask(i) = 0.0;
        //         }
        //     }
        //     for (int k=0;k<Nz;k++){
        //         for (int j=0;j<Ny;j++){
        //             for (int i=0;i<Nx;i++){
        //                 int n=k*Nx*Ny+j*Nx+i;
        //                 if (k==plate_max+1){
        //                     Averages->GradPhiZ(n) = 1.0;
        //                 } 
        //             }
        //         }
        //     }
        // }
        
        Dm->CommunicateMeshHalo(Averages->VFmask);
        
        for (int i = 0; i < N; i++) {
            if (Averages->VFmask(i) > 0.5) LIBB_id[i] = 0;
            else LIBB_id[i] = 2;
        }

        // if (rank==0){
        //     printf("VFMask \n");
        // }
        // for (int k=0;k<Nz;k++){
        //     for (int j=0;j<Ny;j++){
        //         for (int i=0;i<Nx;i++){
        //             int n=k*Nx*Ny+j*Nx+i;
        //             printf("%.2f ",Averages->VFmask(n));
        //         }
        //         printf("\n");
        //     }
        //     printf("\n\n");
        // }
        
        // if (rank==0){
        //     printf("LIBB_id \n");
        // }
        // for (int k=0;k<Nz;k++){
        //     for (int j=0;j<Ny;j++){
        //         for (int i=0;i<Nx;i++){
        //             int n=k*Nx*Ny+j*Nx+i;
        //             printf("%i ",LIBB_id[n]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n\n");
        // }

        /* Check the porosity */
        double porosity = calculate_porosity(LIBB_id, Nx, Ny, Nz, *Dm, rank);
        if (rank == 0) cout << "Actual id porosity=" << porosity << endl;
        
       
    
        
        for (n=0; n<N; n++) Mask->id[n] = 1;
        Mask->CommInit();
        
        if (BoundaryCondition > 0){
            if (Dm->kproc()==0){
                int n;
                int N = Nx*Ny*Nz;
                for (int Slice = 0; Slice < rmin; Slice++) { 
                    for (n=Slice*Nx*Ny; n<(Slice+1)*Nx*Ny; n++){
                        LIBB_id[n] = 2;
                    }
                }
            }
            
            if (Dm->kproc() == nprocz-1){
                int n;
                int N = Nx*Ny*Nz;
                for (int Slice = rmax; Slice < Nz; Slice++) {
                    for (n=Slice*Nx*Ny; n<(Slice+1)*Nx*Ny; n++){
                        LIBB_id[n] = 2;
                    }
                }
            }
        }

        /* May need to modify these fields for pressure bc */
        /* Write signed distance files */
        sprintf(LocalRankFilename,"SignDist.%05i",rank);
        FILE *RTFILE = fopen(LocalRankFilename,"wb");
        if (RTFILE==NULL) ERROR("Error opening file: rt.xxxxx");
        fwrite(rtDistance,RT_DOUBLE,N,RTFILE);
        fclose(RTFILE);
        
        sprintf(LocalRankFilename,"OffsetDist.%05i",rank);
        FILE *OFFSETFILE = fopen(LocalRankFilename,"wb");
        if (OFFSETFILE==NULL) ERROR("Error opening file: offset.xxxxx");
        fwrite(offsetDistance,RT_DOUBLE,N,OFFSETFILE);
        fclose(OFFSETFILE);
        
        for (int i=0; i<Nx*Ny*Nz; i++) Mask->id[i] = LIBB_id[i];
        size_t Np=Mask->PoreCount();
        //...........................................................................
        auto ScaLBL_Comm  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));
        
        size_t Npad=(Np/16 + 2)*16;
        IntArray Map;
        Map.resize(Nx,Ny,Nz);       Map.fill(-2);

        auto neighborList= new int[Nq*Npad];
        auto interpolationList= new int[Nq*Npad];
        auto scalarList= new int[Nq*Npad];
        Np = ScaLBL_Comm->MemoryOptimizedLayoutAA_LIBB(Map,neighborList,interpolationList,scalarList,Mask->id,Np,Nx,Nx*Ny);
        
        MPI_Barrier(comm);
        
        TmpMap=new int[Np];
        for (int k=0; k<Nz; k++){
            for (int j=0; j<Ny; j++){
                for (int i=0; i<Nx; i++){
                    int idx=Map(i,j,k);
                    if (!(idx < 0)) { TmpMap[idx] = k*Nx*Ny+j*Nx+i;
                    }
                }
            }
        }
        
        // check that TmpMap is valid
        for (int idx=0; idx<ScaLBL_Comm->LastExterior(); idx++){
            int n = TmpMap[idx];
            if (n > Nx*Ny*Nz){
                printf("Bad value! idx=%i \n",n);
                TmpMap[idx] = Nx*Ny*Nz-1;
            }
        }
        for (int idx=ScaLBL_Comm->FirstInterior(); idx<ScaLBL_Comm->LastInterior(); idx++){
            int n = TmpMap[idx];
            if (n > Nx*Ny*Nz){
                printf("Bad value! idx=%i \n",n);
                TmpMap[idx] = Nx*Ny*Nz-1;
            }
        }
        
        auto libbqA = new double[18*Np];
        auto libbqBC = new double[18*Np];
        auto libbqD = new double[18*Np];
        
        for (size_t i = 0; i < 18*Np; i++) { libbqA[i] = 0;  libbqBC[i] = 0;  libbqD[i] = 0; }
        
        for (int q = 0; q < 18; q++) {
            for (int n = 0; n < N; n++) {
                if (LIBB_id[n] == 0) {
                    qList[n + q*N] = 0.0;
                }
            }
        }
        for (int q = 0; q < 6; q++) {
            for (int n = 0; n < N; n++) {
                if (qList[n + q*N] == -1.0) {
                    qList[n + q*N] = 1.0;
                }
            }
        }
        for (int q = 6; q < 18; q++) {
            for (int n = 0; n < N; n++) {
                if (qList[n + q*N] == -1.0) {
                    qList[n + q*N] = sqrt(2);
                }
            }
        }
        for (int q = 6; q < 18; q++) {
            for (size_t n = 0; n < N; n++) {
                    qList[n + q*N] /= sqrt(2);
            }
        }
        
        
        MPI_Barrier(comm);
        
        count = LIBB_Init(TmpMap, LIBB_id, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(),
                          ScaLBL_Comm->LastInterior(), N, Np, qList, libbqA, libbqBC, libbqD);
        LIBB_Init(TmpMap, LIBB_id, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(),
                           N, Np, qList, libbqA, libbqBC, libbqD);
        
        
        
        // if (rank==0){
        //     printf("LIBB A \n");
        // }
        // for (int k=0;k<Nz;k++){
        //     for (int j=0;j<Ny;j++){
        //         for (int i=0;i<Nx;i++){
        //             int n=k*Nx*Ny+j*Nx+i;
        //             printf("%.2f ",libbqA[n]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n\n");
        // }

        // if (rank==0){
        //     printf("LIBB BC \n");
        // }
        // for (int k=0;k<Nz;k++){
        //     for (int j=0;j<Ny;j++){
        //         for (int i=0;i<Nx;i++){
        //             int n=k*Nx*Ny+j*Nx+i;
        //             printf("%.2f ",libbqBC[n]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n\n");
        // }

        // if (rank==0){
        //     printf("LIBB D \n");
        // }
        // for (int k=0;k<Nz;k++){
        //     for (int j=0;j<Ny;j++){
        //         for (int i=0;i<Nx;i++){
        //             int n=k*Nx*Ny+j*Nx+i;
        //             printf("%.2f ",libbqD[n]);
        //         }
        //         printf("\n");
        //     }
        //     printf("\n\n");
        // }

        //         if (rank==0){
        //     printf("ns x \n");
        // }
        // for (int k=0;k<Nz;k++){
        //     for (int j=0;j<Ny;j++){
        //         for (int i=0;i<Nx;i++){
        //             int n=k*Nx*Ny+j*Nx+i;
        //             printf("%.2f ",Averages->GradPhiX(n));
        //         }
        //         printf("\n");
        //     }
        //     printf("\n\n");
        // }

        // if (rank==0){
        //     printf("ns y \n");
        // }
        // for (int k=0;k<Nz;k++){
        //     for (int j=0;j<Ny;j++){
        //         for (int i=0;i<Nx;i++){
        //             int n=k*Nx*Ny+j*Nx+i;
        //             printf("%.2f ",Averages->GradPhiY(n));
        //         }
        //         printf("\n");
        //     }
        //     printf("\n\n");
        // }

        // if (rank==0){
        //     printf("ns z \n");
        // }
        // for (int k=0;k<Nz;k++){
        //     for (int j=0;j<Ny;j++){
        //         for (int i=0;i<Nx;i++){
        //             int n=k*Nx*Ny+j*Nx+i;
        //             printf("%.2f ",Averages->GradPhiZ(n));
        //         }
        //         printf("\n");
        //     }
        //     printf("\n\n");
        // }

        sprintf(LocalRankString,"%05d",rank);

        sprintf(LocalRankFilename,"%s%s","libbqA.",LocalRankString);
        FILE *AFILE = fopen(LocalRankFilename,"wb");
        if (AFILE==NULL) ERROR("Error opening file: libbqA.xxxxx");
        fwrite(libbqA,RT_DOUBLE,18*Np,AFILE);
        fclose(AFILE);

        sprintf(LocalRankFilename,"%s%s","libbqBC.",LocalRankString);
        FILE *ACOFILE = fopen(LocalRankFilename,"wb");
        if (ACOFILE==NULL) ERROR("Error opening file: libbqBC.xxxxx");
        fwrite(libbqBC,RT_DOUBLE,18*Np,ACOFILE);
        fclose(ACOFILE);

        sprintf(LocalRankFilename,"%s%s","libbqD.",LocalRankString);
        FILE *DFILE = fopen(LocalRankFilename,"wb");
        if (DFILE==NULL) ERROR("Error opening file: libbqD.xxxxx");
        fwrite(libbqD,RT_DOUBLE,18*Np,DFILE);
        fclose(DFILE);

        for (size_t n = 0; n < N; n++) TemporaryField[n] = Averages->VFmask(n);

        sprintf(LocalRankFilename,"%s%s","VFmask.",LocalRankString);
        FILE *VFFILE = fopen(LocalRankFilename,"wb");
        if (VFFILE==NULL) ERROR("Error opening file: VFmask.xxxxx");
        fwrite(TemporaryField,RT_DOUBLE,N,VFFILE);
        fclose(VFFILE);
    
        for (size_t n = 0; n < N; n++) Averages->ID(n) = LIBB_id[n];
        
        for (size_t n = 0; n < N; n++) TemporaryField[n] = Averages->GradPhiX(n);
        sprintf(LocalRankFilename,"%s%s","ns_x.",LocalRankString);
        FILE *NSX = fopen(LocalRankFilename,"wb");
        if (NSX==NULL) ERROR("Error opening file: NSX.xxxxx");
        fwrite(TemporaryField,RT_DOUBLE,N,NSX);
        fclose(NSX);
        
        for (size_t n = 0; n < N; n++) TemporaryField[n] = Averages->GradPhiY(n);
        sprintf(LocalRankFilename,"%s%s","ns_y.",LocalRankString);
        FILE *NSY = fopen(LocalRankFilename,"wb");
        if (NSY==NULL) ERROR("Error opening file: NSY.xxxxx");
        fwrite(TemporaryField,RT_DOUBLE,N,NSY);
        fclose(NSY);
        
        for (size_t n = 0; n < N; n++) TemporaryField[n] = Averages->GradPhiZ(n);
        sprintf(LocalRankFilename,"%s%s","ns_z.",LocalRankString);
        FILE *NSZ = fopen(LocalRankFilename,"wb");
        if (NSZ==NULL) ERROR("Error opening file: NSZ.xxxxx");
        fwrite(TemporaryField,RT_DOUBLE,N,NSZ);
        fclose(NSZ);
        
        delete[] TemporaryField;
        
        sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
        FILE *DOMAINFILE = fopen(LocalRankFilename,"wb");
        if (DOMAINFILE==NULL) ERROR("Error opening file: DOMAINFILE.xxxxx");
        fwrite(LIBB_id,RT_INT,N,DOMAINFILE);
        fclose(DOMAINFILE);
        
      
        
        /* Clean up memory */
        delete[] LIBB_id;
        delete[] qList;
        delete[] rtDistance;
        
        if (rank == 0) cout << "End of ray tracing pre-processor. " << endl;
    }
    // ****************************************************
    MPI_Barrier(comm);
    MPI_Finalize();
    // ****************************************************
}
