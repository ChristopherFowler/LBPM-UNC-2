/*
 * Pre-processor to parallelize outputs from Blender for use in color simulator
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <exception>      // std::exception
#include <stdexcept>
#include <string>

#include "common/Domain.h"
#include "common/Array.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "common/ScaLBL.h"
#include "analysis/TwoPhase.h"
#include "analysis/runAnalysis.h"



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
                if (id[n] == 2 || id[n] == 1)  count++;
            }
        }
    }
    
    // total Global is the number of nodes in the pore-space
    MPI_Allreduce(&count,&totalGlobal,1,MPI_INT,MPI_SUM,Dm.Comm);
    double volume=double(nprocx*nprocy*nprocz)*double(Nx-2)*double(Ny-2)*double(Nz-2);
    double porosity=double(totalGlobal)/volume;
     if (rank == 0) printf("porosity w/o  halo= %.12f\n",porosity);
    
    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    double dcount = 0;
    double dcountGlobal = 0;
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                for (int p=0;p<8;p++){
                    int n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
                    if (id[n] == 2 || id[n] == 1)  {  dcount += 0.125; }
                }
            }
        }
    }
    MPI_Allreduce(&dcount,&dcountGlobal,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    porosity=double(dcountGlobal)/volume;
     if (rank == 0) printf("porosity with halo= %.12f\n",porosity);
 
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
    
    
    
} // Sets interior to false, keeps halo true

inline void PackID(int *list,
                   int count, char *sendbuf, char *ID) {
    // Fill in the phase ID values from neighboring processors
    // This packs up the values that need to be sent from one processor to another
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        sendbuf[idx] = ID[n];
    }
}

inline void UnpackID(int *list,
                     int count, char *recvbuf, char *ID) {
    // Fill in the phase ID values from neighboring processors
    // This unpacks the values once they have been recieved from neighbors
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        ID[n] = recvbuf[idx];
    }
}

inline void PackID_Double(int *list,
                          int count, double *sendbuf, double *ID) {
    // Fill in the phase ID values from neighboring processors
    // This packs up the values that need to be sent from one processor to another
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        sendbuf[idx] = ID[n];
    }
}

inline void UnpackID_Double(int *list,
                            int count, double *recvbuf, double *ID) {
    // Fill in the phase ID values from neighboring processors
    // This unpacks the values once they have been recieved from neighbors
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        ID[n] = recvbuf[idx];
    }
}

inline void PrintDoubleField(std::string str, double * FIELD, Domain &Dm, int rank, int Nx, int Ny, int Nz) {
    
    if (rank == 0) {
        std::cout << str << std::endl;
        double DVALUE = 0;
        printf("rank=%d\n",rank);
        for (int j=5;j<6;j++){
            for (int k=0;k<Nz;k++){
                for (int i=0;i<Nx;i++){
                    int n=k*Nx*Ny+j*Nx+i;
                    DVALUE = FIELD[n];
                    printf("%.7f ",DVALUE);
                }
                printf("\n");
            }
            printf("\n\n");
        }
    }
//    MPI_Barrier(Dm.Comm);
    if (rank == 1) {
        double DVALUE = 0;
        printf("rank=%d\n",rank);
        for (int j=5;j<6;j++){
            for (int k=0;k<Nz;k++){
                for (int i=0;i<Nx;i++){
                    int n=k*Nx*Ny+j*Nx+i;
                    DVALUE = FIELD[n];
                    printf("%.7f ",DVALUE);
                }
                printf("\n");
            }
            printf("\n\n");
        }
    }
//    MPI_Barrier(Dm.Comm);
    if (rank == 4) {
        double VALUE = 0;
        printf("rank=%d\n",rank);
        for (int j=5;j<6;j++){
            for (int k=0;k<Nz;k++){
                for (int i=0;i<Nx;i++){
                    int n=k*Nx*Ny+j*Nx+i;
                    VALUE = FIELD[n];
                    printf("%.7f ",VALUE);
                }
                printf("\n");
            }
            printf("\n\n");
        }
    }
//    MPI_Barrier(Dm.Comm);
    if (rank == 5) {
        double VALUE = 0;
        printf("rank=%d\n",rank);
        for (int j=5;j<6;j++){
            for (int k=0;k<Nz;k++){
                for (int i=0;i<Nx;i++){
                    int n=k*Nx*Ny+j*Nx+i;
                    VALUE = FIELD[n];
                    printf("%.7f ",VALUE);
                }
                printf("\n");
            }
            printf("\n\n");
        }
    }
//    MPI_Barrier(Dm.Comm);
}

inline void PrintIntegerField(std::string str, char * FIELD, Domain &Dm, int rank, int Nx, int Ny, int Nz)
    {
    if (rank == 0) {
        std::cout << str << std::endl;
        int VALUE = 0;
        printf("rank=%d\n",rank);
        for (int j=5;j<6;j++){
            for (int k=0;k<Nz;k++){
                for (int i=0;i<Nx;i++){
                    int n=k*Nx*Ny+j*Nx+i;
                    VALUE = FIELD[n];
                    printf("%d ",VALUE);
                }
                printf("\n");
            }
            printf("\n\n");
        }
    }
//    MPI_Barrier(Dm.Comm);
    if (rank == 1) {
        int VALUE = 0;
        printf("rank=%d\n",rank);
        for (int j=5;j<6;j++){
            for (int k=0;k<Nz;k++){
                for (int i=0;i<Nx;i++){
                    int n=k*Nx*Ny+j*Nx+i;
                    VALUE = FIELD[n];
                    printf("%d ",VALUE);
                }
                printf("\n");
            }
            printf("\n\n");
        }
    }
//    MPI_Barrier(Dm.Comm);
    if (rank == 4) {
        int VALUE = 0;
        printf("rank=%d\n",rank);
        for (int j=5;j<6;j++){
            for (int k=0;k<Nz;k++){
                for (int i=0;i<Nx;i++){
                    int n=k*Nx*Ny+j*Nx+i;
                    VALUE = FIELD[n];
                    printf("%d ",VALUE);
                }
                printf("\n");
            }
            printf("\n\n");
        }
    }
//    MPI_Barrier(Dm.Comm);
    if (rank == 5) {
        int VALUE = 0;
        printf("rank=%d\n",rank);
        for (int j=5;j<6;j++){
            for (int k=0;k<Nz;k++){
                for (int i=0;i<Nx;i++){
                    int n=k*Nx*Ny+j*Nx+i;
                    VALUE = FIELD[n];
                    printf("%d ",VALUE);
                }
                printf("\n");
            }
            printf("\n\n");
        }
    }
//    MPI_Barrier(Dm.Comm);
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
    int rank, nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);

    {
        //.......................................................................
        // Reading the domain information file
        //.......................................................................

        string filename;
        if (argc > 1){
            filename=argv[1];
        }
        else ERROR("No input database provided\n");

        // read the input database
        auto db = std::make_shared<Database>( filename );
        auto domain_db = db->getDatabase( "Domain" );
        auto analysis_db = db->getDatabase( "Analysis" );
        auto rt_db = db->getDatabase( "RayTrace" );
        
        // Read domain parameters
        auto size = domain_db->getVector<int>( "n" );
        auto size_global = domain_db->getVector<int>( "N" );
        auto nproc = domain_db->getVector<int>( "nproc" );

        string READFILE_FluidPhaseID = "FluidPhaseID.raw";
        string READFILE_ID = "ID.raw";
        string READFILE_SD = "SD.raw";
        string READFILE_SDMC = "SDMC.raw";
        string READFILE_LIBBA = "libbA.raw";
        string READFILE_LIBBBC = "libbBC.raw";
        string READFILE_LIBBD = "libbD.raw";

        char LocalRankFilename[40];
        char LocalRankString[8];

        std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
        std::shared_ptr<Domain> Mask (new Domain(domain_db,comm));

        int n, nx, ny, nz;
        int nprocx, nprocy, nprocz;
        int global_Nx, global_Ny, global_Nz;

        int RANK = Dm->rank();
        nx = size[0];
        ny = size[1];
        nz = size[2];
        nprocx = nproc[0];
        nprocy = nproc[1];
        nprocz = nproc[2];

        global_Nx = size_global[0];
        global_Ny = size_global[1];
        global_Nz = size_global[2];

        int N = (nx+2)*(ny+2)*(nz+2);
        int N_global = global_Nx*global_Ny*global_Nz;
        int nprocs=nprocx*nprocy*nprocz;

        int geometry = rt_db->getScalar<int>( "geometry" );
        int rmin = 1;
        int rmax = nz-1;
        // if (BoundaryCondition > 0){
        //     if ( domain_db->keyExists( "ReservoirMin" ) && kproc == 0){
        //         rmin = domain_db->getScalar<int>( "ReservoirMin" );
        //     }
        //     if ( domain_db->keyExists( "ReservoirMax" ) && kproc == (nprocz - 1)){
        //         rmax = domain_db->getScalar<int>( "ReservoirMax" );
        //         rmax = Nz-1-((Nz-2)*nprocz - rmax);
        //     }
        // }

        int SIZE = global_Nx*global_Ny*global_Nz;

        std::shared_ptr<TwoPhase> Averages;

        char *FluidPhaseSegData = new char[SIZE];
        char *SegData = new char[SIZE];
        double *SignDist = new double[SIZE];
        double *SignDistMC = new double[SIZE];
        double *LibbqA = new double[18*SIZE];
        double *LibbqBC = new double[18*SIZE];
        double *LibbqD = new double[18*SIZE];

        int VALUE = 0;
        double DVALUE = 0;
        
        if (RANK==0){

            //Rank 0 reads the data and distributes to worker processes
            printf("Dimensions: %i x %i x %i \n",global_Nx,global_Ny,global_Nz);

            //Phase
            printf("Read data from %s \n",READFILE_ID.c_str());
            FILE *SEGDAT = fopen(READFILE_ID.c_str(),"rb");
            if (SEGDAT==NULL) ERROR("Error reading id data, null \n");
            int ReadSeg;
            ReadSeg=fread(SegData,1,SIZE,SEGDAT);
            if (ReadSeg != SIZE) printf("Error reading id data, size \n");
            fclose(SEGDAT);
            
            

       //     PrintIntegerField("SegData", SegData, *Dm, rank, global_Nx, global_Ny, global_Nz);
            
            //FluidPhase
            printf("Read data from %s \n",READFILE_FluidPhaseID.c_str());
            FILE *FLUIDPHASESEGDAT = fopen(READFILE_FluidPhaseID.c_str(),"rb");
            if (FLUIDPHASESEGDAT==NULL) ERROR("Error reading fluid phase id data, null \n");
            int ReadFluidPhaseSeg;
            ReadFluidPhaseSeg=fread(FluidPhaseSegData,1,SIZE,FLUIDPHASESEGDAT);
            if (ReadFluidPhaseSeg != SIZE) printf("Error reading fluid phase id data, size \n");
            fclose(FLUIDPHASESEGDAT);
            
        //    PrintIntegerField("FluidPhaseSegData", FluidPhaseSegData, *Dm, rank, global_Nx, global_Ny, global_Nz);
            
//            printf("SegData read:\n");
//            for (int i=5;i<6;i++){
//                for (int j=0;j<global_Ny;j++){
//                    for (int k=0;k<global_Nz;k++){
//                        int n=k*global_Nx*global_Ny+j*global_Nx+i;
//                        VALUE = SegData[n];
//                        printf("%d ",VALUE);
//                    }
//                    printf("\n");
//                }
//                printf("\n\n");
//            }
            
            
            //SignDist
            printf("Read data from %s \n",READFILE_SD.c_str());
            FILE *RayTraceDistance = fopen(READFILE_SD.c_str(),"rb");
            if (RayTraceDistance==NULL) ERROR("Error reading sd data, null \n");
            int ReadSD;
            ReadSD=fread(SignDist,8,SIZE,RayTraceDistance);
            if (ReadSD != SIZE) printf("Error reading sd data, size \n");
            fclose(RayTraceDistance);
            
      //      PrintDoubleField("SignDist", SignDist, *Dm, rank, global_Nx, global_Ny, global_Nz);
            
//            printf("SignDist read:\n");
//            for (int i=5;i<6;i++){
//                for (int j=0;j<global_Ny;j++){
//                    for (int k=0;k<global_Nz;k++){
//                        int n=k*global_Nx*global_Ny+j*global_Nx+i;
//                        DVALUE = SignDist[n];
//                        printf("%.2f ",DVALUE);
//                    }
//                    printf("\n");
//                }
//                printf("\n\n");
//            }

            //SignDistMC
            printf("Read data from %s \n",READFILE_SDMC.c_str());
            FILE *RayTraceDistanceMC = fopen(READFILE_SDMC.c_str(),"rb");
            if (RayTraceDistanceMC==NULL) ERROR("Error reading sd data, null \n");
            int ReadSDMC;
            ReadSDMC=fread(SignDistMC,8,SIZE,RayTraceDistanceMC);
            if (ReadSDMC != SIZE) printf("Error reading sd data, size \n");
            fclose(RayTraceDistanceMC);
            
        //    PrintDoubleField("SignDistMC", SignDistMC, *Dm, rank, global_Nx, global_Ny, global_Nz);

            //LibbA
            printf("Read data from %s \n",READFILE_LIBBA.c_str());
            FILE *AFILE = fopen(READFILE_LIBBA.c_str(),"rb");
            if (AFILE==NULL) ERROR("Error reading libba data, null \n");
            int read_A;
            read_A = fread(LibbqA,8,18*SIZE,AFILE);
            if (read_A != int(18*SIZE)) printf("Error reading libba data, size \n");
            fclose(AFILE);

            //LibbBC
            printf("Read data from %s \n",READFILE_LIBBBC.c_str());
            FILE *BCFILE = fopen(READFILE_LIBBBC.c_str(),"rb");
            if (BCFILE==NULL) ERROR("Error reading libbbc data, null \n");
            int read_bc;
            read_bc = fread(LibbqBC,8,18*SIZE,BCFILE);
            if (read_bc != int(18*SIZE)) printf("Error reading libbbc data, size \n");
            fclose(BCFILE);

            //LibbD
            printf("Read data from %s \n",READFILE_LIBBD.c_str());
            FILE *DFILE = fopen(READFILE_LIBBD.c_str(),"rb");
            if (DFILE==NULL) ERROR("Error reading libbd data, null \n");
            int read_d;
            read_d = fread(LibbqD,8,18*SIZE,DFILE);
            if (read_d != int(18*SIZE)) printf("Error reading libbd data, size \n");
            fclose(DFILE);

        }


       
        
        char * loc_fluidphaseid;
        loc_fluidphaseid = new char [(nx+2)*(ny+2)*(nz+2)];
        char * fluidphaseid;
        fluidphaseid = new char [(nx+2)*(ny+2)*(nz+2)];
    
        char * loc_id;
        loc_id = new char [(nx+2)*(ny+2)*(nz+2)];
        char * id;
        id = new char [(nx+2)*(ny+2)*(nz+2)];
        for (int n=0; n<N; n++) id[n] = 2;
        for (int n=0; n<N; n++) Dm->id[n] = 1;
        Dm->CommInit();
        for (int n=0; n<N; n++) Dm->id[n] = id[n]; // Set wetting phase (id=2)
        
        double * loc_sd;
        loc_sd = new double [(nx+2)*(ny+2)*(nz+2)];
        double * sd;
        sd = new double [(nx+2)*(ny+2)*(nz+2)];

        double * loc_sdmc;
        loc_sdmc = new double [(nx+2)*(ny+2)*(nz+2)];
        double * sdmc;
        sdmc = new double [(nx+2)*(ny+2)*(nz+2)];

        double * loc_libba;
        loc_libba = new double [(nx+2)*(ny+2)*(nz+2)*18];
        double * libba;
        libba = new double [(nx+2)*(ny+2)*(nz+2)*18];

        double * loc_libbbc;
        loc_libbbc = new double [(nx+2)*(ny+2)*(nz+2)*18];
        double * libbbc;
        libbbc = new double [(nx+2)*(ny+2)*(nz+2)*18];

        double * loc_libbd;
        loc_libbd = new double [(nx+2)*(ny+2)*(nz+2)*18];
        double * libbd;
        libbd = new double [(nx+2)*(ny+2)*(nz+2)*18];

        // Set up the sub-domains
        if (RANK==0){
            printf("Distributing subdomains across %i processors \n",nprocs);
            printf("Process grid: %i x %i x %i \n",nprocx,nprocy,nprocz);
            printf("Subdomain size: %i x %i x %i \n",nx,ny,nz);

            for (int kp=0; kp<nprocz; kp++){
                for (int jp=0; jp<nprocy; jp++){
                    for (int ip=0; ip<nprocx; ip++){
                        // rank of the process that gets this subdomain
                        int rnk = kp*nprocx*nprocy + jp*nprocx + ip;
                        // Pack and send the subdomain for rnk
                        for (int k=0;k<nz+2;k++){
                            for (int j=0;j<ny+2;j++){
                                for (int i=0;i<nx+2;i++){
                                    int x = ip*nx + i-1;
                                    int y = jp*ny + j-1;
                                    int z = kp*nz + k-1;
                                    if (x<0)     x=0;
                                    if (!(x<global_Nx))    x=global_Nx-1;
                                    if (y<0)     y=0;
                                    if (!(y<global_Ny))    y=global_Ny-1;
                                    if (z<0)     z=0;
                                    if (!(z<global_Nz))    z=global_Nz-1;
                                    int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
                                    int nglobal = z*global_Nx*global_Ny+y*global_Nx+x;
                                    loc_fluidphaseid[nlocal] = FluidPhaseSegData[nglobal];
                                    loc_id[nlocal] = SegData[nglobal];
                                    loc_sd[nlocal] = SignDist[nglobal];
                                    loc_sdmc[nlocal] = SignDistMC[nglobal];
                                }
                            }
                        }
                        for (int q=0;q<18;q++){
                            for (int k=0;k<nz+2;k++){
                                for (int j=0;j<ny+2;j++){
                                    for (int i=0;i<nx+2;i++){
                                        int x = ip*nx + i-1;
                                        int y = jp*ny + j-1;
                                        int z = kp*nz + k-1;
                                        if (x<0)     x=0;
                                        if (!(x<global_Nx))    x=global_Nx-1;
                                        if (y<0)     y=0;
                                        if (!(y<global_Ny))    y=global_Ny-1;
                                        if (z<0)     z=0;
                                        if (!(z<global_Nz))    z=global_Nz-1;
                                        int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
                                        int nglobal = z*global_Nx*global_Ny+y*global_Nx+x;
                                        nlocal = nlocal + q*N;
                                        nglobal = nglobal + q*N_global;
                                        loc_libba[nlocal]=LibbqA[nglobal];
                                        loc_libbbc[nlocal]=LibbqBC[nglobal];
                                        loc_libbd[nlocal]=LibbqD[nglobal];
                                    }
                                }
                            }
                        }
                        if (rnk==0){
                            for (int k=0;k<nz+2;k++){
                                for (int j=0;j<ny+2;j++){
                                    for (int i=0;i<nx+2;i++){
                                        int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
                                        fluidphaseid[nlocal] = loc_fluidphaseid[nlocal];
                                        id[nlocal] = loc_id[nlocal];
                                        sd[nlocal] = loc_sd[nlocal];
                                        sdmc[nlocal] = loc_sdmc[nlocal];
                                    }
                                }
                            }
                            for (int q=0;q<18;q++){
                                for (int k=0;k<nz+2;k++){
                                    for (int j=0;j<ny+2;j++){
                                        for (int i=0;i<nx+2;i++){
                                            int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
                                            nlocal = nlocal + q*N;
                                            libba[nlocal] = loc_libba[nlocal];
                                            libbbc[nlocal] = loc_libbbc[nlocal];
                                            libbd[nlocal] = loc_libbd[nlocal];
                                        }
                                    }
                                }
                            }
                        }
                        else{
                            //Send data to other ranks
                            MPI_Send(loc_fluidphaseid,N,MPI_CHAR,rnk,15,comm);
                            MPI_Send(loc_id,N,MPI_CHAR,rnk,15,comm);
                            MPI_Send(loc_sd,N,MPI_DOUBLE,rnk,15,comm);
                            MPI_Send(loc_sdmc,N,MPI_DOUBLE,rnk,15,comm);
                            MPI_Send(loc_libba,18*N,MPI_DOUBLE,rnk,15,comm);
                            MPI_Send(loc_libbbc,18*N,MPI_DOUBLE,rnk,15,comm);
                            MPI_Send(loc_libbd,18*N,MPI_DOUBLE,rnk,15,comm);
                        }

                    }
                }
            }

        }
        else{
            // Receive the subdomain from rank 0
            MPI_Recv(fluidphaseid,N,MPI_CHAR,0,15,comm,MPI_STATUS_IGNORE);
            MPI_Recv(id,N,MPI_CHAR,0,15,comm,MPI_STATUS_IGNORE);
            MPI_Recv(sd,N,MPI_DOUBLE,0,15,comm,MPI_STATUS_IGNORE);
            MPI_Recv(sdmc,N,MPI_DOUBLE,0,15,comm,MPI_STATUS_IGNORE);
            MPI_Recv(libba,18*N,MPI_DOUBLE,0,15,comm,MPI_STATUS_IGNORE);
            MPI_Recv(libbbc,18*N,MPI_DOUBLE,0,15,comm,MPI_STATUS_IGNORE);
            MPI_Recv(libbd,18*N,MPI_DOUBLE,0,15,comm,MPI_STATUS_IGNORE);
        }
        
        auto boolfield = new bool[N];  for (n=0; n<N; n++) boolfield[n]=true; fillBoolfield(boolfield, nx+2, ny+2, nz+2);  for (int n = 0; n < N; n++) if (boolfield[n] == true) id[n] = 0;
        std::array<size_t,3> subdivide = { 1, 1, 1 };
        if ( analysis_db->keyExists( "subdivide" ) ) {
            auto tmp = analysis_db->getVector<size_t>( "subdivide" );
            subdivide = { tmp[0], tmp[1], tmp[2] };
        }
        Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm,subdivide) );
        PopulateHalo_Char(id, *Dm, nx+2, ny+2, nz+2, RANK);
        PopulateHalo_Double(sdmc, *Dm, nx+2, ny+2, nz+2, RANK);
        PopulateHalo_Double(sd, *Dm, nx+2, ny+2, nz+2, RANK);
        
       // PrintIntegerField("id", id, *Dm, rank, nx+2, ny+2, nz+2);

        for (size_t n=0; n<N; n++) Averages->SDs(n) = sdmc[n];
        
        Averages->ComputeVolumeFraction(Averages->SDs, Averages->SDs_x, Averages->SDs_y, Averages->SDs_z, Averages->VFmask, true, rmin, rmax, geometry);
        
        
        Dm->CommunicateMeshHalo(Averages->VFmask);
        Dm->CommunicateMeshHalo(Averages->SDs);
        Dm->CommunicateMeshHalo(Averages->SDs_x);
        Dm->CommunicateMeshHalo(Averages->SDs_y);
        Dm->CommunicateMeshHalo(Averages->SDs_z);

        for (int i = 0; i < N; i++) {

            //Change ID based on VFmask
             if (Averages->VFmask(i) > 0.5){
                 id[i] = 0;
             } else{
//                 if (id[i]id[i] = 2;
             }
        }
     //   PrintIntegerField("id", id, *Dm, rank, nx+2, ny+2, nz+2);
        id[0] = id[(nx+2)-1] = id[((ny+2)-1)*(nx+2)] = id[((ny+2)-1)*(nx+2) + (nx+2)-1] = 0;
        id[((nz+2)-1)*(nx+2)*(ny+2)] = id[((nz+2)-1)*(nx+2)*(ny+2)+(nx+2)-1] = id[((nz+2)-1)*(nx+2)*(ny+2)+((ny+2)-1)*(nx+2)] = id[((nz+2)-1)*(nx+2)*(ny+2)+((ny+2)-1)*(nx+2) + (nx+2)-1] = 0;
//        printf("line=%d\n",__LINE__);

        for (int n = 0; n < N; n++) { Averages->SDs(n) = sd[n]; }

//
            double porosity = calculate_porosity(id, nx+2, ny+2, nz+2, *Dm, RANK);


//        int q = 0;
//        for (int i = 0; i < N; i++) {
//            Mask->id[i] = id[i];
//            Averages->ID(i) = id[i];
//            q = 4;
//            Averages->Vel_x(i) = libba[i+q*N];
//            Averages->Vel_y(i) = libbbc[i+q*N];
//            Averages->Vel_z(i) = libbd[i+q*N];
//
//            // q = 5;
//            // Averages->Vel_z(i) = libba[i+q*N];
//            // Averages->Vel_z(i) = libbbc[i+q*N];
//            // Averages->Vel_z(i) = libbd[i+q*N];
//        }

    

        for (int n=0; n<N; n++) Mask->id[n] = 1;
        Mask->CommInit();
        
        for (int i=0; i<N; i++) Mask->id[i] = id[i];
        size_t Np=Mask->PoreCount();
        //Initialize visualization/make compressed data map
        int Nq = 18;
       
        //...........................................................................
        auto ScaLBL_Comm  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));
        
        size_t Npad=(Np/16 + 2)*16;
        IntArray Map;
        Map.resize(nx+2,ny+2,nz+2);       Map.fill(-2);

        auto neighborList= new int[Nq*Npad];
        auto interpolationList= new int[Nq*Npad];
        auto scalarList= new int[Nq*Npad];
        Np = ScaLBL_Comm->MemoryOptimizedLayoutAA_LIBB(Map,neighborList,interpolationList,scalarList,Mask->id,Np,nx+2,(nx+2)*(ny+2));

        const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);

        bool Regular = false;
        runAnalysis analysis( analysis_db, rank_info, ScaLBL_Comm, Dm, Np, Regular, 0.9, Map );

       
        
        //compress LIBB
        double * comp_libbqA = new double[18*Np];
        double * comp_libbqBC = new double[18*Np];
        double * comp_libbqD = new double[18*Np];
        int ijk;
        
        int * TmpMap; TmpMap = new int[Np];
        for (int k=0; k<nz+2; k++){
            for (int j=0; j<ny+2; j++){
                for (int i=0; i<nx+2; i++){
                    int idx=Map(i,j,k);
                    if (!(idx < 0)) { TmpMap[idx] = k*(nx+2)*(ny+2)+j*(nx+2)+i;
                    }
                }
            }
        }
        
//        printf("line=%d\n",__LINE__);
        // check that TmpMap is valid
        for (int idx=0; idx<ScaLBL_Comm->LastExterior(); idx++){
            int n = TmpMap[idx];
            if (n > (nx+2)*(ny+2)*(nz+2)){
                printf("Bad value! idx=%i \n",n);
                TmpMap[idx] = (nx+2)*(ny+2)*(nz+2)-1;
            }
        }
        for (int idx=ScaLBL_Comm->FirstInterior(); idx<ScaLBL_Comm->LastInterior(); idx++){
            int n = TmpMap[idx];
            if (n > (nx+2)*(ny+2)*(nz+2)){
                printf("Bad value! idx=%i \n",n);
                TmpMap[idx] = (nx+2)*(ny+2)*(nz+2)-1;
            }
        }
        
        
        for (int i=0; i<ScaLBL_Comm->LastExterior(); i++){
            ijk = TmpMap[i]; //index in regular layout
            for (int q=0;q<18;q++){
                comp_libbqA[i+q*Np] = libba[ijk+q*N];
                comp_libbqBC[i+q*Np] = libbbc[ijk+q*N];
                comp_libbqD[i+q*Np] = libbd[ijk+q*N];
            }
        }
        for (int i=ScaLBL_Comm->FirstInterior(); i<ScaLBL_Comm->LastInterior(); i++){
            ijk = TmpMap[i]; //index in regular layout
            for (int q=0;q<18;q++){
                comp_libbqA[i+q*Np] = libba[ijk+q*N];
                comp_libbqBC[i+q*Np] = libbbc[ijk+q*N];
                comp_libbqD[i+q*Np] = libbd[ijk+q*N];
            }
        }
        
//        printf("line=%d\n",__LINE__);


        delete[] TmpMap;
     
        
        //Produce visualization
//        analysis.run6(0, *Averages, Np, comp_libbqA);
//        analysis.finish();

        // Write the data
        double * TemporaryField = new double[N];
        for (int n = 0; n < N; n++) TemporaryField[n] = 0.0;

        sprintf(LocalRankFilename,"FluidPhaseID.%05i",RANK);
        FILE *FLUIDPHASEID = fopen(LocalRankFilename,"wb");
        fwrite(fluidphaseid,1,(nx+2)*(ny+2)*(nz+2),FLUIDPHASEID);
        fclose(FLUIDPHASEID);
        
//            printf("FluidPhaseID:\n");
//            for (int k=18;k<19;k++){
//                for (int i=1;i<nx+1;i++){
//                    for (int j=1;j<ny+1;j++){
//
//                            int n=k*(nx+2)*(ny+2)+j*(nx+2)+i;
//                            VALUE = fluidphaseid[n];
//                            printf("%d ",VALUE);
//                        }
//                        printf("\n");
//                    }
//                    printf("\n\n");
//                }
        
        sprintf(LocalRankFilename,"ID.%05i",RANK);
        FILE *ID = fopen(LocalRankFilename,"wb");
        fwrite(id,1,(nx+2)*(ny+2)*(nz+2),ID);
        fclose(ID);
        
        // ID
        if (nprocs > 1)
        {
            sprintf(LocalRankFilename,"ID_p.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(id,1,N,MYFILE);
            fclose(MYFILE);
        }
        else {
            sprintf(LocalRankFilename,"ID_s.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(id,1,N,MYFILE);
            fclose(MYFILE);
        }
        
        
        sprintf(LocalRankFilename,"SignDist.%05i",RANK);
        FILE *SD = fopen(LocalRankFilename,"wb");
        fwrite(sd,8,(nx+2)*(ny+2)*(nz+2),SD);
        fclose(SD);
        
        // SignDist
        if (nprocs > 1)
        {
            sprintf(LocalRankFilename,"SignDist_p.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(sd,8,N,MYFILE);
            fclose(MYFILE);
        }
        else {
            sprintf(LocalRankFilename,"SignDist_s.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(sd,8,N,MYFILE);
            fclose(MYFILE);
        }
        
//        printf("sd:\n");
//        for (int i=5;i<6;i++){
//            for (int j=1;j<nx+1;j++){
//                for (int k=1;k<nz+1;k++){
//                    int n=k*(nx+2)*(ny+2)+j*(nx+2)+i;
//                    DVALUE = sd[n];
//                    printf("%.2f ",DVALUE);
//                }
//                printf("\n");
//            }
//            printf("\n\n");
//        }
        
        
        
        
        

        sprintf(LocalRankString,"%05d",RANK);

        sprintf(LocalRankFilename,"%s%s","libbqA.",LocalRankString);
        FILE *AFILE = fopen(LocalRankFilename,"wb");
        if (AFILE==NULL) ERROR("Error opening file: libbqA.xxxxx");
        fwrite(comp_libbqA,8,18*Np,AFILE);
        fclose(AFILE);

        sprintf(LocalRankFilename,"%s%s","libbqBC.",LocalRankString);
        FILE *BCFILE = fopen(LocalRankFilename,"wb");
        if (BCFILE==NULL) ERROR("Error opening file: libbqBC.xxxxx");
        fwrite(comp_libbqBC,8,18*Np,BCFILE);
        fclose(BCFILE);

        sprintf(LocalRankFilename,"%s%s","libbqD.",LocalRankString);
        FILE *DFILE = fopen(LocalRankFilename,"wb");
        if (DFILE==NULL) ERROR("Error opening file: libbqD.xxxxx");
        fwrite(comp_libbqD,8,18*Np,DFILE);
        fclose(DFILE);

        for (int n = 0; n < N; n++) TemporaryField[n] = Averages->VFmask(n);
        //print tempfield slice
        sprintf(LocalRankFilename,"%s%s","VFmask.",LocalRankString);
        FILE *VFCFILE = fopen(LocalRankFilename,"wb");
        if (VFCFILE==NULL) ERROR("Error opening file: VFmask.xxxxx");
        fwrite(TemporaryField,8,(nx+2)*(ny+2)*(nz+2),VFCFILE);
        fclose(VFCFILE);
        
        // VFmask
        if (nprocs > 1)
        {
            sprintf(LocalRankFilename,"VFmask_p.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(TemporaryField,8,N,MYFILE);
            fclose(MYFILE);
        }
        else {
            sprintf(LocalRankFilename,"VFmask_s.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(TemporaryField,8,N,MYFILE);
            fclose(MYFILE);
        }
        
       
        for (int n = 0; n < N; n++) TemporaryField[n] = Averages->GradPhiX(n);
        sprintf(LocalRankFilename,"%s%s","ns_x.",LocalRankString);
        FILE *NSX = fopen(LocalRankFilename,"wb");
        if (NSX==NULL) ERROR("Error opening file: NSX.xxxxx");
        fwrite(TemporaryField,8,(nx+2)*(ny+2)*(nz+2),NSX);
        fclose(NSX);

        // asdf
        // NormalToSolid_X
        if (nprocs > 1)
        {
            sprintf(LocalRankFilename,"NormalToSolid_X_p.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(TemporaryField,8,N,MYFILE);
            fclose(MYFILE);
        }
        else {
            sprintf(LocalRankFilename,"NormalToSolid_X_s.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(TemporaryField,8,N,MYFILE);
            fclose(MYFILE);
        }
        for (int n = 0; n < N; n++) TemporaryField[n] = Averages->GradPhiY(n);
        sprintf(LocalRankFilename,"%s%s","ns_y.",LocalRankString);
        FILE *NSY = fopen(LocalRankFilename,"wb");
        if (NSY==NULL) ERROR("Error opening file: NSY.xxxxx");
        fwrite(TemporaryField,8,(nx+2)*(ny+2)*(nz+2),NSY);
        fclose(NSY);

        // NormalToSolid_Y
        if (nprocs > 1)
        {
            sprintf(LocalRankFilename,"NormalToSolid_Y_p.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(TemporaryField,8,N,MYFILE);
            fclose(MYFILE);
        }
        else {
            sprintf(LocalRankFilename,"NormalToSolid_Y_s.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(TemporaryField,8,N,MYFILE);
            fclose(MYFILE);
        }
        for (int n = 0; n < N; n++) TemporaryField[n] = Averages->GradPhiZ(n);
        sprintf(LocalRankFilename,"%s%s","ns_z.",LocalRankString);
        FILE *NSZ = fopen(LocalRankFilename,"wb");
        if (NSZ==NULL) ERROR("Error opening file: NSZ.xxxxx");
        fwrite(TemporaryField,8,(nx+2)*(ny+2)*(nz+2),NSZ);
        fclose(NSZ);
        // NormalToSolid_Z
        if (nprocs > 1)
        {
            sprintf(LocalRankFilename,"NormalToSolid_Z_p.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(TemporaryField,8,N,MYFILE);
            fclose(MYFILE);
        }
        else {
            sprintf(LocalRankFilename,"NormalToSolid_Z_s.%05i",rank);
            FILE *MYFILE = fopen(LocalRankFilename,"wb");
            fwrite(TemporaryField,8,N,MYFILE);
            fclose(MYFILE);
        }
        
        delete[] TemporaryField;

        MPI_Barrier(comm);

    }

    MPI_Barrier(comm);
    MPI_Finalize();
    return 0;

}
