/*
 * Pre-processor to generate signed distance function from segmented data
 * segmented data should be stored in a raw binary file as 1-byte integer (type char)
 * will output distance functions for phases
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>


//#include "common/pmmc.h"
#include "common/Domain.h"
#include "common/SpherePack.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "analysis/distance.h"
#include "analysis/analysis.h"
#include "analysis/runAnalysis.h"
#include "analysis/morphology.h"

inline void PackReadID(int *list, int count, char *sendbuf, char *ID){
    // Fill in the phase ID values from neighboring processors
    // This packs up the values that need to be sent from one processor to another
    int idx,n;
    
    for (idx=0; idx<count; idx++){
        n = list[idx];
        sendbuf[idx] = ID[n];
    }
}
//***************************************************************************************

inline void UnpackReadID(int *list, int count, char *recvbuf, char *ID){
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

inline void PopulateHalo_Char(char *id,
                              Domain &Dm, int nx, int ny, int nz, int rank) {
    double GlobalNumber = 1.f;
    double LocalNumber=0.f;
    int Nx = nx;
    int Ny = ny;
    int Nz = nz;
    int i,j,k,n;
    double count,countGlobal,totalGlobal;
    count = 0.f;
    double maxdist=0.f;
    double maxdistGlobal;
    int nprocx=Dm.nprocx();
    int nprocy=Dm.nprocy();
    int nprocz=Dm.nprocz();
    
    
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
    PackReadID(Dm.sendList_x, Dm.sendCount_x ,sendID_x, id);
    PackReadID(Dm.sendList_X, Dm.sendCount_X ,sendID_X, id);
    PackReadID(Dm.sendList_y, Dm.sendCount_y ,sendID_y, id);
    PackReadID(Dm.sendList_Y, Dm.sendCount_Y ,sendID_Y, id);
    PackReadID(Dm.sendList_z, Dm.sendCount_z ,sendID_z, id);
    PackReadID(Dm.sendList_Z, Dm.sendCount_Z ,sendID_Z, id);
    PackReadID(Dm.sendList_xy, Dm.sendCount_xy ,sendID_xy, id);
    PackReadID(Dm.sendList_Xy, Dm.sendCount_Xy ,sendID_Xy, id);
    PackReadID(Dm.sendList_xY, Dm.sendCount_xY ,sendID_xY, id);
    PackReadID(Dm.sendList_XY, Dm.sendCount_XY ,sendID_XY, id);
    PackReadID(Dm.sendList_xz, Dm.sendCount_xz ,sendID_xz, id);
    PackReadID(Dm.sendList_Xz, Dm.sendCount_Xz ,sendID_Xz, id);
    PackReadID(Dm.sendList_xZ, Dm.sendCount_xZ ,sendID_xZ, id);
    PackReadID(Dm.sendList_XZ, Dm.sendCount_XZ ,sendID_XZ, id);
    PackReadID(Dm.sendList_yz, Dm.sendCount_yz ,sendID_yz, id);
    PackReadID(Dm.sendList_Yz, Dm.sendCount_Yz ,sendID_Yz, id);
    PackReadID(Dm.sendList_yZ, Dm.sendCount_yZ ,sendID_yZ, id);
    PackReadID(Dm.sendList_YZ, Dm.sendCount_YZ ,sendID_YZ, id);
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
    UnpackReadID(Dm.recvList_x, Dm.recvCount_x ,recvID_x, id);
    UnpackReadID(Dm.recvList_X, Dm.recvCount_X ,recvID_X, id);
    UnpackReadID(Dm.recvList_y, Dm.recvCount_y ,recvID_y, id);
    UnpackReadID(Dm.recvList_Y, Dm.recvCount_Y ,recvID_Y, id);
    UnpackReadID(Dm.recvList_z, Dm.recvCount_z ,recvID_z, id);
    UnpackReadID(Dm.recvList_Z, Dm.recvCount_Z ,recvID_Z, id);
    UnpackReadID(Dm.recvList_xy, Dm.recvCount_xy ,recvID_xy, id);
    UnpackReadID(Dm.recvList_Xy, Dm.recvCount_Xy ,recvID_Xy, id);
    UnpackReadID(Dm.recvList_xY, Dm.recvCount_xY ,recvID_xY, id);
    UnpackReadID(Dm.recvList_XY, Dm.recvCount_XY ,recvID_XY, id);
    UnpackReadID(Dm.recvList_xz, Dm.recvCount_xz ,recvID_xz, id);
    UnpackReadID(Dm.recvList_Xz, Dm.recvCount_Xz ,recvID_Xz, id);
    UnpackReadID(Dm.recvList_xZ, Dm.recvCount_xZ ,recvID_xZ, id);
    UnpackReadID(Dm.recvList_XZ, Dm.recvCount_XZ ,recvID_XZ, id);
    UnpackReadID(Dm.recvList_yz, Dm.recvCount_yz ,recvID_yz, id);
    UnpackReadID(Dm.recvList_Yz, Dm.recvCount_Yz ,recvID_Yz, id);
    UnpackReadID(Dm.recvList_yZ, Dm.recvCount_yZ ,recvID_yZ, id);
    UnpackReadID(Dm.recvList_YZ, Dm.recvCount_YZ ,recvID_YZ, id);
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
    // if (rank == 0) printf("porosity = %.12f\n",porosity);
    
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
    // if (rank == 0) printf("porosity = %.12f\n",porosity);
    
    return porosity;
}



void lbpm_random_pp(int argc, char **argv,int rank, int nprocs, MPI_Comm comm)
{
    
    string filename = argv[1];
    
    auto db = std::make_shared<Database>( filename );
    auto domain_db = db->getDatabase( "Domain" );
    auto color_db = db->getDatabase( "Color" );
    //bool Restart = color_db->getScalar<bool>( "Restart" );
    
    if (rank == 0) { printf("\nUsing legacy pre-processor: LBPM_RANDOM_PP \n");}
    
    if (rank==0) printf("RUNNING RI! \n");
    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int nprocx, nprocy, nprocz, Nx, Ny, Nz, nspheres;
    double Lx, Ly, Lz;
    int i,j,k,n;
    
    auto Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis

    Nx = Dm->Nx;
    Ny = Dm->Ny;
    Nz = Dm->Nz;
    Lx = Dm->Lx;
    Ly = Dm->Ly;
    Lz = Dm->Lz;
    int iproc = Dm->iproc();
    int jproc = Dm->jproc();
    int kproc = Dm->kproc();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();

    double Saturation = domain_db->getScalar<double>("RandomSaturation");
    int BC = domain_db->getScalar<int>("BC");

    size_t N = size_t(Nx)*size_t(Ny)*size_t(Nz);
    
    /* Initialize domain class communication */
    for (n=0; n<N; n++) Dm->id[n]=1;  Dm->CommInit();  
    //for (n=0; n<N; n++) Dm->id[n]=2; // Set wetting
    
    char LocalRankFilename[40];
    

    char *ReadID;
    ReadID = new char[N];
    
    // Read the ID map from file
    sprintf(LocalRankFilename,"ID.%05i",rank);
    size_t readID;
    FILE *IDFILE = fopen(LocalRankFilename,"rb");
    if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
    readID=fread(ReadID,1,N,IDFILE);
    if (readID != size_t(N)) printf("lbpm_morphdrain_pp: Error reading ID (rank=%i) \n",rank);
    fclose(IDFILE);
    
    double porosity = calculate_porosity(ReadID, Nx, Ny, Nz, *Dm, rank);
    // if (rank == 0) cout << "Read ID file porosity=" << porosity << endl;
    
    double count = 0.f;
    double denom = 0;
    double sum_local = 0;


    int amin = Dm->amin;
    int amax = Dm->amax;
    // int amin = 1;
    // int amax = Nz-1;
    // if (BC > 0){
    //     if ( domain_db->keyExists( "AnalysisMin" ) && Dm->kproc() == 0){
    //         amin = domain_db->getScalar<int>( "AnalysisMin" );
    //     } 
    //     if ( domain_db->keyExists( "AnalysisMax" ) && Dm->kproc() == (nprocz - 1)){
    //         amax = domain_db->getScalar<int>( "AnalysisMax" );
    //         amax = Nz-1-((Nz-2)*nprocz - amax);
    //     }
    // }
    
    
    // int kmin=1; 
    // int kmax=Nz-1;
    // if (BC > 0 && kproc == 0) kmin=5;
    // if (BC > 0 && kproc == (nprocz-1)) kmax=Nz-1-5; 
    // if (BC > 0 && kproc == 0) kmin=126; //5
    // if (BC > 0 && kproc == (nprocz-1)) kmax=1;
    // if (BC > 0 && kproc == (nprocz-2)) kmax=1; 
    // if (BC > 0 && kproc == (nprocz-3)) kmax=Nz-1-135;  //Nz-1-82

    //Determine volume of void space
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                int n = k*Nx*Ny+j*Nx+i;
                if (ReadID[n] > 0) {
                    ReadID[n] = 1;
                } 
            }
        }
    }

    for (int k=amin; k<amax; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                int n = k*Nx*Ny+j*Nx+i;
                if (ReadID[n] > 0)  {
                    count += 1.0;
                }
            }
        }
    }

    sum_local = count; 
    MPI_Allreduce(&sum_local,&denom,1,MPI_DOUBLE,MPI_SUM,comm);
    
    

    int bin, binCount;
    ifstream Dist("BlobSize.in");
    Dist >> binCount;
    int *SizeX, *SizeY, *SizeZ;
    if (binCount == 0){ printf("Warning: BlobSize.in is empty"); }
    SizeX = new int [binCount];
    SizeY = new int [binCount];
    SizeZ = new int [binCount];
    for (bin=0; bin<binCount; bin++){
        Dist >> SizeX[bin];
        Dist >> SizeY[bin];
        Dist >> SizeZ[bin];
    }
    Dist.close();
    
    // Generate the residual NWP
    if (rank==0) printf("Initializing with saturation = %f \n",Saturation);
    //if (rank==0) printf("Current Saturation \n");
    fflush(stdout);
    
    double sum_sat_local = 0;
    double sat = 0.0;
    
    
    int x,y,z;
    int sizeX,sizeY,sizeZ;
    int ii,jj,kk;
    int Number = 0;		// number of features
    while (sat < Saturation) {
        if (rank==0) {
            Number++;
            // Randomly generate a point in the domain
            x = (Nx-2)*nprocx*float(rand())/float(RAND_MAX);
            y = (Ny-2)*nprocy*float(rand())/float(RAND_MAX);
            z = (Nz-2)*nprocz*float(rand())/float(RAND_MAX);
            
            bin = int(floor(binCount*float(rand())/float(RAND_MAX)));
            sizeX = SizeX[bin];
            sizeY = SizeY[bin];
            sizeZ = SizeZ[bin];
        }
        MPI_Bcast(&x,1,MPI_INT,0,comm);
        MPI_Bcast(&y,1,MPI_INT,0,comm);
        MPI_Bcast(&z,1,MPI_INT,0,comm);
        MPI_Bcast(&sizeX,1,MPI_INT,0,comm);
        MPI_Bcast(&sizeY,1,MPI_INT,0,comm);
        MPI_Bcast(&sizeZ,1,MPI_INT,0,comm);
        
        //if (rank==0) printf("Broadcast block at %i,%i,%i \n",x,y,z);
        
        for (k=z;k<z+sizeZ;k++){
            for (j=y;j<y+sizeY;j++){
                for (i=x;i<x+sizeX;i++){
                    
                    // Identify nodes in the domain (periodic BC)
                    ii = i;
                    jj = j;
                    kk = k;
                    
                    if (ii>nprocx*(Nx-2)) ii-=nprocx*(Nx-2);
                    if (jj>nprocy*(Ny-2)) jj-=nprocy*(Ny-2);
                    if (kk>nprocz*(Nz-2)) kk-=nprocz*(Nz-2);
                    
                    // Check if this is in the subdomain
                    if (ii < (iproc+1)*(Nx-2)+1 && jj < (jproc+1)*(Ny-2)+1 && kk < (kproc+1)*(Nz-2)+1 &&
                        ii  > iproc*(Nx-2) && jj > jproc*(Ny-2) && kk > kproc*(Nz-2) ){
                        
                        // Map from global to local coordinates
                        ii -= iproc*(Nx-2);
                        jj -= jproc*(Ny-2);
                        kk -= kproc*(Nz-2);
                        
                        n = kk*Nx*Ny+jj*Nx+ii;

                        // BC > 0? Only add inside analysis region
                        if (ReadID[n] == 1 && kk >= amin && kk < amax) { 
                            ReadID[n] = 2;
                            // if (BC > 0){
                                // if (kproc == 0 && kk >= kmin){
                                //     ReadID[n] = 2;
                                // } else if (kproc == (nprocz-3) && kk < kmax){
                                //     ReadID[n] = 2;
                                // } else if (kproc < (nprocz-3)){
                                //     ReadID[n] = 2;
                                // }
                            // } else{
                            //     ReadID[n] = 2;
                            // }
                            // if (BC > 0 && kproc == 0 && kk >= 20){
                            //     ReadID[n] = 2;
                            // } else if (BC > 0 && kproc == (nprocz-1) && kk < (Nz-20)){
                            //     ReadID[n] = 2;
                            // } else {
                            //     ReadID[n] = 2;
                            // }
                        }
                        
                    }
                }
            }
        }
        
        

        count = 0.f;
        for (int k=amin; k<amax; k++){
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++){
                    int n = k*Nx*Ny+j*Nx+i;
                    if (ReadID[n] == 2)  {count += 1.0;}
                }
            }
        }
        
        sat = 0;
        sum_sat_local = count;
        MPI_Allreduce(&sum_sat_local,&sat,1,MPI_DOUBLE,MPI_SUM,comm);
        
        
        sat = sat/denom;
        //if (rank==0) printf("New count=%i\n",countGlobal);
        if (rank==0) printf("%f\n",sat);
        fflush(stdout);
    }
    
    
    //    for ( k=0;k<2;k++){
    //        for ( j=0;j<2;j++){
    //            for ( i=0;i<2;i++){
    //                n = k*Nx*Ny+j*Nx+i;
    //                if (i == 1 && j == 1 && k == 1) {int m = (Nz-1)*Nx*Ny + (Ny-1)*Nx+(Nz-1);    Dm->id[m] = ReadID[n]; }
    //            }
    //        }
    //    }
    
    //  if (InitialWetting == 1)	FlipID(ReadID,Nx*Ny*Nz);
    
    
    //PopulateHalo_Char(ReadID, *Dm, Nx, Ny, Nz, rank);
    
    //ReadID[0] = ReadID[(Nx-1)*(Ny-1)*(Nz-1)];
    //porosity = calculate_porosity(ReadID, Nx, Ny, Nz, *Dm, rank);
    //if (rank == 0) cout << "Writing new ID files with porosity=" << porosity << endl;
    
    
    if (rank == 0) printf("Writing ID2 file \n");
    sprintf(LocalRankFilename,"ID2.%05i",rank);
    FILE *ID = fopen(LocalRankFilename,"wb");
    fwrite(ReadID,1,N,ID);
    fclose(ID);
    
    if (rank == 0) cout << "End of random saturation pp." << endl << endl;

}




void lbpm_partial_saturation_pp(int argc, char **argv,int rank, int nprocs, MPI_Comm comm)
{
    
    string filename = argv[1];
    
    auto db = std::make_shared<Database>( filename );
    auto domain_db = db->getDatabase( "Domain" );
    auto color_db = db->getDatabase( "Color" );
    auto analysis_db = db->getDatabase( "Analysis" );
    
    
    bool Restart = color_db->getScalar<bool>( "Restart" );
    
    // if (Restart == 0) {
    
    if (rank == 0) { printf("\nUsing legacy pre-processor: LBPM_PARTIAL_SATURATION_PP \n");}
    
    
    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int nprocx, nprocy, nprocz, Nx, Ny, Nz, nspheres;
    double Lx, Ly, Lz;
    int i,j,k,n;
    int BC=0;
    
    
    //  auto db = std::make_shared<Database>( filename );
    //  auto domain_db = db->getDatabase( "Domain" );
    
    auto Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));      // full domain for analysis
    Nx = Dm->Nx;
    Ny = Dm->Ny;
    Nz = Dm->Nz;
    Lx = Dm->Lx;
    Ly = Dm->Ly;
    Lz = Dm->Lz;
    int iproc = Dm->iproc();
    int jproc = Dm->jproc();
    int kproc = Dm->kproc();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();
    nspheres = domain_db->getScalar<int>( "nspheres");
    int InitialWetting = domain_db->getScalar<int>("InitialSaturation");
    double Saturation = domain_db->getScalar<double>("RandomSaturation");
    /* Fact check these print statements... */
    // if (rank==0){ printf("Initializing wetting phase saturation of %f \n",Saturation);
    // if (InitialWetting == 1) {printf("Initial connected phase labeled (1) \n");} else  {printf("Initial connected phase labeled (2) \n");} }
    if (InitialWetting == 1) Saturation=1.0-Saturation;
    int BoundaryCondition=1;
    
//    Nx+=2;
//    Ny+=2;
//    Nz+=2;
    printf("Nx=%d Ny=%d Nz=%d\n",Nx,Ny,Nz);
//    Nx += 2;
//    Nx = Ny = Nz;
    int N = Nx*Ny*Nz;
    
    /* Initialize domain class communication */
    for (n=0; n<N; n++) Dm->id[n]=1;  Dm->CommInit();  for (n=0; n<N; n++) Dm->id[n]=2; // Set wetting
    
    char LocalRankFilename[40];
    

    char *ReadID;
    ReadID = new char[N];
    
    // Read the signed distance from file
    sprintf(LocalRankFilename,"ID.%05i",rank);
    FILE *DIST = fopen(LocalRankFilename,"rb");
    size_t ReadSignDist;
    ReadSignDist=fread(ReadID,1,N,DIST);
    if (ReadSignDist != size_t(N)) printf("lbpm_random_pp: Error reading ID file (rank=%i)\n",rank);
    fclose(DIST);
    
//    double porosity = calculate_porosity(ReadID, Nx, Ny, Nz, *Dm, rank);
//    if (rank == 0) cout << "Read ID file porosity=" << porosity << endl;
//
//    double count = 0;
//    double swcount = 0;
//    double totalGlobal = 0;
//    double swGlobal = 0;
//
//    double denom = 0;
//    double sum_local = 0;
//
//   // for (int i = 0; i < Nx*Ny*Nz; i++) if (ReadID[n]==)
//
//    /* Update the Dm->id */
//    for (int k=1; k<Nz-1; k++){
//        for (int j=1; j<Ny-1; j++){
//            for (int i=1; i<Nx-1; i++){
//                n = k*Nx*Ny+j*Nx+i;
//                if (ReadID[n] > 0)  {
//			count = count + 1.0;
//			ReadID[n] = 2;
//		}
//		if (ReadID[n] == 2) {swcount = swcount + 1.0;}
//
//            }
//        }
//    }
//
//    MPI_Allreduce(&count,&totalGlobal,1,MPI_DOUBLE,MPI_SUM,comm);
//    MPI_Allreduce(&swcount,&swGlobal,1,MPI_DOUBLE,MPI_SUM,comm);
//    double sat = swGlobal/totalGlobal;
//    if (Dm->rank()==0){
//	    printf("count pre: %f \n",count);
//            printf("swcount pre: %f \n",swcount);
//	    printf("SW Pre: %f \n",sat);
//    }
//    //printf("rank Nx Ny Nz: %i %i %i %i \n",Dm->rank(),Nx,Ny,Nz);
//    //printf("rank iproc jproc kproc: %i %i %i %i \n",Dm->rank(),iproc,jproc,kproc);
//
//
//    //double center_k = 2.*double(Nz)/3.-8 + double((Nz-2));
//    //double center_ij = double(Nx)/2.;
//    printf("woof\n");
//
//    //double lower_bd = 0.3*nprocz*(Nz-2);
//    //double upper_bd = 0.7*nprocz*(Nz-2);
//    //int lowerxbd = int(0.3*double(nprocx*(Nx-2)));
//    //int upperxbd = int(0.7*double(nprocx*(Nx-2)));
////std::cout << "lower bounds:" << lowerxbd << " upper bounds:" << upperxbd << std::endl;
//
//   // for (int n=0; n<N; n++){ if (ReadID[n] == 2) ReadID[n] = 1; }
//    //double beta = 0.8;
//    //double betadelta = 1.0-beta;
//    int extra_counter = 0;
//    for (int k=0; k<Nz; k++){
//        for (int j=0; j<Ny; j++){
//            for (int i=0; i<Nx; i++){
//                size_t n = k*Nx*Ny+j*Nx+i;
//
//                // //write a cube
//                // if (ReadID[n] > 0){
//                //     if ( ((i+(Nx-2)*iproc) < 40) && ((i+(Nx-2)*iproc) > 20) ) {
//                //         if ( ((j+(Ny-2)*jproc) < 40) && ((j+(Ny-2)*jproc) > 20) ) {
//                //             if ( ((k+(Nz-2)*kproc) < 40)  && ((k+(Nz-2)*kproc) > 20) ) {
//                //                 ReadID[n] = 1;
//                //                 extra_counter = extra_counter + 1;
//                //             }
//                //         }
//                //     }
//                // }
//
//                // //write hemisphere ORIGINAL for akai solid sphere
//                // if (ReadID[n] > 0 && (k+(Nz-2)*kproc) > 30){
//                //    if ( ((i+(Nx-2)*iproc) - 50)*((i+(Nx-2)*iproc) - 50) + ((j+(Ny-2)*jproc) - 50)*((j+(Ny-2)*jproc) - 50) +
//                //    ((k+(Nz-2)*kproc) - 30)*((k+(Nz-2)*kproc) - 30) <= 28*28   ){
//                //        ReadID[n] = 1;
//                //    }
//                // }
//
//                // // //write hemisphere excel match
//                // if (ReadID[n] > 0 && (k+(Nz-2)*kproc) > 30){
//                //    if ( ((i+(Nx-2)*iproc) - 50)*((i+(Nx-2)*iproc) - 50) + ((j+(Ny-2)*jproc) - 50)*((j+(Ny-2)*jproc) - 50) +
//                //    ((k+(Nz-2)*kproc) - 30)*((k+(Nz-2)*kproc) - 30) <= 28*28   ){
//                //        ReadID[n] = 1;
//                //        //extra_counter = extra_counter + 1;
//                //    }
//                // }
//
//                // //write a cube for akai plate
////                if (ReadID[n] > 0){
////                    if ( ((i+(Nx-2)*iproc) < 35) && ((i+(Nx-2)*iproc) > 15) ) {
////                        if ( ((j+(Ny-2)*jproc) < 35) && ((j+(Ny-2)*jproc) > 15) ) {
////                            if ( ((k+(Nz-2)*kproc) < 20)  && ((k+(Nz-2)*kproc) > 0) ) {
////                                ReadID[n] = 1;
////                                extra_counter = extra_counter + 1;
////                            }
////                        }
////                    }
////                }
//
//                // 36^3 case
////                {
////                    if ( ((i+(Nx-2)*iproc) > 9) && ((i+(Nx-2)*iproc) < 30) ) {
////                        if ( ((j+(Ny-2)*jproc) > 9) && ((j+(Ny-2)*jproc) < 30) ) {
////                            if ( ((k+(Nz-2)*kproc) > 9)  && ((k+(Nz-2)*kproc) < 30) ) {
////                                if (ReadID[n] == 2) {ReadID[n] = 1;  }
////                              //  extra_counter = extra_counter + 1;
////                            }
////                        }
////                    }
////                }
//
//                // 12^3 case
//                {
//                    if ( ((i+(Nx-2)*iproc) > 3) && ((i+(Nx-2)*iproc) < 10) ) {
//                        if ( ((j+(Ny-2)*jproc) > 3) && ((j+(Ny-2)*jproc) < 10) ) {
//                            if ( ((k+(Nz-2)*kproc) >= 1)  && ((k+(Nz-2)*kproc) <= 12) ) {
//                                if (ReadID[n] == 2) {ReadID[n] = 1;  }
//                              //  extra_counter = extra_counter + 1;
//                            }
//                        }
//                    }
//                }
//
//            }
//        }
//    }
//    //sum_local = 1.0*double(count); MPI_Allreduce(&sum_local,&denom,1,MPI_DOUBLE,MPI_SUM,comm);
//
//
//   //PopulateHalo_Char(ReadID, *Dm, Nx, Ny, Nz, rank);
//
//
//    //porosity = calculate_porosity(ReadID, Nx, Ny, Nz, *Dm, rank);
//    //if (rank == 0) cout << "Writing new ID files with porosity=" << porosity << endl;
//
//
//
//    count = 0;
//    swcount = 0;
//    totalGlobal = 0;
//    swGlobal = 0;
//
//   // if (Dm->rank()==0){
//   // 	printf("count,swcount,totalGlobal,swGlobal: %f %f %f %f \n",count,swcount,totalGlobal,swGlobal);
//   // }
//
//    for (int k=1; k<Nz-1; k++){
//        for (int j=1; j<Ny-1; j++){
//            for (int i=1; i<Nx-1; i++){
//                n = k*Nx*Ny+j*Nx+i;
//                if (ReadID[n] > 0)  {
//                    count += 1;
//                }
//                if (ReadID[n] == 2) {swcount +=1;}
//
//            }
//        }
//    }
//
//    MPI_Allreduce(&count,&totalGlobal,1,MPI_DOUBLE,MPI_SUM,comm);
//    MPI_Allreduce(&swcount,&swGlobal,1,MPI_DOUBLE,MPI_SUM,comm);
//    sat = swGlobal/totalGlobal;
//    if (Dm->rank()==0){
//        printf("count post: %f \n",count);
//            printf("swcount post: %f \n",swcount);
//        printf("SW Post: %f \n",sat);
//    }
//

    
    //sprintf(LocalRankFilename,"%s%s/ID.%05i","/mnt/bb/",usr,rank);
    sprintf(LocalRankFilename,"ID2.%05i",rank);
    FILE *ID = fopen(LocalRankFilename,"wb");
    fwrite(ReadID,1,N,ID);
    fclose(ID);
    
   // if (rank == 0 ) cout << "End of random saturation pp." << endl << endl;
}


void lbpm_morphopen_pp(int argc, char **argv,int rank, int nprocs, MPI_Comm comm)
{
    {
        //.......................................................................
        // Reading the domain information file
        //.......................................................................
        int nprocx, nprocy, nprocz;
        
        int i,j,k,n;
        
        //  char fluidValue,solidValue;
        
        
        
        
        char LocalRankFilename[40];
        
        string filename;
        double Rcrit_new, SW;
        if (argc > 1){
            filename=argv[1];
            Rcrit_new=0.f;
            //SW=strtod(argv[2],NULL);
        }
        else ERROR("No input database provided\n");
        // read the input database
        auto db = std::make_shared<Database>( filename );
        auto domain_db = db->getDatabase( "Domain" );
        
        // Read domain parameters
        auto L = domain_db->getVector<double>( "L" );
        auto size = domain_db->getVector<int>( "n" );
        auto nproc = domain_db->getVector<int>( "nproc" );
        auto ReadValues = domain_db->getVector<char>( "ReadValues" );
        auto WriteValues = domain_db->getVector<char>( "WriteValues" );
        SW = domain_db->getScalar<double>("RandomSaturation");
        
        if (rank == 0) { printf("\nUsing legacy pre-processor: LBPM_MORPHOPEN_PP \n");}
        
        if (rank==0)    printf("Target saturation %f \n",SW);
        
        
        
        std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
        
        
        int nx = Dm->Nx; int ny = Dm->Ny;  int nz = Dm->Nz;
        nprocx = Dm->nprocx();  nprocy = Dm->nprocy();  nprocz = Dm->nprocz();
        size_t N = size_t(nx)*size_t(ny)*size_t(nz);
        
        for (n=0; n<N; n++) Dm->id[n]=1;
        Dm->CommInit();
        
        
        
        
        char *ReadID;
        ReadID = new char[N];
       
        char *SDsID;
        SDsID = new char[N];
        
        // Read the signed distance from file
        sprintf(LocalRankFilename,"ID.%05i",rank);
        FILE *DIST = fopen(LocalRankFilename,"rb");
        size_t ReadSignDist;
        ReadSignDist=fread(ReadID,1,N,DIST);
        if (ReadSignDist != size_t(N)) printf("lbpm_random_pp: Error reading ID file (rank=%i)\n",rank);
        fclose(DIST);
        
        for (int n = 0; n < nx*ny*nz; n++) SDsID[n] = 1;
        
        
        
        // ReadID is the original domain.
        
        double porosity = calculate_porosity(ReadID, nx, ny, nz, *Dm, rank); // For permeability test of BCC
        // std::cout << "Before Porosity=" << porosity << std::endl;
        
        // Generate the signed distance map
        // Initialize the domain and communication
        Array<char> id_solid(nx,ny,nz);
        DoubleArray SignDist(nx,ny,nz);
        
        
                int xdim,ydim,zdim;
                    xdim=Dm->Nx-2;
                    ydim=Dm->Ny-2;
                    zdim=Dm->Nz-2;
                fillHalo<double> fillData(Dm->Comm, Dm->rank_info,{xdim,ydim,zdim},{1,1,1},0,1);
        
        
                    // Solve for the position of the solid phase
                    for (int k=0;k<nz;k++){
                        for (int j=0;j<ny;j++){
                            for (int i=0;i<nx;i++){
                                int n = k*nx*ny+j*nx+i;
                                // Initialize the solid phase
                                if (ReadID[n] > 0)    id_solid(i,j,k) = 1;
                                else             id_solid(i,j,k) = 0;
                            }
                        }
                    }
                    // Initialize the signed distance function
                    for (int k=0;k<nz;k++){
                        for (int j=0;j<ny;j++){
                            for (int i=0;i<nx;i++){
                                int n = k*nx*ny+j*nx+i;
                                // Initialize distance to +/- 1
                                SignDist(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
                            }
                        }
                    }
        
                    if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
                    CalcDist(SignDist,id_solid,*Dm);
        
        
        
        Dm->CommunicateMeshHalo(SignDist);
        
        for (int n = 0; n < nx*ny*nz; n++) SignDist(n) +=2;
        
        
        // Read the signed distance from file
//        sprintf(LocalRankFilename,"SignDist.%05i",rank);
//        FILE *SDS = fopen(LocalRankFilename,"rb");
//        
//        ReadSignDist=fread(SignDist.data(),8,N,SDS);
//        if (ReadSignDist != size_t(N)) printf("lbpm_morphdrain_pp: Error reading signed distance function (rank=%i)\n",rank);
//        fclose(SDS);
        
        // fillData.fill(SignDist);

        int kproc = Dm->kproc();

        MPI_Barrier(comm);
        double count,countGlobal,totalGlobal;
        count = 0.f;
        double maxdist=-200.f;
        double maxdistGlobal;
        for (int k=1; k<nz-1; k++){
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++){
                    n = k*nx*ny+j*nx+i;
                    // extract distance for critical radius
                    if ( SignDist(i,j,k) > maxdist) maxdist=SignDist(i,j,k);
                }
            }
        }
        for (int k=0; k<nz; k++){
            for (int j=0; j<ny; j++){
                for (int i=0; i<nx; i++){
                    n = k*nx*ny+j*nx+i;
                    if (SignDist(i,j,k) < 0) {
                        // don't do anything
                    }
                    else{
                        // initially saturated with wetting phase
                        SDsID[n] = 2;
                        // if (kproc < 1){
                        //     SDsID[n] = 1;
                        // } 
                        // if (kproc == 1 && k < 35){
                        //     SDsID[n] = 1;
                        // }
                        count+=1.0;
                    }
                    // don't let halo be the maximum dist
                    if ( SignDist(i,j,k) > maxdist) SignDist(i,j,k) = maxdist;
                }
            }
        }
       
        
      //  for (int n = 0; n < nx*ny*nz; n++) SignDist(n) -=3;
        
        MPI_Barrier(comm);
        // total Global is the number of nodes in the pore-space
        MPI_Allreduce(&count,&totalGlobal,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&maxdist,&maxdistGlobal,1,MPI_DOUBLE,MPI_MAX,comm);
        double volume=double(nprocx*nprocy*nprocz)*double(nx-2)*double(ny-2)*double(nz-2);
        porosity=totalGlobal/volume;
        if (rank==0) printf("Media Porosity: %f \n",porosity);
        if (rank==0) printf("Maximum pore size: %f \n",maxdistGlobal);
        
        
        
        // Generate the NWP configuration
        //if (rank==0) printf("Initializing morphological distribution with critical radius %f \n", Rcrit);
        if (rank==0) printf("Performing morphological opening with target saturation %f \n", SW);
        //    GenerateResidual(id,nx,ny,nz,Saturation);
        
        // Communication buffers
        char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
        char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
        char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
        char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
        char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
        char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
        // send buffers
        sendID_x = new char [Dm->sendCount_x];
        sendID_y = new char [Dm->sendCount_y];
        sendID_z = new char [Dm->sendCount_z];
        sendID_X = new char [Dm->sendCount_X];
        sendID_Y = new char [Dm->sendCount_Y];
        sendID_Z = new char [Dm->sendCount_Z];
        sendID_xy = new char [Dm->sendCount_xy];
        sendID_yz = new char [Dm->sendCount_yz];
        sendID_xz = new char [Dm->sendCount_xz];
        sendID_Xy = new char [Dm->sendCount_Xy];
        sendID_Yz = new char [Dm->sendCount_Yz];
        sendID_xZ = new char [Dm->sendCount_xZ];
        sendID_xY = new char [Dm->sendCount_xY];
        sendID_yZ = new char [Dm->sendCount_yZ];
        sendID_Xz = new char [Dm->sendCount_Xz];
        sendID_XY = new char [Dm->sendCount_XY];
        sendID_YZ = new char [Dm->sendCount_YZ];
        sendID_XZ = new char [Dm->sendCount_XZ];
        //......................................................................................
        // recv buffers
        recvID_x = new char [Dm->recvCount_x];
        recvID_y = new char [Dm->recvCount_y];
        recvID_z = new char [Dm->recvCount_z];
        recvID_X = new char [Dm->recvCount_X];
        recvID_Y = new char [Dm->recvCount_Y];
        recvID_Z = new char [Dm->recvCount_Z];
        recvID_xy = new char [Dm->recvCount_xy];
        recvID_yz = new char [Dm->recvCount_yz];
        recvID_xz = new char [Dm->recvCount_xz];
        recvID_Xy = new char [Dm->recvCount_Xy];
        recvID_xZ = new char [Dm->recvCount_xZ];
        recvID_xY = new char [Dm->recvCount_xY];
        recvID_yZ = new char [Dm->recvCount_yZ];
        recvID_Yz = new char [Dm->recvCount_Yz];
        recvID_Xz = new char [Dm->recvCount_Xz];
        recvID_XY = new char [Dm->recvCount_XY];
        recvID_YZ = new char [Dm->recvCount_YZ];
        recvID_XZ = new char [Dm->recvCount_XZ];
        //......................................................................................
        int sendtag,recvtag;
        sendtag = recvtag = 7;
        
        
        int ii,jj,kk;
        int Nx = nx;
        int Ny = ny;
        int Nz = nz;
        
        // std::cout << "nx=" << nx << " ny=" << ny << " nz=" << nz << " N=" << N <<  std::endl;
        
        double sw_old=1.0;
        double sw_new=1.0;
        double sw_diff_old = 1.0;
        double sw_diff_new = 1.0;
        
        // Increase the critical radius until the target saturation is met
        double deltaR=0.05; // amount to change the radius in voxel units
        double Rcrit_old;
        
        double GlobalNumber = 1.f;
        int imin,jmin,kmin,imax,jmax,kmax;
        
        Rcrit_new = maxdistGlobal;
        //if (argc>2){
        //    Rcrit_new = strtod(argv[2],NULL);
        //    if (rank==0) printf("Max. distance =%f, Initial critical radius = %f \n",maxdistGlobal,Rcrit_new);
        //}
        while (sw_new > SW)
        {
            sw_diff_old = sw_diff_new;
            sw_old = sw_new;
            Rcrit_old = Rcrit_new;
            Rcrit_new -= deltaR*Rcrit_old;
            int Window=round(Rcrit_new);
            if (Window == 0) Window = 1; // If Window = 0 at the begining, after the following process will have sw=1.0
            // and sw<Sw will be immediately broken
            double LocalNumber=0.f;
            for(k=0; k<Nz; k++){
                for(j=0; j<Ny; j++){
                    for(i=0; i<Nx; i++){
                        n = k*nx*ny + j*nx+i;
                        if (SignDist(i,j,k) > Rcrit_new){
                            // loop over the window and update
                            imin=max(1,i-Window);
                            jmin=max(1,j-Window);
                            kmin=max(1,k-Window);
                            imax=min(Nx-1,i+Window);
                            jmax=min(Ny-1,j+Window);
                            kmax=min(Nz-1,k+Window);
                            for (kk=kmin; kk<kmax; kk++){
                                for (jj=jmin; jj<jmax; jj++){
                                    for (ii=imin; ii<imax; ii++){
                                        int nn = kk*nx*ny+jj*nx+ii;
                                        double dsq = double((ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k));
                                        if (SDsID[nn] == 2 && dsq <= Rcrit_new*Rcrit_new){
                                            LocalNumber+=1.0;
                                            SDsID[nn]=1;
                                        }
                                    }
                                }
                            }
                            
                        }
                        // move on
                    }
                }
            }
            
            
            
            
            
            DoubleArray TempDomain; TempDomain.resize(Nx,Ny,Nz); TempDomain.fill(2);
            
            for (int k=1; k<Nz-1; k++)
            for (int j=1; j<Ny-1; j++)
            for (int i=1; i<Nx-1; i++){
                size_t n = k*Nx*Ny+j*Nx+i;
                TempDomain(n) = (double)SDsID[n];
            }
            
            Dm->CommunicateMeshHalo(TempDomain);
            
            for (int n = 0; n < Nx*Ny*Nz; n++) SDsID[n] = (char)TempDomain(n);
            
            for (int n = 0; n < Nx*Ny*Nz; n++) {
                if (ReadID[n] > 0) ReadID[n] = SDsID[n];
                
                
            }
            
            
            MPI_Allreduce(&LocalNumber,&GlobalNumber,1,MPI_DOUBLE,MPI_SUM,comm);
            
            count = 0.f;
            for (int k=1; k<Nz-1; k++){
                for (int j=1; j<Ny-1; j++){
                    for (int i=1; i<Nx-1; i++){
                        n=k*Nx*Ny+j*Nx+i;
                        if (ReadID[n] == 2){
                            count+=1.0;
                        }
                    }
                }
            }
            MPI_Allreduce(&count,&countGlobal,1,MPI_DOUBLE,MPI_SUM,comm);
            sw_new = countGlobal/totalGlobal;
            sw_diff_new = abs(sw_new-SW);
            if (rank==0){
                printf("     %f ",sw_new);
                printf("     %f\n",Rcrit_new);
            }
        }
        
        if (sw_diff_new<sw_diff_old){
            if (rank==0){
                printf("Final saturation=%f\n",sw_new);
                printf("Final critical radius=%f\n",Rcrit_new);
            }
        }
        else{
            if (rank==0){
                printf("Final saturation=%f\n",sw_old);
                printf("Final critical radius=%f\n",Rcrit_old);
            }
        }
        
        // DoubleArray TempDomain; TempDomain.resize(Nx,Ny,Nz); TempDomain.fill(2);
        
        // for (int k=1; k<Nz-1; k++)
        // for (int j=1; j<Ny-1; j++)
        // for (int i=1; i<Nx-1; i++){
        //     size_t n = k*Nx*Ny+j*Nx+i;
        //     TempDomain(n) = (double)SDsID[n];
        // }
        
        // Dm->CommunicateMeshHalo(TempDomain);
        
        // for (int n = 0; n < Nx*Ny*Nz; n++) SDsID[n] = (char)TempDomain(n);
        
        // for (int n = 0; n < Nx*Ny*Nz; n++) {
        //     if (ReadID[n] > 0) ReadID[n] = SDsID[n];
            
            
        // }
        
        if (rank==0) printf("Writing ID2 file \n");
        sprintf(LocalRankFilename,"ID2.%05i",rank);
        FILE *IDFILE = fopen(LocalRankFilename,"wb");
        fwrite(ReadID,1,N,IDFILE);
        fclose(IDFILE);
        
       // porosity = calculate_porosity(ReadID, Nx, Ny, Nz, *Dm, rank); // For permeability test of BCC
        //std::cout << "Porosity=" << porosity << std::endl;
    } 
}
void lbpm_morphdrain_pp(int argc, char **argv,int rank, int nprocs, MPI_Comm comm)
{
    {
        //.......................................................................
        // Reading the domain information file
        //.......................................................................
        int nprocx, nprocy, nprocz;
        int i,j,k,n;

        int MAXTIME=1000;
        int READ_FROM_BLOCK=0;

        char LocalRankFilename[40];
        char LocalRankFoldername[40];

        char LocalRankFilenameSub_format[] = "SW%03f/ID2.%05i"; 
        char LocalRankFilenameSub[sizeof(LocalRankFilenameSub_format)+20];
        
        string filename;
        double Rcrit_new, SW;
        if (argc > 1){
            filename=argv[1];
            Rcrit_new=0.f;
        }
        else ERROR("No input database provided\n");

        // read the input database
        auto db = std::make_shared<Database>( filename );
        auto domain_db = db->getDatabase( "Domain" );
        
        // Read domain parameters
        auto L = domain_db->getVector<double>( "L" );
        auto size = domain_db->getVector<int>( "n" );
        auto nproc = domain_db->getVector<int>( "nproc" );
        auto ReadValues = domain_db->getVector<char>( "ReadValues" );
        auto WriteValues = domain_db->getVector<char>( "WriteValues" );
        auto BC = domain_db->getScalar<int>("BC");
        auto READFILE = domain_db->getScalar<std::string>( "Filename" );
	    int nspheres = domain_db->getScalar<int>( "nspheres");

        bool ScanningCurve = false;
        bool FromRaw = false;
        double deltaR = 0.05;

        if (domain_db->keyExists( "ScanningCurve" )){
            ScanningCurve = domain_db->getScalar<bool>( "ScanningCurve" );
        }
        if (domain_db->keyExists( "FromRaw" )){
            FromRaw = domain_db->getScalar<bool>( "FromRaw" );
        }
        if (domain_db->keyExists( "Rcrit" )){
            Rcrit_new = domain_db->getScalar<double>("Rcrit");
        }
        if (domain_db->keyExists( "deltaR")){
            deltaR = domain_db->getScalar<double>("deltaR");
        } 

        int target_count = 0;
        double *SW_List;

        //Read in single sw target from RandomSaturation or multiple from targets.in
        //If reading multiple targets, create directories to store results
        if (domain_db->keyExists( "RandomSaturation" )){
            SW = domain_db->getScalar<double>("RandomSaturation");
            if (rank == 0) printf("Single Saturation Target: %f \n",SW);
            target_count = 1;
            SW_List = new double[target_count];
            SW_List[0] = SW;
        } else {
            double value;
            FILE * pFile;
            pFile = fopen ("targets.in","r");
            char line[ 100 ];
            while (fgets(line,100,pFile)) {
                target_count += 1;
            }
            fclose (pFile);

            pFile = fopen ("targets.in","r");
            SW_List = new double[target_count];
            target_count = 0;
            while (fgets(line,100,pFile)) {
                value = strtod(line,NULL);
                SW_List[target_count] = value;
                target_count += 1;
            }
            fclose (pFile);

            if (rank == 0) {
                printf("Quantity of sw targets: %i \n",target_count);
                sort(SW_List,SW_List+target_count,greater<double>{});
                for (int i=0;i<target_count;i++){
                    sprintf(LocalRankFoldername,"SW%03f",SW_List[i]);
                    mkdir(LocalRankFoldername,S_IRWXU|S_IRGRP);
                }
                printf("\n");
            }
        }

        if (SW_List[0] == 1.0 && BC > 0 && rank == 0){
            printf("Warning: Initializing to 1.0 with BC > 0 may produce nonwetting bubbles\n");
        }
        
        std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
        std::shared_ptr<Domain> Mask (new Domain(domain_db,comm));
        
        int rank = Dm->rank();

        int nx = Dm->Nx; int ny = Dm->Ny;  int nz = Dm->Nz;
        nprocx = Dm->nprocx();  nprocy = Dm->nprocy();  nprocz = Dm->nprocz();
        size_t N = size_t(nx)*size_t(ny)*size_t(nz);

        int amin = Dm->amin;
        int amax = Dm->amax;

        for (n=0; n<N; n++) Dm->id[n]=1;
        Dm->CommInit();
      
        char *ReadID;
        ReadID = new char[N];
       
        char *SDsID;
        SDsID = new char[N];

        //double *SDs; SDs = new double[N];
        
        Array<char> id_solid(nx,ny,nz);
        DoubleArray SignDist(nx,ny,nz);
        DoubleArray phase(nx,ny,nz);
        IntArray phase_label(nx,ny,nz);

        if (FromRaw == false){
        // Read the ID from file
            if (rank==0) printf("Read ID from RT Output \n");
            sprintf(LocalRankFilename,"ID.%05i",rank);
            FILE *IDFILE = fopen(LocalRankFilename,"rb");
            size_t readID;
            if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
            readID=fread(ReadID,1,N,IDFILE);
            if (readID != size_t(N)) printf("lbpm_morphdrain_pp: Error reading ID (rank=%i) \n",rank);
            fclose(IDFILE);
        } else {
            if (rank==0) printf("Read ID from raw File \n");
            Mask->Decomp(READFILE,db);
            Mask->CommInit();
            for (int k=0; k<nz; k++){
                for (int j=0; j<ny; j++){
                    for (int i=0; i<nx; i++){
                        n=k*nx*ny+j*nx+i;
                        ReadID[n] = Mask->id[n];
                    }
                }
            }
        }

	// Print ID
        //for (int k=50; k<51; k++){
        //    printf("%i \n",k);
        //    for (int j=0; j<ny; j++){
        //        for (int i=0; i<nx; i++){
        //            n = k*nx*ny+j*nx+i;
        //            printf("%i ",ReadID[n]);
        //        }
        //        printf("\n");
        //    }
        //    printf("\n \n");
        //}


        //Read SignDist from file
        // sprintf(LocalRankFilename,"SignDist.%05i",rank);
        // FILE *SDist = fopen(LocalRankFilename,"rb");
        // size_t ReadSignDist;
        // ReadSignDist=fread(SDs,8,N,SDist);
        // if (ReadSignDist != size_t(N)) printf("lbpm_morphdrain_pp: Error reading SDist \n");
        // fclose(SDist);

        //Calculate SignDist from ID field
        for (int k=0;k<nz;k++){
            for (int j=0;j<ny;j++){
                for (int i=0;i<nx;i++){
                    int n = k*nx*ny+j*nx+i;
                    // Initialize the solid phase
                    if (ReadID[n] > 0){
                        id_solid(i,j,k) = 1;
                    }
                    else	    
                        id_solid(i,j,k) = 0;
                }
            }
        }
        // Initialize the signed distance function
        for (int k=0;k<nz;k++){
            for (int j=0;j<ny;j++){
                for (int i=0;i<nx;i++){
                    // Initialize distance to +/- 1
                    SignDist(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
                }
            }
        }

        CalcDist(SignDist,id_solid,*Dm);
        MPI_Barrier(comm);

        for (int n = 0; n < nx*ny*nz; n++) SDsID[n] = 1;
        
        MPI_Barrier(comm);
       
        // Initialize the domain and communication
        int xdim,ydim,zdim;
        xdim=Dm->Nx-2;
        ydim=Dm->Ny-2;
        zdim=Dm->Nz-2;
        fillHalo<double> fillData(comm, Dm->rank_info,{xdim,ydim,zdim},{1,1,1},0,1);

        MPI_Barrier(comm);

        double count,countGlobal,totalGlobal;
        double swCount,swGlobal;
        swCount = 0.f;
        count = 0.f;
        double maxdist=-200.f;
        double maxdistGlobal;
        int SolidCount = 0;
	int SolidCountTop = 0;
        int zSolid = 100;
	int zSolidTop = 0;
        int zSolidGlobal = 0;
	int zSolidGlobalTop = 0;

        for (int k=0; k<nz; k++){
            for (int j=0; j<ny; j++){
                for (int i=0; i<nx; i++){
                    n = k*nx*ny+j*nx+i;
                    if ( SignDist(i,j,k) > maxdist) maxdist=SignDist(i,j,k);
                    if (BC > 0 && Dm->kproc() == (nprocz - 1) && k > (nz-2)){
                        SignDist(i,j,k) = 0.0;
                    }
                    if (ReadID[n] > 0){
                        if (FromRaw == false){
                            // initially saturated with wetting phase
                            SDsID[n] = 2;
                        } else if(FromRaw==true && ReadID[n] == 2){
                            SDsID[n] = ReadID[n];
                        }
                        // count+=1.0;
                        // if (SDsID[n] == 2){
                        //     swCount+=1.0;
                        // }
                    }
                }
            }
        }

	// Print ID
        //for (int k=50; k<51; k++){
        //    printf("%i \n",k);
        //    for (int j=0; j<ny; j++){
        //        for (int i=0; i<nx; i++){
        //            n = k*nx*ny+j*nx+i;
        //            printf("%i ",SDsID[n]);
        //        }
        //        printf("\n");
        //    }
        //    printf("\n \n");
        //}


        for (int k=amin; k<amax; k++){
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++){
                    n = k*nx*ny+j*nx+i;
                    if (ReadID[n] > 0){
                        count+=1.0;
                        if (SDsID[n] == 2){
                            swCount+=1.0;
                        }
                    }
                }
            }
        }        

        //Find the lowest part of the lowest solid in z
        if (BC > 0){
            for (int k=1; k<nz-1; k++){
                for (int j=1; j<ny-1; j++){
                    for (int i=1; i<nx-1; i++){
                        n = k*nx*ny+j*nx+i;
                        if (ReadID[n] == 0 && SolidCount == 0 && Dm->kproc() == 0){
                            SolidCount = 1;
                            zSolid = k;
                        }
                    }
                }
            }
	    for (int k=(nz-2); k>0; k--){
		for (int j=1; j<ny-1; j++){
                    for (int i=1; i<nx-1; i++){
                        n = k*nx*ny+j*nx+i;
			if (ReadID[n] == 0 && SolidCountTop == 0 && Dm->kproc() == (nprocz-1)){
                            SolidCountTop = 1;
                            zSolidTop = k;
                        }
		    }
		}
	    }
        }

        MPI_Allreduce(&zSolid,&zSolidGlobal,1,MPI_INT,MPI_MIN,comm);
	MPI_Allreduce(&zSolidTop,&zSolidGlobalTop,1,MPI_INT,MPI_MAX,comm);
        if (zSolidGlobal == 100 || zSolidGlobal < 2) zSolidGlobal = 2;
        if (rank==0){
		printf("Lower bound of lowest solid: %i \n",zSolidGlobal);
		printf("Upper bound of highest solid: %i \n",zSolidGlobalTop);
		printf("amin %i \n", amin);
		printf("amax %i \n", amax);
		if (amin != zSolidGlobal){
			printf("Warning: analysis starts at voxel %i and solid starts at voxel %i \n",amin,zSolidGlobal);
		}
	}

        //if (BC > 0 && Dm->kproc() == 0 && nspheres > 0){
        if (BC > 0 && Dm->kproc() == 0){
	    for (int k=1; k<zSolidGlobal; k++){
                for (int j=1; j<ny-1; j++){
                    for (int i=1; i<nx-1; i++){
                        n = k*nx*ny+j*nx+i;
                        //if (k < zSolidGlobal){
                        SDsID[n] = 1; //add n phase to build off 
                        //}
                    }
                }
            }
            //for (int k=zSolidGlobal; k<(zSolidGlobal+3); k++){
            //    for (int j=1; j<ny-1; j++){
            //        for (int i=1; i<nx-1; i++){
            //            n = k*nx*ny+j*nx+i;
            //            if (ReadID[n] > 0){
            //                SDsID[n] = 1; //add n phase to build off 
            //            }
            //        }
            //    }               
            //}
        }

        PopulateHalo_Char(SDsID, *Dm, nx, ny, nz, rank);

        // for (int k=0; k<nz; k++){
        //     for (int j=0; j<ny; j++){
        //         for (int i=0; i<nx; i++){
        //             n = k*nx*ny+j*nx+i;
        //             SignDist(i,j,k) = (double)SDs[n];
        //         }
        //     }
        // }

        MPI_Barrier(comm);
        // total Global is the number of nodes in the pore-space
        MPI_Allreduce(&count,&totalGlobal,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&maxdist,&maxdistGlobal,1,MPI_DOUBLE,MPI_MAX,comm);
        MPI_Allreduce(&swCount,&swGlobal,1,MPI_DOUBLE,MPI_SUM,comm);
        double volume=double(nprocx*nprocy*nprocz)*double(nx-2)*double(ny-2)*double(nz-2);
        //double porosity=totalGlobal/volume;
        //if (rank==0) printf("Full domain porosity: %f \n",porosity);
        if (rank==0) {
            if (Rcrit_new > 0.0){
                printf("Maximum pore size: %f \n",Rcrit_new);
            } else {
                printf("Maximum pore size: %f \n",maxdistGlobal);
            }
            fflush(stdout);
        }
        Dm->CommInit();
        int iproc = Dm->iproc();
        int jproc = Dm->jproc();
        int kproc = Dm->kproc();
        
        //......................................................................................

        int sendtag,recvtag;
        sendtag = recvtag = 7;
        
        int ii,jj,kk;
        int imin,jmin,kmin,imax,jmax,kmax;
       
        //double sw_new=1.0;
        double sw_new = swGlobal/totalGlobal;
        double sw_diff_new = 1.0;
        double sum;
        double sum_sat;

        // Increase the critical radius until the target saturation is met
        double Rcrit_old;
        if (Rcrit_new <= 0.0){
            Rcrit_new = maxdistGlobal;
        }

        // Print ID
        fflush(stdout);
        if (rank == 0){
            // printf("rank: %i \n",rank);
            // for (int k=10; k<11; k++){
            //     for (int j=0; j<ny; j++){
            //         for (int i=0; i<nx; i++){
            //             n = k*nx*ny+j*nx+i;
            //             printf("%i ",SDsID[n]);
            //         }
            //         printf("\n");
            //     }
            //     printf("\n \n");
            // }
            printf("rank: %i \n",rank);
            for (int k=10; k<11; k++){
                for (int j=0; j<ny; j++){
                    for (int i=0; i<nx; i++){
                        n = k*nx*ny+j*nx+i;
                        printf("%.2f ",SignDist(i,j,k));
                    }
                    printf("\n");
                }
                printf("\n \n");
            }
        }



    double max_inlet = 0.f;
    double max_inlet_global = 0.f;
    if (kproc == 0){
        for (int k=amin; k<(amin+1); k++){
            for (int j=0; j<ny; j++){
                for (int i=0; i<nx; i++){
                    if (SignDist(i,j,k) > max_inlet){
                        max_inlet = SignDist(i,j,k);
                    }
                }
            }
        }
    }

    MPI_Allreduce(&max_inlet,&max_inlet_global,1,MPI_DOUBLE,MPI_MAX,comm);
	if (rank == 0) printf("Recommend min %f inlet \n",max_inlet);
	
	MPI_Barrier(comm);

        FILE *DRAIN = fopen("morphdrain.csv","w");
        fprintf(DRAIN,"sw radius\n");	

        fflush(stdout);
        if (rank == 1){
            printf("rank: %i \n",rank);
            for (int k=10; k<11; k++){
                for (int j=0; j<ny; j++){
                    for (int i=0; i<nx; i++){
                        n = k*nx*ny+j*nx+i;
                        printf("%.2f ",SignDist(i,j,k));
                    }
                    printf("\n");
                }
                printf("\n \n");
            }
        }

       
        for (int aa=0;aa<target_count;aa++){

            SW = SW_List[aa];
            if (rank==0) {
                printf("Running morphological drainage to %f \n",SW);
            }

            fflush(stdout);

            while (sw_new > SW && Rcrit_new > 0.5)
            {
                // if (rank==0) {printf("Rcrit new: %.2f \n",Rcrit_new);}
                Rcrit_old = Rcrit_new;
                Rcrit_new -= deltaR*Rcrit_old;
                int Window=round(Rcrit_new)+1;
                if (Rcrit_new < 15.44 && rank == 0){ 
                    printf("Rcrit_new: %.2f \n",Rcrit_new);
                    printf("Window: %i \n",Window);
                }
                if (Window == 0) Window = 1; // If Window = 0 at the begining, after the following process will have sw=1.0
		        // and sw<Sw will be immediately broken
                double LocalNumber=0.f;
                // for(k=1; k<nz-1; k++){
                //     for(j=1; j<ny-1; j++){
                //         for(i=1; i<nx-1; i++){
                //             n = k*nx*ny + j*nx+i;
                //             if (SignDist(i,j,k) > Rcrit_new){
                //                 // if (kproc==0 && k == 3) printf("ijk: %i %i %i \n",i,j,k);
                //                 // loop over the window and update
                //                 imin=max(0,i-Window);
                //                 jmin=max(0,j-Window);
                //                 kmin=max(0,k-Window);
                //                 imax=min(nx-1,i+Window);
                //                 jmax=min(ny-1,j+Window);
                //                 kmax=min(nz-1,k+Window);
                //                 for (kk=kmin; kk<kmax; kk++){
                //                     for (jj=jmin; jj<jmax; jj++){
                //                         for (ii=imin; ii<imax; ii++){
                //                             int nn = kk*nx*ny+jj*nx+ii;
                //                             double dsq = double((ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k));
                //                             if (SDsID[nn] == 2 && dsq <= (Rcrit_new+1)*(Rcrit_new+1)){
                //                                 LocalNumber+=1.0;
                //                                 SDsID[nn]=1;
                //                             }
                //                         }
                //                     }
                //                 }
                //             }
                //         }
                //     }
                // }

                //Find a point? mark it as "solid"
                for(k=1; k<nz-1; k++){
                    for(j=1; j<ny-1; j++){
                        for(i=1; i<nx-1; i++){
                            if (SignDist(i,j,k) > Rcrit_new){
                                //mark point as 1
                                //id_solid(i,j,k) = 1;
                            } 
                        }
                    }
                }

                for (int k=0;k<nz;k++){
                    for (int j=0;j<ny;j++){
                        for (int i=0;i<nx;i++){
                            //if point not equal to 1, make it 0

                            // Initialize distance to +/- 1
                            //SignDist(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
                        }
                    }
                }

                //How far are global points from marked points
                //CalcDist(SignDist,id_solid,*Dm);

                for(k=1; k<nz-1; k++){
                    for(j=1; j<ny-1; j++){
                        for(i=1; i<nx-1; i++){
                            //if dist to marked pt is smaller than rcrit:
                            //update phase
                        }
                    }
                }


                //......................................................................................

                // Print ID
                // if (Rcrit_new < 15.44 && kproc == 0){
                //     printf("rank: %i \n",rank);
                //     for (int k=10; k<11; k++){
                //         for (int j=0; j<ny; j++){
                //             for (int i=0; i<nx; i++){
                //                 n = k*nx*ny+j*nx+i;
                //                 printf("%i ",SDsID[n]);
                //             }
                //             printf("\n");
                //         }
                //         printf("\n \n");
                //     }
                // }
                // fflush(stdout);

                PopulateHalo_Char(SDsID, *Dm, nx, ny, nz, rank);

                for (int k=0; k<nz; k++){
                    for (int j=0; j<ny; j++){
                        for (int i=0; i<nx; i++){
                            n=k*nx*ny+j*nx+i;
                            if (SDsID[n] == 1){
                                phase(i,j,k) = 1.0;
                            } else {
                                phase(i,j,k) = -1.0;
                            }
                        }
                    }
                }

                // Find connected nonwetting phase
                if (rank == 0) printf("Pre Reduction "); ////Blobs: %i \n (from ComputeBlobalBlobIDs)
                BlobIDstruct new_index;
                double vF=0.0; double vS=0.0;
                ComputeGlobalBlobIDs(nx-2,ny-2,nz-2,Dm->rank_info,phase,SignDist,vF,vS,phase_label,Dm->Comm);
                MPI_Barrier(comm);  

                //If nonperiodic BC, find the blob number of the nonwetting phase inlet
                //If periodic, 1 is the largest blob 
                int target_blob_label = 0;
                int target_blob_label_x = 0;
                int phase_count = 0;

                if (BC > 0){
                    if (kproc == 0){
                        for (int k=1; k<2; k++){
                            for (int j=0; j<(ny); j++){
                                for (int i=0; i<(nx); i++){
                                    n=k*nx*ny+j*nx+i;
                                    if (phase_count == 0 && SDsID[n] == 1 && phase_label(i,j,k) >= 0){
                                        target_blob_label_x = phase_label(i,j,k);
                                        phase_count = phase_count + 1;
                                    } else if (phase_count > 0 && phase_label(i,j,k) != target_blob_label_x && phase_label(i,j,k) > -1){
                                        printf("%i \n",phase_label(i,j,k));
                                    }
                                }
                            }
                        }
                    }
                } else{
                    target_blob_label = 1;
                }

                MPI_Allreduce(&target_blob_label_x,&target_blob_label,1,MPI_INT,MPI_MAX,Dm->Comm);

                //Remove unconnected blobs
                if (BC > 0 && kproc == 0){
                    for (int k=zSolidGlobal+2; k<nz; k++){
                        for (int j=0; j<ny; j++){
                            for (int i=0; i<nx; i++){
                                n=k*nx*ny+j*nx+i;
                                if (SDsID[n] == 1 && phase_label(i,j,k) != target_blob_label && phase_label(i,j,k) != -1){
                                    SDsID[n] = 2;
                                }
                            }
                        }
                    }  
                } else {
                    for (int k=0; k<nz; k++){
                        for (int j=0; j<ny; j++){
                            for (int i=0; i<nx; i++){
                                n=k*nx*ny+j*nx+i;
                                if (SDsID[n] == 1 && phase_label(i,j,k) != target_blob_label && phase_label(i,j,k) != -1){
                                    SDsID[n] = 2;
                                }
                            }
                        }
                    }      
                }

                PopulateHalo_Char(SDsID, *Dm, nx, ny, nz, rank);

                // // Print ID
                // if (Rcrit_new < 15.44 && kproc == 0){
                //     printf("rank: %i \n",rank);
                //     for (int k=10; k<11; k++){
                //         for (int j=0; j<ny; j++){
                //             for (int i=0; i<nx; i++){
                //                 n = k*nx*ny+j*nx+i;
                //                 printf("%i ",SDsID[n]);
                //             }
                //             printf("\n");
                //         }
                //         printf("\n \n");
                //     }
                // }
                // fflush(stdout);

                //2=wetting,1=nonwetting,0=solid,MD=replacing 2's with 1's
                for (int n = 0; n < nx*ny*nz; n++) {
                    if (ScanningCurve == false && ReadID[n] > 0){
                        ReadID[n] = SDsID[n];
                    } else if (ScanningCurve == true && ReadID[n] == 2){
                        ReadID[n] = SDsID[n];
                    }
                }
                
                //PopulateHalo_Char(SDsID, *Dm, nx, ny, nz, rank);        

                sum_sat = 0.f;
                sum = 0.f;
                count = 0.f;
                double porecount = 0.f;

                for (int k=amin; k<amax; k++){
                    for (int j=1; j<ny-1; j++){
                        for (int i=1; i<nx-1; i++){
                            n = k*nx*ny+j*nx+i;
                            if (ReadID[n] == 2)  {count += 1.0;}
                            if (ReadID[n] > 0)  {porecount += 1.0;}
                        }
                    }
                }


                MPI_Allreduce(&porecount,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
                MPI_Allreduce(&count,&sum_sat,1,MPI_DOUBLE,MPI_SUM,comm);
                sw_new = sum_sat/sum;

                if (rank==0){
                    fprintf(DRAIN,"%f ",sw_new);
                    fprintf(DRAIN,"%f\n",Rcrit_new);
                    printf("     %f ",sw_new);
                    printf("     %f\n",Rcrit_new);
                    fflush(stdout);
                }
            }
            
            //2=wetting,1=nonwetting,0=solid,MD=replacing 2's with 1's
            for (int n = 0; n < nx*ny*nz; n++) {
                if (ScanningCurve == false && ReadID[n] > 0){
                    ReadID[n] = SDsID[n];
                } else if (ScanningCurve == true && ReadID[n] == 2){
                    ReadID[n] = SDsID[n];
                }
            }


	                    // Print ID
                //for (int k=5; k<8; k++){
                //        printf("%i \n",k);
                //        for (int j=0; j<ny; j++){
                //                for (int i=0; i<nx; i++){
                //                        n = k*nx*ny+j*nx+i;
                //                        printf("%i ",ReadID[n]);
                //                }
                //                printf("\n");
                //        }
                //        printf("\n \n");
                //}


            // if (rank==0){
            //     printf("Final saturation=%f\n",sw_new);
            //     printf("Final critical radius=%f\n",Rcrit_new);
            // }
            
            for (int k=0;k<nz;k++){
                for (int j=0;j<ny;j++){ 
                    for (int i=0;i<nx;i++){
                        int n = k*nx*ny+j*nx+i;
                        Mask->id[n] = ReadID[n];
                        // if (ReadID[n] > 2){
                        //     if (rank==0){
                        //         printf ("Error C %i \n",ReadID[n]);
                        //     } 
                        // } else if (ReadID[n] < 0){
                        //     if (rank==0){
                        //         printf ("Error D %i \n",ReadID[n]);
                        //     } 
                        // }
                    }
                }
            }
            MPI_Barrier(comm);

            auto filename2 = "SW." + to_string(SW) + ".morphdrain.raw";
            if (rank==0) printf("Writing file to: %s \n", filename2.data() );
            Mask->AggregateLabels( filename2 );

            if (target_count > 1){
                if (rank==0) printf("Saving result %f to folder %f\n",sw_new,SW);
                snprintf(LocalRankFilenameSub,sizeof(LocalRankFilenameSub),LocalRankFilenameSub_format,SW,rank);
                FILE *IDSub = fopen(LocalRankFilenameSub,"wb"); // w: write new file, b: in binary
                fwrite(ReadID,1,N,IDSub);
                fclose(IDSub);
            } else {
                sprintf(LocalRankFilename,"ID2.%05i",rank);
                FILE *ID = fopen(LocalRankFilename,"wb");
                fwrite(ReadID,1,N,ID);
                fclose(ID);
            }

            if (rank ==0) {
                printf("Ending SW %f Target\n",SW);
            }

        }
    } 
}

void lbpm_morphimb_pp(int argc, char **argv,int rank, int nprocs, MPI_Comm comm)
{
    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int n, nprocx, nprocy, nprocz;

    char LocalRankFilename[40];
    char LocalRankFoldername[40];

    char LocalRankFilenameSub_format[] = "SW%03f/ID2.%05i"; 
    char LocalRankFilenameSub[sizeof(LocalRankFilenameSub_format)+20];

    string filename;
    double SW;
    if (argc > 1){
        filename=argv[1];
    }
    else ERROR("No input database provided\n");
    // read the input database 
    auto db = std::make_shared<Database>( filename );
    auto domain_db = db->getDatabase( "Domain" );

    // Read domain parameters
    auto size = domain_db->getVector<int>( "n" );
    auto nproc = domain_db->getVector<int>( "nproc" );
    auto ReadValues = domain_db->getVector<int>( "ReadValues" );
    auto WriteValues = domain_db->getVector<int>( "WriteValues" );
    auto BC = domain_db->getScalar<int>("BC");
    auto READFILE = domain_db->getScalar<std::string>( "Filename" );

    int target_count = 0;
    double *SW_List;
    double sw_new;

    double Rcrit_new = 0.f;
    double deltaR = 0.05;

    if (domain_db->keyExists( "Rcrit" )){
        Rcrit_new = domain_db->getScalar<double>("Rcrit");
    }
    if (domain_db->keyExists( "deltaR")){
        deltaR = domain_db->getScalar<double>("deltaR");
    } 

    //Read in single sw target from RandomSaturation or multiple from targets.in
    //If reading multiple targets, create directories to store results

    if (domain_db->keyExists( "RandomSaturation" )){
        SW = domain_db->getScalar<double>("RandomSaturation");
        if (rank == 0) printf("Single Saturation Target: %f \n",SW);
    } else {
        double value;
        FILE * pFile;
        pFile = fopen ("targets.in","r");
        char line[ 100 ];
        while (fgets(line,100,pFile)) {
            target_count += 1;
        }
        fclose (pFile);

        pFile = fopen ("targets.in","r");
        SW_List = new double[target_count];
        target_count = 0;
        while (fgets(line,100,pFile)) {
            value = strtod(line,NULL);
            SW_List[target_count] = value;
            target_count += 1;
        }
        fclose (pFile);

        if (rank == 0) {
            printf("Quantity of sw targets: %i \n",target_count);
            sort(SW_List,SW_List+target_count);
            for (int i=0;i<target_count;i++){
                sprintf(LocalRankFoldername,"SW%03f",SW_List[i]);
                mkdir(LocalRankFoldername,S_IRWXU|S_IRGRP);
            }
            printf("\n");
        }

    }
    
    if (target_count == 0) {
        target_count = 1;
        SW_List = new double[target_count];
        SW_List[0] = SW;
    }

    auto filename2 = "SW." + to_string(SW) + ".morphimb.raw";

    std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
    std::shared_ptr<Domain> Mask (new Domain(domain_db,comm));
    
    int nx = Dm->Nx; int ny = Dm->Ny;  int nz = Dm->Nz;
    nprocx = Dm->nprocx();  nprocy = Dm->nprocy();  nprocz = Dm->nprocz();
    size_t N = size_t(nx)*size_t(ny)*size_t(nz);

    int amin = Dm->amin;
    int amax = Dm->amax;
    // int amin = 1;
    // int amax = nz-1;
    // if (BC > 0){
    //     if ( domain_db->keyExists( "AnalysisMin" ) && Dm->kproc() == 0){
    //         amin = domain_db->getScalar<int>( "AnalysisMin" );
    //     } 
    //     if ( domain_db->keyExists( "AnalysisMax" ) && Dm->kproc() == (nprocz - 1)){
    //         amax = domain_db->getScalar<int>( "AnalysisMax" );
    //         amax = nz-1-((nz-2)*nprocz - amax);
    //     }
    // }

    for (n=0; n<N; n++) Dm->id[n]=1;
    Dm->CommInit();

    char *id;
    id = new char [N];
    char *id_connected;
    id_connected = new char [N];

    Array<char> id_solid(nx,ny,nz);
    DoubleArray SignDist(nx,ny,nz);
    DoubleArray phase(nx,ny,nz);
    IntArray phase_label(nx,ny,nz);

    double sum_sat = 0.f;
    double sum = 0.f;
    double count = 0.f;
    double porecount = 0.f;

    for (int aa=0;aa<target_count;aa++){
        
        if (target_count > 1) SW = SW_List[aa];
        if (aa > 0){
            filename2 = "SW." + to_string(SW_List[aa-1]) + ".morphimb.raw";
            READFILE = filename2;
        }

        //if (rank==0) printf("Performing morphological imbibition with target saturation %f \n", SW);

        Mask->Decomp(READFILE,db);
        Mask->CommInit();

       for (int k=0; k<nz; k++){
            for (int j=0; j<ny; j++){
                for (int i=0; i<nx; i++){
                    n=k*nx*ny+j*nx+i;
                    id[n] = Mask->id[n];
                    if (BC > 0 && Dm->kproc() == 0 && k < 5){
                        id[n] = 1; //add strip of n phase to work from
                    }
                }
            }
       }
       
        PopulateHalo_Char(id, *Dm, nx, ny, nz, rank);


        //Initial Saturation
        sum_sat = 0.f;
        sum = 0.f;
        count = 0.f;
        porecount = 0.f;

        for (int k=amin; k<amax; k++){
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++){
                    n = k*nx*ny+j*nx+i;
                    if (id[n] == 2)  {count += 1.0;}
                    if (id[n] > 0)  {porecount += 1.0;}
                }
            }
        }


        MPI_Allreduce(&porecount,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&count,&sum_sat,1,MPI_DOUBLE,MPI_SUM,comm);
        sw_new = sum_sat/sum;

        if (rank==0) printf("Performing morphological imbibition from saturation %f to target saturation %f from file %s \n", sw_new, SW, READFILE.c_str());
        
        // Generate the signed distance map
        // Initialize the domain and communication
        for (int k=0; k<nz; k++){
            for (int j=0; j<ny; j++){
                for (int i=0; i<nx; i++){
                    n=k*nx*ny+j*nx+i;
                    //id[n] = Mask->id[n];
                    // if (BC > 0 && Dm->kproc() == 0 && k < 5){
                    //     id[n] = 1; //add strip of n phase to work from
                    // }
                    if (id[n] == 1){
                        phase(i,j,k) = 1.0;
                    }
                    else
                        phase(i,j,k) = -1.0;
                }
            }
        }
        
        // Solve for the position of the solid phase
        for (int k=0;k<nz;k++){
            for (int j=0;j<ny;j++){
                for (int i=0;i<nx;i++){
                    int n = k*nx*ny+j*nx+i;
                    // Initialize the solid phase
                    if (id[n] > 0){
                        id_solid(i,j,k) = 1;
                    }
                    else	    
                        id_solid(i,j,k) = 0;
                }
            }
        }
        // Initialize the signed distance function
        for (int k=0;k<nz;k++){
            for (int j=0;j<ny;j++){
                for (int i=0;i<nx;i++){
                    // Initialize distance to +/- 1
                    SignDist(i,j,k) = 2.0*double(id_solid(i,j,k))-1.0;
                }
            }
        }

        CalcDist(SignDist,id_solid,*Dm);
        MPI_Barrier(comm);

        //Set SD at zmax to zero to prevent wrap around effects
        if (BC > 0 && Dm->kproc() == (nprocz - 1)){
            for (int k=(nz-1); k<nz; k++){
                for (int j=0; j<ny; j++){
                    for (int i=0; i<nx; i++){
                        n = k*nx*ny+j*nx+i;
                        if ( SignDist(i,j,k) > 0.0 ){
                            SignDist(i,j,k) = 0.0; 
                        }
                    }
                }
            }
        }
        MPI_Barrier(comm);

        // Extract only the connected part of NWP
        BlobIDstruct new_index;
        double vF=0.0; double vS=0.0;
        if (rank == 0) printf("Blob count pre morph imb: ");
        ComputeGlobalBlobIDs(nx-2,ny-2,nz-2,Dm->rank_info,phase,SignDist,vF,vS,phase_label,Dm->Comm);
        MPI_Barrier(Dm->Comm);

        int target_blob_label = 0;
        int target_blob_label_x = 0;
        int phase_count = 0;
        int phase_count_total;

        if (BC > 0){
            if (Dm->kproc() == 0 && Dm->iproc()==0 && Dm->jproc()==0){
                for (int k=1; k<2; k++){
                    for (int j=0; j<(ny); j++){
                        for (int i=0; i<(nx); i++){
                            n=k*nx*ny+j*nx+i;
                            if (phase_count == 0 && id[n] == 1 && phase_label(i,j,k) >= 0){
                                target_blob_label_x = phase_label(i,j,k);
                                phase_count = phase_count + 1;
                            }
                        }
                    }
                }
            }
        } else{
            target_blob_label = 1;  
        }

        MPI_Allreduce(&target_blob_label_x,&target_blob_label,1,MPI_INT,MPI_SUM,Dm->Comm);
        MPI_Allreduce(&phase_count,&phase_count_total,1,MPI_INT,MPI_SUM,Dm->Comm);  


        int count_connected=0;
        int count_porespace=0;
        int count_water=0;
        for (int k=amin; k<amax; k++){
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++){
                    n=k*nx*ny+j*nx+i;
                    // only apply opening to connected component 
                    if ( phase_label(i,j,k) == target_blob_label){
                        count_connected++;
                    }
                    if (id[n] > 0){
                        count_porespace++;
                    }
                    if (id[n] == 2){
                        count_water++;
                    }
                }
            }
        }
        count_connected=sumReduce( Dm->Comm, count_connected);
        count_porespace=sumReduce( Dm->Comm, count_porespace);
        count_water=sumReduce( Dm->Comm, count_water);
        
        for (int k=0; k<nz; k++){
            for (int j=0; j<ny; j++){
                for (int i=0; i<nx; i++){
                    n=k*nx*ny+j*nx+i;
                    // only apply opening to connected component 
                    if ( phase_label(i,j,k) == target_blob_label){
                        id_solid(i,j,k) = 1;
                        id_connected[n] = 2;
                        id[n] = 2;
                    }
                    else{
                        id_solid(i,j,k) = 0;
                        id_connected[n] = 0;
                    }
                }
            }
        }

        CalcDist(SignDist,id_solid,*Dm);
        //Put SD at zmax back to zero to prevent wrap around effects
        if (BC > 0 && Dm->kproc() == (nprocz - 1)){
            for (int k=(nz-2); k<nz; k++){
                for (int j=0; j<ny; j++){
                    for (int i=0; i<nx; i++){
                        n = k*nx*ny+j*nx+i;
                        if ( SignDist(i,j,k) > 0.0 ){
                            SignDist(i,j,k) = 0.0; 
                        }
                    }
                }
            }
        }

        int unconnected_nw = count_porespace - count_water - count_connected;
        int target_volume = (1-SW)*count_porespace - unconnected_nw;
        //double St = (SW*count_porespace - count_water)/count_porespace; 
        double St = (SW*count_porespace - count_water);

        if (St > double(count_connected)){
            St = 1.0;
        }

        if(rank==0){
            printf("Nonwetting volume connected available: %i \n",count_connected);
            //printf("Nonwetting target connected volume: %i \n",target_volume);
            //printf("Total porespace volume: %i \n",count_porespace);
            //printf("Count water: %i \n",count_water);
            //printf("St %f \n",St);
        }  
        
        //if (rank==0) printf("Here 1A \n");

        fflush(stdout);
        
        char water=2;
        char notwater=1;
        // Run the morphological opening
        if (phase_count_total > 0 && St > 0.0) {
            if (rank==0) printf("St: %f \n",St);
            MorphOpen(SignDist, id_connected, id, Dm, St, water, notwater, amin, amax, deltaR, Rcrit_new, count_connected);
        } else {
            if(rank==0) printf("Initial condition satisfies condition for saturation target or imbibition cannot proceed further \n");
        }
        
        // re-label 
        for (int k=0; k<nz; k++){
            for (int j=0; j<ny; j++){
                for (int i=0; i<nx; i++){
                    n=k*nx*ny+j*nx+i;
                    // only apply opening to connected component 
                    if ( id_connected[n] == 1){
                        id[n] = 1;
                    }
                }
            }
        }
            

        double sum_sat = 0.f;
        double sum = 0.f;
        double count = 0.f;
        double porecount = 0.f;

        //for (int k=1; k<nz-1; k++){
        for (int k=amin; k<amax; k++){
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++){
                    n = k*nx*ny+j*nx+i;
                    if (id[n] == 2)  {count += 1.0;}
                    if (id[n] > 0)  {porecount += 1.0;}
                }
            }
        }


        MPI_Allreduce(&porecount,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&count,&sum_sat,1,MPI_DOUBLE,MPI_SUM,comm);
        sw_new = sum_sat/sum;

        if(rank==0) printf("Final saturation: %f \n", sw_new);

        //if (rank==0) printf("Writing ID file \n");

        //write the geometry to a single file
        for (int k=0;k<nz;k++){
            for (int j=0;j<ny;j++){ 
                for (int i=0;i<nx;i++){
                    int n = k*nx*ny+j*nx+i;
                    Mask->id[n] = id[n];
                }
            }
        }
        MPI_Barrier(comm);       

        auto filename2 = "SW." + to_string(SW) + ".morphimb.raw";
        //if (rank==0) printf("Writing file to: %s \n", filename2.c_str());
        Mask->AggregateLabels(filename2);


        if (target_count > 1){
            if (rank==0) printf("Saving result %f to folder %f\n",sw_new,SW);
            snprintf(LocalRankFilenameSub,sizeof(LocalRankFilenameSub),LocalRankFilenameSub_format,SW,rank);
            FILE *IDSub = fopen(LocalRankFilenameSub,"wb"); // w: write new file, b: in binary
            fwrite(id,1,N,IDSub);
            fclose(IDSub);
        } else {
            sprintf(LocalRankFilename,"ID2.%05i",rank);
            FILE *ID = fopen(LocalRankFilename,"wb");
            fwrite(id,1,N,ID);
            fclose(ID);
        }

        if (rank ==0) {
            printf("Ending SW %f Target\n",SW);
        }

        //MPI_Barrier(comm);
    }
}
