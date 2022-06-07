/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
Implementation of color lattice boltzmann model
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/Communication.h"
#include "analysis/TwoPhase.h"
#include "analysis/runAnalysis.h"
#include "common/MPI_Helpers.h"
#include "ProfilerApp.h"
#include "threadpool/thread_pool.h"

class ScaLBL_ColorModel{
public:
    ScaLBL_ColorModel(int RANK, int NP, MPI_Comm COMM);
    ~ScaLBL_ColorModel();
    
    // functions in they should be run
    void ReadParams(string filename);
    void ReadParams(std::shared_ptr<Database> db0);
    void SetDomain();
    void ReadInput();
    bool ReadInputAndValidate();
    void Create();
    bool GridIndependence();
    void Initialize();
    void Run(string filename);
    void WriteDebug();
    
    bool legacySDs;
    bool Restart,pBC;
    int timestep,timestepMax;
    int BoundaryCondition;
    double tauA,tauB,rhoA,rhoB,alpha,beta;
    double Fx,Fy,Fz,flux;
    double din,dout,inletA,inletB,outletA,outletB;
    double input_angle;
    int WBCFlag;
    
    int Nx,Ny,Nz,N,Np,Ni,Nsb;
    int rank,nprocx,nprocy,nprocz,nprocs;
    double Lx,Ly,Lz;
    int method;

    std::shared_ptr<Domain> Dm;   // this domain is for analysis
    std::shared_ptr<Domain> Mask; // this domain is for lbm
    std::shared_ptr<Domain> blank;
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm;
    std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm_Regular;
    std::shared_ptr<TwoPhase> Averages;
    
    // input database
    std::shared_ptr<Database> db;
    std::shared_ptr<Database> domain_db;
    std::shared_ptr<Database> color_db;
    std::shared_ptr<Database> analysis_db;

    IntArray InactiveTemporaryMap;
    IntArray TemporaryMap;
    IntArray Map;
    IntArray ActiveMap;
    char *id;
    int *NeighborList;
    int *InterpolationList;

    
    double IIprimevwndnn;
    
    char * MaskDomain;
    double * LIBBqA;
    double *LIBBqBC;
    double *LIBBqD;
    double * NormalToSolid_X;
    double * NormalToSolid_Y;
    double * NormalToSolid_Z;
    double * GradSDsX;
    double * GradSDsY;
    double * GradSDsZ;
    double * CField;
    
    double * gradphix;
    double * gradphiy;
    double * gradphiz;
    
    double * GradPhiX;
    double * GradPhiY;
    double * GradPhiZ;
    
    double * GradPhiX2;
    double * GradPhiY2;
    double * GradPhiZ2;
    
    int * InactiveScalarList;
    int * SBScalarList;
    int * ActiveScalarList;
    
    int * dvcInactiveMap;
    int * dvcSBMap;
    
    int * activeMap;
    int * inactiveMap;
    int * sbMap;
    
    double * LIBBqA_SBB;
    double *LIBBqBC_SBB;
    double *LIBBqD_SBB;
    int *dvcMap;
    double *fq, *Aq, *Bq, *Aq2, *Bq2, *Qdist;
    double * StencilCoefs;
    double * fq2;
    double * savedfq;
    double * savedAq, * savedBq;
    double * vfmask;
    double * VFmask;
    double * DenA, *DenB, *DenA2, *DenB2;
    double * dvcSignDist;
    double * SignDist;
    double *Den, *Phi, *Phi2;
    double *ColorGrad;
    double *Velocity;
    double *Velocity2;
    
    double *Velx;
    double * Vely;
    double * Velz;
    double * Velx2;
    double * Vely2;
    double * Velz2;
    double *Pressure;
    double *TimeAveragedVelocity;
    int * Xi;
    double * force;
    
private:
    MPI_Comm comm;
    
    int dist_mem_size;
    int neighborSize;
    // filenames
    char LocalRankString[8];
    char LocalRankFilename[40];
    char LocalRestartFile[40];
   
    //int rank,nprocs;
    void LoadParams(std::shared_ptr<Database> db0);
    void AssignComponentLabels(double *phase);
    double MorphInit(const double beta, const double morph_delta);
        
};

