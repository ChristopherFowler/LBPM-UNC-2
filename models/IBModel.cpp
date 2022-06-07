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
 color lattice boltzmann model
 */
#include "models/IBModel.h"
#include "analysis/distance.h"
#include <vector>
#include <string>


inline void MeanFilter(Array<double> & tmp, DoubleArray &Mesh) {}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsometimes-uninitialized"
#pragma GCC diagnostic ignored "-Wreorder"

ScaLBL_ColorModel::ScaLBL_ColorModel(int RANK, int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tauA(0),tauB(0),rhoA(0),rhoB(0),alpha(0),beta(0),
Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),inletA(0),inletB(0),outletA(0),outletB(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),comm(COMM)
{
    
}
#pragma GCC diagnostic pop
ScaLBL_ColorModel::~ScaLBL_ColorModel(){}


void ScaLBL_ColorModel::ReadParams(string filename) {}
void ScaLBL_ColorModel::SetDomain() {}

bool ScaLBL_ColorModel::GridIndependence() { return false; }

bool ScaLBL_ColorModel::ReadInputAndValidate() { return false; }

void ScaLBL_ColorModel::ReadInput() {}

void ScaLBL_ColorModel::AssignComponentLabels(double *phase) {}


void ScaLBL_ColorModel::Create() {}

/********************************************************
 * AssignComponentLabels                                 *
 ********************************************************/

void ScaLBL_ColorModel::Initialize() {}

std::string getLastLine(std::ifstream& in)
{
    std::string line;
    while (in >> std::ws && std::getline(in, line)) // skip empty lines
        ;

    return line;
}

void ScaLBL_ColorModel::Run() {}

double ScaLBL_ColorModel::MorphInit(const double beta, const double morph_delta){
    const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
    
    double vF = 0.f;
    double vS = 0.f;
    
    DoubleArray phase(Nx,Ny,Nz);
    IntArray phase_label(Nx,Ny,Nz);;
    DoubleArray phase_distance(Nx,Ny,Nz);
    Array<char> phase_id(Nx,Ny,Nz);
    
    // Basic algorithm to
    // 1. Copy phase field to CPU
    ScaLBL_CopyToHost(phase.data(), Phi, N*sizeof(double));
    
    double count,count_global,volume_initial,volume_final;
    count = 0.f;
    for (int k=0; k<Nz; k++){
        for (int j=0; j<Ny; j++){
            for (int i=0; i<Nx; i++){
                if (phase(i,j,k) > 0.f && Averages->SDs(i,j,k) > 0.f) count+=1.f;
            }
        }
    }
    MPI_Allreduce(&count,&count_global,1,MPI_DOUBLE,MPI_SUM,comm);
    volume_initial = count_global;
    
    // 2. Identify connected components of phase field -> phase_label
    BlobIDstruct new_index;
    ComputeGlobalBlobIDs(Nx-2,Ny-2,Nz-2,rank_info,phase,Averages->SDs,vF,vS,phase_label,comm);
    MPI_Barrier(comm);
    // only operate on component "0"
    for (int k=0; k<Nz; k++){
        for (int j=0; j<Ny; j++){
            for (int i=0; i<Nx; i++){
                int label = phase_label(i,j,k);
                if (label == 0 )     phase_id(i,j,k) = 0;
                else 		     phase_id(i,j,k) = 1;
            }
        }
    }
    // 3. Generate a distance map to the largest object -> phase_distance
    CalcDist(phase_distance,phase_id,*Dm);
    
    double temp,value;
    double factor=0.5/beta;
    for (int k=0; k<Nz; k++){
        for (int j=0; j<Ny; j++){
            for (int i=0; i<Nx; i++){
                if (phase_distance(i,j,k) < 3.f ){
                    value = phase(i,j,k);
                    if (value > 1.f)   value=1.f;
                    if (value < -1.f)  value=-1.f;
                    // temp -- distance based on analytical form McClure, Prins et al, Comp. Phys. Comm.
                    temp = -factor*log((1.0+value)/(1.0-value));
                    /// use this approximation close to the object
                    if (fabs(value) < 0.8 && Averages->SDs(i,j,k) > 1.f ){
                        phase_distance(i,j,k) = temp;
                    }
                }
            }
        }
    }
    
    // 4. Apply erosion / dilation operation to phase_distance
    for (int k=0; k<Nz; k++){
        for (int j=0; j<Ny; j++){
            for (int i=0; i<Nx; i++){
                double walldist=Averages->SDs(i,j,k);
                double wallweight = 1.f / (1+exp(-5.f*(walldist-1.f)));
                phase_distance(i,j,k) -= wallweight*morph_delta;
            }
        }
    }
    
    // 5. Update phase indicator field based on new distnace
    for (int k=0; k<Nz; k++){
        for (int j=0; j<Ny; j++){
            for (int i=0; i<Nx; i++){
               
                double d = phase_distance(i,j,k);
                if (Averages->SDs(i,j,k) > 0.f){
                    // only update phase field in immediate proximity of largest component
                    if (d < 3.f){
                        phase(i,j,k) = (2.f*(exp(-2.f*beta*d))/(1.f+exp(-2.f*beta*d))-1.f);
                    }
                }
            }
        }
    }
    
    count = 0.f;
    for (int k=0; k<Nz; k++){
        for (int j=0; j<Ny; j++){
            for (int i=0; i<Nx; i++){
                if (phase(i,j,k) > 0.f && Averages->SDs(i,j,k) > 0.f) count+=1.f;
            }
        }
    }
    MPI_Allreduce(&count,&count_global,1,MPI_DOUBLE,MPI_SUM,comm);
    volume_final=count_global;
    
    double delta_volume = (volume_final-volume_initial);
    if (rank == 0)  printf("MorphInit: change fluid volume fraction by %f \n", delta_volume/volume_initial);
    
    // 6. copy back to the device
    //if (rank==0)  printf("MorphInit: copy data  back to device\n");
    ScaLBL_CopyToDevice(Phi,phase.data(),N*sizeof(double));
    
    // 7. Re-initialize phase field and density
    ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0, ScaLBL_Comm->LastExterior(), Np);
    ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np);
    if (BoundaryCondition >0 ){
        if (Dm->kproc()==0){
            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
            ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
        }
        if (Dm->kproc() == nprocz-1){
            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-1);
            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-2);
            ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-3);
        }
    }
    return delta_volume;
}

void ScaLBL_ColorModel::WriteDebug() {}
