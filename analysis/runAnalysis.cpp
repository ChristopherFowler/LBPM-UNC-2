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
#include "analysis/runAnalysis.h"
#include "analysis/analysis.h"
#include "common/Array.h"
#include "common/Communication.h"
#include "common/MPI_Helpers.h"
#include "common/ScaLBL.h"


#include "IO/MeshDatabase.h"
#include "threadpool/thread_pool.h"

#include "ProfilerApp.h"


AnalysisType& operator |=(AnalysisType &lhs, AnalysisType rhs)
{
    lhs = static_cast<AnalysisType>(
                                    static_cast<std::underlying_type<AnalysisType>::type>(lhs) |
                                    static_cast<std::underlying_type<AnalysisType>::type>(rhs)
                                    );
    return lhs;
}
bool matches( AnalysisType x, AnalysisType y )
{
    return ( static_cast<std::underlying_type<AnalysisType>::type>(x) &
            static_cast<std::underlying_type<AnalysisType>::type>(y) ) != 0;
}


template<class TYPE>
void DeleteArray( const TYPE *p )
{
    delete [] p;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsometimes-uninitialized"
#pragma GCC diagnostic ignored "-Wreorder"
// Helper class to write the restart file from a seperate thread
class WriteRestartWorkItem: public ThreadPool::WorkItemRet<void>
{
public:
    WriteRestartWorkItem( const char* filename_, std::shared_ptr<double> cDenA_,std::shared_ptr<double> cDenB_, std::shared_ptr<double> cfq_,  std::shared_ptr<double> cPhase_,
                         std::shared_ptr<double> cVx_,
                         std::shared_ptr<double> cVy_,
                         std::shared_ptr<double> cVz_,
                         int Np_, int N_, double Fz_,
                         std::shared_ptr<double> cGPhiX_,
                         std::shared_ptr<double> cGPhiY_,
                         std::shared_ptr<double> cGPhiZ_,
                         std::shared_ptr<double> cCField_):
    filename(filename_), cDenA(cDenA_),cDenB(cDenB_), cfq(cfq_), cPhase(cPhase_), cVx(cVx_), cVy(cVy_), cVz(cVz_), Np(Np_), N(N_), Fz(Fz_), cGPhiX(cGPhiX_), cGPhiY(cGPhiY_), cGPhiZ(cGPhiZ_),cCField(cCField_)  {}
    virtual void run() {
        PROFILE_START("Save Checkpoint",1);
        double value;
        ofstream File(filename,ios::binary);
        
        for (size_t n=0; n<N; n++){
            value = cDenA.get()[n];
            File.write((char*) &value, sizeof(value));
        }
        for (size_t n=0; n<N; n++){
            value = cDenB.get()[n];
            File.write((char*) &value, sizeof(value));
        }
        for (size_t n=0; n<N; n++) {
            value = cPhase.get()[n];
            File.write((char*) &value, sizeof(value));
        }
        for (int n=0; n<Np; n++){
            for (int q=0; q<19; q++){
                size_t element = q*Np+n;
                value = cfq.get()[element];
                File.write((char*) &value, sizeof(value));
            }
        }
        for (size_t n=0; n<N; n++) {
            value = cVx.get()[n];
            File.write((char*) &value, sizeof(value));
        }
        for (size_t n=0; n<N; n++) {
            value = cVy.get()[n];
            File.write((char*) &value, sizeof(value));
        }
        for (size_t n=0; n<N; n++) {
            value = cVz.get()[n];
            File.write((char*) &value, sizeof(value));
        }
        
        value = Fz;
        File.write((char*) &value, sizeof(value));
        
        for (size_t n=0; n<N; n++) {
            value = cGPhiX.get()[n];
            File.write((char*) &value, sizeof(value));
        }
        for (size_t n=0; n<N; n++) {
            value = cGPhiY.get()[n];
            File.write((char*) &value, sizeof(value));
        }
        for (size_t n=0; n<N; n++) {
            value = cGPhiZ.get()[n];
            File.write((char*) &value, sizeof(value));
        }
        
        for (size_t n=0; n<N; n++) {
            value = cCField.get()[n];
            File.write((char*) &value, sizeof(value));
        }
       
        File.close();
        PROFILE_STOP("Save Checkpoint",1);
    };
private:
    WriteRestartWorkItem();
    const char* filename;
    std::shared_ptr<double> cDenA,cDenB,cfq,cPhase,cVx,cVy,cVz,cGPhiX,cGPhiY,cGPhiZ,cCField;
    const int Np;
    const int N;
    const double Fz;
};
#pragma GCC diagnostic pop


// Helper class to compute the blob ids
static const std::string id_map_filename = "lbpm_id_map.txt";
class BlobIdentificationWorkItem1: public ThreadPool::WorkItemRet<void>
{
public:
    /*BlobIdentificationWorkItem1( int timestep_, int Nx_, int Ny_, int Nz_, const RankInfoStruct& rank_info_,
                                std::shared_ptr<const DoubleArray> phase_, const DoubleArray& dist_,
                                BlobIDstruct last_id_, BlobIDstruct new_index_, BlobIDstruct new_id_, BlobIDList new_list_, runAnalysis::commWrapper&& comm_ ):*/
    /*timestep(timestep_), Nx(Nx_), Ny(Ny_), Nz(Nz_), rank_info(rank_info_),
    phase(phase_), dist(dist_), last_id(last_id_), new_index(new_index_), new_id(new_id_), new_list(new_list_), comm(std::move(comm_))*/
    BlobIdentificationWorkItem1( int timestep_, int Nx_, int Ny_, int Nz_, const RankInfoStruct& rank_info_,
            const DoubleArray& SDN_, const DoubleArray& dist_,
            BlobIDstruct last_id_, BlobIDstruct new_index_, BlobIDstruct new_id_, BlobIDList new_list_, runAnalysis::commWrapper&& comm_):
                timestep(timestep_), Nx(Nx_), Ny(Ny_), Nz(Nz_), rank_info(rank_info_),
                SDN(SDN_), dist(dist_), last_id(last_id_), new_index(new_index_), new_id(new_id_), new_list(new_list_), comm(std::move(comm_))
    {
    }
    ~BlobIdentificationWorkItem1() { }
    virtual void run() {
        // Compute the global blob id and compare to the previous version
        PROFILE_START("Identify blobs",1);
        double vF = 0.0;
        double vS = -1.0; // one voxel buffer region around solid
        IntArray& ids = new_index->second;
        //new_index->first = ComputeGlobalBlobIDs(Nx-2,Ny-2,Nz-2,rank_info,*phase,dist,vF,vS,ids,comm.comm);
        new_index->first = ComputeGlobalBlobIDs(Nx-2,Ny-2,Nz-2,rank_info,SDN,dist,vF,vS,ids,comm.comm);
        PROFILE_STOP("Identify blobs",1);
    }
private:
    BlobIdentificationWorkItem1();
    int timestep;
    int Nx, Ny, Nz;
    const RankInfoStruct& rank_info;
    //std::shared_ptr<const DoubleArray> phase;
    const DoubleArray& SDN;
    const DoubleArray& dist;
    BlobIDstruct last_id, new_index, new_id;
    BlobIDList new_list;
    runAnalysis::commWrapper comm;
};
class BlobIdentificationWorkItem2: public ThreadPool::WorkItemRet<void>
{
public:
    /*BlobIdentificationWorkItem2( int timestep_, int Nx_, int Ny_, int Nz_, const RankInfoStruct& rank_info_,
                                std::shared_ptr<const DoubleArray> phase_, const DoubleArray& dist_,
                                BlobIDstruct last_id_, BlobIDstruct new_index_, BlobIDstruct new_id_, BlobIDList new_list_ , runAnalysis::commWrapper&& comm_ ):
    timestep(timestep_), Nx(Nx_), Ny(Ny_), Nz(Nz_), rank_info(rank_info_),
    phase(phase_), dist(dist_), last_id(last_id_), new_index(new_index_), new_id(new_id_), new_list(new_list_), comm(std::move(comm_))*/
    BlobIdentificationWorkItem2( int timestep_, int Nx_, int Ny_, int Nz_, const RankInfoStruct& rank_info_,
                                const DoubleArray& SDN_, const DoubleArray& dist_,
                                BlobIDstruct last_id_, BlobIDstruct new_index_, BlobIDstruct new_id_, BlobIDList new_list_ , runAnalysis::commWrapper&& comm_):
    timestep(timestep_), Nx(Nx_), Ny(Ny_), Nz(Nz_), rank_info(rank_info_),
    SDN(SDN_), dist(dist_), last_id(last_id_), new_index(new_index_), new_id(new_id_), new_list(new_list_), comm(std::move(comm_))
    {
    }
    ~BlobIdentificationWorkItem2() { }
    virtual void run() {
        // Compute the global blob id and compare to the previous version
        PROFILE_START("Identify blobs maps",1);
        const IntArray& ids = new_index->second;
        static int max_id = -1;
        new_id->first = new_index->first;
        new_id->second = new_index->second;
        if ( last_id.get()!=NULL ) {
            // Compute the timestep-timestep map
            const IntArray& old_ids = last_id->second;
            ID_map_struct map = computeIDMap(Nx,Ny,Nz,old_ids,ids,comm.comm);
            // Renumber the current timestep's ids
            getNewIDs(map,max_id,*new_list);
            renumberIDs(*new_list,new_id->second);
            writeIDMap(map,timestep,id_map_filename);
        } else {
            max_id = -1;
            ID_map_struct map(new_id->first);
            getNewIDs(map,max_id,*new_list);
            writeIDMap(map,timestep,id_map_filename);
        }
        PROFILE_STOP("Identify blobs maps",1);
    }
private:
    BlobIdentificationWorkItem2();
    int timestep;
    int Nx, Ny, Nz;
    const RankInfoStruct& rank_info;
    //std::shared_ptr<const DoubleArray> phase;
    const DoubleArray& SDN;
    const DoubleArray& dist;
    BlobIDstruct last_id, new_index, new_id;
    BlobIDList new_list;
    runAnalysis::commWrapper comm;
};


// Helper class to write the vis file from a thread
class WriteVisWorkItem: public ThreadPool::WorkItemRet<void>
{
public:
    WriteVisWorkItem( int timestep_, std::vector<IO::MeshDataStruct>& visData_,
                     TwoPhase& Avgerages_, fillHalo<double>& fillData_, runAnalysis::commWrapper&& comm_ ):
    timestep(timestep_), visData(visData_), Averages(Avgerages_), fillData(fillData_), comm(std::move(comm_))
    {
    }
    ~WriteVisWorkItem() { }
    virtual void run() {
        PROFILE_START("Save Vis",1);

//         ASSERT(visData[0].vars[0]->name=="SDn");
         ASSERT(visData[0].vars[0]->name=="Pressure");
         ASSERT(visData[0].vars[1]->name=="Velocity_x");
         ASSERT(visData[0].vars[2]->name=="Velocity_y");
         ASSERT(visData[0].vars[3]->name=="Velocity_z");
         ASSERT(visData[0].vars[4]->name=="ID");
//         ASSERT(visData[0].vars[6]->name=="VFmask");
//
//         Array<double>& PhaseData = visData[0].vars[0]->data;
         Array<double>& PressData = visData[0].vars[0]->data;
         Array<double>& VelxData = visData[0].vars[1]->data;
         Array<double>& VelyData = visData[0].vars[2]->data;
         Array<double>& VelzData = visData[0].vars[3]->data;
         Array<double>& IDData = visData[0].vars[4]->data;
//         Array<double>& VFmaskData = visData[0].vars[6]->data;
//
//         fillData.copy(Averages.SDn,PhaseData);
         fillData.copy(Averages.Press,PressData);
         fillData.copy(Averages.Vel_x,VelxData);
         fillData.copy(Averages.Vel_y,VelyData);
         fillData.copy(Averages.Vel_z,VelzData);
         fillData.copy(Averages.ID,IDData);
//         fillData.copy(Averages.VFmask,VFmaskData); //vfmask?

        ASSERT(visData[0].vars[5]->name=="SDn");
        ASSERT(visData[0].vars[6]->name=="SDs");
        ASSERT(visData[0].vars[7]->name=="VFmask");
        ASSERT(visData[0].vars[8]->name=="SDsx");
        ASSERT(visData[0].vars[9]->name=="SDsy");
        ASSERT(visData[0].vars[10]->name=="SDsz");
        ASSERT(visData[0].vars[11]->name=="libba");
        ASSERT(visData[0].vars[12]->name=="libbbc");
        ASSERT(visData[0].vars[13]->name=="libbd");
        ASSERT(visData[0].vars[14]->name=="complibba");

        Array<double>& SDnData = visData[0].vars[5]->data;
        Array<double>& SDsData = visData[0].vars[6]->data;
        Array<double>& VFmaskData = visData[0].vars[7]->data;
        Array<double>& SDsxData = visData[0].vars[8]->data;
        Array<double>& SDsyData = visData[0].vars[9]->data;
        Array<double>& SDszData = visData[0].vars[10]->data;
        Array<double>& libbaData = visData[0].vars[11]->data;
        Array<double>& libbbcData = visData[0].vars[12]->data;
        Array<double>& libbdData = visData[0].vars[13]->data;
        Array<double>& complibbaData = visData[0].vars[14]->data;

        fillData.copy(Averages.SDn,SDnData);
        fillData.copy(Averages.SDs,SDsData);
        fillData.copy(Averages.VFmask,VFmaskData);
        fillData.copy(Averages.SDs_x,SDsxData);
        fillData.copy(Averages.SDs_y,SDsyData);
        fillData.copy(Averages.SDs_z,SDszData);
        fillData.copy(Averages.Vel_x,libbaData);
        fillData.copy(Averages.Vel_y,libbbcData);
        fillData.copy(Averages.Vel_z,libbdData);
        fillData.copy(Averages.Press,complibbaData);
        
        ASSERT(visData[0].vars[15]->name=="DenA");
        ASSERT(visData[0].vars[16]->name=="DenB");
        
        Array<double>& DenAData = visData[0].vars[15]->data;
        Array<double>& DenBData = visData[0].vars[16]->data;
        
        fillData.copy(Averages.DensityA,DenAData);
        fillData.copy(Averages.DensityB,DenBData);

        
        
        IO::writeData( timestep, visData, comm.comm );

        char CurrentIDFilename[40];
        sprintf(CurrentIDFilename,"id_t%d.raw",timestep);
        Averages.AggregateLabels(CurrentIDFilename);

        PROFILE_STOP("Save Vis",1);
    };
private:
    WriteVisWorkItem();
    int timestep;
    std::vector<IO::MeshDataStruct>& visData;
    TwoPhase& Averages;
    fillHalo<double>& fillData;
    runAnalysis::commWrapper comm;
};


// Helper class to run the analysis from within a thread
// Note: Averages will be modified after the constructor is called
class AnalysisWorkItem: public ThreadPool::WorkItemRet<void>
{
public:
    AnalysisWorkItem( AnalysisType type_, int timestep_, TwoPhase& Averages_,
                     BlobIDstruct ids, BlobIDList id_list_, double beta_ ):
    type(type_), timestep(timestep_), Averages(Averages_),
    blob_ids(ids), id_list(id_list_), beta(beta_) { }
    ~AnalysisWorkItem() { }
    virtual void run() {
        
        if ( matches(type,AnalysisType::ComputeAverages_tminus) ) {
            Averages.Initialize();
            Averages.ColorToSignedDistance(1.0,Averages.Phase,Averages.SDn_tminus);
            Averages.ColorToSignedDistance(1.0,Averages.SDn_tminus,Averages.SDn);
            Averages.time_flag = -1;
            Averages.UpdateMeshValues();
            Averages.ComputeDelPhi();
            Averages.ComputeEwn();
        }
        
        if ( matches(type,AnalysisType::ComputeAverages_tcenter) ) {
            Averages.ColorToSignedDistance(1.0,Averages.Phase,Averages.SDn);
        }
        
        if ( matches(type,AnalysisType::ComputeAverages_tplus) ) {
            Averages.ColorToSignedDistance(1.0,Averages.Phase,Averages.SDn_tplus);
            Averages.time_flag = 0;
            Averages.UpdateMeshValues();
            Averages.ComputeDelPhi();
            Averages.Speed();
            Averages.ComputeLocal();
            Averages.Reduce();
            Averages.PrintAll(timestep);
        }

    }
private:
    AnalysisWorkItem();
    AnalysisType type;
    int timestep;
    TwoPhase& Averages;
    BlobIDstruct blob_ids;
    BlobIDList id_list;
    double beta;
    bool geometry;
    bool ComputeAccurateNWPhaseVolume;
};

class ComponentWorkItem: public ThreadPool::WorkItemRet<void>
{
public:
    ComponentWorkItem( AnalysisType type_, int timestep_, TwoPhase& Averages_,
            BlobIDstruct ids, BlobIDList id_list_, double beta_ ):
                type(type_), timestep(timestep_), Averages(Averages_),
                blob_ids(ids), id_list(id_list_), beta(beta_) { }
    ~ComponentWorkItem() { }
    virtual void run() {
        Averages.NumberComponents_NWP = blob_ids->first;
        Averages.Label_NWP = blob_ids->second;
        Averages.Label_NWP_map = *id_list;
        Averages.NumberComponents_WP = 1;
        Averages.Label_WP.fill(0.0);
        if ( matches(type,AnalysisType::ComputeComponents) ) {
            PROFILE_START("Compute comp",1);

            Averages.Initialize();
            Averages.ColorToSignedDistance(1.0,Averages.Phase,Averages.SDn);
            Averages.time_flag = -1;
            Averages.UpdateMeshValues();
            Averages.ComponentAverages();
            Averages.SortBlobs();
            Averages.PrintComponents(timestep);

            PROFILE_STOP("Compute comp",1);
        }
    }
private:
    ComponentWorkItem();
    AnalysisType type;
    int timestep;
    TwoPhase& Averages;
    BlobIDstruct blob_ids;
    BlobIDList id_list;
    double beta;
    bool geometry;
    bool ComputeAccurateNWPhaseVolume;
};



/******************************************************************
 *  MPI comm wrapper for use with analysis                         *
 ******************************************************************/
runAnalysis::commWrapper::commWrapper( int tag_, MPI_Comm comm_, runAnalysis* analysis_ ):
comm(comm_),
tag(tag_),
analysis(analysis_)
{
}
runAnalysis::commWrapper::commWrapper( commWrapper &&rhs ):
comm(rhs.comm),
tag(rhs.tag),
analysis(rhs.analysis)
{
    rhs.tag = -1;
}
runAnalysis::commWrapper::~commWrapper()
{
    if ( tag == -1 )
        return;
    MPI_Barrier( comm );
    analysis->d_comm_used[tag] = false;
}
runAnalysis::commWrapper runAnalysis::getComm( )
{
    // Get a tag from root
    int tag = -1;
    if ( d_rank == 0 ) {
        for (int i=0; i<1024; i++) {
            if ( !d_comm_used[i] ) {
                tag = i;
                break;
            }
        }
        if ( tag == -1 )
            ERROR("Unable to get comm");
    }
    MPI_Bcast( &tag, 1, MPI_INT, 0, d_comm );
    d_comm_used[tag] = true;
    if ( d_comms[tag] == MPI_COMM_NULL )
        MPI_Comm_dup( MPI_COMM_WORLD, &d_comms[tag] );
    return commWrapper(tag,d_comms[tag],this);
}


/******************************************************************
 *  Constructor/Destructors                                        *
 ******************************************************************/
runAnalysis::runAnalysis( std::shared_ptr<Database> db,
                         const RankInfoStruct& rank_info, std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm, std::shared_ptr <Domain> Dm,
                         int Np, bool Regular, double beta, IntArray Map ):
d_Np( Np ),
d_beta( beta ),
d_regular ( Regular),
d_rank_info( rank_info ),
d_Map( Map ),
d_fillData(Dm->Comm,Dm->rank_info,{Dm->Nx-2,Dm->Ny-2,Dm->Nz-2},{1,1,1},0,1),
d_ScaLBL_Comm( ScaLBL_Comm)
{
    
    // Ids of work items to use for dependencies
    ThreadPool::thread_id_t d_wait_blobID;
    ThreadPool::thread_id_t d_wait_analysis;
    ThreadPool::thread_id_t d_wait_vis;
    ThreadPool::thread_id_t d_wait_restart;
    ThreadPool::thread_id_t d_wait_component;
    
    char rankString[20];
    sprintf(rankString,"%05d",Dm->rank());
    d_N[0] = Dm->Nx;
    d_N[1] = Dm->Ny;
    d_N[2] = Dm->Nz;
    d_rank = Dm->rank();
    d_restart_interval = db->getScalar<int>( "restart_interval" );
    d_analysis_interval = db->getScalar<int>( "analysis_interval" );
    d_blobid_interval = db->getScalar<int>( "blobid_interval" );
    d_visualization_interval = db->getScalar<int>( "visualization_interval" );
    d_component_interval = db->getScalar<int>( "component_interval");
    auto restart_file = db->getScalar<std::string>( "restart_file" );
    d_restartFile = restart_file + "." + rankString;
    d_rank = MPI_WORLD_RANK();
    writeIDMap(ID_map_struct(),0,id_map_filename);
    // Initialize IO for silo
    IO::initialize("","silo","false");
    d_meshData.resize(1);
    d_meshData[0].meshName = "domain";
    d_meshData[0].mesh = std::make_shared<IO::DomainMesh>( Dm->rank_info,Dm->Nx-2,Dm->Ny-2,Dm->Nz-2,Dm->Lx,Dm->Ly,Dm->Lz );
    
    Point P,A,B,C;

    auto PressVar = std::make_shared<IO::Variable>();
    auto VxVar = std::make_shared<IO::Variable>();
    auto VyVar = std::make_shared<IO::Variable>();
    auto VzVar = std::make_shared<IO::Variable>();
    auto IDVar = std::make_shared<IO::Variable>();
    auto SDnVar = std::make_shared<IO::Variable>();
    auto SDsVar = std::make_shared<IO::Variable>();
    auto VFmaskVar = std::make_shared<IO::Variable>();
    auto SDsxVar = std::make_shared<IO::Variable>();
    auto SDsyVar = std::make_shared<IO::Variable>();
    auto SDszVar = std::make_shared<IO::Variable>();
    auto libbaVar = std::make_shared<IO::Variable>();
    auto libbbcVar = std::make_shared<IO::Variable>();
    auto libbdVar = std::make_shared<IO::Variable>();
    auto complibbaVar = std::make_shared<IO::Variable>();
    auto DenAVar = std::make_shared<IO::Variable>();
    auto DenBVar = std::make_shared<IO::Variable>();
  
    
     PressVar->name = "Pressure";
     PressVar->type = IO::VariableType::VolumeVariable;
     PressVar->dim = 1;
     PressVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
     d_meshData[0].vars.push_back(PressVar);

     VxVar->name = "Velocity_x";
     VxVar->type = IO::VariableType::VolumeVariable;
     VxVar->dim = 1;
     VxVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
     d_meshData[0].vars.push_back(VxVar);
     VyVar->name = "Velocity_y";
     VyVar->type = IO::VariableType::VolumeVariable;
     VyVar->dim = 1;
     VyVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
     d_meshData[0].vars.push_back(VyVar);
     VzVar->name = "Velocity_z";
     VzVar->type = IO::VariableType::VolumeVariable;
     VzVar->dim = 1;
     VzVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
     d_meshData[0].vars.push_back(VzVar);

     IDVar->name = "ID";
     IDVar->type = IO::VariableType::VolumeVariable;
     IDVar->dim = 1;
     IDVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
     d_meshData[0].vars.push_back(IDVar);
    
    
    SDnVar->name = "SDn";
    SDnVar->type = IO::VariableType::VolumeVariable;
    SDnVar->dim = 1;
    SDnVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(SDnVar);

    SDsVar->name = "SDs";
    SDsVar->type = IO::VariableType::VolumeVariable;
    SDsVar->dim = 1;
    SDsVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(SDsVar);

    VFmaskVar->name = "VFmask";
    VFmaskVar->type = IO::VariableType::VolumeVariable;
    VFmaskVar->dim = 1;
    VFmaskVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(VFmaskVar);

    SDsxVar->name = "SDsx";
    SDsxVar->type = IO::VariableType::VolumeVariable;
    SDsxVar->dim = 1;
    SDsxVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(SDsxVar);

    SDsyVar->name = "SDsy";
    SDsyVar->type = IO::VariableType::VolumeVariable;
    SDsyVar->dim = 1;
    SDsyVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(SDsyVar);

    SDszVar->name = "SDsz";
    SDszVar->type = IO::VariableType::VolumeVariable;
    SDszVar->dim = 1;
    SDszVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(SDszVar);

    libbaVar->name = "libba";
    libbaVar->type = IO::VariableType::VolumeVariable;
    libbaVar->dim = 1;
    libbaVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(libbaVar);

    libbbcVar->name = "libbbc";
    libbbcVar->type = IO::VariableType::VolumeVariable;
    libbbcVar->dim = 1;
    libbbcVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(libbbcVar);

    libbdVar->name = "libbd";
    libbdVar->type = IO::VariableType::VolumeVariable;
    libbdVar->dim = 1;
    libbdVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(libbdVar);

    complibbaVar->name = "complibba";
    complibbaVar->type = IO::VariableType::VolumeVariable;
    complibbaVar->dim = 1;
    complibbaVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(complibbaVar);
    
    DenAVar->name = "DenA";
    DenAVar->type = IO::VariableType::VolumeVariable;
    DenAVar->dim = 1;
    DenAVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(DenAVar);
    
    DenBVar->name = "DenB";
    DenBVar->type = IO::VariableType::VolumeVariable;
    DenBVar->dim = 1;
    DenBVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
    d_meshData[0].vars.push_back(DenBVar);

//     auto PhaseVar = std::make_shared<IO::Variable>();
//     auto PressVar = std::make_shared<IO::Variable>();
//     auto VxVar = std::make_shared<IO::Variable>();
//     auto VyVar = std::make_shared<IO::Variable>();
//     auto VzVar = std::make_shared<IO::Variable>();
//     auto IDVar = std::make_shared<IO::Variable>();
//     auto VFmaskVar = std::make_shared<IO::Variable>();
//
//     PhaseVar->name = "SDn";
//     PhaseVar->type = IO::VariableType::VolumeVariable;
//     PhaseVar->dim = 1;
//     PhaseVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
//     d_meshData[0].vars.push_back(PhaseVar);
//

//
//     VFmaskVar->name = "VFmask";
//     VFmaskVar->type = IO::VariableType::VolumeVariable;
//     VFmaskVar->dim = 1;
//     VFmaskVar->data.resize(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2);
//     d_meshData[0].vars.push_back(VFmaskVar);


    // Initialize the comms
    MPI_Comm_dup(MPI_COMM_WORLD,&d_comm);
    for (int i=0; i<1024; i++) {
        d_comms[i] = MPI_COMM_NULL;
        d_comm_used[i] = false;
    }
    // Initialize the threads
    int N_threads = db->getWithDefault<int>( "N_threads", 4 );
    auto method = db->getWithDefault<std::string>( "load_balance", "default" );
    createThreads( method, N_threads );
}
runAnalysis::~runAnalysis( )
{
    // Finish processing analysis
    finish();
    // Clear internal data
    MPI_Comm_free( &d_comm );
    for (int i=0; i<1024; i++) {
        if ( d_comms[i] != MPI_COMM_NULL )
            MPI_Comm_free(&d_comms[i]);
    }
}
void runAnalysis::finish( )
{
    PROFILE_START("finish");
    // Wait for the work items to finish
    d_tpool.wait_pool_finished();
    // Clear the wait ids
    d_wait_blobID.reset();
    d_wait_analysis.reset();
    d_wait_vis.reset();
    d_wait_restart.reset();
    d_wait_component.reset();
    // Syncronize
    MPI_Barrier( d_comm );
    PROFILE_STOP("finish");
}


/******************************************************************
 *  Set the thread affinities                                      *
 ******************************************************************/
void print( const std::vector<int>& ids )
{
    if ( ids.empty() )
        return;
    printf("%i",ids[0]);
    for (size_t i=1; i<ids.size(); i++)
        printf(", %i",ids[i]);
    printf("\n");
}
void runAnalysis::createThreads( const std::string& method, int N_threads )
{
    // Check if we are not using analysis threads
    if ( method == "none" )
        return;
    // Check if we have thread support
    int thread_support;
    MPI_Query_thread( &thread_support );
    if ( thread_support < MPI_THREAD_MULTIPLE ) {
        std::cerr << "Warning: Failed to start MPI with necessary thread support, thread support will be disabled" << std::endl;
        return;
    }
    // Create the threads
    const auto cores = d_tpool.getProcessAffinity();
    if ( cores.empty() ) {
        d_tpool.setNumThreads( N_threads );
    } else if ( method == "default" ) {
        d_tpool.setNumThreads( N_threads );
    } else if ( method == "independent" ) {
        int N = cores.size() - 1;
        d_tpool.setNumThreads( N );
        d_tpool.setThreadAffinity( { cores[0] } );
        for ( int i=0; i<N; i++)
            d_tpool.setThreadAffinity( i, { cores[i+1] } );
    }
    // Print the current affinities
    if ( d_rank == 0 ) {
        printf("Affinities - rank 0:\n");
        printf("Main: ");
        print(d_tpool.getProcessAffinity());
        for (int i=0; i<d_tpool.getNumThreads(); i++) {
            printf("Thread %i: ",i+1);
            print(d_tpool.getThreadAffinity(i));
        }
    }
}


/******************************************************************
 *  Check which analysis we want to perform                        *
 ******************************************************************/

AnalysisType runAnalysis::computeAnalysisType( int timestep )
{
    AnalysisType type = AnalysisType::AnalyzeNone;
    if ( timestep%d_analysis_interval + 6 == d_analysis_interval ) {
      
        type |= AnalysisType::ComputeAverages_tminus;
    }

    if ( timestep%d_analysis_interval + 4 == d_analysis_interval ) {
        type |= AnalysisType::ComputeAverages_tcenter;
    }
    if ( timestep%d_analysis_interval + 2 ==  d_analysis_interval) {
  
        type |= AnalysisType::ComputeAverages_tplus;
    }
    
    if ( timestep%d_component_interval + 4 == d_component_interval ) {
        // Copy the averages to the CPU and identify blobs
        //type |= AnalysisType::CopySimState;
        type |= AnalysisType::IdentifyBlobs;
        type |= AnalysisType::ComputeComponents;
    }
    
    
    if (timestep%d_restart_interval + 4 == d_restart_interval) {
    
        type |= AnalysisType::CreateRestart;
    }
    if (timestep%d_visualization_interval + 2 == d_visualization_interval ) {
      
        type |= AnalysisType::WriteVis;
        type |= AnalysisType::CopySimState;
     
    }
    return type;
}


/******************************************************************
 *  Run the analysis                                               *
 ******************************************************************/
void runAnalysis::run( int timestep, TwoPhase& Averages, const double *Phi,
                      double *Pressure, double *Velocity, double *fq, double *Den)
{
    //
}

void runAnalysis::run3( int timestep, TwoPhase& Averages, const double *Phi,
                       double *Pressure, double *Velocity, double *fq, double *Aq, double *Bq, double *Den, int Np, double * TimeAveragedVelocity, double Fx, double Fy, double Fz) {}
/******************************************************************
 *  Run the analysis                                               *
 ******************************************************************/
void runAnalysis::run5( int timestep, TwoPhase& Averages,
                       const double *Phi, double *Pressure, double *Velx, double *Vely,
                       double *Velz, double *fq,
                       double *GradPhiX, double *GradPhiY, double *GradPhiZ, double * CField, double *DenA, double *DenB,
                       int Np, double Fx, double Fy, double Fz)
{
    size_t N = d_N[0]*d_N[1]*d_N[2];
    
    // Check which analysis steps we need to perform
    auto type = computeAnalysisType( abs(timestep) );
    if ( type == AnalysisType::AnalyzeNone )
        return;
    
    // Check how may queued items we have
    if ( d_tpool.N_queued() > 20 ) {
        std::cerr << "Analysis queue is getting behind, waiting ...\n";
        finish();
    }
    
    // Copy the appropriate variables to the host (so we can spawn new threads)
    ScaLBL_DeviceBarrier();
    PROFILE_START("Copy data to host",1);
    std::shared_ptr<DoubleArray> phase;

    if ( timestep%d_analysis_interval + 6 == d_analysis_interval ) {
        ScaLBL_CopyToHost(Averages.Phase.data(),Phi,N*sizeof(double));
        d_ScaLBL_Comm->RegularLayout(d_Map,&Pressure[0],Averages.Press);
     //   std::cout << "-4 timestep=" << timestep << std::endl;
    }
    
    if ( timestep%d_analysis_interval + 4 == d_analysis_interval ) {
        Averages.Fx = Fx;
        Averages.Fy = Fy;
        Averages.Fz = Fz;
        
        ScaLBL_CopyToHost(Averages.Phase.data(),Phi,N*sizeof(double));
        ScaLBL_CopyToHost(Averages.DensityA.data(),DenA,N*sizeof(double));
        ScaLBL_CopyToHost(Averages.DensityB.data(),DenB,N*sizeof(double));
        ScaLBL_CopyToHost(Averages.GradPhiX.data(),&GradPhiX[0],N*sizeof(double));
        ScaLBL_CopyToHost(Averages.GradPhiY.data(),&GradPhiY[0],N*sizeof(double));
        ScaLBL_CopyToHost(Averages.GradPhiZ.data(),&GradPhiZ[0],N*sizeof(double));
        ScaLBL_CopyToHost(Averages.CField.data(),&CField[0],N*sizeof(double));
    
        ScaLBL_CopyToHost(Averages.Vel_x.data(),&Velx[0],N*sizeof(double));
        ScaLBL_CopyToHost(Averages.Vel_y.data(),&Vely[0],N*sizeof(double));
        ScaLBL_CopyToHost(Averages.Vel_z.data(),&Velz[0],N*sizeof(double));
        d_ScaLBL_Comm->RegularLayout(d_Map,&Pressure[0],Averages.Press);
      //  std::cout << "-2 timestep=" << timestep << std::endl;
    }

    if ( timestep%d_analysis_interval + 2 == d_analysis_interval ) {
        ScaLBL_CopyToHost(Averages.Phase.data(),Phi,N*sizeof(double));
      //  std::cout << "0 timestep=" << timestep << std::endl;
          
    }

    //Spawn threads to do blob identification work
    //if ( matches(type,AnalysisType::IdentifyBlobs) ) {
    if ( timestep%d_component_interval + 4 == d_component_interval ) {
        if ( d_rank == 0 )  { printf("IdentifyBlobs being called at timestep %i \n",timestep);  }

                phase = std::shared_ptr<DoubleArray>(new DoubleArray(d_N[0],d_N[1],d_N[2]));
        
                ScaLBL_CopyToHost(phase->data(),Phi,N*sizeof(double));
                Averages.Dm->CommunicateMeshHalo(*phase);

                Averages.ColorToSignedDistance(1.0,*phase,Averages.SDn);

                Averages.UpdateMeshValues();

                BlobIDstruct new_index(new std::pair<int,IntArray>(0,IntArray()));
                BlobIDstruct new_ids(new std::pair<int,IntArray>(0,IntArray()));
                BlobIDList new_list(new std::vector<BlobIDType>());
                auto work1 = new BlobIdentificationWorkItem1(timestep,d_N[0],d_N[1],d_N[2],d_rank_info,
                    Averages.SDn,Averages.SDs,d_last_ids,new_index,new_ids,new_list,getComm());
                auto work2 = new BlobIdentificationWorkItem2(timestep,d_N[0],d_N[1],d_N[2],d_rank_info,
                    Averages.SDn,Averages.SDs,d_last_ids,new_index,new_ids,new_list,getComm());
                work1->add_dependency(d_wait_blobID);
                work2->add_dependency(d_tpool.add_work(work1));
                d_wait_blobID = d_tpool.add_work(work2);
                d_last_index = new_index;
                d_last_ids = new_ids;
                d_last_id_map = new_list;
                
    }    
    
    
    std::shared_ptr<double> cfq,cDenA,cDenB,cPhase,cVx,cVy,cVz,cGPhiX,cGPhiY,cGPhiZ,cCField;
    if (timestep%d_restart_interval + 4 == d_restart_interval){
        // Copy restart data to the CPU
        cDenA = std::shared_ptr<double>(new double[N],DeleteArray<double>);
        cDenB = std::shared_ptr<double>(new double[N],DeleteArray<double>);
        cfq = std::shared_ptr<double>(new double[19*d_Np],DeleteArray<double>);
        cPhase = std::shared_ptr<double>(new double[N],DeleteArray<double>);
        
        cVx = std::shared_ptr<double>(new double[N],DeleteArray<double>);
        cVy = std::shared_ptr<double>(new double[N],DeleteArray<double>);
        cVz = std::shared_ptr<double>(new double[N],DeleteArray<double>);
        
        cGPhiX = std::shared_ptr<double>(new double[N],DeleteArray<double>);
        cGPhiY = std::shared_ptr<double>(new double[N],DeleteArray<double>);
        cGPhiZ = std::shared_ptr<double>(new double[N],DeleteArray<double>);
        cCField = std::shared_ptr<double>(new double[N],DeleteArray<double>);
        
        ScaLBL_CopyToHost(cfq.get(),fq,19*d_Np*sizeof(double));
        ScaLBL_CopyToHost(cDenA.get(),DenA,N*sizeof(double));
        ScaLBL_CopyToHost(cDenB.get(),DenB,N*sizeof(double));
        ScaLBL_CopyToHost(cPhase.get(),Phi,N*sizeof(double));
        
        ScaLBL_CopyToHost(cVx.get(),Velx,N*sizeof(double));
        ScaLBL_CopyToHost(cVy.get(),Vely,N*sizeof(double));
        ScaLBL_CopyToHost(cVz.get(),Velz,N*sizeof(double));
        
        ScaLBL_CopyToHost(cGPhiX.get(),GradPhiX,N*sizeof(double));
        ScaLBL_CopyToHost(cGPhiY.get(),GradPhiY,N*sizeof(double));
        ScaLBL_CopyToHost(cGPhiZ.get(),GradPhiZ,N*sizeof(double));
        ScaLBL_CopyToHost(cCField.get(),CField,N*sizeof(double));
        
    }
    PROFILE_STOP("Copy data to host",1);
    
    // Spawn threads to do blob identification work
    // if ( matches(type,AnalysisType::IdentifyBlobs) ) {}
    
    if ( timestep%d_analysis_interval + 6 == d_analysis_interval ) {
        auto work = new AnalysisWorkItem(type,timestep,Averages,d_last_index,d_last_id_map,d_beta);
        work->add_dependency(d_wait_analysis);
        d_wait_analysis = d_tpool.add_work(work);
    }
    
    if ( timestep%d_analysis_interval + 4 == d_analysis_interval ) {
        auto work = new AnalysisWorkItem(type,timestep,Averages,d_last_index,d_last_id_map,d_beta);
        work->add_dependency(d_wait_analysis);
        d_wait_analysis = d_tpool.add_work(work);
    }

    if ( timestep%d_component_interval + 4 == d_component_interval ) {
        auto work = new ComponentWorkItem(type,timestep,Averages,d_last_index,d_last_id_map,d_beta);
        work->add_dependency(d_wait_blobID);
        //work->add_dependency(d_wait_analysis);
        work->add_dependency(d_wait_component);
        d_wait_component = d_tpool.add_work(work);
    }
    
    
    if ( timestep%d_analysis_interval + 2 == d_analysis_interval ) {
        auto work = new AnalysisWorkItem(type,timestep,Averages,d_last_index,d_last_id_map,d_beta);
        work->add_dependency(d_wait_analysis);
        work->add_dependency(d_wait_vis);
        d_wait_analysis = d_tpool.add_work(work);
    }
 
    if (timestep%d_visualization_interval + 2 == d_visualization_interval) {
        // Write the vis files
        auto work = new WriteVisWorkItem( timestep, d_meshData, Averages, d_fillData, getComm() );
        //  work->add_dependency(d_wait_blobID);
        work->add_dependency(d_wait_analysis);
        work->add_dependency(d_wait_vis);
        d_wait_vis = d_tpool.add_work(work);
    }

    if (timestep%d_restart_interval + 4 == d_restart_interval){
        if (d_rank==0) {
            FILE *Rst = fopen("Restart.txt","w");
            fprintf(Rst,"%i\n",timestep);
            fclose(Rst);
        }
        // Write the restart file (using a seperate thread)
        auto work = new WriteRestartWorkItem(d_restartFile.c_str(),cDenA,cDenB,cfq,cPhase,cVx,cVy,cVz,Np,N,Fz,cGPhiX,cGPhiY,cGPhiZ,cCField);
        work->add_dependency(d_wait_restart);
        d_wait_restart = d_tpool.add_work(work);
    }

}

void runAnalysis::run6( int timestep, TwoPhase& Averages, int Np, double *libb)
{
    size_t N = d_N[0]*d_N[1]*d_N[2];

    // Copy the appropriate variables to the host (so we can spawn new threads)
    ScaLBL_DeviceBarrier();
    PROFILE_START("Copy data to host",1);

    d_ScaLBL_Comm->RegularLayout(d_Map,&libb[0+4*Np],Averages.Press);
//    ScaLBL_CopyToHost(Averages.SDn.data(),&libb[0],N*sizeof(double));

    // Write the vis files
    auto work = new WriteVisWorkItem( timestep, d_meshData, Averages, d_fillData, getComm() );
    work->add_dependency(d_wait_vis);
    d_wait_vis = d_tpool.add_work(work);

}

void runAnalysis::run7( int timestep, TwoPhase& Averages, int Np, double *libb, double * Velx, double * Vely, double * Velz, double * DenA, double * DenB)
{
    size_t N = d_N[0]*d_N[1]*d_N[2];

    // Copy the appropriate variables to the host (so we can spawn new threads)
    ScaLBL_DeviceBarrier();
    PROFILE_START("Copy data to host",1);

//    d_ScaLBL_Comm->RegularLayout(d_Map,&libb[0+4*Np],Averages.Press);
    ScaLBL_CopyToHost(Averages.SDn.data(),&libb[0],N*sizeof(double));
    ScaLBL_CopyToHost(Averages.Vel_x.data(),&Velx[0],N*sizeof(double));
    ScaLBL_CopyToHost(Averages.Vel_y.data(),&Vely[0],N*sizeof(double));
    ScaLBL_CopyToHost(Averages.Vel_z.data(),&Velz[0],N*sizeof(double));
    ScaLBL_CopyToHost(Averages.DensityA.data(),&DenA[0],N*sizeof(double));
    ScaLBL_CopyToHost(Averages.DensityB.data(),&DenB[0],N*sizeof(double));
    // Write the vis files
    auto work = new WriteVisWorkItem( timestep, d_meshData, Averages, d_fillData, getComm() );
    work->add_dependency(d_wait_vis);
    d_wait_vis = d_tpool.add_work(work);

}



void runAnalysis::RayTraceRun( int timestep, TwoPhase& Averages, int Np)
{
   // int N = d_N[0]*d_N[1]*d_N[2];
    
    // Check which analysis steps we need to perform
    auto type = computeAnalysisType( timestep );
    if ( type == AnalysisType::AnalyzeNone )
        return;
    
    // Check how may queued items we have
    if ( d_tpool.N_queued() > 20 ) {
        std::cerr << "Analysis queue is getting behind, waiting ...\n";
        finish();
    }
    PROFILE_START("run");
    
    // Copy the appropriate variables to the host (so we can spawn new threads)
    ScaLBL_DeviceBarrier();
    PROFILE_START("Copy data to host",1);
    std::shared_ptr<DoubleArray> phase;
    
    if ( timestep%d_analysis_interval + 2 == d_analysis_interval ) {
        auto work = new AnalysisWorkItem(type,timestep,Averages,d_last_index,d_last_id_map,d_beta);
        //  work->add_dependency(d_wait_blobID);
        work->add_dependency(d_wait_analysis);
        work->add_dependency(d_wait_vis);     // Make sure we are done using analysis before modifying
        d_wait_analysis = d_tpool.add_work(work);
    }
    
    if (timestep%d_visualization_interval + 2 ==  d_visualization_interval){
        // Write the vis files
        auto work = new WriteVisWorkItem( timestep, d_meshData, Averages, d_fillData, getComm() );
        //  work->add_dependency(d_wait_blobID);
        work->add_dependency(d_wait_analysis);
        work->add_dependency(d_wait_vis);
        d_wait_vis = d_tpool.add_work(work);
    }
    PROFILE_STOP("run");
}
