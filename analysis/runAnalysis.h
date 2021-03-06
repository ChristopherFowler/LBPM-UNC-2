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
#ifndef RunAnalysis_H_INC
#define RunAnalysis_H_INC

#include "analysis/analysis.h"
#include "analysis/TwoPhase.h"
#include "common/Communication.h"
#include "common/ScaLBL.h"
#include "threadpool/thread_pool.h"


typedef std::shared_ptr<std::pair<int,IntArray>> BlobIDstruct;
typedef std::shared_ptr<std::vector<BlobIDType>> BlobIDList;


// Types of analysis
enum class AnalysisType : uint64_t { AnalyzeNone=0, IdentifyBlobs=0x01, CopySimState=0x02, ComputeAverages_tcenter=0x04,
 CreateRestart=0x8, WriteVis=0x10, ComputeAverages_tplus=0x20, ComputeAverages_tminus=0x40, ComputeComponents=0x80 };


//! Class to run the analysis in multiple threads
class runAnalysis
{
public:

    //! Constructor
    runAnalysis( std::shared_ptr<Database> db, const RankInfoStruct& rank_info,
            std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm, std::shared_ptr <Domain> dm, int Np, bool Regular, double beta, IntArray Map );

    //! Destructor
    ~runAnalysis();

    //! Run the next analysis
    void run( int timestep, TwoPhase &Averages, const double *Phi,
        double *Pressure, double *Velocity, double *fq, double *Den );
    //! Run the next analysis
    void run3( int timestep, TwoPhase& Averages, const double *Phi,
              double *Pressure, double *Velocity, double *fq, double *Aq, double *Bq, double *Den, int Np, double * TimeAveragedVelocity,double Fx, double Fy, double Fz);
    
    void run4( int timestep, TwoPhase& Averages,
              const double *Phi, double *Pressure, double *Velocity,
              double *fq, double *Aq, double *Bq,
              double *GradPhiX, double * GradPhiY, double * GradPhiZ, double *DenA, double *DenB,
              int Np, double Fx, double Fy, double Fz);
    
    void run5( int timestep, TwoPhase& Averages,
              const double *Phi, double *Pressure, double *Velx, double *Vely,
              double *Velz, double *fq, 
              double *GradPhiX, double *GradPhiY, double *GradPhiZ, double * CFeild, double *DenA, double *DenB,
              int Np, double Fx, double Fy, double Fz);

    void run6( int timestep, TwoPhase& Averages, int Np, double *libb);
    void run7( int timestep, TwoPhase& Averages, int Np, double *libb, double *Velx, double *Vely, double *Velz, double * DenA, double * DenB);

    void RayTraceRun( int timestep, TwoPhase& Averages, int Np);

    //! Finish all active analysis
    void finish();

    /*!
     *  \brief    Set the affinities
     *  \details  This function will create the analysis threads and set the affinity
     *      of this thread and all analysis threads.  If MPI_THREAD_MULTIPLE is not
     *      enabled, the analysis threads will be disabled and the analysis will run in the current thread.
     * @param[in] method    Method used to control the affinities:
     *                      none - Don't use threads (runs all analysis in the current thread)
     *                      default - Create the specified number of threads, but don't load balance
     *                      independent - Create the necessary number of threads to fill all cpus,
     *                                and set the affinities based on the current process such
     *                                that all threads run on independent cores
     * @param[in] N_threads Number of threads, only used by some of the methods
     */
    void createThreads( const std::string& method = "default", int N_threads = 4 );


private:

    runAnalysis();

    // Determine the analysis to perform
    AnalysisType computeAnalysisType( int timestep );

public:

    class commWrapper
    {
      public:
        MPI_Comm comm;
        int tag;
        runAnalysis *analysis;
        commWrapper( int tag, MPI_Comm comm, runAnalysis *analysis );
        commWrapper( ) = delete;
        commWrapper( const commWrapper &rhs ) = delete;
        commWrapper& operator=( const commWrapper &rhs ) = delete;
        commWrapper( commWrapper &&rhs );
        ~commWrapper();
    };

    // Get a comm (not thread safe)
    commWrapper getComm( );

private:

    int d_N[3];
    int d_Np;
    int d_rank;
    int d_restart_interval, d_analysis_interval, d_blobid_interval, d_visualization_interval,d_component_interval;
    double d_beta;
    bool d_regular;
    ThreadPool d_tpool;
    RankInfoStruct d_rank_info;
    IntArray d_Map;
    BlobIDstruct d_last_ids;
    BlobIDstruct d_last_index;
    BlobIDList d_last_id_map;
    std::vector<IO::MeshDataStruct> d_meshData;
    fillHalo<double> d_fillData;
    std::string d_restartFile;
    MPI_Comm d_comm;
    MPI_Comm d_comms[1024];
    volatile bool d_comm_used[1024];
    std::shared_ptr<ScaLBL_Communicator> d_ScaLBL_Comm;

    // Ids of work items to use for dependencies
    ThreadPool::thread_id_t d_wait_blobID;
    ThreadPool::thread_id_t d_wait_analysis;
    ThreadPool::thread_id_t d_wait_vis;
    ThreadPool::thread_id_t d_wait_restart;
    ThreadPool::thread_id_t d_wait_component;

    // Friends
    friend commWrapper::~commWrapper();

};

#endif
