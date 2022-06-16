

 #define WBC 1


#include "models/ColorModel.h"
#include "analysis/distance.h"
#include <vector>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

int CM_DOUBLE = 8;
int CM_INT = 1;


void ScaLBL_ColorModel::AssignComponentLabels(double *phase) {}

std::string getLastLine(std::ifstream& in){
    std::string line;
    while (in >> std::ws && std::getline(in, line));
    
    return line;
}

inline double analytical_phase(double i, double j, double k, double cx, double cy, double cz, double beta, double radius) {}
bool ScaLBL_ColorModel::GridIndependence() { return true; }
bool ScaLBL_ColorModel::ReadInputAndValidate(){ return 1; } // Not used
inline void ComputeFiniteDifferenceStencilCoefficients(double *qDistances,double *stencilCoefs, int N) {}
inline void save_3Ddouble_parallel(int rank, double* field, int Nx, int Ny, int Nz, int depth, std::string s) {}
inline void save_3Dchar_parallel(int rank, char* field, int Nx, int Ny, int Nz, int depth, std::string s) {}
inline double calculate_porosity(char * id,
                                 int Nx, int Ny, int Nz, Domain &Dm, int rank) {}



ScaLBL_ColorModel::ScaLBL_ColorModel(int RANK,
                                     int NP, MPI_Comm COMM):
rank(RANK), nprocs(NP), Restart(0),timestep(0),timestepMax(0),tauA(0),tauB(0),rhoA(0),rhoB(0),alpha(0),beta(0),
Fx(0),Fy(0),Fz(0),flux(0),din(0),dout(0),inletA(0),inletB(0),outletA(0),outletB(0),
Nx(0),Ny(0),Nz(0),N(0),Np(0),Ni(0),Nsb(0),nprocx(0),nprocy(0),nprocz(0),BoundaryCondition(0),Lx(0),Ly(0),Lz(0),WBCFlag(0),comm(COMM){}

ScaLBL_ColorModel::~ScaLBL_ColorModel(){}

void ScaLBL_ColorModel::ReadParams(string filename){
    // read the input database
    db = std::make_shared<Database>( filename );
    domain_db = db->getDatabase( "Domain" );
    color_db = db->getDatabase( "Color" );
    analysis_db = db->getDatabase( "Analysis" );
    
    // Color Model parameters
    timestepMax = color_db->getScalar<int>( "timestepMax" );
    tauA = color_db->getScalar<double>( "tauA" );
    tauB = color_db->getScalar<double>( "tauB" );
    rhoA = color_db->getScalar<double>( "rhoA" );
    rhoB = color_db->getScalar<double>( "rhoB" );
    Fx = color_db->getVector<double>( "F" )[0];
    Fy = color_db->getVector<double>( "F" )[1];
    Fz = color_db->getVector<double>( "F" )[2];
    alpha = color_db->getScalar<double>( "alpha" );
    beta = color_db->getScalar<double>( "beta" );
    Restart = color_db->getScalar<bool>( "Restart" );
    din = color_db->getScalar<double>( "din" );
    dout = color_db->getScalar<double>( "dout" );
    flux = color_db->getScalar<double>( "flux" );
    method = color_db->getScalar<int>("method"); // method = 1=BGK, 2=MRT, 3=LIBB
    input_angle = color_db->getScalar<double>( "contact_angle" );
    inletA=1.f;
    inletB=0.f;
    outletA=0.f;
    outletB=1.f;
    
    offsetdistance = domain_db->getScalar<double>("offsetDistance");
    
    // Read domain parameters
    auto L = domain_db->getVector<double>( "L" );
    auto size = domain_db->getVector<int>( "n" );
    auto nproc = domain_db->getVector<int>( "nproc" );
    BoundaryCondition = domain_db->getScalar<int>( "BC" );
    legacySDs = domain_db->getScalar<bool>("legacySDs");
    Nx = size[0];
    Ny = size[1];
    Nz = size[2];
    Lx = L[0];
    Ly = L[1];
    Lz = L[2];
    nprocx = nproc[0];
    nprocy = nproc[1];
    nprocz = nproc[2];

    int WBCFlag = 0;
    if ( domain_db->keyExists( "WBCFlag" ) ) {
        WBCFlag = domain_db->getScalar<int>( "WBCFlag" ); //0: Original, 1: Situational off (add if statement)
    }
    if (rank == 0) printf("WBCFlag: %i \n",WBCFlag);

    string susr = domain_db->getScalar<string>( "usr" );
   // int fs = domain_db->getScalar<int>( "filesystem" );
    
    if (rank == 0) {
        printf("INPUT PARAMETERS:\n");
        //printf("subdomain size: Nx=%d Ny=%d Nz=%d\n",Nx,Ny,Nz);
        printf("nprocx=%d nprocy=%d nprocz=%d\n", nprocx, nprocy, nprocz);
        printf("Lx=%-5.5E Ly=%-5.5E Lz=%-5.5E\n",Lx,Ly,Lz);
        printf("BC=%d\n",BoundaryCondition);
        printf("tauA=%-5.5E tauB=%-5.5E\n",tauA,tauB); //rhoA=%-5.5E rhoB=%-5.5E\n",tauA,tauB,rhoA,rhoB);
        printf("Fx=%-5.5E Fy=%-5.5E Fz=%-5.5E\n",Fx,Fy,Fz);
        
        printf("contact_angle=%-5.5E\n",input_angle);
       /*
        printf("din=%-5.5E dout=%-5.5E flux=%-5.5E\n",din,dout,flux);
        printf("filesystem=%d usr=%s\n",fs,susr.c_str());
         */
        printf("timestepMax=%d\n",timestepMax);
        printf("\n");
        printf("LBM METHOD:\n");

        if (method == 7 ) cout << "NON-AA MRT-LIBB" << endl << endl;
    }
    if (BoundaryCondition==4) flux *= rhoA;
}





void ScaLBL_ColorModel::SetDomain(){
    
    Dm  = std::shared_ptr<Domain>(new Domain(domain_db,comm));
    Mask  = std::shared_ptr<Domain>(new Domain(domain_db,comm));
    Nx+=2; Ny+=2; Nz += 2;
    N = Nx*Ny*Nz;
    id = new char [N];
    for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = 1;
    std::array<size_t,3> subdivide = { 1, 1, 1 };
    if ( analysis_db->keyExists( "subdivide" ) ) {
        auto tmp = analysis_db->getVector<size_t>( "subdivide" );
        subdivide = { tmp[0], tmp[1], tmp[2] };
    }
    Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm,subdivide) );
    MPI_Barrier(comm);
    
    Dm->CommInit();
    MPI_Barrier(comm);
    rank = Dm->rank();
}







void ScaLBL_ColorModel::ReadInput(){

  
    size_t readID,read_file;
    Mask->ReadIDs(db);
    for (int i=0; i<Nx*Ny*Nz; i++) Averages->ID(i) = Mask->id[i];
    for (int i=0; i<Nx*Ny*Nz; i++) id[i] = Mask->id[i];  // save what was read
    
        if (rank == 0) cout << "Using Ray-traced signed distance" << endl;
        char LocalRankFilename[40];
        //sprintf(LocalRankFilename,"OffsetDist.%05i",rank);
        sprintf(LocalRankFilename,"SignDist.%05i",rank);
        SignDist = new double[Nx*Ny*Nz];
        FILE *RayTraceDistance = fopen(LocalRankFilename,"rb");
        readID=fread(SignDist,8,N,RayTraceDistance);
        if (readID != size_t(N)) printf("lbpm_segmented_pp: Error reading RayTraceDistance \n");
        fclose(RayTraceDistance);

        for (size_t i=0; i < Nx*Ny*Nz; i++) Averages->SDs(i) = SignDist[i];
        Dm->CommunicateMeshHalo(Averages->SDs);
        for (size_t n = 0; n < N; n++) SignDist[n] = Averages->SDs(n);
        
        gradphix  = new double[N];
        gradphiy  = new double[N];
        gradphiz  = new double[N];
        
        sprintf(LocalRankFilename,"ns_x.%05i",rank);
        FILE *NSX = fopen(LocalRankFilename,"rb");
        read_file = fread(gradphix,8,N,NSX);
        if (read_file != size_t(N)) printf("Error reading NSX \n");
        fclose(NSX);
       
        sprintf(LocalRankFilename,"ns_y.%05i",rank);
        FILE *NSY = fopen(LocalRankFilename,"rb");
        read_file = fread(gradphiy,8,N,NSY);
        if (read_file != size_t(N)) printf("Error reading NSY \n");
        fclose(NSY);
        
        sprintf(LocalRankFilename,"ns_z.%05i",rank);
        FILE *NSZ = fopen(LocalRankFilename,"rb");
        read_file = fread(gradphiz,8,N,NSZ);
        if (read_file != size_t(N)) printf("Error reading NSZ \n");
        fclose(NSZ);
        
        vfmask = new double[N];
        sprintf(LocalRankFilename,"VFmask.%05i",rank);
        FILE *VFFILE = fopen(LocalRankFilename,"rb");
        read_file = fread(vfmask,8,N,VFFILE);
        if (read_file != size_t(N)) printf("Error reading VFFILE \n");
        fclose(VFFILE);

        for (size_t n = 0; n < N; n++) Averages->VFmask(n) = vfmask[n];
        Dm->CommunicateMeshHalo(Averages->VFmask);
        for (size_t n = 0; n < N; n++) vfmask[n] = Averages->VFmask(n);
   
  
}




static inline void fgetl( char * str,
                         int num, FILE * stream ) {
    char* ptr = fgets( str, num, stream );
    if ( 0 ) {char *temp = (char *)&ptr; temp++;}
}

double ReadSpherePacking() {
    FILE *fid = fopen("pack.out","rb");
    char line[100];
    int count = 0;
    double *List_rad;  List_rad = new double[4];
    fgets(line,100,fid);
    char* line2 = line;
    List_rad[0]= strtod(line2,&line2);
    List_rad[1] = strtod(line2,&line2);
    List_rad[2]= strtod(line2,&line2);
    List_rad[3] = strtod(line2,&line2);

    return List_rad[3];
}


void ScaLBL_ColorModel::Create() {
    blank = std::shared_ptr<Domain>(new Domain(domain_db,comm));
    for (int i = 0; i < Nx*Ny*Nz; i++) blank->id[i] = 1;
    blank->CommInit();
    
    //.........................................................
    // Initialize communication structures in averaging domain
    for (int i=0; i<Nx*Ny*Nz; i++) Dm->id[i] = Mask->id[i];
    Mask->CommInit();
    Np=Mask->PoreCount();
  
    //...........................................................................
    ScaLBL_Comm  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Mask));
    ScaLBL_Comm_Regular  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(blank));
    
    int Npad=(Np/16 + 2)*16;
    Map.resize(Nx,Ny,Nz);       Map.fill(-2);
    /* Fill in the qDistances arrays */
    int Nq = 18;

    auto neighborList= new int[Nq*Npad];
    auto interpolationList= new int[Nq*Npad];
    auto scalarList= new int[Nq*Npad];
    
    //for (int i= 0; i < N; i++) (i) = Mask->id[i];
    
    Np = ScaLBL_Comm->MemoryOptimizedLayoutAA_LIBB(Map,neighborList,interpolationList,scalarList,Mask->id,Np,Nx,Nx*Ny);
    MPI_Barrier(comm);
    
    //......................device distributions.................................
    dist_mem_size = Np*sizeof(double);
    neighborSize=Nq*(Np*sizeof(int));
    int libbSize=Nq*(Np*sizeof(double));
    int signSize = N*sizeof(double);
   // int qsize = 18*N*sizeof(double);
    //...........................................................................
       ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
    ScaLBL_AllocateDeviceMemory((void **) &InterpolationList, neighborSize);

    ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);

    ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &fq2, 19*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &savedfq, 19*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &savedAq, 1*sizeof(double));//19*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &savedBq, 1*sizeof(double));//19*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &Aq, 1*sizeof(double));//19*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &Bq, 1*sizeof(double));//19*dist_mem_size);

//    ScaLBL_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
    ScaLBL_AllocateDeviceMemory((void **) &DenA,  N*sizeof(double));
    ScaLBL_AllocateDeviceMemory((void **) &DenB,  N*sizeof(double));
    ScaLBL_AllocateDeviceMemory((void **) &DenA2, N*sizeof(double));
    ScaLBL_AllocateDeviceMemory((void **) &DenB2, N*sizeof(double));
    ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &Phi2, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &Pressure, sizeof(double)*Np);
    
    ScaLBL_AllocateDeviceMemory((void **) &Velx, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &Vely, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &Velz, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &Velx2, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &Vely2, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &Velz2, sizeof(double)*N);
    
    ScaLBL_AllocateDeviceMemory((void **) &GradSDsX, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &GradSDsY, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &GradSDsZ, sizeof(double)*N);
    
    ScaLBL_AllocateDeviceMemory((void **) &GradPhiX, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &GradPhiY, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &GradPhiZ, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &CField, sizeof(double)*N);

    ScaLBL_AllocateDeviceMemory((void **) &MaskDomain, sizeof(char)*N);
    ScaLBL_AllocateDeviceMemory((void **) &LIBBqA, libbSize);
    ScaLBL_AllocateDeviceMemory((void **) &LIBBqBC,libbSize);
    ScaLBL_AllocateDeviceMemory((void **) &LIBBqD, libbSize);
    ScaLBL_AllocateDeviceMemory((void **) &VFmask, sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &NormalToSolid_X,   sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &NormalToSolid_Y,   sizeof(double)*N);
    ScaLBL_AllocateDeviceMemory((void **) &NormalToSolid_Z,   sizeof(double)*N);
    
    ScaLBL_CopyToDevice(NormalToSolid_X,   gradphix,          sizeof(double)*N);
    ScaLBL_CopyToDevice(NormalToSolid_Y,   gradphiy,          sizeof(double)*N);
    ScaLBL_CopyToDevice(NormalToSolid_Z,   gradphiz,          sizeof(double)*N);
    ScaLBL_CopyToDevice(VFmask,            vfmask,            sizeof(double)*N);
    fflush(stdout);
    
    ScaLBL_Comm_Regular->SendHalo(VFmask);
ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
    ScaLBL_Comm_Regular->RecvHalo(VFmask);
    ScaLBL_DeviceBarrier();    MPI_Barrier(comm);
	
    ScaLBL_Comm_Regular->SendHalo(NormalToSolid_X);
ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
    ScaLBL_Comm_Regular->RecvHalo(NormalToSolid_X);

    ScaLBL_Comm_Regular->SendHalo(NormalToSolid_Y);
ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
    ScaLBL_Comm_Regular->RecvHalo(NormalToSolid_Y);

    ScaLBL_Comm_Regular->SendHalo(NormalToSolid_Z);
ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
    ScaLBL_Comm_Regular->RecvHalo(NormalToSolid_Z);

   double *tmpfq; tmpfq = new double[19*Np];
        for (int a = 0; a < 19*Np; a++) tmpfq[a] = 0;

        ScaLBL_CopyToDevice(fq, tmpfq, 19*sizeof(double)*Np);
        ScaLBL_CopyToDevice(fq2, tmpfq, 19*sizeof(double)*Np);
        ScaLBL_CopyToDevice(savedfq, tmpfq, 19*sizeof(double)*Np);

delete[] tmpfq; 
    
    fflush(stdout);
    
    activeMap = new int[Np];
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                int idx=Map(i,j,k);
                if (!(idx < 0))
                    activeMap[idx] = k*Nx*Ny+j*Nx+i;
            }
        }
    }
    // check that TmpMap is valid
    for (int idx=0; idx<ScaLBL_Comm->LastExterior(); idx++){
        int n = activeMap[idx];
        if (n > Nx*Ny*Nz){
            printf("Bad value! idx=%i \n",n);
            activeMap[idx] = Nx*Ny*Nz-1;
        }
    }
    for (int idx=ScaLBL_Comm->FirstInterior(); idx<ScaLBL_Comm->LastInterior(); idx++){
        int n = activeMap[idx];
        if (n > Nx*Ny*Nz){
            printf("Bad value! idx=%i \n",n);
            activeMap[idx] = Nx*Ny*Nz-1;
        }
    }
    ScaLBL_CopyToDevice(dvcMap, activeMap, sizeof(int)*Np);
    ScaLBL_DeviceBarrier();
    
    TemporaryMap.resize(Nx,Ny,Nz); TemporaryMap.fill(-2);
    
    InactiveTemporaryMap.resize(Nx,Ny,Nz); InactiveTemporaryMap.fill(-2);
    Ni = 0;

    Ni = ScaLBL_Comm->MemoryOptimizedInactiveLayout(InactiveTemporaryMap, Mask->id, vfmask, Ni);
    ScaLBL_AllocateDeviceMemory((void **) &dvcInactiveMap, sizeof(int)*Ni);

    DoubleArray TempDomain; TempDomain.resize(Nx,Ny,Nz);
    TempDomain.fill(2);
    for (int k=1; k<Nz-1; k++)
    for (int j=1; j<Ny-1; j++)
    for (int i=1; i<Nx-1; i++){
        size_t n = k*Nx*Ny+j*Nx+i;
        TempDomain(n) = (double)(Mask->id[n]);
    }
    Dm->CommunicateMeshHalo(TempDomain);
    for (int n = 0; n < Nx*Ny*Nz; n++) Mask->id[n] = (char)TempDomain(n);
    
    inactiveMap=new int[Ni];
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                int idx=InactiveTemporaryMap(i,j,k);
                if (!(idx < 0)) inactiveMap[idx] = k*Nx*Ny+j*Nx+i;
            }
        }
    }

    for (int idx=0; idx<ScaLBL_Comm->LastInactiveExterior(); idx++){
        int n = inactiveMap[idx];
        if (n > Nx*Ny*Nz){
            printf("Inactive Site Bad value! idx=%i \n",n);
            inactiveMap[idx] = Nx*Ny*Nz-1;
        }
    }
    for (int idx=ScaLBL_Comm->FirstInactiveInterior(); idx<ScaLBL_Comm->LastInactiveInterior(); idx++){
        int n = inactiveMap[idx];
        if (n > Nx*Ny*Nz){
            printf("Inactive Site Bad value! idx=%i \n",n);
            inactiveMap[idx] = Nx*Ny*Nz-1;
        }
    }

    ScaLBL_CopyToDevice(dvcInactiveMap, inactiveMap, sizeof(int)*Ni);
    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    

    TemporaryMap.fill(-2);
    Nsb = 0;

    Nsb = ScaLBL_Comm->MemoryOptimizedSBLayout(TemporaryMap, Mask->id, vfmask, Nsb);
    ScaLBL_AllocateDeviceMemory((void **) &dvcSBMap, sizeof(int)*Nsb);
    
    TempDomain.fill(2);
    for (int k=1; k<Nz-1; k++)
    for (int j=1; j<Ny-1; j++)
    for (int i=1; i<Nx-1; i++){
        size_t n = k*Nx*Ny+j*Nx+i;
        TempDomain(n) = (double)(Mask->id[n]);
    }
    
    Dm->CommunicateMeshHalo(TempDomain);
    
    for (int n = 0; n < Nx*Ny*Nz; n++) Mask->id[n] = (char)TempDomain(n);

    sbMap = new int[Nsb];
    for (int k=1; k<Nz-1; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                int idx=TemporaryMap(i,j,k);
                if (!(idx < 0))
                    sbMap[idx] = k*Nx*Ny+j*Nx+i;
            }
        }
    }

    for (int idx=0; idx<ScaLBL_Comm->LastSBExterior(); idx++){
        int n = sbMap[idx];
        if (n > Nx*Ny*Nz){
            printf("Bad value! idx=%i \n",n);
            sbMap[idx] = Nx*Ny*Nz-1;
        }
    }
    for (int idx=ScaLBL_Comm->FirstSBInterior(); idx<ScaLBL_Comm->LastSBInterior(); idx++){
        int n = sbMap[idx];
        if (n > Nx*Ny*Nz){
            printf("Bad value! idx=%i \n",n);
            sbMap[idx] = Nx*Ny*Nz-1;
        }
    }
    ScaLBL_CopyToDevice(dvcSBMap, sbMap, sizeof(int)*Nsb);
    ScaLBL_DeviceBarrier();
    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
    ScaLBL_CopyToDevice(InterpolationList, interpolationList, neighborSize);

    size_t read_file;
    auto libbqA= new double[Nq*Np];
    sprintf(LocalRankFilename,"libbqA.%05i",rank);
    FILE *AFILE = fopen(LocalRankFilename,"rb");
    read_file = fread(libbqA,8,Nq*Np,AFILE);
    if (read_file != size_t(Nq*Np)) printf("Error reading AFILE \n");
    fclose(AFILE);
    auto libbqBC= new double[Nq*Np];
    sprintf(LocalRankFilename,"libbqBC.%05i",rank);
    FILE *BCFILE = fopen(LocalRankFilename,"rb");
    read_file = fread(libbqBC,8,Nq*Np,BCFILE);
    if (read_file != size_t(Nq*Np)) printf("Error reading BCFILE \n");
    fclose(BCFILE);
    auto libbqD= new double[Nq*Np];
    sprintf(LocalRankFilename,"libbqD.%05i",rank);
    FILE *DFILE = fopen(LocalRankFilename,"rb");
    read_file = fread(libbqD,8,Nq*Np,DFILE);
    if (read_file != size_t(Nq*Np)) printf("Error reading DFILE \n");
    fclose(DFILE);
    
    
    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    ScaLBL_CopyToDevice(LIBBqA, libbqA, libbSize);
    ScaLBL_CopyToDevice(LIBBqBC, libbqBC, libbSize);
    ScaLBL_CopyToDevice(LIBBqD, libbqD, libbSize);
    
    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

    //for (int i = 0; i < N; i++) Averages->ID(i) = Mask->id[i];

    delete[] gradphix;
    delete[] gradphiy;
    delete[] gradphiz;
}

void ScaLBL_ColorModel::Initialize() {
    
   

    ScaLBL_D3Q19_Init(fq, Np);
    ScaLBL_D3Q19_Init(savedfq,Np);
    ScaLBL_D3Q19_Init(fq2,Np);    


    if (Restart == false) {
        
        
        size_t read_file;
        auto fluidphaseID= new char[N];
        sprintf(LocalRankFilename,"FluidPhaseID.%05i",rank);
        FILE *FlUIDPHASEIDFILE = fopen(LocalRankFilename,"rb");
        read_file = fread(fluidphaseID,1,N,FlUIDPHASEIDFILE);
        if (read_file != size_t(N)) printf("Error reading FlUIDPHASEIDFILE \n");
        fclose(FlUIDPHASEIDFILE);
        
        
        auto PhaseLabel = new double[N];
        auto DenALabel = new double[N];
        auto DenBLabel = new double[N];
        for (int k=1; k<Nz-1; k++)
        for (int j=1; j<Ny-1; j++)
        for (int i=1; i<Nx-1; i++) {
            int n = i + j*Nx + k*Nx*Ny;
            if (Mask->id[n] == 0) { PhaseLabel[n] = -1.0; DenALabel[n] = 0; DenBLabel[n] = 1;}
            if (Mask->id[n] == 1) { PhaseLabel[n] = -1.0; DenALabel[n] = 0; DenBLabel[n] = 1;}
            if (Mask->id[n] == 2) { PhaseLabel[n] = -1.0; DenALabel[n] = 0; DenBLabel[n] = 1;}
            if (Mask->id[n] == 3) { PhaseLabel[n] = -1.0; DenALabel[n] = 0; DenBLabel[n] = 1;}
            if (Mask->id[n] == 4) { PhaseLabel[n] = -1.0; DenALabel[n] = 0; DenBLabel[n] = 1;}
        }
        
        
        {
            int n; int VALUE;
        for (int k=1; k<Nz-1; k++)
        for (int j=1; j<Ny-1; j++)
        for (int i=1; i<Nx-1; i++) {
            n = i + j*Nx + k*Nx*Ny;
            VALUE = fluidphaseID[n];
            if (VALUE == 1) { PhaseLabel[n] = 1;  DenALabel[n] = 1; DenBLabel[n] = 0; }
            if (VALUE == 2) { PhaseLabel[n] = -1; DenALabel[n] = 0; DenBLabel[n] = 1; }
        }
        }
//    {
//        int n; int VALUE; double DVALUE;
//        printf("FluidPhaseID:\n");
//        for (int k=18;k<19;k++){
//            for (int i=1;i<Nx-1;i++){
//                for (int j=1;j<Ny-1;j++){
//                    n=k*(Nx)*(Ny)+j*(Nx)+i;
//                    // VALUE = fluidphaseID[n];
//                    DVALUE = PhaseLabel[n];
//                    // printf("%d ",VALUE);
//                    printf("%.0f ",DVALUE);
//                }
//                printf("\n");
//            }
//            printf("\n\n");
//        }
//    }
//
//
//                if ( i > 3 && i <= 10 ) {
//                    if ( j > 3 && j <= 10 ) {
//                        if ( k >= 1  && k <= 12 ) {
//
//                            if (PhaseLabel[n] == -1) {PhaseLabel[n] = 1; DenALabel[n] = 1; DenBLabel[n] = 0;  }
//                          //  extra_counter = extra_counter + 1;
//                        }
//                    }
//                }
//
//
//        }



        for (size_t n = 0; n < N; n++) Averages->Phase(n)  = PhaseLabel[n];
        for (size_t n = 0; n < N; n++) Averages->ID(n)  = fluidphaseID[n];
        ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));
        ScaLBL_CopyToDevice(Phi2, PhaseLabel, N*sizeof(double));
        ScaLBL_CopyToDevice(DenA, DenALabel, N*sizeof(double));
        ScaLBL_CopyToDevice(DenB, DenBLabel, N*sizeof(double));
        ScaLBL_CopyToDevice(DenA2, DenALabel, N*sizeof(double));
        ScaLBL_CopyToDevice(DenB2, DenBLabel, N*sizeof(double));    
    
	delete[] PhaseLabel;
        delete[] DenALabel;
        delete[] DenBLabel;

        ScaLBL_Comm_Regular->SendHalo(Phi);
        ScaLBL_Comm_Regular->RecvHalo(Phi);

        ScaLBL_Comm_Regular->SendHalo(DenA);
        ScaLBL_Comm_Regular->RecvHalo(DenA);

        ScaLBL_Comm_Regular->SendHalo(DenB);
        ScaLBL_Comm_Regular->RecvHalo(DenB);

	ScaLBL_Comm_Regular->SendHalo(DenA2);
        ScaLBL_Comm_Regular->RecvHalo(DenA2);

        ScaLBL_Comm_Regular->SendHalo(DenB2);
        ScaLBL_Comm_Regular->RecvHalo(DenB2);
        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    }





    /* Overwrite things if restart is true */
    if (Restart == true) {

        if (rank==0){
            printf("Reading restart file! \n");
            ifstream restart("Restart.txt");
            if (restart.is_open()){
                restart  >> timestep;
                printf("Restarting from timestep %i \n",timestep);
            }
            else{
                printf("WARNING:No Restart.txt file, setting timestep=0 \n");
                timestep=0;
            }
        }

        MPI_Bcast(&timestep,1,MPI_INT,0,comm);
        // Read in the restart file to CPU buffers
        double *cPhi, *cfq, *cDenA, *cDenB, *cVx, *cVy, *cVz, *cGPhiX, *cGPhiY, *cGPhiZ, *cCField;
        cPhi = new double[N];
        cDenA = new double[N];
        cDenB = new double[N];
        cfq = new double[19*Np];
        cVx = new double[N];
        cVy = new double[N];
        cVz = new double[N];
        cGPhiX = new double[N];
        cGPhiY = new double[N];
        cGPhiZ = new double[N];
        cCField = new double[N];

        char LocalRestartFile[40];
        sprintf(LocalRestartFile,"Restart.%05i",rank);
        ifstream File(LocalRestartFile,ios::binary);
        double value,va,vb,vphase,vx,vy,vz,fz;

        for (int n=0; n < N; n++){
            File.read((char*) &va, sizeof(va));
            cDenA[n] = va;
        }

        for (int n=0; n < N; n++){
            File.read((char*) &vb, sizeof(vb));
            cDenB[n] = vb;
        }


        for (int n=0; n<N; n++){
            File.read((char*) &vphase, sizeof(value));
            cPhi[n] = vphase;
        }
        for (int n=0; n<Np; n++){
            // Read the distributions
            for (int q=0; q<19; q++){
                File.read((char*) &value, sizeof(value));
                cfq[q*Np+n] = value;
            }
        }
        for (int n=0; n < N; n++){
            File.read((char*) &vx, sizeof(vx));
            cVx[n] = vx;
        }
        for (int n=0; n < N; n++){
            File.read((char*) &vy, sizeof(vy));
            cVy[n] = vy;
        }
        for (int n=0; n < N; n++){
            File.read((char*) &vz, sizeof(vz));
            cVz[n] = vz;
        }

        File.read((char*) &fz, sizeof(fz));

        Fz = fz;

        for (int n=0; n < N; n++){
            File.read((char*) &vx, sizeof(vx));
            cGPhiX[n] = vx;
        }
        for (int n=0; n < N; n++){
            File.read((char*) &vy, sizeof(vy));
            cGPhiY[n] = vy;
        }
        for (int n=0; n < N; n++){
            File.read((char*) &vz, sizeof(vz));
            cGPhiZ[n] = vz;
        }
        for (int n=0; n < N; n++){
            File.read((char*) &value, sizeof(value));
            cCField[n] = value;
        }




        File.close();


        // Copy the restart data to the GPU
        ScaLBL_CopyToDevice(DenA2,cDenA,N*sizeof(double));
        ScaLBL_CopyToDevice(DenB2,cDenB,N*sizeof(double));
        ScaLBL_CopyToDevice(fq,cfq,19*Np*sizeof(double));
        ScaLBL_CopyToDevice(fq2,cfq,19*Np*sizeof(double));
        ScaLBL_CopyToDevice(savedfq,cfq,19*Np*sizeof(double));
        ScaLBL_CopyToDevice(Phi,cPhi,N*sizeof(double));
        ScaLBL_CopyToDevice(Velx2,cVx,N*sizeof(double));
        ScaLBL_CopyToDevice(Vely2,cVy,N*sizeof(double));
        ScaLBL_CopyToDevice(Velz2,cVz,N*sizeof(double));
        ScaLBL_CopyToDevice(GradPhiX,cGPhiX,N*sizeof(double));
        ScaLBL_CopyToDevice(GradPhiY,cGPhiY,N*sizeof(double));
        ScaLBL_CopyToDevice(GradPhiZ,cGPhiZ,N*sizeof(double));
        ScaLBL_CopyToDevice(CField,cCField,N*sizeof(double));
        ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
        delete[] cDenA;
        delete[] cDenB;
        delete[] cfq;
        delete[] cPhi;
        delete[] cVx;
        delete[] cVy;
        delete[] cVz;
        delete[] cGPhiX;
        delete[] cGPhiY;
        delete[] cGPhiZ;
        delete[] cCField;
    }

    MPI_Bcast(&Fz,1,MPI_DOUBLE,0,comm);
    MPI_Bcast(&Fy,1,MPI_DOUBLE,0,comm);
    MPI_Bcast(&Fx,1,MPI_DOUBLE,0,comm);
    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
 

    if (BoundaryCondition >0 ){
            if (Dm->kproc()==0){
                ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
                ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
                ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
                ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,3);
                ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,4);
            }
            if (Dm->kproc() == nprocz-1){
                ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-1);
                ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-2);
                ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-3);
                ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-4);
                ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-5);
            }
        } // BC > 0
}


void ScaLBL_ColorModel::Run(string filename){

    printf("COLORMODEL: offsetdistance=%f\n",offsetdistance);
    Averages->offset_distance = offsetdistance;
  //   for (int i = 0; i < 18*Np; i++) { LIBBqA[i] = 0; LIBBqBC[i] = 1; LIBBqD[i] = 0; }

    db = std::make_shared<Database>( filename );
    domain_db = db->getDatabase( "Domain" );

    int WBCFlag = 0;
    if ( domain_db->keyExists( "WBCFlag" ) ) {
        WBCFlag = domain_db->getScalar<int>( "WBCFlag" ); //0: Original, 1: Situational off (add if statement)
    }

    int nprocs=nprocx*nprocy*nprocz;
    const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
    
    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
    
    double starttime,stoptime,cputime;
    ScaLBL_DeviceBarrier();
    MPI_Barrier(comm);
   
    
  //  double Ca_previous = 0.0;
    double Ca, muA, momentum_z, relError;
    Ca = muA = momentum_z = 0.0;
    bool Regular = false;
    
    
    runAnalysis analysis( analysis_db, rank_info, ScaLBL_Comm, Dm, Np, Regular, beta, Map );
    
   
    
  
    
//    printf("Before run5\n");
//    while (timestep < 300) {
//        timestep++;
//        analysis.run6(0, *Averages, Np, LIBBqBC);
//    analysis.finish();
//    }
    
    ScaLBL_Comm_Regular->SendHalo(VFmask);
    ScaLBL_Comm_Regular->RecvHalo(VFmask);

    MPI_Barrier(comm);  ScaLBL_DeviceBarrier();
    ScaLBL_Comm_Regular->SendHalo(Phi);
    ScaLBL_Comm_Regular->RecvHalo(Phi);

    ScaLBL_Comm_Regular->SendHalo(DenA);
    ScaLBL_Comm_Regular->RecvHalo(DenA);
    ScaLBL_Comm_Regular->SendHalo(DenB);
    ScaLBL_Comm_Regular->RecvHalo(DenB);

    ScaLBL_Comm_Regular->SendHalo(DenA2);
    ScaLBL_Comm_Regular->RecvHalo(DenA2);
    ScaLBL_Comm_Regular->SendHalo(DenB2);
    ScaLBL_Comm_Regular->RecvHalo(DenB2);

    ScaLBL_AllocateDeviceMemory((void **) &ActiveScalarList, 18*sizeof(int)*Np);
    auto activeScalarList = new int[18*Np];
    ScaLBL_Comm->CreateActiveMap(Map, Mask->id, Np, Nx, Nx*Ny, activeScalarList);
    ScaLBL_CopyToDevice(ActiveScalarList, activeScalarList, 18*sizeof(int)*Np);
    delete[] activeScalarList;



#ifdef WBC
    ScaLBL_CopyToDevice(MaskDomain, Mask->id, sizeof(char)*N);

    save_scalar(Phi ,Phi2 ,N);
    InitExtrapolatePhaseFieldInactive(dvcInactiveMap, MaskDomain, Phi , Phi2 , ScaLBL_Comm->FirstInactiveInterior(), ScaLBL_Comm->LastInactiveInterior(), Nx, Nx*Ny,Ni);
    InitExtrapolatePhaseFieldInactive(dvcInactiveMap, MaskDomain, Phi , Phi2 , 0, ScaLBL_Comm->LastInactiveExterior(), Nx, Nx*Ny,Ni);
    save_scalar(Phi2 ,Phi ,N);
    
   

    save_scalar(Phi ,Phi2 ,N);
    InitExtrapolateScalarField(dvcSBMap, MaskDomain, Phi , Phi2 , ScaLBL_Comm->FirstSBInterior(), ScaLBL_Comm->LastSBInterior(), Nsb, Nx, Nx*Ny);
    InitExtrapolateScalarField(dvcSBMap, MaskDomain, Phi , Phi2 , 0, ScaLBL_Comm->LastSBExterior(), Nsb, Nx, Nx*Ny);
    save_scalar(Phi2 ,Phi ,N);

    ScaLBL_AllocateDeviceMemory((void **) &InactiveScalarList, 18*sizeof(int)*Ni);
    auto inactiveScalarList = new int[18*Ni];
    ScaLBL_Comm->CreateInactiveMap(InactiveTemporaryMap, Mask->id, Ni, Nx, Nx*Ny, inactiveScalarList);
    ScaLBL_CopyToDevice(InactiveScalarList, inactiveScalarList, 18*sizeof(int)*Ni);
    delete[] inactiveScalarList;
    ScaLBL_DeviceBarrier();  MPI_Barrier(comm);

    ScaLBL_AllocateDeviceMemory((void **) &SBScalarList, 18*sizeof(int)*Nsb);
    auto sbScalarList = new int[18*Nsb];
    ScaLBL_Comm->CreateSBMap(TemporaryMap, Mask->id, Nsb, Nx, Nx*Ny, sbScalarList);
    ScaLBL_CopyToDevice(SBScalarList, sbScalarList, 18*sizeof(int)*Nsb);
    delete[] sbScalarList;

 
// double *tmp; tmp = new double[N];
//	for (int a = 0; a < N; a++) tmp[a] = 0;
//
//
//	ScaLBL_CopyToDevice(GradPhiX, tmp, sizeof(double)*N);
//	ScaLBL_CopyToDevice(GradPhiY, tmp, sizeof(double)*N);
//	ScaLBL_CopyToDevice(GradPhiZ, tmp, sizeof(double)*N);
//	ScaLBL_CopyToDevice(Velx, tmp, sizeof(double)*N);
//	ScaLBL_CopyToDevice(Vely, tmp, sizeof(double)*N);
//	ScaLBL_CopyToDevice(Velz, tmp, sizeof(double)*N);
//	ScaLBL_CopyToDevice(Velx2, tmp, sizeof(double)*N);
//	ScaLBL_CopyToDevice(Vely2, tmp, sizeof(double)*N);
//	ScaLBL_CopyToDevice(Velz2, tmp, sizeof(double)*N);
//	ScaLBL_CopyToDevice(CField, tmp, sizeof(double)*N);
//
//	delete[] tmp;



    
#endif
    starttime = MPI_Wtime();

    if (Restart == false) {
        save_scalar(DenA,DenA2,N);
        save_scalar(DenB,DenB2,N);
    }

    Averages->alpha = alpha;
    Averages->average_mu = (0.5*(tauB+tauA) -0.5)/3.;
    Averages->sphere_diameter = ReadSpherePacking();
    Averages->sphere_diameter *= double(Nz-2);

    if (Averages->sphere_diameter == 0) Averages->sphere_diameter = double(Nz-2)/3.;
    
    
//    int VALUE = 0;
//    double DVALUE = 0;
//    printf("VFmask:\n");
//    for (int i=5;i<6;i++){
//        for (int j=1;j<Nx-1;j++){
//            for (int k=1;k<Nz-1;k++){
//                int n=k*(Nx)*(Ny)+j*(Nx)+i;
//                DVALUE = Averages->VFmask(n);
//                printf("%.2f ",DVALUE);
//            }
//            printf("\n");
//        }
//        printf("\n\n");
//    }
//
   

//    while(timestep < timestepMax) {
//
//        timestep++;
//        analysis.run5(timestep, *Averages, Phi, Pressure, Velx2, Vely2, Velz2, fq, GradPhiX, GradPhiY, GradPhiZ, CField, DenA2, DenB2,Np,Fx,Fy,Fz);
//
////        analysis.run6(0, *Averages, Np, Phi);
//
////        analysis.run7(0, *Averages, Np, Phi,Velx,Vely,Velz,DenA,DenB);
//        analysis.finish();
//
//
//    }
//    timestep+=10000;
  

    while(timestep < timestepMax) {
        ScaLBL_Comm_Regular->SendHaloMany(Phi,DenA2,DenB2,GradPhiX,GradPhiY,GradPhiZ,CField,Velx2,Vely2,Velz2);
        ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
        ScaLBL_Comm_Regular->RecvHaloMany(Phi,DenA2,DenB2,GradPhiX,GradPhiY,GradPhiZ,CField,Velx2,Vely2,Velz2);

        save_scalar(DenA2,DenA,N);
        save_scalar(DenB2,DenB,N);
        Inactive_Color_LIBB(InactiveScalarList, dvcInactiveMap, DenA2 , DenB2, DenA , DenB, Phi, Velx2, Vely2, Velz2, beta,  Nx, Nx*Ny, ScaLBL_Comm->FirstInactiveInterior(), ScaLBL_Comm->LastInactiveInterior(), Ni, N, GradPhiX, GradPhiY, GradPhiZ, CField);
        Inactive_Color_LIBB(InactiveScalarList, dvcInactiveMap, DenA2, DenB2, DenA, DenB, Phi, Velx2, Vely2, Velz2, beta,  Nx, Nx*Ny, 0, ScaLBL_Comm->LastInactiveExterior(), Ni, N, GradPhiX, GradPhiY, GradPhiZ, CField);

        Inactive_Color_LIBB(SBScalarList, dvcSBMap, DenA2 , DenB2, DenA , DenB, Phi, Velx2, Vely2, Velz2, beta,  Nx, Nx*Ny, ScaLBL_Comm->FirstSBInterior(), ScaLBL_Comm->LastSBInterior(), Nsb, N, GradPhiX, GradPhiY, GradPhiZ, CField);
        Inactive_Color_LIBB(SBScalarList, dvcSBMap, DenA2, DenB2, DenA, DenB, Phi, Velx2, Vely2, Velz2, beta,  Nx, Nx*Ny, 0, ScaLBL_Comm->LastSBExterior(), Nsb, N, GradPhiX, GradPhiY, GradPhiZ, CField);
        save_scalar(DenA,DenA2,N);
        save_scalar(DenB,DenB2,N);

        ComputeGradPhi(input_angle, dvcMap, Phi, GradPhiX, GradPhiY, GradPhiZ, CField, NormalToSolid_X, NormalToSolid_Y, NormalToSolid_Z, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(),Np, WBCFlag);
        ComputeGradPhi(input_angle, dvcMap, Phi, GradPhiX, GradPhiY, GradPhiZ, CField, NormalToSolid_X, NormalToSolid_Y, NormalToSolid_Z, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(),Np, WBCFlag);
        ComputeGradPhi(input_angle, dvcInactiveMap, Phi, GradPhiX, GradPhiY, GradPhiZ, CField, NormalToSolid_X, NormalToSolid_Y, NormalToSolid_Z, Nx, Nx*Ny, ScaLBL_Comm->FirstInactiveInterior(), ScaLBL_Comm->LastInactiveInterior(),Ni, WBCFlag);
        ComputeGradPhi(input_angle, dvcInactiveMap, Phi, GradPhiX, GradPhiY, GradPhiZ, CField, NormalToSolid_X, NormalToSolid_Y, NormalToSolid_Z, Nx, Nx*Ny, 0, ScaLBL_Comm->LastInactiveExterior(),Ni, WBCFlag);




        ScaLBL_Comm->SendD3Q19AA(fq);
        ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
        ScaLBL_Comm->RecvD3Q19AA(fq);
        save_state(fq,savedfq,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(),19,Np);
        save_state(fq,savedfq,0, ScaLBL_Comm->LastExterior(),19,Np);
        ScaLBL_D3Q19_Color_LIBB(ActiveScalarList, InterpolationList, NeighborList,
                                dvcMap, fq, fq2, savedfq, Aq, Bq, DenA2, DenB2, DenA, DenB, Phi, Velx , Vely, Velz, Velx2 , Vely2, Velz2 , Pressure, rhoA, rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, N, LIBBqA, LIBBqBC, LIBBqD, GradPhiX, GradPhiY, GradPhiZ, CField);


        if (BoundaryCondition == 3){
            ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }
        if (BoundaryCondition == 4){
            din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
            ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
        }

        ScaLBL_DeviceBarrier();
        ScaLBL_D3Q19_Color_LIBB(ActiveScalarList, InterpolationList, NeighborList,
                                dvcMap, fq, fq2, savedfq, Aq, Bq, DenA2, DenB2, DenA, DenB, Phi, Velx, Vely, Velz, Velx2, Vely2, Velz2, Pressure, rhoA, rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np, N, LIBBqA, LIBBqBC, LIBBqD, GradPhiX, GradPhiY, GradPhiZ, CField);

        if (BoundaryCondition > 0) {
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, DenA, DenB, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, DenA, DenB, outletA, outletB);
           }
        if (BoundaryCondition > 0) {
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, DenA2, DenB2, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, DenA2, DenB2, outletA, outletB);
           }


        timestep++;

        ScaLBL_Comm_Regular->SendHaloMany(Phi,DenA,DenB,GradPhiX,GradPhiY,GradPhiZ,CField,Velx,Vely,Velz);
        ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
        ScaLBL_Comm_Regular->RecvHaloMany(Phi,DenA,DenB,GradPhiX,GradPhiY,GradPhiZ,CField,Velx,Vely,Velz);

        save_scalar(DenA,DenA2,N);
        save_scalar(DenB,DenB2,N);
        Inactive_Color_LIBB(InactiveScalarList, dvcInactiveMap, DenA , DenB, DenA2 , DenB2, Phi, Velx2, Vely2, Velz2, beta,  Nx, Nx*Ny, ScaLBL_Comm->FirstInactiveInterior(), ScaLBL_Comm->LastInactiveInterior(), Ni, N, GradPhiX, GradPhiY, GradPhiZ, CField);
        Inactive_Color_LIBB(InactiveScalarList, dvcInactiveMap, DenA, DenB, DenA2, DenB2, Phi, Velx2, Vely2, Velz2, beta,  Nx, Nx*Ny, 0, ScaLBL_Comm->LastInactiveExterior(), Ni, N, GradPhiX, GradPhiY, GradPhiZ, CField);
        Inactive_Color_LIBB(SBScalarList, dvcSBMap, DenA , DenB, DenA2 , DenB2, Phi, Velx2, Vely2, Velz2, beta,  Nx, Nx*Ny, ScaLBL_Comm->FirstSBInterior(), ScaLBL_Comm->LastSBInterior(), Nsb, N, GradPhiX, GradPhiY, GradPhiZ, CField);
        Inactive_Color_LIBB(SBScalarList, dvcSBMap, DenA, DenB, DenA2, DenB2, Phi, Velx2, Vely2, Velz2, beta,  Nx, Nx*Ny, 0, ScaLBL_Comm->LastSBExterior(), Nsb, N, GradPhiX, GradPhiY, GradPhiZ, CField);
        save_scalar(DenA2,DenA,N);
        save_scalar(DenB2,DenB,N);

        ComputeGradPhi(input_angle, dvcMap, Phi, GradPhiX, GradPhiY, GradPhiZ, CField, NormalToSolid_X, NormalToSolid_Y, NormalToSolid_Z, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(),Np, WBCFlag);
        ComputeGradPhi(input_angle, dvcMap, Phi, GradPhiX, GradPhiY, GradPhiZ, CField, NormalToSolid_X, NormalToSolid_Y, NormalToSolid_Z, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(),Np, WBCFlag);
        ComputeGradPhi(input_angle, dvcInactiveMap, Phi, GradPhiX, GradPhiY, GradPhiZ, CField, NormalToSolid_X, NormalToSolid_Y, NormalToSolid_Z, Nx, Nx*Ny, ScaLBL_Comm->FirstInactiveInterior(), ScaLBL_Comm->LastInactiveInterior(),Ni, WBCFlag);
        ComputeGradPhi(input_angle, dvcInactiveMap, Phi, GradPhiX, GradPhiY, GradPhiZ, CField, NormalToSolid_X, NormalToSolid_Y, NormalToSolid_Z, Nx, Nx*Ny, 0, ScaLBL_Comm->LastInactiveExterior(),Ni, WBCFlag);




        ScaLBL_Comm->SendD3Q19AA(fq2);
        ScaLBL_DeviceBarrier();  MPI_Barrier(comm);
        ScaLBL_Comm->RecvD3Q19AA(fq2);
        save_state(fq2,savedfq,ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(),19,Np);
        save_state(fq2,savedfq,0, ScaLBL_Comm->LastExterior(),19,Np);
        ScaLBL_D3Q19_Color_LIBB(ActiveScalarList, InterpolationList, NeighborList,
                                dvcMap, fq2, fq, savedfq, Aq, Bq, DenA , DenB , DenA2 , DenB2 , Phi, Velx2 , Vely2, Velz2, Velx, Vely, Velz, Pressure, rhoA, rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm->FirstInterior(), ScaLBL_Comm->LastInterior(), Np, N, LIBBqA, LIBBqBC, LIBBqD, GradPhiX, GradPhiY, GradPhiZ, CField);
        ScaLBL_DeviceBarrier();

        if (BoundaryCondition == 3){
               ScaLBL_Comm->D3Q19_Pressure_BC_z(NeighborList, fq2, din, timestep);
               ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq2, dout, timestep);
           }
       if (BoundaryCondition == 4){
           din = ScaLBL_Comm->D3Q19_Flux_BC_z(NeighborList, fq2, flux, timestep);
           ScaLBL_Comm->D3Q19_Pressure_BC_Z(NeighborList, fq2, dout, timestep);
       }

        ScaLBL_D3Q19_Color_LIBB(ActiveScalarList, InterpolationList, NeighborList,
                                dvcMap, fq2, fq, savedfq, Aq, Bq, DenA, DenB, DenA2, DenB2, Phi, Velx2, Vely2, Velz2, Velx, Vely, Velz, Pressure, rhoA, rhoB, tauA, tauB, alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm->LastExterior(), Np, N, LIBBqA, LIBBqBC, LIBBqD, GradPhiX, GradPhiY, GradPhiZ, CField);

        if (BoundaryCondition > 0){
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, DenA, DenB, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, DenA, DenB, outletA, outletB);
        }
        if (BoundaryCondition > 0){
            ScaLBL_Comm->Color_BC_z(dvcMap, Phi, DenA2, DenB2, inletA, inletB);
            ScaLBL_Comm->Color_BC_Z(dvcMap, Phi, DenA2, DenB2, outletA, outletB);
        }

        timestep++;
        analysis.run5(timestep, *Averages, Phi, Pressure, Velx2, Vely2, Velz2, fq, GradPhiX, GradPhiY, GradPhiZ, CField, DenA2, DenB2,Np,Fx,Fy,Fz);
        analysis.finish();

    }

    ScaLBL_DeviceBarrier();
    MPI_Barrier(comm);
    stoptime = MPI_Wtime();
    
    if (rank == 0) cout << "Finished lattice Boltzmann iterations." << endl;
    if (rank==0) cout << std::setprecision(3);
    
    analysis.finish();
    PROFILE_STOP("Loop");
    PROFILE_SAVE("lbpm_color_simulator",1);
    //************************************************************************
    
    // Compute the walltime per timestep
    cputime = (stoptime - starttime)/timestep;
    // Performance obtained from each node
    double MLUPS = double(Np+Ni+Nsb)/cputime/1000000;
    if (rank==0) {
        printf("\nPERFORMANCE METRICS:\n");
   printf("CPU time = %f \n", cputime);
    printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
    MLUPS *= nprocs;
   printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
   printf("\n");
   cout << "End of lattice Boltzmann simulator." << endl;
    }
    // End of the lattice Boltzmann simulator
} // LBM iteration function

double ScaLBL_ColorModel::MorphInit(const double beta,
                                    const double morph_delta){
    return 0;
}

void ScaLBL_ColorModel::WriteDebug(){}
