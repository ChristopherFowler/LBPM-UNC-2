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

#include "common/Domain.h"
#include "common/Array.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "common/ScaLBL.h"
#include "analysis/TwoPhase.h"
#include "analysis/runAnalysis.h"

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
            
            //FluidPhase
            printf("Read data from %s \n",READFILE_FluidPhaseID.c_str());
            FILE *FLUIDPHASESEGDAT = fopen(READFILE_FluidPhaseID.c_str(),"rb");
            if (FLUIDPHASESEGDAT==NULL) ERROR("Error reading fluid phase id data, null \n");
            int ReadFluidPhaseSeg;
            ReadFluidPhaseSeg=fread(FluidPhaseSegData,1,SIZE,FLUIDPHASESEGDAT);
            if (ReadFluidPhaseSeg != SIZE) printf("Error reading fluid phase id data, size \n");
            fclose(FLUIDPHASESEGDAT);
            
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
									if (x<0) 	x=0;
									if (!(x<global_Nx))	x=global_Nx-1;
									if (y<0) 	y=0;
									if (!(y<global_Ny))	y=global_Ny-1;
									if (z<0) 	z=0;
									if (!(z<global_Nz))	z=global_Nz-1;
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
										if (x<0) 	x=0;
										if (!(x<global_Nx))	x=global_Nx-1;
										if (y<0) 	y=0;
										if (!(y<global_Ny))	y=global_Ny-1;
										if (z<0) 	z=0;
										if (!(z<global_Nz))	z=global_Nz-1;
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

       
     

		std::array<size_t,3> subdivide = { 1, 1, 1 };

		if ( analysis_db->keyExists( "subdivide" ) ) { 
			auto tmp = analysis_db->getVector<size_t>( "subdivide" );  
			subdivide = { tmp[0], tmp[1], tmp[2] }; 
		}

		Averages = std::shared_ptr<TwoPhase> ( new TwoPhase(Dm,subdivide) );

        
        
		for (int k=0;k<nz+2;k++){
			for (int j=0;j<ny+2;j++){
				for (int i=0;i<nx+2;i++){
					n = k*(nx+2)*(ny+2) + j*(nx+2) + i;
					Averages->SDs(n) = sdmc[n];
				}
			}
		}
        
//        printf("sdmc used to compute VFmask:\n");
//        for (int i=5;i<6;i++){
//            for (int j=1;j<nx+1;j++){
//                for (int k=1;k<nz+1;k++){
//                    int n=k*(nx+2)*(ny+2)+j*(nx+2)+i;
//                    DVALUE = sdmc[n];
//                    printf("%.2f ",DVALUE);
//                }
//                printf("\n");
//            }
//            printf("\n\n");
//        }

        
        Dm->CommInit();

        int xdim,ydim,zdim;
        xdim=Dm->Nx-2;
        ydim=Dm->Ny-2;
        zdim=Dm->Nz-2;
        fillHalo<double> fillData(Dm->Comm, Dm->rank_info,{xdim,ydim,zdim},{1,1,1},0,1);

      
        //Averages->SDs.fill(-0.00001);
      //  Averages->SDn.fill(-1);
      //
        for (int i = 0; i < N; i++) {

            //Change ID based on VFmask
            // if (Averages->VFmask(i) > 0.5){
            //     id[i] = 0;
            // } else{
            //     id[i] = 2;
            // }

            Mask->id[i] = id[i];
            if (id[i] == 0){
                Averages->SDn(i) = 1.0;
            } else {
                Averages->SDn(i) = 1.0;
            }
        }
        
        for (int k=1; k<nz+1; k++){
            for (int j=1; j<ny+1; j++){
                for (int i=1; i<nx+1; i++){
                    size_t n = k*(nx+2)*(ny+2)+j*(nx+2)+i;
                    if (Averages->SDn(n) == -1) {
                        Averages->SDn(n) = 1;
//                        if (k == 1) {
//                            Averages->SDn(n) = 1;
//                        }
//                        
//                        
//                        if (k == 3) {
//                            Averages->SDn(n) = 1;
//                        }
//                        
//                        if (k == 6) {
//                            Averages->SDn(n) = 1;
//                        }
                    }
                }
            }
        }
        
        fillData.fill(Averages->SDn);
        
        
        Averages->ComputeVolumeFraction(Averages->SDs, Averages->SDs_x, Averages->SDs_y, Averages->SDs_z, Averages->VFmask, true, rmin, rmax, geometry);
        
        
//        printf("Averages->GradPhiZ:\n");
//        for (int i=5;i<6;i++){
//            for (int j=1;j<nx+1;j++){
//                for (int k=1;k<nz+1;k++){
//                    int n=k*(nx+2)*(ny+2)+j*(nx+2)+i;
//                    DVALUE = Averages->GradPhiZ(n);
//                    printf("%.2f ",DVALUE);
//                }
//                printf("\n");
//            }
//            printf("\n\n");
//        }

     
        fillData.fill(Averages->SDs);
        fillData.fill(Averages->SDs_x);
        fillData.fill(Averages->SDs_y);
        fillData.fill(Averages->SDs_z);
        fillData.fill(Averages->VFmask);

        

		for (int k=0;k<nz+2;k++){
			for (int j=0;j<ny+2;j++){
				for (int i=0;i<nx+2;i++){
					n = k*(nx+2)*(ny+2) + j*(nx+2) + i;
					Averages->SDs(n) = sd[n];
				}
			}
		}

//                printf("sdmc for MC:\n");
//                for (int i=5;i<6;i++){
//                    for (int j=1;j<nx+1;j++){
//                        for (int k=1;k<nz+1;k++){
//                            int n=k*(nx+2)*(ny+2)+j*(nx+2)+i;
//                            DVALUE = sdmc[n];
//                            printf("%.2f ",DVALUE);
//                        }
//                        printf("\n");
//                    }
//                    printf("\n\n");
//                }

		
//        printf("VFmask:\n");
//        for (int i=5;i<6;i++){
//            for (int j=1;j<nx+1;j++){
//                for (int k=1;k<nz+1;k++){
//                    int n=k*(nx+2)*(ny+2)+j*(nx+2)+i;
//                    DVALUE = Averages->VFmask(n);
//                    printf("%.2f ",DVALUE);
//                }
//                printf("\n");
//            }
//            printf("\n\n");
//        }
//
//        printf("id:\n");
//        for (int i=5;i<6;i++){
//            for (int j=1;j<nx+1;j++){
//                for (int k=1;k<nz+1;k++){
//                    int n=k*(nx+2)*(ny+2)+j*(nx+2)+i;
//                    VALUE = id[n];
//                    printf("%d ",VALUE);
//                }
//                printf("\n");
//            }
//            printf("\n\n");
//        }
//

		int q = 0;
		for (int i = 0; i < N; i++) {
			Mask->id[i] = id[i];
            Averages->ID(i) = id[i];
			q = 4;
			Averages->Vel_x(i) = libba[i+q*N];
			Averages->Vel_y(i) = libbbc[i+q*N];
			Averages->Vel_z(i) = libbd[i+q*N];

			// q = 5;
			// Averages->Vel_z(i) = libba[i+q*N];
			// Averages->Vel_z(i) = libbbc[i+q*N];
			// Averages->Vel_z(i) = libbd[i+q*N];
		}

	

		//Initialize visualization/make compressed data map
		int Nq = 18;
        int Np=Mask->PoreCount();
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
        
        


        delete[] TmpMap;
     
        
		//Produce visualization
		analysis.run6(0, *Averages, Np, comp_libbqA); 
    	analysis.finish();

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

		sprintf(LocalRankFilename,"SignDist.%05i",RANK);
		FILE *SD = fopen(LocalRankFilename,"wb");
		fwrite(sd,8,(nx+2)*(ny+2)*(nz+2),SD);
		fclose(SD);
        
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

        for (int n = 0; n < N; n++) TemporaryField[n] = Averages->GradPhiX(n);
        sprintf(LocalRankFilename,"%s%s","ns_x.",LocalRankString);
        FILE *NSX = fopen(LocalRankFilename,"wb");
        if (NSX==NULL) ERROR("Error opening file: NSX.xxxxx");
        fwrite(TemporaryField,8,(nx+2)*(ny+2)*(nz+2),NSX);
        fclose(NSX);

        for (int n = 0; n < N; n++) TemporaryField[n] = Averages->GradPhiY(n);
        sprintf(LocalRankFilename,"%s%s","ns_y.",LocalRankString);
        FILE *NSY = fopen(LocalRankFilename,"wb");
        if (NSY==NULL) ERROR("Error opening file: NSY.xxxxx");
        fwrite(TemporaryField,8,(nx+2)*(ny+2)*(nz+2),NSY);
        fclose(NSY);

        for (int n = 0; n < N; n++) TemporaryField[n] = Averages->GradPhiZ(n);
        sprintf(LocalRankFilename,"%s%s","ns_z.",LocalRankString);
        FILE *NSZ = fopen(LocalRankFilename,"wb");
        if (NSZ==NULL) ERROR("Error opening file: NSZ.xxxxx");
        fwrite(TemporaryField,8,(nx+2)*(ny+2)*(nz+2),NSZ);
        fclose(NSZ);

		delete[] TemporaryField;

		MPI_Barrier(comm);

	}

	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;

}
