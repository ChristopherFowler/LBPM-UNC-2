/*
 * Pre-processor to generate signed distance function from segmented data
 * segmented data should be stored in a raw binary file as 1-byte integer (type char)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common/Domain.h"
#include "common/SpherePack.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "analysis/distance.h"



int main(int argc, char **argv)
{
    int rank, nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);

    {

        int i,j,k,n;
        int nx, ny, nz;

        char LocalRankFilename[40];
        
        string filename;
        if (argc > 1){
            filename=argv[1];
        }
        else ERROR("No input database provided\n");

        // read the input database
        auto db = std::make_shared<Database>( filename );
        auto domain_db = db->getDatabase( "Domain" );
        
        // Read domain parameters
        auto size = domain_db->getVector<int>( "n" );
        auto READFILE = domain_db->getScalar<std::string>( "Filename" );

        std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
        std::shared_ptr<Domain> Mask (new Domain(domain_db,comm));

        nx = size[0]; ny = size[1]; nz = size[2];
        nx+=2; ny+=2; nz+=2;
        int N = nx*ny*nz;

        char *ReadID;
        ReadID = new char[N];
        
        Array<char> id_solid(nx,ny,nz);
        DoubleArray SignDist(nx,ny,nz);

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

        auto SD_new = new double[N];
        for (int k=0;k<nz;k++){
            for (int j=0;j<ny;j++){
                for (int i=0;i<nx;i++){
                    int n = k*nx*ny+j*nx+i;
                    SD_new[n] = SignDist(i,j,k);
                }
            }
        }

		sprintf(LocalRankFilename,"SignDist.%05i",rank);
		FILE *DIST = fopen(LocalRankFilename,"wb");
		fwrite(SD_new,8,N,DIST);
		fclose(DIST);

	}

	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;

}
