#include <analysis/morphology.h>
// Implementation of morphological opening routine


inline void PackID(int *list, int count, char *sendbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This packs up the values that need to be sent from one processor to another
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = ID[n];
	}
}
//***************************************************************************************

inline void UnpackID(int *list, int count, char *recvbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This unpacks the values once they have been recieved from neighbors
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		ID[n] = recvbuf[idx];
	}
}

double MorphOpen(DoubleArray &SignDist, char *id, char *id_original, std::shared_ptr<Domain> Dm, double VoidFraction, char ErodeLabel, char NewLabel, int amin, int amax,
	double deltaR, double Rcrit_new, int count_connected){
	// SignDist is the distance to the object that you want to constaing the morphological opening
	// VoidFraction is the the empty space where the object inst
	// id is a labeled map
	// Dm contains information about the domain structure

	int nx = Dm->Nx;
	int ny = Dm->Ny;
	int nz = Dm->Nz;
	int nprocx = Dm->nprocx();
	int nprocy = Dm->nprocy();
	int nprocz = Dm->nprocz();
	int rank = Dm->rank();

	int n;
	double final_void_fraction;
	double count,countGlobal,totalGlobal,sw_new,sum_sat;
	count = 0.f;
	double maxdist=-200.f;
	double maxdistGlobal;

	for (int k=amin; k<amax; k++){
		for (int j=1; j<ny-1; j++){
			for (int i=1; i<nx-1; i++){
				n = k*nx*ny+j*nx+i;
				if ( id_original[n] > 0){
					count += 1.0;
					//id[n]  = 2;
				}
			}
		}
	}

	for (int k=1; k<nz-1; k++){
		for (int j=1; j<ny-1; j++){
			for (int i=1; i<nx-1; i++){
				n = k*nx*ny+j*nx+i;
				// extract maximum distance for critical radius
				if ( SignDist(i,j,k) > maxdist) maxdist=SignDist(i,j,k);
			}
		}
	}
	MPI_Barrier(Dm->Comm);


	// total Global is the number of nodes in the pore-space
	MPI_Allreduce(&count,&totalGlobal,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	MPI_Allreduce(&maxdist,&maxdistGlobal,1,MPI_DOUBLE,MPI_MAX,Dm->Comm);
	double volume=double(nprocx*nprocy*nprocz)*double(nx-2)*double(ny-2)*double(nz-2);
	double volume_fraction=totalGlobal/volume;
	//if (rank==0) printf("Volume fraction for morphological opening: %f \n",volume_fraction);
	//if (rank==0) printf("Maximum pore size: %f \n",maxdistGlobal);
	final_void_fraction = volume_fraction; //initialize

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
	int imin,jmin,kmin,imax,jmax,kmax;
	int Nx = nx;
	int Ny = ny;
	int Nz = nz;

	double void_fraction_old=1.0;
	double void_fraction_new=2.0*VoidFraction; 
	double void_fraction_diff_old = 1.0;
	double void_fraction_diff_new = 1.0;
	if (ErodeLabel == 1){
		VoidFraction = 1.0 - VoidFraction;
	}

	// Increase the critical radius until the target saturation is met
	//double deltaR=0.05; // amount to change the radius in voxel units
	double Rcrit_old= maxdistGlobal;
	//double Rcrit_new = maxdistGlobal;


	if (rank==0) {
		printf("VoidFraction: %f \n",VoidFraction);
		printf("voidfractionnew: %f \n",void_fraction_new);
	}
	//while (void_fraction_new > VoidFraction)
	while (void_fraction_new > VoidFraction && Rcrit_new > 0.5)
	{
		void_fraction_diff_old = void_fraction_diff_new;
		void_fraction_old = void_fraction_new;
		Rcrit_old = Rcrit_new;
		Rcrit_new -= deltaR*Rcrit_old;
		int Window=round(Rcrit_new);
		if (Window == 0) Window = 1; // If Window = 0 at the begining, after the following process will have sw=1.0
		// and sw<Sw will be immediately broken
		double LocalNumber=0.f;
		for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
				for(int i=0; i<Nx; i++){
					n = k*nx*ny + j*nx+i;
					if (SignDist(i,j,k) > Rcrit_new){
						// loop over the window and update
						imin=max(0,i-Window);
						jmin=max(0,j-Window);
						kmin=max(0,k-Window);
						imax=min(Nx,i+Window);
						jmax=min(Ny,j+Window);
						kmax=min(Nz,k+Window);
						for (kk=kmin; kk<kmax; kk++){
							for (jj=jmin; jj<jmax; jj++){
								for (ii=imin; ii<imax; ii++){
									int nn = kk*nx*ny+jj*nx+ii;
									double dsq = double((ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k));
									if (id[nn] == ErodeLabel && dsq <= Rcrit_new*Rcrit_new){
										if (Dm->kproc() == 0 && k < 5){
										    //do nothing
										} else{
											LocalNumber+=1.0;
											id[nn]=NewLabel;
										}
									}
								}
							}
						}

					}
					// move on
				}
			}
		}
		// Pack and send the updated ID values
		PackID(Dm->sendList_x, Dm->sendCount_x ,sendID_x, id);
		PackID(Dm->sendList_X, Dm->sendCount_X ,sendID_X, id);
		PackID(Dm->sendList_y, Dm->sendCount_y ,sendID_y, id);
		PackID(Dm->sendList_Y, Dm->sendCount_Y ,sendID_Y, id);
		PackID(Dm->sendList_z, Dm->sendCount_z ,sendID_z, id);
		PackID(Dm->sendList_Z, Dm->sendCount_Z ,sendID_Z, id);
		PackID(Dm->sendList_xy, Dm->sendCount_xy ,sendID_xy, id);
		PackID(Dm->sendList_Xy, Dm->sendCount_Xy ,sendID_Xy, id);
		PackID(Dm->sendList_xY, Dm->sendCount_xY ,sendID_xY, id);
		PackID(Dm->sendList_XY, Dm->sendCount_XY ,sendID_XY, id);
		PackID(Dm->sendList_xz, Dm->sendCount_xz ,sendID_xz, id);
		PackID(Dm->sendList_Xz, Dm->sendCount_Xz ,sendID_Xz, id);
		PackID(Dm->sendList_xZ, Dm->sendCount_xZ ,sendID_xZ, id);
		PackID(Dm->sendList_XZ, Dm->sendCount_XZ ,sendID_XZ, id);
		PackID(Dm->sendList_yz, Dm->sendCount_yz ,sendID_yz, id);
		PackID(Dm->sendList_Yz, Dm->sendCount_Yz ,sendID_Yz, id);
		PackID(Dm->sendList_yZ, Dm->sendCount_yZ ,sendID_yZ, id);
		PackID(Dm->sendList_YZ, Dm->sendCount_YZ ,sendID_YZ, id);
		//......................................................................................
		MPI_Sendrecv(sendID_x,Dm->sendCount_x,MPI_CHAR,Dm->rank_x(),sendtag,
				recvID_X,Dm->recvCount_X,MPI_CHAR,Dm->rank_X(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_X,Dm->sendCount_X,MPI_CHAR,Dm->rank_X(),sendtag,
				recvID_x,Dm->recvCount_x,MPI_CHAR,Dm->rank_x(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_y,Dm->sendCount_y,MPI_CHAR,Dm->rank_y(),sendtag,
				recvID_Y,Dm->recvCount_Y,MPI_CHAR,Dm->rank_Y(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Y,Dm->sendCount_Y,MPI_CHAR,Dm->rank_Y(),sendtag,
				recvID_y,Dm->recvCount_y,MPI_CHAR,Dm->rank_y(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_z,Dm->sendCount_z,MPI_CHAR,Dm->rank_z(),sendtag,
				recvID_Z,Dm->recvCount_Z,MPI_CHAR,Dm->rank_Z(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Z,Dm->sendCount_Z,MPI_CHAR,Dm->rank_Z(),sendtag,
				recvID_z,Dm->recvCount_z,MPI_CHAR,Dm->rank_z(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xy,Dm->sendCount_xy,MPI_CHAR,Dm->rank_xy(),sendtag,
				recvID_XY,Dm->recvCount_XY,MPI_CHAR,Dm->rank_XY(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_XY,Dm->sendCount_XY,MPI_CHAR,Dm->rank_XY(),sendtag,
				recvID_xy,Dm->recvCount_xy,MPI_CHAR,Dm->rank_xy(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Xy,Dm->sendCount_Xy,MPI_CHAR,Dm->rank_Xy(),sendtag,
				recvID_xY,Dm->recvCount_xY,MPI_CHAR,Dm->rank_xY(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xY,Dm->sendCount_xY,MPI_CHAR,Dm->rank_xY(),sendtag,
				recvID_Xy,Dm->recvCount_Xy,MPI_CHAR,Dm->rank_Xy(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xz,Dm->sendCount_xz,MPI_CHAR,Dm->rank_xz(),sendtag,
				recvID_XZ,Dm->recvCount_XZ,MPI_CHAR,Dm->rank_XZ(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_XZ,Dm->sendCount_XZ,MPI_CHAR,Dm->rank_XZ(),sendtag,
				recvID_xz,Dm->recvCount_xz,MPI_CHAR,Dm->rank_xz(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Xz,Dm->sendCount_Xz,MPI_CHAR,Dm->rank_Xz(),sendtag,
				recvID_xZ,Dm->recvCount_xZ,MPI_CHAR,Dm->rank_xZ(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xZ,Dm->sendCount_xZ,MPI_CHAR,Dm->rank_xZ(),sendtag,
				recvID_Xz,Dm->recvCount_Xz,MPI_CHAR,Dm->rank_Xz(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_yz,Dm->sendCount_yz,MPI_CHAR,Dm->rank_yz(),sendtag,
				recvID_YZ,Dm->recvCount_YZ,MPI_CHAR,Dm->rank_YZ(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_YZ,Dm->sendCount_YZ,MPI_CHAR,Dm->rank_YZ(),sendtag,
				recvID_yz,Dm->recvCount_yz,MPI_CHAR,Dm->rank_yz(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Yz,Dm->sendCount_Yz,MPI_CHAR,Dm->rank_Yz(),sendtag,
				recvID_yZ,Dm->recvCount_yZ,MPI_CHAR,Dm->rank_yZ(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_yZ,Dm->sendCount_yZ,MPI_CHAR,Dm->rank_yZ(),sendtag,
				recvID_Yz,Dm->recvCount_Yz,MPI_CHAR,Dm->rank_Yz(),recvtag,Dm->Comm,MPI_STATUS_IGNORE);
		//......................................................................................
		UnpackID(Dm->recvList_x, Dm->recvCount_x ,recvID_x, id);
		UnpackID(Dm->recvList_X, Dm->recvCount_X ,recvID_X, id);
		UnpackID(Dm->recvList_y, Dm->recvCount_y ,recvID_y, id);
		UnpackID(Dm->recvList_Y, Dm->recvCount_Y ,recvID_Y, id);
		UnpackID(Dm->recvList_z, Dm->recvCount_z ,recvID_z, id);
		UnpackID(Dm->recvList_Z, Dm->recvCount_Z ,recvID_Z, id);
		UnpackID(Dm->recvList_xy, Dm->recvCount_xy ,recvID_xy, id);
		UnpackID(Dm->recvList_Xy, Dm->recvCount_Xy ,recvID_Xy, id);
		UnpackID(Dm->recvList_xY, Dm->recvCount_xY ,recvID_xY, id);
		UnpackID(Dm->recvList_XY, Dm->recvCount_XY ,recvID_XY, id);
		UnpackID(Dm->recvList_xz, Dm->recvCount_xz ,recvID_xz, id);
		UnpackID(Dm->recvList_Xz, Dm->recvCount_Xz ,recvID_Xz, id);
		UnpackID(Dm->recvList_xZ, Dm->recvCount_xZ ,recvID_xZ, id);
		UnpackID(Dm->recvList_XZ, Dm->recvCount_XZ ,recvID_XZ, id);
		UnpackID(Dm->recvList_yz, Dm->recvCount_yz ,recvID_yz, id);
		UnpackID(Dm->recvList_Yz, Dm->recvCount_Yz ,recvID_Yz, id);
		UnpackID(Dm->recvList_yZ, Dm->recvCount_yZ ,recvID_yZ, id);
		UnpackID(Dm->recvList_YZ, Dm->recvCount_YZ ,recvID_YZ, id);
		//......................................................................................

		//MPI_Allreduce(&LocalNumber,&GlobalNumber,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

		//ErodeLabel = 2, NewLabel = 1
		count = 0.f;
		for (int k=amin; k<amax; k++){
			for (int j=1; j<Ny-1; j++){
				for (int i=1; i<Nx-1; i++){
					n=k*Nx*Ny+j*Nx+i;
					if (id[n] == ErodeLabel){
						count+=1.0;
					}
				}
			}
		}

		MPI_Allreduce(&count,&countGlobal,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
		//void_fraction_new = countGlobal/totalGlobal;
		void_fraction_new = countGlobal;
		void_fraction_diff_new = abs(void_fraction_new-VoidFraction);

		count = 0.f;
		//porecount = 0.f;
		for (int k=amin; k<amax; k++){
            for (int j=1; j<Ny-1; j++){
                for (int i=1; i<Nx-1; i++){
                    n = k*Nx*Ny+j*Nx+i;
                    if ( id_original[n] == 1 || id[n] == NewLabel){
						count+=1.0;
					}
                }
            }
        }


        //MPI_Allreduce(&porecount,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&count,&sum_sat,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
        sw_new = 1.0-(sum_sat/totalGlobal);

		if (Dm->rank()==0){
			//printf("     %f ",sw_new);
			printf("voidfractionnew: %f voidfraction: %f \n",void_fraction_new,VoidFraction);
			//printf("     %f ",void_fraction_new);
			//printf("     %f\n",Rcrit_new);
		}
	}

	if (void_fraction_diff_new<void_fraction_diff_old){
		final_void_fraction=void_fraction_new;
		if (rank==0){
		// 	printf("Final void fraction =%f\n",void_fraction_new);
		 	printf("Final critical radius=%f\n",Rcrit_new);
		}
	}
	else{
		final_void_fraction=void_fraction_old;
		// if (rank==0){
		// 	printf("Final void fraction=%f\n",void_fraction_old);
		// 	printf("Final critical radius=%f\n",Rcrit_old);
		// }
	}
	return final_void_fraction;
}
/*
double morph_open()
{
	fillHalo<char> fillChar(Dm->Comm,Dm->rank_info,{Nx-2,Ny-2,Nz-2},{1,1,1},0,1);
	MPI_Allreduce(&LocalNumber,&GlobalNumber,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	count = 0.f;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				n=k*Nx*Ny+j*Nx+i;
				if (id[n] == 2){
					count+=1.0;
				}
			}
		}
	}
	MPI_Allreduce(&count,&countGlobal,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	return countGlobal;
}
*/
