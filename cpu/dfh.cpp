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
#include <math.h>
#include <stdio.h>

extern "C" void ScaLBL_Gradient_Unpack(double weight, double Cqx, double Cqy, double Cqz, 
		int *list, int start, int count, double *recvbuf, double *phi, double *grad, int N){

}

extern "C" void ScaLBL_DFH_Init(double *Phi, double *Den, double *Aq, double *Bq, int start, int finish, int Np){

}


// LBM based on density functional hydrodynamics
extern "C" void ScaLBL_D3Q19_AAeven_DFH(int *neighborList, double *dist, double *Aq, double *Bq, double *Den, double *Phi,
		double *Gradient, double *SolidForce, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int start, int finish, int Np){

	
}

extern "C" void ScaLBL_D3Q19_AAodd_DFH(int *neighborList, double *dist, double *Aq, double *Bq, double *Den, 
		double *Phi, double *Gradient, double *SolidForce, double rhoA, double rhoB, double tauA, double tauB, double alpha, double beta,
		double Fx, double Fy, double Fz, int start, int finish, int Np){
	
	
}

extern "C" void ScaLBL_D3Q7_AAodd_DFH(int *neighborList, double *Aq, double *Bq, 
			double *Den, double *Phi, int start, int finish, int Np){


}

extern "C" void ScaLBL_D3Q7_AAeven_DFH(double *Aq, double *Bq, double *Den, double *Phi, 
			int start, int finish, int Np){

}

extern "C" void ScaLBL_D3Q19_Gradient_DFH(int *neighborList, double *Phi, double *ColorGrad, int start, int finish, int Np){
	

}
