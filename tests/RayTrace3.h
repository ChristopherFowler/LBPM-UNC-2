#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <array>
#include <math.h>
#include "common/Array.h"

#include "models/ColorModel.h"
#include "analysis/TwoPhase.h"

using namespace std;

class RayTrace {
public:
    RayTrace() {
        ray_origin = new double[3];
        ray_direction = new double[3];
        geometry = 0;
        
        /* s: geometry = 2 */
        plane_s = new double[3];
        /* */
        plane_pp1 = new double[3];
        plane_pp2 = new double[3];
        
        
        plane_st1 = new double[3];
        plane_st2 = new double[3];
        plane_st3 = new double[3];
        plane_st4 = new double[3];
        
        
        
        /* sp: geometry = 1 */
        sphere_center = new double[3];
        sphere_radius = 0;
        o_minus_c = new double[3];
    }
    
    /* Setters and Getters */
    void setSphereRadius(double & rad) {
        sphere_radius = rad;
    }
    void setRayOrigin(double *  tmp) {
        ray_origin[0] = tmp[0];
        ray_origin[1] = tmp[1];
        ray_origin[2] = tmp[2];
    }
    void setSphereCenter(double *  tmp) {
        sphere_center[0] = tmp[0];
        sphere_center[1] = tmp[1];
        sphere_center[2] = tmp[2];
    }
    double * getRayOrigin() {
        return ray_origin;
    }
    void setRayDirection(double *  tmp) {
        ray_direction[0] = tmp[0];
        ray_direction[1] = tmp[1];
        ray_direction[2] = tmp[2];
    }
    double * getRayDirection_sp() { return ray_direction; }
    double getRayTraceDistance_sp() {
        return ray_trace_sp(ray_origin, sphere_center, ray_direction, sphere_radius);
    }
    
    void setPlane_s(double A, double B, double C, double D) {
        plane_s[0] = A;
        plane_s[1] = B;
        plane_s[2] = C;
        D_s = D;
    }
    
    
    void setPlane_pp1(double A, double B, double C, double D) {
        plane_pp1[0] = A;
        plane_pp1[1] = B;
        plane_pp1[2] = C;
        D_pp1 = D;
    }
    
    void setPlane_pp2(double A, double B, double C, double D) {
        plane_pp2[0] = A;
        plane_pp2[1] = B;
        plane_pp2[2] = C;
        D_pp2 = D;
    }
    
    
    
    void setPlane_st1(double A, double B, double C, double D) {
        plane_st1[0] = A;
        plane_st1[1] = B;
        plane_st1[2] = C;
        D_st1 = D;
    }
    void setPlane_st2(double A, double B, double C, double D) {
        plane_st2[0] = A;
        plane_st2[1] = B;
        plane_st2[2] = C;
        D_st2 = D;
    }
    void setPlane_st3(double A, double B, double C, double D) {
        plane_st3[0] = A;
        plane_st3[1] = B;
        plane_st3[2] = C;
        D_st3 = D;
    }
    void setPlane_st4(double A, double B, double C, double D) {
        plane_st4[0] = A;
        plane_st4[1] = B;
        plane_st4[2] = C;
        D_st4 = D;
    }
    
    
    
    double * getRayDirection_s() { return ray_direction; }
    double getRayTraceDistance_s() {
        return ray_trace_s(ray_origin, plane_s, D_s);
    }
    
    double getRayTraceDistance_s_q() {
        return ray_trace_s_q(ray_origin, plane_s, D_s);
    }
    
    double * getRayDirection_pp() { return ray_direction; }
    double getRayTraceDistance_pp() {
        return ray_trace_pp(ray_origin, plane_pp1, D_pp1, plane_pp2, D_pp2);
    }
    
    double getRayTraceDistance_pp_q() {
           return ray_trace_pp_q(ray_origin, ray_direction, plane_pp1, D_pp1, plane_pp2, D_pp2);
       }
    
    
    double * getRayDirection_st() { return ray_direction; }
    double getRayTraceDistance_st() {
        return ray_trace_st(ray_origin,
                            plane_st1, D_st1,
                            plane_st2, D_st2,
                            plane_st3, D_st3,
                            plane_st4, D_st4);
    }
    
    double getRayTraceDistance_st_q() {
        return ray_trace_st_q(ray_origin,
                              plane_st1, D_st1,
                              plane_st2, D_st2,
                              plane_st3, D_st3,
                              plane_st4, D_st4);
    }
    
    double * getRayDirection_ct() { return ray_direction; }
    double getRayTraceDistance_ct() {
        return ray_trace_ct(ray_origin, sphere_center, ray_direction, sphere_radius);
    }
    
    
    
    
private:
    // 1 periodic sphere pack, 2 slab, 3 parallel plates, 4 square tube, 5 capillary tube
    int geometry;
    
    /* Operators needed to compute the ray_trace function */
    double dot(double *  vec1, double *  vec2);
    double mag(double *  vec1 /* Get || . ||  of this vector */);
    double magSq(double *  vec1 /* Get || . ||^2 of this vector */);
    void subtract_vectors(double *  vec1 /* ray_origin */, double *  vec2 /* cell_center */,double *  vec3 /* Return vector */);
    
    
    
    /* Signatures */
    double dist_sp(double & term1 /* u_dot_o_minus_c */, double & term2 /* under_sqrt */);
    
    
    
    double dist_ct(double & term1 /* "u.(o-c)" */,
                   double & term2 /* discriminant: "( u.(o-c) )^2 - ( ||o-c||^2 - r^2 )" */);
    
    
    double ray_trace_sp(double *  vec1 /* ray_origin */, double *  vec2 /* sphere center */, double *  vec3 /*unit_vector */, double rad /* sphere_radius */);
    
    double ray_trace_s(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* plane translation */);
    double ray_trace_s_q(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* plane translation */);
    
    double ray_trace_pp(double *  vec1 /* ray_origin */, double *  vec2 /* plane_pp1 */, double D_pp1, double *  vec3 /* plane_pp2 */, double D_pp2);
    double ray_trace_pp_q(double *  lnaught /* ray_origin */, double * l /* ray direction */,
                          double *  plane_pp1 /* plane_pp1 */, double D_pp1 /* distance from origin- need to generalize this for phi,theta!=0 */,
                          double *  plane_pp2 /* plane_pp2 */, double D_pp2);
    
    double ray_trace_st(double *  vec1 /* ray_origin */,
                        double *  vec2 /* plane_pp1 */, double D_st1,
                        double *  vec3 /* plane_pp2 */, double D_st2,
                        double *  vec4 /* plane_st3 */, double D_st3,
                        double *  vec5 /* plane_st4 */, double D_st4);
    
    double ray_trace_st_q(double *  vec1 /* ray_origin */,
                          double *  vec2 /* plane_pp1 */, double D_st1,
                          double *  vec3 /* plane_pp2 */, double D_st2,
                          double *  vec4 /* plane_st3 */, double D_st3,
                          double *  vec5 /* plane_st4 */, double D_st4);
    
    double ray_trace_ct(double *  vec1 /* ray_origin */, double *  vec2 /* sphere center */, double *  vec3 /*unit_vector */, double rad /* sphere_radius */);
    
private:
    
    /* general ray information */
    double * ray_origin;
    double * ray_direction;
    
    /* sphere pack */
    double * sphere_center;
    double sphere_radius;
    
    /* slab */
    double * plane_s;
    double D_s;
    
    /* parallel plates */
    double * plane_pp1, * plane_pp2;
    double D_pp1,D_pp2;
    
    /* square tube */
    double * plane_st1, * plane_st2, * plane_st3, * plane_st4;
    double D_st1,D_st2, D_st3, D_st4;
    
    /* capillary tube */
    // NYI
    
    double * o_minus_c;
    double u_dot_o_minus_c;
    double dp,dp1,dp2,dp3,dp4;
    double under_sqrt;
    
    double A,B,C,D;
};

double RayTrace::dot(double *  vec1, double *  vec2) {
    return (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]);
}

double RayTrace::mag(double *  vec1 /* Get || . ||  of this vector */) {
    return sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
}

double RayTrace::magSq(double *  vec1 /* Get || . ||^2 of this vector */) {
    return (vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
}

void RayTrace::subtract_vectors(double *  vec1 /* ray_origin */, double *  vec2 /* cell_center */,
                                double *  vec3 /* Return vector */) {
    vec3[0] = vec1[0] - vec2[0];
    vec3[1] = vec1[1] - vec2[1];
    vec3[2] = vec1[2] - vec2[2];
}


double RayTrace::dist_sp(double & term1 /* "u.(o-c)" */,
                         double & term2 /* discriminant: "( u.(o-c) )^2 - ( ||o-c||^2 - r^2 )" */) {
    // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    
    double d1,d2,pos;
    d1 = 0;
    d2 = 0;
    pos = 0;
    /* No solution exists */
    if (term2 < 0) pos = 0;
    
    /* One solution exists */
    if (term2 == 0) {
        d1 = -term1;
        pos = d1;
    }
    
    /* Two solutions exist */
    if (term2 > 0) {
        d1 = -term1 + sqrt(term2);
        d2 = -term1 - sqrt(term2);
        if (d2 >= d1) std::swap(d2,d1);
        if (d2 < 0) d2 = d1;
        pos = d2;
    }
    
    return pos;
}



double RayTrace::dist_ct(double & term1 /* "u.(o-c)" */,
                         double & term2 /* discriminant: "( u.(o-c) )^2 - ( ||o-c||^2 - r^2 )" */) { return 0.0; }





double RayTrace::ray_trace_sp(double *  vec1 /* ray_origin */,
                              double *  vec2 /*sphere_center */,
                              double *  vec3 /* ray direction -> will become unit vector */,
                              double rad /* analytical sphere radius */) {
    // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    // Note: using "u" for unit vector instead of "l"
    
    /* Ray tracing work */
    subtract_vectors(vec1 /* ray origin "o" */,
                     vec2 /* sphere center "c" */,
                     o_minus_c /* Return: "o-c"*/);
    
    /* Normalize "unit vector" */
    double tmp = mag(vec3);
    if (tmp == 0.0) tmp = 1.0; // Needed for domain (LIBB_id) generation
    vec3[0] /= tmp; vec3[1] /= tmp; vec3[2] /= tmp;
    dp = dot(vec3 /* unit vector */ ,o_minus_c); // Return: "u.(o-c)"
    under_sqrt = dp*dp - (magSq(o_minus_c) - rad*rad); // Return the discriminant: "( u.(o-c) )^2 - ( ||o-c||^2 - r^2 )"
    /* Get the solution */
    return dist_sp(dp, under_sqrt);
}

double RayTrace::ray_trace_s(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* displacement */) {
    
    double tmp = mag(vec2); // Denominator of distance
    if (tmp == 0) tmp = 1.0;
    
    
    dp = dot(vec2, vec1);
    // cout << "dp=" << dp << endl;
    
    return (-(dp - D_s)/tmp);
}


double RayTrace::ray_trace_s_q(double *  vec1 /* ray_origin */, double *  vec2 /* plane_s */, double D_s /* displacement */)
{
    double tmp = mag(vec2); if (tmp == 0) tmp = 1.0;
    subtract_vectors(vec2 /* ray origin "o" */, vec1 /* distance from world origin point "c" */, o_minus_c /* Return: "l0-p0"*/);
    o_minus_c[0] += 0; o_minus_c[1] += 0; o_minus_c[2] += 0;
    
    double * normal; normal = new double[3];
    normal[0] = vec2[0]/=tmp;
    normal[1] = vec2[1]/=tmp;
    normal[2] = vec2[2]/=tmp;
    
    tmp = dot(ray_direction,normal); if (tmp == 0) tmp = 1.0;
    dp = dot(o_minus_c, normal) / tmp;
    
    delete[] normal;
    
    // if (dp < 1.5) cout << "dp=" << dp << endl;
    
    return (-(dp-(D_s)));
}





//double RayTrace::ray_trace_st(double *  vec1 /* ray_origin */, double *  vec2 /* sphere center */, double *  vec3 /*unit_vector */, double rad /* sphere_radius */) { return 0.0; }
//
//double RayTrace::ray_trace_ct(double *  vec1 /* ray_origin */, double *  vec2 /* sphere center */, double *  vec3 /*unit_vector */, double rad /* sphere_radius */) { return 0.0; }


static inline void fgetl( char * str,
                         int num, FILE * stream ) {
    char* ptr = fgets( str, num, stream );
    if ( 0 ) {char *temp = (char *)&ptr; temp++;}
}

void ReadSpherePacking(int nspheres,
                       double *List_cx, double *List_cy, double *List_cz, double *List_rad) {
    FILE *fid = fopen("pack.out","rb");
    char line[100];
    int count = 0;
    while ( !feof(fid) && fgets(line,100,fid)!=NULL ) {
        char* line2 = line;
        List_cx[count] = strtod(line2,&line2);
        List_cy[count] = strtod(line2,&line2);
        List_cz[count] = strtod(line2,&line2);
        List_rad[count] = strtod(line2,&line2);
        count++;
    }
    INSIST( count==nspheres, "Specified number of spheres is probably incorrect!" );
    // .............................................................
}

void SetIdAndCount(char * id,
                   size_t & count, double & hx, int & Nx, int & Ny, int & Nz, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz,int rank) {} // Not used

void ComputeLIBBqDistances_sp(char * LIBB_id,
                              double * qDistances, double * signedDistance, double * offsetDistance, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, Domain &Dm, int BoundaryCondition,
                              int rmin, int rmax) {
    
    
    
    size_t m;
   // size_t count = 0;
   // size_t c = 0;
    
  //  double volume = double((Nx-2)*(Ny-2)*(Nz-2)*nprocx*nprocy*nprocz);
    
    RayTrace rt;
    double * sphere_center;  sphere_center = new double[3];
    double * ray_origin;     ray_origin = new double[3];
    double * ray_direction;  ray_direction = new double[3];
    double sphere_radius;
    
    
    
    
  //  size_t cross_count = 0;
    double *ccx,*ccy,*ccz,*r;
    ccx = new double[nspheres];
    ccy = new double[nspheres];
    ccz = new double[nspheres];
    r = new double[nspheres];
    ReadSpherePacking(nspheres,
                      ccx,ccy,ccz,r);
    double cx,cy,cz,rad;
    
    /* Parallel registers */
    int imin,imax,jmin,jmax,kmin,kmax;
    double min_x,min_y,min_z;
    
   // double x,y,z;
    
    double dat = 0;
 //   double a = 0;
    
    //cout << "Lx=" << Lx << endl;
    double hx = Lx/double((Nx-2)*nprocx);
    
    double sr = r[0];
    if (rank == 0) cout << "All sphere radii are " << sr << " lattice units" << endl;
    
    //double analytical_porosity = (1.0-double(2)* (4./3.*3.14159265359 * sr*sr*sr))/1.0;
    //if (rank == 0) printf("Analytical porosity=%.15f (valid for equally sized spheres)\n",analytical_porosity);
    //............................................
    // Get maximum and minimum for this domain
    // Halo is included !
    min_x = double(iproc*(Nx-2))*hx;
    min_y = double(jproc*(Ny-2))*hx;
    min_z = double(kproc*(Nz-2))*hx;
    
    /* GENERATE CORRECT ID */
    for (int p=0;p<nspheres;p++){
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];

        /* converting sphere center to LBM units */
        sphere_center[0] = cx/hx+0.5;
        sphere_center[1] = cy/hx+0.5;
        sphere_center[2] = cz/hx+0.5;
        sphere_radius = rad/hx;

        /* Check for range of a particular sphere */
        imin = int ((cx-rad)/hx)-12;
        imax = int ((cx+rad)/hx)+12;
        jmin = int ((cy-rad)/hx)-12;
        jmax = int ((cy+rad)/hx)+12;
        kmin = int ((cz-rad)/hx)-12;
        kmax = int ((cz+rad)/hx)+12;

        if (imin<0)     imin = 0;
        if (imin>Nx)    imin = Nx;
        if (imax<0)     imax = 0;
        if (imax>Nx)    imax = Nx;
        if (jmin<0)     jmin = 0;
        if (jmin>Ny)    jmin = Ny;
        if (jmax<0)     jmax = 0;
        if (jmax>Ny)    jmax = Ny;
        if (kmin<0)     kmin = 0;
        if (kmin>Nz)    kmin = Nz;
        if (kmax<0)     kmax = 0;
        if (kmax>Nz)    kmax = Nz;

        rt.setSphereCenter(sphere_center);
        rt.setSphereRadius(sphere_radius);

        // Loop over the domain for this sphere (may be null)
        for (int i=imin;i<imax;i++){
            for (int j=jmin;j<jmax;j++){
                for (int k=kmin;k<kmax;k++){
                    m = k*Nx*Ny+j*Nx+i;
                    /* ray origin: LBM units */
                    ray_origin[0] = double(i);
                    ray_origin[1] = double(j);
                    ray_origin[2] = double(k);

                    rt.setRayOrigin(ray_origin);
                    ray_direction[0] = 0.0;
                    ray_direction[1] = 0.0;
                    ray_direction[2] = 0.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat >= 0.5) {
                       
                        LIBB_id[m] = 0;
                        
                    }

                }
            }
        }
    }

    // if (rank==0){
	//     printf("LIBB_id \n");
    // }
    // for (int k=0;k<Nz;k++){
	//     for (int j=0;j<Ny;j++){
	// 	    for (int i=0;i<Nx;i++){
	// 		    int n=k*Nx*Ny+j*Nx+i;
	// 		    printf("%i ",LIBB_id[n]);
	// 	    }
	// 	    printf("\n");
	//     }
	//     printf("\n\n");
    // }




    if (BoundaryCondition > 0){
        if (Dm.kproc() == 0){
            int n;
            int N = Nx*Ny*Nz;
            for (int Slice = 0; Slice < rmin; Slice++) {
                for (n=Slice*Nx*Ny; n<(Slice+1)*Nx*Ny; n++){
                    LIBB_id[n] = 2;
                }
            }
        }
        
        if (Dm.kproc() == nprocz-1){
            int n;
            int N = Nx*Ny*Nz;
            for (int Slice = rmax; Slice < Nz; Slice++) {
                for (n=Slice*Nx*Ny; n<(Slice+1)*Nx*Ny; n++){
                    LIBB_id[n] = 2;
                }
            }
        }
    }
    
    
    for (int p=0;p<nspheres;p++){
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];

        // printf("cx cy cz r: %.2f %.2f %.2f %.2f\n",cx,cy,cz,rad);
        // printf("ccx ccy ccz rad: %.2f %.2f %.2f %.2f\n",ccx[p],ccy[p],ccz[p],r[p]);
        // printf("min x y z: %.2f %.2f %.2f \n",min_x,min_y,min_z);

        /* converting sphere center to LBM units */
        sphere_center[0] = cx/hx+0.5+0.5; // Extra 0.5 needed for tri mesh symmetry with the porous medium
        sphere_center[1] = cy/hx+0.5+0.5;
        sphere_center[2] = cz/hx+0.5+0.5;
        sphere_radius = rad/hx;

        // printf("sphere x y z r: %.2f %.2f %.2f %.2f\n",sphere_center[0],sphere_center[1],sphere_center[2],sphere_radius);

        for (int i=1;i<Nx-1;i++){
            for (int j=1;j<Ny-1;j++){
                for (int k=1;k<Nz-1;k++){
                    m = k*Nx*Ny+j*Nx+i;
                    double dist = sqrt( (double(i)-sphere_center[0])*(double(i)-sphere_center[0]) + (double(j)-sphere_center[1])*(double(j)-sphere_center[1]) + (double(k)-sphere_center[2])*(double(k)-sphere_center[2]) ) - sphere_radius ;
                    // if (i==4 && j==4 && k==4){
                    //     printf("dist at 4,4,4: %.2f \n",dist);
                    // }
                    if (dist < signedDistance[m]) signedDistance[m] = dist;
                }
            }
        }
    } // n spheres loop
    
    
    for (int p=0;p<nspheres;p++){
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];
        
        /* converting sphere center to LBM units */
        sphere_center[0] = cx/hx+0.5; // Extra 0.5 needed for tri mesh symmetry with the porous medium
        sphere_center[1] = cy/hx+0.5;
        sphere_center[2] = cz/hx+0.5;
        sphere_radius = rad/hx;
        for (int i=1;i<Nx-1;i++){
            for (int j=1;j<Ny-1;j++){
                for (int k=1;k<Nz-1;k++){
                    m = k*Nx*Ny+j*Nx+i;
                    double dist = sqrt( (double(i)-sphere_center[0])*(double(i)-sphere_center[0]) + (double(j)-sphere_center[1])*(double(j)-sphere_center[1]) + (double(k)-sphere_center[2])*(double(k)-sphere_center[2]) ) - sphere_radius ;
                    if (dist < offsetDistance[m]) offsetDistance[m] = dist;
                }
            }
        }
    } // n spheres loop
    

    // if (rank==0){
    //         printf("signedDistance \n");
    // }
    // for (int k=0;k<Nz;k++){
    //         for (int j=0;j<Ny;j++){
    //                 for (int i=0;i<Nx;i++){
    //                         int n=k*Nx*Ny+j*Nx+i;
    //                         printf("%.2f ",signedDistance[n]);
    //                 }
    //                 printf("\n");
    //         }
    //         printf("\n\n");
    // }

    // if (rank==0){
    //         printf("offsetDistance \n");
    // }
    // for (int k=0;k<Nz;k++){
    //         for (int j=0;j<Ny;j++){
    //                 for (int i=0;i<Nx;i++){
    //                         int n=k*Nx*Ny+j*Nx+i;
    //                         printf("%.2f ",offsetDistance[n]);
    //                 }
    //                 printf("\n");
    //         }
    //         printf("\n\n");
    // }


    
    if (BoundaryCondition > 0){
        if (Dm.kproc() == 0){
            int n;
            int N = Nx*Ny*Nz;
            for (int Slice = 0; Slice < rmin; Slice++) {
                for (n=Slice*Nx*Ny; n<(Slice+1)*Nx*Ny; n++){
                    offsetDistance[n] = rmin-1;
                    signedDistance[n] = rmin-1;
                }
            }
        }
        
        if (Dm.kproc() == nprocz-1){
            int n;
            int N = Nx*Ny*Nz;
            for (int Slice = rmax; Slice < Nz; Slice++) {
                for (n=Slice*Nx*Ny; n<(Slice+1)*Nx*Ny; n++){
                    offsetDistance[n] = Nz-rmax-1;
                    signedDistance[n] = Nz-rmax-1;
                }
            }
        }
    }
    
//    for (int i=1;i<Nx-1;i++)
//    for (int j=1;j<Ny-1;j++)
//    for (int k=1;k<Nz-1;k++){
//        m = k*Nx*Ny+j*Nx+i;
//        if (offsetDistance[m] < 0) LIBB_id[m] = 0;
//    }
    
    
    // .........Loop over the spheres.............
    for (int p=0;p<nspheres;p++){
        // Get the sphere from the list, map to local min
        cx = ccx[p] - min_x;
        cy = ccy[p] - min_y;
        cz = ccz[p] - min_z;
        rad = r[p];
        
        /* converting sphere center to LBM units */
        sphere_center[0] = cx/hx+0.5;
        sphere_center[1] = cy/hx+0.5;
        sphere_center[2] = cz/hx+0.5;
        sphere_radius = rad/hx;
        
        /* Check for range of a particular sphere */
        imin = int ((cx-rad)/hx)-2;
        imax = int ((cx+rad)/hx)+3;
        jmin = int ((cy-rad)/hx)-2;
        jmax = int ((cy+rad)/hx)+3;
        kmin = int ((cz-rad)/hx)-2;
        kmax = int ((cz+rad)/hx)+3;
        if (imin<0)     imin = 0;
        if (imin>Nx)    imin = Nx;
        if (imax<0)     imax = 0;
        if (imax>Nx)    imax = Nx;
        if (jmin<0)     jmin = 0;
        if (jmin>Ny)    jmin = Ny;
        if (jmax<0)     jmax = 0;
        if (jmax>Ny)    jmax = Ny;
        if (kmin<0)     kmin = 0;
        if (kmin>Nz)    kmin = Nz;
        if (kmax<0)     kmax = 0;
        if (kmax>Nz)    kmax = Nz;
        
        rt.setSphereCenter(sphere_center);
        rt.setSphereRadius(sphere_radius);
        
        
        // Loop over the domain for this sphere (may be null)
        for (int i=imin;i<imax;i++){
            for (int j=jmin;j<jmax;j++){
                for (int k=kmin;k<kmax;k++){
                    m = k*Nx*Ny+j*Nx+i;
                    
                    /* ray origin: LBM units */
                    ray_origin[0] = double(i);
                    ray_origin[1] = double(j);
                    ray_origin[2] = double(k);
                    rt.setRayOrigin(ray_origin);
                    
                    /* Generate ray trace distance map */
                    // 1
                    ray_direction[0] = 1.0;
                    ray_direction[1] = 0.0;
                    ray_direction[2] = 0.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= 1 && dat > 0) qDistances[m] = dat;
                    
                    // 2
                    ray_direction[0] = -1.0;
                    ray_direction[1] = 0.0;
                    ray_direction[2] = 0.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= 1 && dat > 0) qDistances[m+Nx*Ny*Nz] = dat;
                    
                    // 3
                    ray_direction[0] = 0.0;
                    ray_direction[1] = 1.0;
                    ray_direction[2] = 0.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= 1 && dat > 0) qDistances[m+2*Nx*Ny*Nz] = dat;
                    
                    // 4
                    ray_direction[0] = 0.0;
                    ray_direction[1] = -1.0;
                    ray_direction[2] = 0.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= 1 && dat > 0) qDistances[m+3*Nx*Ny*Nz] = dat;
                    
                    // 5
                    ray_direction[0] = 0.0;
                    ray_direction[1] = 0.0;
                    ray_direction[2] = 1.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= 1 && dat > 0) qDistances[m+4*Nx*Ny*Nz] = dat;
                    
                    // 6
                    ray_direction[0] = 0.0;
                    ray_direction[1] = 0.0;
                    ray_direction[2] = -1.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= 1 && dat > 0) qDistances[m+5*Nx*Ny*Nz] = dat;
                    
                    // 7
                    ray_direction[0] = 1.0;
                    ray_direction[1] = 1.0;
                    ray_direction[2] = 0.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+6*Nx*Ny*Nz] = dat;
                    
                    // 8
                    ray_direction[0] = -1.0;
                    ray_direction[1] = -1.0;
                    ray_direction[2] = 0.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+7*Nx*Ny*Nz] = dat;
                    
                    // 9
                    ray_direction[0] = 1.0;
                    ray_direction[1] = -1.0;
                    ray_direction[2] = 0.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+8*Nx*Ny*Nz] = dat;
                    
                    // 10
                    ray_direction[0] = -1.0;
                    ray_direction[1] = 1.0;
                    ray_direction[2] = 0.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+9*Nx*Ny*Nz] = dat;
                    
                    // 11
                    ray_direction[0] = 1.0;
                    ray_direction[1] = 0.0;
                    ray_direction[2] = 1.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+10*Nx*Ny*Nz] = dat;
                    
                    // 12
                    ray_direction[0] = -1.0;
                    ray_direction[1] = 0.0;
                    ray_direction[2] = -1.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+11*Nx*Ny*Nz] = dat;
                    
                    // 13
                    ray_direction[0] = 1.0;
                    ray_direction[1] = 0.0;
                    ray_direction[2] = -1.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+12*Nx*Ny*Nz] = dat;
                    
                    // 14
                    ray_direction[0] = -1;
                    ray_direction[1] = 0.0;
                    ray_direction[2] = 1;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+13*Nx*Ny*Nz] = dat;
                    
                    // 15
                    ray_direction[0] = 0.0;
                    ray_direction[1] = 1;
                    ray_direction[2] = 1;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+14*Nx*Ny*Nz] = dat;
                    
                    // 16
                    ray_direction[0] = 0.0;
                    ray_direction[1] = -1;
                    ray_direction[2] = -1;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    dat = abs(dat);
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+15*Nx*Ny*Nz] = dat;
                    
                    // 17
                    ray_direction[0] = 0.0;
                    ray_direction[1] = 1.0;
                    ray_direction[2] = -1.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+16*Nx*Ny*Nz] = dat;
                    
                    // 18
                    ray_direction[0] = 0.0;
                    ray_direction[1] = -1.0;
                    ray_direction[2] = 1.0;
                    rt.setRayDirection(ray_direction);
                    dat = rt.getRayTraceDistance_sp();
                    if (dat <= sqrt(2) && dat > 0) qDistances[m+17*Nx*Ny*Nz] = dat;
                    // cout << "dat=" << dat << endl;
                }
            }
        }
    }
    
    
    // if (rank==0){
    //         printf("qDistances \n");
    // }
    // for (int k=0;k<Nz;k++){
    //         for (int j=0;j<Ny;j++){
    //                 for (int i=0;i<Nx;i++){
    //                         int n=k*Nx*Ny+j*Nx+i;
    //                         printf("%.2f ",qDistances[n]);
    //                 }
    //                 printf("\n");
    //         }
    //         printf("\n\n");
    // }


    
}

void ComputeLIBBqDistances_s(char * LIBB_id, double * qDistances, double * signedDistance, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, 
    int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, Domain &Dm, int plate_min, int plate_max) {

    size_t N = Nx*Ny*Nz;
    double dist = 0;
    auto signedDistance_Invert = new double[N];  for (size_t n = 0; n < N; n++) signedDistance_Invert[n] = 100.0;
    
    for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            for (int k=0;k<Nz;k++){
                int n = k*Nx*Ny+j*Nx+i;
                if (k>=plate_min){
                    if (k<=plate_max){
                        LIBB_id[n] = 0;
                        qDistances[n] = 1.0;
                        signedDistance[n] = 0;
                    }
                    if (k>plate_max) {
                        signedDistance_Invert[n] = 0;
                    }
                } else {
                    signedDistance_Invert[n] = 0;
                }
            }
        }
    }

    for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            for (int k=0;k<Nz;k++){
                int n = k*Nx*Ny+j*Nx+i;
                if (k<plate_max+1){

                    for (int ii=0;ii<Nx;ii++){
                        for (int jj=0;jj<Ny;jj++){
                            for (int kk=0;kk<plate_min;kk++){
                                //Calculate SD
                                int m = kk*Nx*Ny+jj*Nx+ii;
                                dist = sqrt( (ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k) );
                                if (dist < signedDistance[m]){
                                    signedDistance[m] = dist;
                                }
                            }
                            for (int kk=(plate_max+1);kk<Nz;kk++){
                                //Calculate SD
                                int m = kk*Nx*Ny+jj*Nx+ii;
                                dist = sqrt( (ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k) );
                                if (dist < signedDistance[m]){
                                    signedDistance[m] = dist;
                                }
                            }
                        }
                    }

                } 
                if(k>plate_max){

                    for (int ii=0;ii<Nx;ii++){
                        for (int jj=0;jj<Ny;jj++){
                            for (int kk=plate_min;kk<plate_max+1;kk++){
                                //Calculate SD
                                int m = kk*Nx*Ny+jj*Nx+ii;
                                dist = sqrt( (ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k) );
                                if (dist < signedDistance_Invert[m]){
                                    signedDistance_Invert[m] = dist;
                                }
                            }
                        }
                    }                    


                }
            }
        }
    }



    // if (rank==0){
	//     printf("LIBB_id \n");
    // }
    // for (int k=0;k<Nz;k++){
	//     for (int j=0;j<Ny;j++){
	// 	    for (int i=0;i<Nx;i++){
	// 		    int n=k*Nx*Ny+j*Nx+i;
	// 		    printf("%i ",LIBB_id[n]);
	// 	    }
	// 	    printf("\n");
	//     }
	//     printf("\n\n");
    // }

    // if (rank==0){
    //         printf("signedDistance \n");
    // }
    // for (int k=0;k<Nz;k++){
    //         for (int j=0;j<Ny;j++){
    //                 for (int i=0;i<Nx;i++){
    //                         int n=k*Nx*Ny+j*Nx+i;
    //                         printf("%.2f ",signedDistance[n]);
    //                 }
    //                 printf("\n");
    //         }
    //         printf("\n\n");
    // }

    // if (rank==0){
    //         printf("signedDistance Invert\n");
    // }
    // for (int k=0;k<Nz;k++){
    //         for (int j=0;j<Ny;j++){
    //                 for (int i=0;i<Nx;i++){
    //                         int n=k*Nx*Ny+j*Nx+i;
    //                         printf("%.2f ",signedDistance_Invert[n]);
    //                 }
    //                 printf("\n");
    //         }
    //         printf("\n\n");
    // }


    for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            for (int k=0;k<Nz;k++){
                int m = k*Nx*Ny+j*Nx+i;
                signedDistance[m] = signedDistance[m] - signedDistance_Invert[m];
            }
        }
    }    

    
}

void ComputeLIBBqDistances_pp(char * LIBB_id, double * qDistances, double * signedDistance, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, 
    int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, Domain &Dm, int plate_min, int plate_max) {

    size_t N = Nx*Ny*Nz;
    double dist = 0;
    auto signedDistance_Invert = new double[N];  for (size_t n = 0; n < N; n++) signedDistance_Invert[n] = 100.0;

    int top_plate_max = Nz-plate_min-1;
    int top_plate_min = Nz-plate_max-1;

    printf("Top Plate: %i to %i voxels \n",top_plate_min,top_plate_max);
    printf("Bot Plate: %i to %i voxels \n",plate_min,plate_max);

    for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            for (int k=0;k<Nz;k++){
                int n = k*Nx*Ny+j*Nx+i;
                if (k>=plate_min){
                    if (k<=plate_max){
                        LIBB_id[n] = 0;
                        qDistances[n] = 1.0;
                        signedDistance[n] = 0;
                    } else if (k < top_plate_min) {
                        signedDistance_Invert[n] = 0;
                    } else if (k >= top_plate_min && k <= top_plate_max ){
                        LIBB_id[n] = 0;
                        qDistances[n] = 1.0;
                        signedDistance[n] = 0;
                    } else {
                        signedDistance_Invert[n] = 0;
                    }
                } else {
                    signedDistance_Invert[n] = 0;
                }
            }
        }
    }

    // if (rank==0){
	//     printf("LIBB_id \n");
    // }
    // for (int k=0;k<Nz;k++){
	//     for (int j=0;j<Ny;j++){
	// 	    for (int i=0;i<Nx;i++){
	// 		    int n=k*Nx*Ny+j*Nx+i;
	// 		    printf("%i ",LIBB_id[n]);
	// 	    }
	// 	    printf("\n");
	//     }
	//     printf("\n\n");
    // }

    if (rank==0){
            printf("signedDistance Invert\n");
    }
    for (int k=0;k<Nz;k++){
            for (int j=0;j<Ny;j++){
                    for (int i=0;i<Nx;i++){
                            int n=k*Nx*Ny+j*Nx+i;
                            printf("%.2f ",signedDistance_Invert[n]);
                    }
                    printf("\n");
            }
            printf("\n\n");
    }

    bool inside_plate = false;
    double count = 0.0;

    for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            for (int kk=plate_min;kk<plate_max+1;kk++){
                int m = kk*Nx*Ny+j*Nx+i;
                signedDistance_Invert[m] = double(kk+1);
                // if (i==0 && j ==0 && kk == 0){
                //     printf("%.2f \n",signedDistance_Invert[m]);
                // }
            }
            for (int kk=plate_max;kk>(plate_min-1);kk--){
                int m = kk*Nx*Ny+j*Nx+i;
                dist = double(plate_max-kk+1);
                if (dist < signedDistance_Invert[m]){
                    signedDistance_Invert[m] = dist;
                } 
                // if (i==0 && j ==0 && kk == 0){
                //     printf("%.2f \n",signedDistance_Invert[m]);
                // }
            }
            for (int kk=(plate_max+1);kk<top_plate_min;kk++){
                //Calculate SD
                int m = kk*Nx*Ny+j*Nx+i;
                signedDistance[m] = double(kk-plate_max);
            }
            for (int kk=(top_plate_min-1);kk>plate_max;kk--){
                //Calculate SD
                int m = kk*Nx*Ny+j*Nx+i;
                dist = double(top_plate_min-kk);
                if (dist < signedDistance[m]){
                    signedDistance[m] = dist;
                }
            }

            // for (int k=0;k<Nz;k++){
            //     int n = k*Nx*Ny+j*Nx+i;

            //     //if inside plate, calculate outside plate TRUE
            //     if (k >= plate_min && k <= plate_max){
            //         inside_plate = true;
            //     } 
            //     // else if (k >= top_plate_min && k <= top_plate_max){
            //     //     inside_plate = true;
            //     // } 


            //     if (inside_plate == true){

            //         for (int ii=0;ii<Nx;ii++){
            //             for (int jj=0;jj<Ny;jj++){
            //                 for (int kk=(plate_max+1);kk<top_plate_min;kk++){
            //                     //Calculate SD
            //                     int m = kk*Nx*Ny+jj*Nx+ii;
            //                     signedDistance[m] = double(kk-plate_max);
            //                 }
            //                 for (int kk=(top_plate_min-1);kk>plate_max;kk--){
            //                     //Calculate SD
            //                     int m = kk*Nx*Ny+jj*Nx+ii;
            //                     dist = double(top_plate_min-kk);
            //                     if (dist < signedDistance[m]){
            //                         signedDistance[m] = dist;
            //                     }
            //                 }
            //             }
            //         }

            //     } else {
            //         // for (int ii=0;ii<Nx;ii++){
            //         //     for (int jj=0;jj<Ny;jj++){
            //         //         for (int kk=plate_min;kk<plate_max+1;kk++){
            //         //             //Calculate SD
            //         //             int m = kk*Nx*Ny+jj*Nx+ii;
            //         //             dist = sqrt( (ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k) );
            //         //             if (dist < signedDistance_Invert[m]){
            //         //                 signedDistance_Invert[m] = dist;
            //         //             }
            //         //         }
            //         //         for (int kk=top_plate_min;kk<top_plate_max+1;kk++){
            //         //             //Calculate SD
            //         //             int m = kk*Nx*Ny+jj*Nx+ii;
            //         //             dist = sqrt( (ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k) );
            //         //             if (dist < signedDistance_Invert[m]){
            //         //                 signedDistance_Invert[m] = dist;
            //         //             }
            //         //         }
            //         //     }
            //         // }   
            //     }
            //}
        }
    }



    if (rank==0){
	    printf("LIBB_id \n");
    }
    for (int k=0;k<Nz;k++){
	    for (int j=0;j<Ny;j++){
		    for (int i=0;i<Nx;i++){
			    int n=k*Nx*Ny+j*Nx+i;
			    printf("%i ",LIBB_id[n]);
		    }
		    printf("\n");
	    }
	    printf("\n\n");
    }

    if (rank==0){
            printf("signedDistance \n");
    }
    for (int k=0;k<Nz;k++){
            for (int j=0;j<Ny;j++){
                    for (int i=0;i<Nx;i++){
                            int n=k*Nx*Ny+j*Nx+i;
                            printf("%.2f ",signedDistance[n]);
                    }
                    printf("\n");
            }
            printf("\n\n");
    }

    if (rank==0){
            printf("signedDistance Invert\n");
    }
    for (int k=0;k<Nz;k++){
            for (int j=0;j<Ny;j++){
                    for (int i=0;i<Nx;i++){
                            int n=k*Nx*Ny+j*Nx+i;
                            printf("%.2f ",signedDistance_Invert[n]);
                    }
                    printf("\n");
            }
            printf("\n\n");
    }


    for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            for (int k=0;k<Nz;k++){
                int m = k*Nx*Ny+j*Nx+i;
                signedDistance[m] = signedDistance[m] - signedDistance_Invert[m];
            }
        }
    }    

    
}




double RayTrace::ray_trace_pp(double *  vec1 /* ray_origin */,
                              double *  vec2 /* plane_pp1 */,
                              double D_pp1,
                              double *  vec3 /* plane_pp2 */,
                              double D_pp2) {
    
    double tmp1 = mag(vec2); if (tmp1 == 0) tmp1 = 1.0;
    double tmp2 = mag(vec3); if (tmp2 == 0) tmp2 = 1.0;
    
    dp1 = dot(vec1, vec2); dp2 = dot(vec1, vec3);
    
    double distance1 = -(dp1 - D_pp1)/tmp1; double distance2 = -(dp2 + D_pp2)/tmp2;
    
    return std::max(distance1,distance2);
}


/*  Ray-Plane Intersection  https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
*/
double RayTrace::ray_trace_pp_q(double *  lnaught /* ray_origin */, double * l /* ray direction */,
                      double *  plane_pp1 /* plane_pp1 */, double D_pp1 /* distance from origin- need to generalize this for phi,theta!=0 */,
                      double *  plane_pp2 /* plane_pp2 */, double D_pp2) {
    
    /* ptrs */
    double * pnaught; pnaught = new double[3];
     double * normal; normal = new double[3];
    
    /* registers */
    double p0x,p0y,p0z;
    double tmp;
    double denominator = 1;
    double numerator = 1;
    
    /* get normalized ray direction */
    tmp = mag(l); if (tmp == 0) tmp = 1.0;
    l[0] /= tmp; l[1] /= tmp;  l[2] /= tmp;
    
    tmp = mag(plane_pp1); if (tmp == 0) tmp = 1.0; // used if plane is rotated
    /* normal vector to the plane pp1 */
    normal[0] = plane_pp1[0]/tmp;
    normal[1] = plane_pp1[1]/tmp;
    normal[2] = plane_pp1[2]/tmp;
    

    /* distance point from world origin (needs to be parallelized!) */
    p0x = 0;  p0y = 0;  p0z = D_pp1;
    pnaught[0] = p0x;  pnaught[1] = p0y;  pnaught[2] = p0z;
    
    /* l0 - p0 */
    subtract_vectors(pnaught /* ray origin "lo" */,
                     lnaught /* distance from world origin point "" */,
                     o_minus_c /* Return: "l0-p0"*/);
    
    /* -(p0-l0).n */
    numerator = dot(o_minus_c,normal); if (tmp == 0) tmp = 1.0;
    /* -(p0-l0).n / (l.n) */
    denominator = dot(l,normal);
    if (denominator == 0.0) denominator = 1.0;
    dp1 = numerator / denominator;  // Result for plane 1
    if (dp1 < 0 ) dp1 = 9;
    

    tmp = mag(plane_pp2); if (tmp == 0) tmp = 1.0; // used if plane is rotated
    /* normal vector to the plane pp1 */
    normal[0] = plane_pp2[0]/tmp;
    normal[1] = plane_pp2[1]/tmp;
    normal[2] = plane_pp2[2]/tmp;
    

    /* distance point from world origin (needs to be parallelized!) */
    p0x = 0;  p0y = 0;  p0z = D_pp2;
    pnaught[0] = p0x;  pnaught[1] = p0y;  pnaught[2] = p0z;
    
    /* l0 - p0 */
    subtract_vectors(pnaught /* ray origin "lo" */,
                     lnaught /* distance from world origin point "" */,
                     o_minus_c /* Return: "l0-p0"*/);


    
    /* -(p0-l0).n */
    numerator = dot(o_minus_c,normal); if (tmp == 0) tmp = 1.0;
    /* -(p0-l0).n / (l.n) */
    denominator = dot(l,normal);
    if (denominator == 0.0) denominator = 1.0;
    dp2 = numerator / denominator;  // Result for plane 1
    if (dp2 < 0 ) dp2 = 9;
    
    
    
    // std::cout << "dp1=" << dp1 << " dp2=" << dp2 << std::endl;
  //  double lg = 1000;
    delete[] pnaught;
    delete[] normal;
    return std::min(dp2,dp1);
}






double RayTrace::ray_trace_st(double *  vec1 /* ray_origin */,
                              double *  vec2 /* plane_pp1 */, double D_st1,
                              double *  vec3 /* plane_pp2 */, double D_st2,
                              double *  vec4 /* plane_st3 */, double D_st3,
                              double *  vec5 /* plane_st4 */, double D_st4) {}

double RayTrace::ray_trace_st_q(double *  vec1 /* ray_origin */,
                                double *  vec2 /* plane_st1 */, double D_st1,
                                double *  vec3 /* plane_st2 */, double D_st2,
                                double *  vec4 /* plane_st3 */, double D_st3,
                                double *  vec5 /* plane_st4 */, double D_st4) {}


void ComputeLIBBqDistances_st(char * LIBB_id,
                              double * qDistances, double * signedDistance, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, Domain &Dm) {}

void ComputeLIBBqDistances_ct(char * LIBB_id,
                              double * qDistances, double * signedDistance, int & Nx, int & Ny, int & Nz, double & Lx, int & nspheres, int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz, int rank, Domain &Dm) {
    // NYI
}

int poreCount(char * id,
              int Nx, int Ny, int Nz) {
    size_t count = 0;
    int dat;
    size_t m;
    for (int k=1;k<Nz-1;k++){ // printf("k=%d\n",k);
        for (int j=1;j<Ny-1;j++){
            for (int i=1;i<Nx-1;i++){
                m = i + j*Nx + k*Nx*Ny;
                dat = id[m];
                //if (dat > 0) count++;
                count++;
            }
        }
    }
    return count;
}


void ComputeHWBBqDistances(char * sid, double * qDistances, int & Nx, int & Ny, int & Nz, size_t & count) {
    size_t n,nn,nn2,N;
    int i,j,k;
    N = Nx*Ny*Nz;
    
    size_t c = 0; // For compressed qDistance array
    /*
     You can pre-compute some q's based upon the neighbors, so
     go get your periodic BC code and run through all pairs and populate qDistance
     if both left and right neighbors are solid
     */
    for (n=0; n<N; n++){  // Maybe this is not correct - not sure yet.
        if (!(sid[n] > 0)) {
            k = n/(Nx*Ny);
            j = (n-Nx*Ny*k)/Nx;
            i = n-Nx*Ny*k-Nx*j;
            
            
            nn = n-1;                               // neighbor index (get convention)
            if (i-1<0)        nn += Nx;            // periodic BC along the x-boundary
            //........................................................................
            nn2 = n+1;                              // neighbor index (get convention)
            if (!(i+1<Nx))    nn2 -= Nx;            // periodic BC along the x-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+0*count] = 0.5;  qDistances[c+1*count] = 0.5;  }
            
            
            nn = n-Nx;                            // neighbor index (get convention)
            if (j-1<0)        nn += Nx*Ny;        // Perioidic BC along the y-boundary
            //........................................................................
            nn2 = n+Nx;                            // neighbor index (get convention)
            if (!(j+1<Ny))    nn2 -= Nx*Ny;        // Perioidic BC along the y-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+2*count] = 0.5;  qDistances[c+3*count] = 0.5;  }
            
            
            nn = n-Nx*Ny;                        // neighbor index (get convention)
            if (k-1<0)        nn += Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            //........................................................................
            nn2 = n+Nx*Ny;                        // neighbor index (get convention)
            if (!(k+1<Nz))    nn2 -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+4*count] = 0.5;  qDistances[c+5*count] = 0.5;  }
            
            
            nn = n-Nx-1;                        // neighbor index (get convention)
            if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
            if (j-1<0)            nn += Nx*Ny;    // Perioidic BC along the y-boundary
            //........................................................................
            nn2 = n+Nx+1;                        // neighbor index (get convention)
            if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
            if (!(j+1<Ny))        nn2 -= Nx*Ny;    // Perioidic BC along the y-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+6*count] = 0.5;  qDistances[c+7*count] = 0.5;  }
            
            
            nn = n+Nx-1;                        // neighbor index (get convention)
            if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
            if (!(j+1<Ny))        nn -= Nx*Ny;    // Perioidic BC along the y-boundary
            //........................................................................
            nn2 = n-Nx+1;                        // neighbor index (get convention)
            if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
            if (j-1<0)            nn2 += Nx*Ny;    // Perioidic BC along the y-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+8*count] = 0.5;  qDistances[c+9*count] = 0.5;  }
            
            
            nn = n-Nx*Ny-1;                        // neighbor index (get convention)
            if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
            if (k-1<0)            nn += Nx*Ny*Nz;    // Perioidic BC along the z-boundary
            //........................................................................
            nn2 = n+Nx*Ny+1;                        // neighbor index (get convention)
            if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
            if (!(k+1<Nz))        nn2 -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+10*count] = 0.5;  qDistances[c+11*count] = 0.5;  }
            
            
            nn = n+Nx*Ny-1;                        // neighbor index (get convention)
            if (i-1<0)            nn += Nx;        // periodic BC along the x-boundary
            if (!(k+1<Nz))        nn -= Nx*Ny*Nz;    // Perioidic BC along the z-boundary
            //........................................................................
            nn2 = n-Nx*Ny+1;                        // neighbor index (get convention)
            if (!(i+1<Nx))        nn -= Nx;        // periodic BC along the x-boundary
            if (k-1<0)            nn2 += Nx*Ny*Nz;    // Perioidic BC along the z-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+12*count] = 0.5;  qDistances[c+13*count] = 0.5;  }
            
            
            nn = n-Nx*Ny-Nx;                    // neighbor index (get convention)
            if (j-1<0)        nn += Nx*Ny;        // Perioidic BC along the y-boundary
            if (k-1<0)        nn += Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            //........................................................................
            nn2 = n+Nx*Ny+Nx;                    // neighbor index (get convention)
            if (!(j+1<Ny))    nn -= Nx*Ny;        // Perioidic BC along the y-boundary
            if (!(k+1<Nz))    nn2 -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+14*count] = 0.5;  qDistances[c+15*count] = 0.5;  }
            
            
            nn = n+Nx*Ny-Nx;                    // neighbor index (get convention)
            if (j-1<0)        nn += Nx*Ny;        // Perioidic BC along the y-boundary
            if (!(k+1<Nz))    nn -= Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            //........................................................................
            nn2 = n-Nx*Ny+Nx;                    // neighbor index (get convention)
            if (!(j+1<Ny))    nn -= Nx*Ny;        // Perioidic BC along the y-boundary
            if (k-1<0)        nn2 += Nx*Ny*Nz;        // Perioidic BC along the z-boundary
            if (sid[nn] == 0 && sid[nn2] == 0) {  qDistances[c+16*count] = 0.5;  qDistances[c+17*count] = 0.5;  }
            
            
            
        }
      c++; // For compressed qDistance array
    }
   
    
}


inline void save_3Ddouble_parallel(int rank,
                                   double* field, int Nx, int Ny, int Nz, int depth, std::string s) {
    //assert(depth >= 0);
    char cstr[s.size() + 1];
    s.copy(cstr, s.size() + 1);
    cstr[s.size()] = '\0';
    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(LocalRankString,".%05d",rank);
    sprintf(LocalRankFilename,"%s%s.txt",cstr,LocalRankString);
    FILE *f2 = fopen(LocalRankFilename, "w");
    
    double var = 0;
    size_t c = 0;
    for (int z =depth; z < Nz-depth; z++) {
        for (int y =depth; y < Ny-depth; y++) {
            for (int x=depth; x < Nx-depth; x++) {
            //    int n = x + Nx*y + Nx*Ny*z;
                
                var = field[c]; c++;
                fprintf(f2, "%d %d %d %.2f\n",x, y, z, var);
                
            } // i // iproc //    printf("\n");
        } // j // jproc //  printf("\n");
    }  // k // kproc
    fclose(f2);
    
    // printf("LINE=%d\n",__LINE__);
    //    delete[] rankArrayNoPM;
}


