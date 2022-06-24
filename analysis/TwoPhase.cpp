#include "analysis/TwoPhase.h"

#include <stdint.h>
#include <limits.h>

#if SIZE_MAX == UCHAR_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#error "- -"
#endif


#include <math.h>
#include "analysis/pmmc.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "common/Array.hpp"
#include "analysis/analysis.h"

#include "shared_ptr.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"

#include <iterator>
#include <vector>

#include "analysis/distance.h"

#define BLOB_AVG_COUNT 40

// Array access for averages defined by the following
#define VOL 0
#define TRIMVOL 1
#define PRS 2
#define AWN 3
#define ANS 4
#define LWNS 5
#define JWN 6
#define KWN 7
#define CWNS 8
#define KNWNS 9
#define KGWNS 10
#define VX 11
#define VY 12
#define VZ 13
#define VSQ 14
#define VWNX 15
#define VWNY 16
#define VWNZ 17
#define VNSX 18
#define VNSY 19
#define VNSZ 20
#define GWNXX 21
#define GWNYY 22
#define GWNZZ 23
#define GWNXY 24
#define GWNXZ 25
#define GWNYZ 26
#define TRAWN 27
#define TRJWN 28
#define CMX 29
#define CMY 30
#define CMZ 31
#define EULER 32
#define INTCURV 33
#define GXX 34
#define GYY 35
#define GZZ 36
#define GXY 37
#define GXZ 38
#define GYZ 39
using std::cout; using std::ofstream;
using std::endl; using std::string;

#define PI 3.14159265359

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsometimes-uninitialized"
#pragma GCC diagnostic ignored "-Wreorder"
TwoPhase::TwoPhase(std::shared_ptr <Domain> dm, std::array<size_t,3> sub):
n_nw_pts(0), n_ns_pts(0), n_ws_pts(0), n_nws_pts(0), n_local_sol_pts(0), n_local_nws_pts(0),
n_nw_tris(0), n_ns_tris(0), n_ws_tris(0), n_nws_seg(0), n_local_sol_tris(0),
nc(0), kstart(0), kfinish(0), fluid_isovalue(0), solid_isovalue(0),    Volume(0),nodeVolume(0),subNodeVolume(0),
TIMELOG_local(nullptr), TIMELOG_global(nullptr), NWPLOG(nullptr), WPLOG(nullptr), GEOMETRY_local(nullptr), GEOMETRY_global(nullptr),
Dm(dm), NumberComponents_WP(0), NumberComponents_NWP(0), trimdist(0),
porosity(0), volume(0), poreVol(0), awn(0), ans(0), aws(0), lwns(0), aw(0), an(0), Xn(0), pwn(0), wp_volume(0), nwp_volume(0),nwSecondCC(0),
As(0), vol_w(0), vol_n(0), vol_w_tminus(0), vol_n_tminus(0), vwndnw(0), vwndnw_global(0), Jws(0), Kws(0), euler(0), euler_global(0), eulerW(0), eulerW_global(0),
pn(0), pw(0), pn_global(0), pw_global(0), pn_tminus(0), pw_tminus(0), pn_tminus_global(0), pw_tminus_global(0), swn(0), swn_global(0), vol_w_global(0), vol_n_global(0), vol_w_tminus_global(0), vol_n_tminus_global(0),
awn_global(0), ans_tcenter_global(0),aws_global(0), lwns_global(0), efawns(0), efawns_count(0), efawns_global(0), JnsS(0), JnsS_global(0), KnsS(0),KnsS_global(0),
JwsS(0),JwsS_global(0), KwsS(0),KwsS_global(0),eulerS(0),eulerS_global(0),
Jwn(0), Jwn_global(0),  Kwn(0), Kwn_global(0), KNwns(0), KNwns_global(0), JwnW_global(0), KwnW_global(0), JwnW(0), KwnW(0),nwndvwns(0),time_flag(-1),an_old(0),an_new(0),
KGwns(0), KGwns_global(0), trawn(0), trawn_global(0), trJwn(0), trJwn_global(0),
trRwn(0), trRwn_global(0), nwp_volume_global(0), wp_volume_global(0), Kws_global(0), Jws_global(0), Kns_global(0), Jns_global(0),deldel_IIprimevwndnn(0),
As_global(0), wwndnw_global(0), wwnsdnwn_global(0), Jwnwwndnw_global(0), dEs(0), dAwn(0), dAns(0), LEFT_an(0), RIGHT_an(0),  left_wwnJnx(0),IIprimevwndnn(0),
left_wwnJny(0), awn_tcenter(0), centroid_ewn(0),centroid_Jwn(0),centroid_Knvwndnn(0),
left_wwnJnz(0), IIprimezzwwnz(0),deldelddIIprimezznnz(0),nnz_ns(0),
right_wwnJnx(0),gradzgradzIIprimezzvwndnn(0),deld_nnz(0), gradzgradz_an(0),nnz_left(0),nnz_right(0), nnz_center(0),ewnwwny(0),deld_ewnwwn2(0), left_gradz(0),Jnvwndnn(0),wwnnn(0),nw_speed(0),
right_gradz(0),wwnz_tcenter(0),wwnx_tcenter(0),wwny_tcenter(0),vwnz(0),vwnx(0),vwny(0),nnnnxz(0),nnnnyz(0),nnnnzz(0),vnsdnn_tcenter(0),gradpz(0),pot_nnw(0),ans_tcenter(0),
right_wwnJny(0),del_IIprimevwndnn(0),HOvwnz(0),vwndnn_tcenter(0), mesh_euler(0),delpdnn(0),Jnnnz(0),deldIIprimez(0),deld_JwnIIprimeMINUSnnnnz(0),deldelddIIprimennz(0),nzmagvwn(0),
right_wwnJnz(0),nw_nnnnxz(0), nw_nnnnyz(0),nw_nnnnzz(0), gradz_an(0),grady_an(0),gradx_an(0),deld_wwnJn(0),deld_ewnwwn(0),deldelddIIprimezzvwndnn(0),vwndnnnnx(0), vwndnnnny(0), vwndnnnnz(0),JnJnnnz(0),
Jns_tcenter(0),Kns_tcenter(0),grad_potential(0),TwoKnnz(0),deld_wwnKn(0),ens(0),ews(0), vnsx(0), vnsy(0), vnsz(0), nnnn_nszz(0), nnnn_nsyz(0), nnnn_nsxz(0),nnsz(0),vwsz(0), n_wsz(0),nwndvwnsnnz(0),offset_distance(0),
nwnwwndnnz(0),grav_nw(0),Jwn_tcenter_global(0), Kwn_tcenter_global(0),awn_tcenter_global(0),awn_tminus_global(0),Jwn_tminus_global(0),an_tminus_global(0),aw_tminus_global(0),wwnz_global(0),wwny_global(0),wwnx_global(0),nGaussianCurvature(0),
nGaussianCurvature_global(0),integrated_cpx(0),integrated_cpy(0),integrated_cpz(0),Jnx(0),Jny(0),Jnz(0),gradz_an_tminus(0),gradz_an_tplus(0), d_gradz_an_dt(0),awn_tplus(0),awn_tminus(0),an_tminus(0),aw_tminus(0),Jwn_tminus(0),deriv(0),
subdivide( sub )
{
    Nx=dm->Nx; Ny=dm->Ny; Nz=dm->Nz;
    npx = Dm->nprocx();
    npy = Dm->nprocy();
    npz = Dm->nprocz();

    int amin = Dm->amin;
    int amax = Dm->amax;
    nodeVolume = 1.0*(double(amax-amin)*double(Nx-2)*double(Ny-2));
    MPI_Allreduce(&nodeVolume,&Volume,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    VolumeFraction.resize(Nx,Ny,Nz); VolumeFraction.fill(0);
    
    Null3.resize(3);
    Null9.resize(9);
    ewnwwn.resize(3); ewnwwn.fill(0);
    wwnJn.resize(3);  wwnJn.fill(0);
    vwndnnIIprime.resize(9); vwndnnIIprime.fill(0);
    
    dendt_.resize(Nx,Ny,Nz);  dendt_.fill(0);
    
    PhaseVolume.resize(Nx,Ny,Nz); PhaseVolume.fill(0);
    PhaseVolume_old.resize(Nx,Ny,Nz); PhaseVolume_old.fill(0);
    
    temp_tplus.resize(Nx,Ny,Nz);   temp_tplus.fill(0);
    
    zero.resize(Nx,Ny,Nz);         zero.fill(0);
    VFmask.resize(Nx,Ny,Nz);       VFmask.fill(0);
    
    DensityA.resize(Nx,Ny,Nz);     DensityA.fill(0);
    DensityB.resize(Nx,Ny,Nz);     DensityB.fill(0);
    
    PhiField.resize(Nx,Ny,Nz);     PhiField.fill(0);
    ID.resize(Nx,Ny,Nz);           ID.fill(0);
    Label_WP.resize(Nx,Ny,Nz);      Label_WP.fill(0);
    Label_NWP.resize(Nx,Ny,Nz);     Label_NWP.fill(0);
    NxArray.resize(Nx,Ny,Nz);      NxArray.fill(0);
    NyArray.resize(Nx,Ny,Nz);      NyArray.fill(0);
    NzArray.resize(Nx,Ny,Nz);      NzArray.fill(0);
    
    nwPress.resize(Nx,Ny,Nz);   nwPress.fill(0);
    
    Jns_.resize(Nx,Ny,Nz);   Jns_.fill(0);
    Kns_.resize(Nx,Ny,Nz);   Kns_.fill(0);
    
    nnz_.resize(Nx,Ny,Nz);   nnz_.fill(0);
    Knwwnx.resize(Nx,Ny,Nz);    Knwwnx.fill(0);
    Knwwny.resize(Nx,Ny,Nz);    Knwwny.fill(0);
    Knwwnz.resize(Nx,Ny,Nz);    Knwwnz.fill(0);
    
    ewn_.resize(Nx,Ny,Nz);  ewn_.fill(0);
    Jwn_.resize(Nx,Ny,Nz);  Jwn_.fill(0);
    
    speed_.resize(Nx,Ny,Nz); speed_.fill(0);
    nw_speed_.resize(Nx,Ny,Nz); nw_speed_.fill(0);
    maggradphi.resize(Nx,Ny,Nz);  maggradphi.fill(0);
    
    Kwn_.resize(Nx,Ny,Nz); Kwn_.fill(0);
    
    PhaseField.resize(Nx,Ny,Nz);  PhaseField.fill(0);
    PhaseField_old.resize(Nx,Ny,Nz);  PhaseField_old.fill(0);
    SDn.resize(Nx,Ny,Nz);           SDn.fill(0);
    SDw.resize(Nx,Ny,Nz);           SDw.fill(0);
    
    SDn_tplus.resize(Nx,Ny,Nz);       SDn_tplus.fill(0);
    SDn_tminus.resize(Nx,Ny,Nz);      SDn_tminus.fill(0);
    
    A.resize(Nx,Ny,Nz); A.fill(0);
    Aco.resize(Nx,Ny,Nz);  Aco.fill(0);
    D.resize(Nx,Ny,Nz);  D.fill(0);
    Dco.resize(Nx,Ny,Nz);  Dco.fill(0);
    E.resize(Nx,Ny,Nz);  E.fill(0);
    
    A2.resize(Nx,Ny,Nz); A2.fill(0);
    Aco2.resize(Nx,Ny,Nz);  Aco2.fill(0);
    D2.resize(Nx,Ny,Nz);  D2.fill(0);
    Dco2.resize(Nx,Ny,Nz);  Dco2.fill(0);
    E2.resize(Nx,Ny,Nz);  E2.fill(0);
    
    wwnnn_.resize(Nx,Ny,Nz); wwnnn_.fill(0);
    
    ewnwwnx_.resize(Nx,Ny,Nz);      ewnwwnx_.fill(0);
    ewnwwny_.resize(Nx,Ny,Nz);      ewnwwny_.fill(0);
    ewnwwnz_.resize(Nx,Ny,Nz);      ewnwwnz_.fill(0);

    Jnwwnx.resize(Nx,Ny,Nz);       Jnwwnx.fill(0);
    Jnwwny.resize(Nx,Ny,Nz);       Jnwwny.fill(0);
    Jnwwnz.resize(Nx,Ny,Nz);       Jnwwnz.fill(0);
    
    SDs.resize(Nx,Ny,Nz);           SDs.fill(0);
    Phase.resize(Nx,Ny,Nz);         Phase.fill(0);
    Press.resize(Nx,Ny,Nz);         Press.fill(0);
    dPdt.resize(Nx,Ny,Nz);          dPdt.fill(0);
    
    SDs_x.resize(Nx,Ny,Nz);         SDs_x.fill(0);      // Gradient of the signed distance
    SDs_y.resize(Nx,Ny,Nz);         SDs_y.fill(0);
    SDs_z.resize(Nx,Ny,Nz);         SDs_z.fill(0);
    SDn_x.resize(Nx,Ny,Nz);         SDn_x.fill(0);      // Gradient of the signed distance
    SDn_y.resize(Nx,Ny,Nz);         SDn_y.fill(0);
    SDn_z.resize(Nx,Ny,Nz);         SDn_z.fill(0);
    
    Delpvwndnnz.resize(Nx,Ny,Nz); Delpvwndnnz.fill(0);
    
    SDw_x.resize(Nx,Ny,Nz);         SDw_x.fill(0);      // Gradient of the signed distance
    SDw_y.resize(Nx,Ny,Nz);         SDw_y.fill(0);
    SDw_z.resize(Nx,Ny,Nz);         SDw_z.fill(0);
    
    vwndnn_.resize(Nx,Ny,Nz);      vwndnn_.fill(0);
    
    normalw_x.resize(Nx,Ny,Nz);         normalw_x.fill(0);      // Gradient of the signed distance
    normalw_y.resize(Nx,Ny,Nz);         normalw_y.fill(0);      // Gradient of the signed distance
    normalw_z.resize(Nx,Ny,Nz);         normalw_z.fill(0);      // Gradient of the signed distance
    
    CField.resize(Nx,Ny,Nz);  CField.fill(0);
    
    nw_nnnnzz_.resize(Nx,Ny,Nz);    nw_nnnnzz_.fill(0);
    vwnz_.resize(Nx,Ny,Nz);         vwnz_.fill(0);
  
    Phasemc.resize(Nx,Ny,Nz);  Phasemc.fill(0);
    
    MagVel.resize(Nx*Ny*Nz);   std::fill(MagVel.begin(), MagVel.end(), 0);
    
    
    DelPhi.resize(Nx,Ny,Nz);        DelPhi.fill(0);
    Phase_tplus.resize(Nx,Ny,Nz);   Phase_tplus.fill(0);
    Phase_tminus.resize(Nx,Ny,Nz);  Phase_tminus.fill(0);
    Vel_x.resize(Nx,Ny,Nz);         Vel_x.fill(0);        // Gradient of the phase indicator field
    Vel_y.resize(Nx,Ny,Nz);         Vel_y.fill(0);
    Vel_z.resize(Nx,Ny,Nz);         Vel_z.fill(0);
    BlobID.resize(Nx,Ny,Nz);        BlobID.fill(0);
    
    GradPhiX.resize(Nx,Ny,Nz);  GradPhiX.fill(0);
    GradPhiY.resize(Nx,Ny,Nz);  GradPhiY.fill(0);
    GradPhiZ.resize(Nx,Ny,Nz);  GradPhiZ.fill(0);
    
    //.........................................
    // Allocate cube storage space
    CubeValues.resize(2,2,2);
    nw_tris.resize(3,20);
    ns_tris.resize(3,20);
    ws_tris.resize(3,20);
    nws_seg.resize(2,20);
    local_sol_tris.resize(3,18);
    nw_pts=DTMutableList<Point>(20);
    ns_pts=DTMutableList<Point>(20);
    ws_pts=DTMutableList<Point>(20);
    nws_pts=DTMutableList<Point>(20);
    local_nws_pts=DTMutableList<Point>(20);
    local_sol_pts=DTMutableList<Point>(20);
    tmp=DTMutableList<Point>(20);
    //.........................................
    Values.resize(20);
    DistanceValues.resize(20);
    KGwns_values.resize(20);
    KNwns_values.resize(20);
    InterfaceSpeed.resize(20);
    NormalVector.resize(60);
    //.........................................
    vn.resize(3);
    vw.resize(3);
    
    nw.resize(3);
    vwn.resize(3);
    vws.resize(3);
    vns.resize(3);
    vwns.resize(3);
    Gwn.resize(6);
    Gwn_tminus.resize(6);
    Gns.resize(6);
    Gws.resize(6);
    Gwns.resize(6);
    vn_global.resize(3);
    vw_global.resize(3);
    
    
    nw_global.resize(3);
    vwn_global.resize(3);
    vws_global.resize(3);
    vns_global.resize(3);
    vwns_global.resize(3);
    Gwn_global.resize(6);
    Gns_global.resize(6);
    Gws_global.resize(6);
    Gwns_global.resize(6);
    Gwn_tminus_global.resize(6);
    // Allocate subdomain variables
    
    Eta_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Eta_sub1.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    trimdist_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    swn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    porosity_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    volume_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    poreVol_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    awn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    tmp_tmp.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    ans_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    aws_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    lwns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    wp_volume_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    nwp_volume_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    As_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    vol_w_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    vol_n_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    pn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    pw_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    efawns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    euler_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    eulerW_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Kn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Jn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    An_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    
    Jwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Kwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    JwnW_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    KwnW_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Jws_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Kws_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Jns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Kns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    
    JwsS_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    KwsS_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    JnsS_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    KnsS_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    
    eulerS_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    
    KNwns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    KGwns_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    trawn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    trJwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    trRwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    wwndnw_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    wwnsdnwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    Jwnwwndnw_sub.resize( { subdivide[0], subdivide[1], subdivide[2] } );
    
    
    nw_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    vn_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    vw_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    
    vws_sub.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    vns_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    vwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    vwns_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 3 } );
    Gwn_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    Gns_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    Gws_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    Gwns_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    pwn_sub.resize(  { subdivide[0], subdivide[1], subdivide[2], 6 } );
    Xn_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    aw_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    an_sub.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    
    an_sub1.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    an_sub2.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    an_sub3.resize( { subdivide[0], subdivide[1], subdivide[2], 6 } );
    
    qDistances.resize(Nx*Ny*Nz*18); qDistances.fill(-9);
    
    
    
    // w^wn
    wwn.resize(3);
    wwn_global.resize(3);
    wwn_sub.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    
    
    nwOwnO.resize(3);
    nwOwnO_global.resize(3);
    nwOwnO_sub.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    nwOwnO2.resize(3);
    nwOwnO_global2.resize(3);
    nwOwnO_sub2.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    nwOwsO.resize(3);
    nwOwsO_global.resize(3);
    nwOwsO_sub.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    nwOwsO2.resize(3);
    nwOwsO_global2.resize(3);
    nwOwsO_sub2.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    nwOnsO2.resize(3);
    nwOnsO_global2.resize(3);
    nwOnsO_sub2.resize({ subdivide[0], subdivide[1], subdivide[2], 3 } );
    
    
    // Create global files
    if (Dm->rank()==0){
        TIMELOG_global = fopen("timelog.tcat","a+");
        if (fseek(TIMELOG_global,0,SEEK_SET) == fseek(TIMELOG_global,0,SEEK_CUR)) {
            // If timelog is empty, write a short header to list the averages
            const char header[] = "time volume porosity "
            "s^w p^w* p^n* eps^{wn}* eps^{ns}* eps^{ws}* J_n^{wn}* K_n^{wn}* l^{wns}* c^{wns} kappa_N^{wns}* kappa_G^{wns}* "
            "v_x^w* v_y^w* v_z^w* v_x^n* v_y^n* v_z^n* v_x^{wn}* v_y^{wn}* v_z^{wn}* v_x^{wns}* v_y^{wns}* v_z^{wns}* "
            "G_{xx}^{wn}* G_{yy}^{wn}* G_{zz}^{wn}* G_{xy}^{wn}* G_{xz}^{wn}* G_{yz}^{wn}* "
            "G_{xx}^{ws}* G_{yy}^{ws}* G_{zz}^{ws}* G_{xy}^{ws}* G_{xz}^{ws}* G_{yz}^{ws}* "
            "G_{xx}^{ns}* G_{yy}^{ns}* G_{zz}^{ns}* G_{xy}^{ns}* G_{xz}^{ns}* G_{yz}^{ns}* "
            "p_n^{wn}* X^n* eps^w eps^n "
            "G_{xx}^{wns}* G_{yy}^{wns}* G_{zz}^{wns}* G_{xy}^{wns}* G_{xz}^{wns}* G_{yz}^{wns}* "
            "v_x^{ws}* v_y^{ws}* v_z^{ws}* v_x^{ns}* v_y^{ns}* v_z^{ns}* "
            "wwnx wwny wwnz Ca Re Jns Kns p^w_minus p^n_minus eps^w_minus eps^n_minus eps*{wn}_minus J_n^{wn}_minus G_{xx}^{wn}_minus";
            fprintf(TIMELOG_global,"%s\n",header);
        }
        
        /*NWPLOG = fopen("components.NWP.tcat","a+");
        fprintf(NWPLOG,"time label vol p^n awn ans J_n^{wn} K_n^{wn} lwns c^{wns} kappa_N^{wns} kappa_G^{wns} ");
        fprintf(NWPLOG,"v_x^n v_y^n v_z^n v_x^{wn} v_y^{wn} v_z^{wn} v_x^{wns} v_y^{wns} v_z^{wns} ");
        fprintf(NWPLOG,"G_{xx}^{wn} G_{yy}^{wn} G_{zz}^{wn} G_{xy}^{wn} G_{xz}^{wn} G_{yz}^{wn} ");
        fprintf(NWPLOG,"G_{xx}^{ns} G_{yy}^{ns} G_{zz}^{ns} G_{xy}^{ns} G_{xz}^{ns} G_{yz}^{ns} ");
        fprintf(NWPLOG,"Euler Cx Cy Cz vsq \n");
        
        WPLOG = fopen("components.WP.tcat","a+");
        fprintf(WPLOG,"time label vol p^w awn aws J_n^{wn} K_n^{wn} lwns c^{wns} kappa_N^{wns} kappa_G^{wns} ");
        fprintf(WPLOG,"v_x^w v_y^w v_z^w v_x^{wn} v_y^{wn} v_z^{wn} v_x^{wns} v_y^{wns} v_z^{wns} ");
        fprintf(WPLOG,"G_{xx}^{wn} G_{yy}^{wn} G_{zz}^{wn} G_{xy}^{wn} G_{xz}^{wn} G_{yz}^{wn} ");
        fprintf(WPLOG,"G_{xx}^{ws} G_{yy}^{ws} G_{zz}^{ws} G_{xy}^{ws} G_{xz}^{ws} G_{yz}^{ws} ");
        fprintf(WPLOG,"Cx Cy Cz vsq \n");*/

        NWPLOG = fopen("components.NWP.tcat","a+");
        fprintf(NWPLOG,"time label vol p^n awn ans v_x^n v_y^n v_z^n v_x^{wn} v_y^{wn} v_z^{wn} v_x^{ns} v_y^{ns} v_z^{ns} \n");
        
        // WPLOG = fopen("components.WP.tcat","a+");
        // fprintf(WPLOG,"time label vol \n");     
    }
    
    // Create the local files
    
    if ( subdivide[0] * subdivide[1] * subdivide[2] > 1 ) {
        char tmp[64];
        sprintf(tmp,"%05d",Dm->rank());
        std::string rankStr( tmp );
        std::string filename = "timelog.tcat." + rankStr;
        TIMELOG_local = fopen(filename.c_str(),"a+");
        
        
        
        const char header[] = "time volume porosity "
        "s^w p^w p^n eps^{wn} eps^{ns} eps^{ws} J_n^{wn} K_n^{wn} l^{wns} c^{wns} kappa_N^{wns} kappa_G^{wns} "
        "v_x^w v_y^w v_z^w v_x^n v_y^n v_z^n v_x^{wn} v_y^{wn} v_z^{wn} v_x^{wns} v_y^{wns} v_z^{wns} "
        "G_{xx}^{wn} G_{yy}^{wn} G_{zz}^{wn} G_{xy}^{wn} G_{xz}^{wn} G_{yz}^{wn} "
        "G_{xx}^{ws} G_{yy}^{ws} G_{zz}^{ws} G_{xy}^{ws} G_{xz}^{ws} G_{yz}^{ws} "
        "G_{xx}^{ns} G_{yy}^{ns} G_{zz}^{ns} G_{xy}^{ns} G_{xz}^{ns} G_{yz}^{ns} "
        "p_n^{wn} X^n eps^w eps^n "
        "G_{xx}^{wns} G_{yy}^{wns} G_{zz}^{wns} G_{xy}^{wns} G_{xz}^{wns} G_{yz}^{wns} "
        "v_x^{ws} v_y^{ws} v_z^{ws} v_x^{ns} v_y^{ns} v_z^{ns} wwnx wwny wwnz "
        "Ca Re Jns Kns pw_top pw_bot pn_top pn_bot eps^w_minus eps^n_minus eps*{wn}_minus";
        
        fprintf(TIMELOG_local,"%s\n",header);
       
        }
    
    Initialize();
    }
#pragma GCC diagnostic pop

// Destructor
TwoPhase::~TwoPhase()
{
    //delete [] TempID;
    if ( NWPLOG!=nullptr ) { fclose(NWPLOG); }
    if ( WPLOG!=nullptr ) { fclose(WPLOG); }
    if ( TIMELOG_local!=nullptr ) { fclose(TIMELOG_local); }
  
    
}

void TwoPhase::CommunicateField(int time_flag) {
    if (time_flag == -1) Dm->CommunicateMeshHalo(SDn_tminus);
    if (time_flag == 0) Dm->CommunicateMeshHalo(SDn);
    if (time_flag == 1) Dm->CommunicateMeshHalo(SDn_tplus);
}

void TwoPhase::ColorToSignedDistance(double Beta, DoubleArray &ColorData, DoubleArray &DistData)
{
    Dm->CommunicateMeshHalo(ColorData);
    for (int k=1; (size_t)k<Nz-1; k++){
        for (int j=1; (size_t)j<Ny-1; j++){
            for (int i=1; (size_t)i<Nx-1; i++){
                DistData(i,j,k) = Beta*ColorData(i,j,k);
            }
        }
    }
    Dm->CommunicateMeshHalo(DistData);
}



void TwoPhase::ComputeDelPhi()
{
    Dm->CommunicateMeshHalo( SDn);
    double nx,ny,nz;
    int strideY = Nx;
    int strideZ = Nx*Ny;
    int nn;
    int ijk;
    int i,j,k;
    double m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
    for (k=1; (size_t)k<Nz-1; k++){
        for (j=1; (size_t)j<Ny-1; j++){
            for (i=1; (size_t)i<Nx-1; i++){
                ijk = i + j*Nx + k*Nx*Ny;
                
                nn=ijk+1; m1= SDn(nn);
                nn=ijk-1; m2= SDn(nn);
                nn=ijk+strideY; m3= SDn(nn);
                nn=ijk-strideY; m4= SDn(nn);
                nn=ijk+strideZ; m5= SDn(nn);
                nn=ijk-strideZ; m6= SDn(nn);

           
                nn=ijk+1+strideY; m7= SDn(nn);
                nn=ijk-1-strideY; m8= SDn(nn);
                nn=ijk+1-strideY; m9= SDn(nn);
                nn=ijk-1+strideY; m10= SDn(nn);
                nn=ijk+1+strideZ; m11= SDn(nn);
                nn=ijk-1-strideZ; m12= SDn(nn);
                nn=ijk+1-strideZ; m13= SDn(nn);
                nn=ijk-1+strideZ; m14= SDn(nn);
                nn=ijk+strideY+strideZ; m15= SDn(nn);
                nn=ijk-strideY-strideZ; m16= SDn(nn);
                nn=ijk+strideY-strideZ; m17= SDn(nn);
                nn=ijk-strideY+strideZ; m18= SDn(nn);


                    // (2,4)
                nx = m1/6. - m10/12. + m11/12. - m12/12. + m13/12. - m14/12. - m2/6. + m7/12. - m8/12. + m9/12.;
                ny = m10/12. + m15/12. - m16/12. + m17/12. - m18/12. + m3/6. - m4/6. + m7/12. - m8/12. - m9/12.;
                nz = m11/12. - m12/12. - m13/12. + m14/12. + m15/12. - m16/12. - m17/12. + m18/12. + m5/6. - m6/6.;
             
               
                
                DelPhi(ijk) = sqrt(nx*nx + ny*ny + nz*nz);
            }
        }
    }
    Dm->CommunicateMeshHalo(DelPhi);
}

void TwoPhase::ComputeSDw() {
    
    int i,j,k,n;
    for (k=0; (size_t)k<Nz; k++){
        for (j=0; (size_t)j<Ny; j++){
            for (i=0; (size_t)i<Nx; i++){
                n = k*Nx*Ny+j*Nx+i;
                if (!(Dm->id[n] > 0)){
                    SDw(i,j,k) = -1;
                }
                else {
                    SDw(i,j,k) = -SDn(i,j,k);
                }
            }
        }
    }
}



void TwoPhase::Initialize()
{
    
    pw_lhs = 0;
    pw_rhs = 0;
    pn_lhs = 0;
    pn_lhs = 0;
    pw_lhs_count = 0;
    pw_rhs_count = 0;
    pn_lhs_count = 0;
    pn_lhs_count = 0;
    
    MAXZ = Nz-1;
    MINZ = 1;
    an2 = 0; an_old = 0; an_new = 0;
    awn_tplus = awn_tminus = Jwn_tminus = 0;
    nGaussianCurvature = nGaussianCurvature_global = 0;
    Kns_global = Jns_global = 0;
    ans_tcenter = ans_tcenter_global = 0;
    Jwn_tcenter_global = Kwn_tcenter_global = 0;
    awn_tcenter_global = awn_tminus_global = wwnz_global = wwny_global = wwnx_global = 0;
    dendt = 0;
    nw_speed = 0;
    nw_nnnnxz = nw_nnnnyz = nw_nnnnzz = 0;
    Ca_global = 0;
    Re_global = 0;

    wwnnn = 0;
  
    

    
  
    awn2 = 0;
    awn2_global = 0;
    Jwn2 = Kwn2 = Jwn2_global = Kwn2_global = 0;
    local_nw_phase_volume = 0.0;
    global_nw_phase_volume = 0.0;
    mass_local = 0.0;
    mass_global = 0.0;
    porosity = 0.0;
    porosity_sub.fill(0);
    volume = 0.0;
    volume_sub.fill(0);
    porosity_global = 0;
    volume_global = 0;
 
    
    Jns_tcenter = Kns_tcenter = 0;
    deld_wwndnn = 0;
    nwndvwns = 0;
    // tri mesh volume from nw Comp
    aw_sub.fill(0);
    aw = 0;
    aw_global = 0;
    aw_tminus = 0;
    aw_tminus_global = 0;
    
    Jwn_tminus_global = 0;
    
    gradpz = 0;
    
    // tri mesh volume from SignedVolume function
    an_sub.fill(0);
    an = 0;
    an_global = 0;
    an_tminus = 0;
    an_tminus_global = 0;
    
    nwEuler_global = nsEuler_global = 0;
    
 
    an_sub3.fill(0);
    SDn.fill(0);
    SDn_tminus.fill(0);
    SDn_tplus.fill(0);
    nw_speed_.fill(0);
    speed_.fill(0);
    dendt_.fill(0);
    ewnwwny_.fill(0);
    ewnwwnz_.fill(0);
    
    vnsdnn_tcenter = 0;
    
    efawns_count = 0;
    efawns = 0.0;
    
    
    vn.fill(0);
    vn_global.fill(0);
    
    vol_w_sub.fill( 0 );
    vol_n_sub.fill( 0 );
    vol_w = vol_n =0.0;
    vol_w_tminus = vol_n_tminus = 0.0;
    vol_w_global = 0;
    vol_n_global = 0;
    vol_w_tminus_global = vol_n_tminus_global = 0.0;
    
   
    
    pn_sub.fill( 0.0 );
    pw_sub.fill( 0.0 );
    pn = pw = pwn = 0.0;
    pn_tminus = pw_tminus = 0.0;
    pn_global = pw_global = 0.0;
    pn_tminus_global = pw_tminus_global = 0.0;
    
    
    //
    Eta_sub.fill(0);
    Eta = Eta_global = 0;
    total_volume_global = 0;
    
    
    Eta_sub1.fill(0);

    an_sub.fill(0);
    
    
    wsGlobal_dS = 0;
    nsGlobal_dS =0;
    nwGlobal_dS =0;
    
    wsGlobal_CCs =0;
    nsGlobal_CCs =0;
    nwGlobal_CCs =0;
    
  
    
    nEuler = wEuler = sEuler = 0.0;
    nwEuler = nsEuler = wsEuler = 0.0;
    nGlobalEuler = sGlobalEuler = wGlobalEuler = 0.0;
    
    trimdist=0.0;
    fluid_isovalue = 0.0;
    solid_isovalue= 0.0;
    // Initialize the averaged quantities
    porosity = 0.0;
    volume = 0.0;
    awn = aws = ans = lwns = 0.0;
    nwp_volume = wp_volume = 0.0;
    As = 0.0;
    pn = pw = pwn = 0.0;
    vw.fill( 0 );
    vn.fill( 0 );
    vwn.fill( 0 );
    vws.fill( 0 );
    vns.fill( 0 );
    vwns.fill( 0 );
    Gwn.fill( 0 );
    Gwn_tminus.fill( 0 );
    Gws.fill( 0 );
    Gns.fill( 0 );
    Gwns.fill( 0 );
    vol_w = vol_n =0.0;
    KGwns = KNwns = 0.0;
    Jwn = Kwn = efawns = efawns_count = efawns_global_count = 0.0;
    Jws = Kws = Jns = Kns = 0;
    trJwn = trawn = trRwn = 0.0;
    euler = Jn = An = Kn = eulerW = 0.0;
    wwndnw = 0.0; wwnsdnwn = 0.0; Jwnwwndnw=0.0;
    vwndnw = 0.0;
    // Initialize subdomain variables
    trimdist_sub.fill( 0 );
    porosity_sub.fill( 0 );
    volume_sub.fill( 0 );
    porosity_global = 0.0;
    volume_global = 0.0;
    poreVol_sub.fill( 0 );
    awn_sub.fill( 0 );
    tmp_tmp.fill( 0 );
    ans_sub.fill( 0 );
    aws_sub.fill( 0 );
    lwns_sub.fill( 0 );
    wp_volume_sub.fill( 0 );
    nwp_volume_sub.fill( 0 );
    As_sub.fill( 0 );
    vol_w_sub.fill( 0 );
    vol_n_sub.fill( 0 );
    pn_sub.fill( 0 );
    pw_sub.fill( 0 );
    pwn_sub.fill ( 0 );
    efawns_sub.fill( 0 );
    euler_sub.fill( 0 );
    euler = 0;
    euler_global = 0;
    eulerW_sub.fill( 0 );
    eulerW = 0;
    eulerW_global = 0;
    Kn_sub.fill( 0 );
    Jn_sub.fill( 0 );
    An_sub.fill( 0 );
    
    
    Jwn_sub.fill( 0 );
    Kwn_sub.fill( 0 );
    Jws_sub.fill( 0 );
    Kws_sub.fill( 0 );
    Jns_sub.fill( 0 );
    Kns_sub.fill( 0 );


    KwnW_sub.fill( 0 );
    JwnW_sub.fill( 0 );
    JwnW = 0;
    KwnW = 0;
    JwnW_global = 0;
    KwnW_global = 0;
    
    
    JwsS_sub.fill( 0 );
    KwsS_sub.fill( 0 );
    JnsS_sub.fill( 0 );
    KnsS_sub.fill( 0 );
    eulerS_sub.fill( 0 );
    eulerS = 0;
    eulerS_global = 0;
    
    KNwns_sub.fill( 0 );
    KGwns_sub.fill( 0 );
    trawn_sub.fill( 0 );
    trJwn_sub.fill( 0 );
    trRwn_sub.fill( 0 );
    wwndnw_sub.fill( 0 );
    wwnsdnwn_sub.fill( 0 );
    Jwnwwndnw_sub.fill( 0 );
    vn_sub.fill( 0 );
    vw_sub.fill( 0 );
    vws_sub.fill( 0 );
    vns_sub.fill( 0 );
    vwn_sub.fill( 0 );
    vwns_sub.fill( 0 );
    Gwn_sub.fill( 0 );
    Gns_sub.fill( 0 );
    Gws_sub.fill( 0 );
    Gwns_sub.fill( 0 );

    wwn.fill( 0 );
    wwn_global.fill( 0 );
    wwn_sub.fill( 0 );
    
    nwOwnO_sub.fill( 0 );
    nwOwnO.fill( 0 );
    nwOwnO_global.fill( 0 );
    
    nwOwnO_sub2.fill( 0 );
    nwOwnO2.fill( 0 );
    nwOwnO_global2.fill( 0 );
    
    nwOwsO_sub2.fill( 0 );
    nwOwsO2.fill( 0 );
    nwOwsO_global2.fill( 0 );
    
    nwOnsO_sub2.fill( 0 );
    nwOnsO2.fill( 0 );
    nwOnsO_global2.fill( 0 );
}


void TwoPhase::UpdateSolid()
{
    Dm->CommunicateMeshHalo(SDs);
    //...........................................................................
    // Gradient of the Signed Distance function
    //...........................................................................
    pmmc_MeshGradient(SDs,SDs_x,SDs_y,SDs_z,Nx,Ny,Nz);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_z);
    //...........................................................................
}

inline void MeanFilterX(DoubleArray &Mesh /*in*/, DoubleArray &Out, IntArray & ID){
        for (int k=1; k<(int)Mesh.size(2)-1; k++){
                for (int j=1; j<(int)Mesh.size(1)-1; j++){
                        for (int i=1; i<(int)Mesh.size(0)-1; i++){
                                double sum;
                            if (ID(i,j,k) == 0) {
                                sum=Mesh(i+1,j,k)+Mesh(i-1,j,k)+Mesh(i,j+1,k)+Mesh(i,j-1,k)+
                                                +Mesh(i,j,k+1)+Mesh(i,j,k-1)
                                +0.5*(Mesh(i+1,j+1,k) +Mesh(i+1,j-1,k) +Mesh(i-1,j+1,k) +Mesh(i-1,j-1,k)
                                +Mesh(i+1,j,k+1) +Mesh(i+1,j,k-1) +Mesh(i-1,j,k+1) +Mesh(i-1,j,k-1)
                                
                                +Mesh(i,j+1,k+1) +Mesh(i,j-1,k+1) +Mesh(i,j+1,k-1) +Mesh(i,j-1,k-1));
                            
                            
                                Out(i,j,k) = sum/18.0;
                            }
                        }
                }
        }
}



void TwoPhase::UpdateMeshValues()
{
 
    MPI_Bcast(&time_flag,1,MPI_INT,0,Dm->Comm);
    
    Dm->CommunicateMeshHalo(SDs);
    
    pmmc_MeshGradient(SDs,SDs_x,SDs_y,SDs_z,Nx,Ny,Nz);
    
    // Gradient of the phase indicator field
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_z);
    //...........................................................................
    
    if (time_flag > -1){
        Dm->CommunicateMeshHalo(SDn_tminus);
            
        pmmc_MeshGradient( SDn_tminus, SDn_x, SDn_y, SDn_z,Nx,Ny,Nz);
        //...........................................................................
        // Gradient of the phase indicator field
        //...........................................................................
        Dm->CommunicateMeshHalo( SDn_x);
        //...........................................................................
        Dm->CommunicateMeshHalo( SDn_y);
        //...........................................................................
        Dm->CommunicateMeshHalo( SDn_z); 
    }

      
    Dm->CommunicateMeshHalo(SDn);
        
    Dm->CommunicateMeshHalo(Vel_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(Vel_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(Vel_z);
    //...........................................................................
    
    Dm->CommunicateMeshHalo(Press);
    //...........................................................................
   
    
   
    Dm->CommunicateMeshHalo(SDn);
    
//    MeanFilterX(SDn , Phase , ID);
//    MeanFilterX(Phase , SDn , ID);
//    MeanFilterX(SDn , Phase , ID);
//    MeanFilterX(Phase , SDn , ID);
//    MeanFilterX(SDn , Phase , ID);
//    MeanFilterX(Phase , SDn , ID);
//    MeanFilterX(SDn , Phase , ID);
//    MeanFilterX(Phase , SDn , ID);
//    MeanFilterX(SDn , Phase , ID);
//    MeanFilterX(Phase , SDn , ID);
//    MeanFilterX(SDn , Phase , ID);
//    MeanFilterX(Phase , SDn , ID);
//    MeanFilterX(SDn , Phase , ID);
//    MeanFilterX(Phase , SDn , ID);
//    MeanFilterX(SDn , Phase , ID);
//    MeanFilterX(Phase , SDn , ID);
    
    Dm->CommunicateMeshHalo(SDn);
        
    pmmc_MeshGradient( SDn, SDn_x, SDn_y, SDn_z,Nx,Ny,Nz);

    Dm->CommunicateMeshHalo( SDn_x);

    Dm->CommunicateMeshHalo( SDn_y);

    Dm->CommunicateMeshHalo( SDn_z);
   
    if (time_flag > -1){
        pmmc_MeshGradient( SDn_tplus, SDn_x, SDn_y, SDn_z,Nx,Ny,Nz);
        //...........................................................................
        // Gradient of the phase indicator field
        //...........................................................................
        Dm->CommunicateMeshHalo( SDn_x);
        //...........................................................................
        Dm->CommunicateMeshHalo( SDn_y);
        //...........................................................................
        Dm->CommunicateMeshHalo( SDn_z);
    }  
    
}








void TwoPhase::ComputePhiDistance(DoubleArray &Field) {
    
    
    Array<char> id_solid(Nx,Ny,Nz);
    
    double factor,temp;
    factor=0.5/0.9;
    // Initialize to -1,1 (segmentation)
    for (int k=0; (size_t)k<Nz; k++){
        for (int j=0; (size_t)j<Ny; j++){
            for (int i=0; (size_t)i<Nx; i++){
                int n = k*Nx*Ny+j*Nx+i;
                double value = Field(n);
                if (value > 0)    id_solid(n) = 1;
                else              id_solid(n) = 0;
                temp = factor*log((1.0+value)/(1.0-value));
                if (value > 0.8) Field(n) = 2.94*factor;
                else if (value < -0.8) Field(n) = -2.94*factor;
                else Field(n) = temp;
            }
        }
    }
    
    
    CalcDist(Field,id_solid,*Dm);
    
}

void TwoPhase::SavePhaseIntoOldPhase(DoubleArray & NewField, DoubleArray & Field ) {
    
    int ijk;
    int strideY = Nx;
    int strideZ = Nx*Ny;
    for (int k=1; (size_t)k < Nz-1; k++) {
        for (int j=1; (size_t)j < Ny-1; j++) {
            for (int i=1; (size_t)i < Nx-1; i++) {
                ijk = i + j * strideY + k * strideZ;
                NewField(ijk) = Field(ijk);
            }
        }
    }
    
}


void TwoPhase::ComputeVolumeFraction(DoubleArray & Field, DoubleArray & Field_x, DoubleArray& Field_y, DoubleArray& Field_z, DoubleArray& VolumeFraction, bool UpdateNormalToSolid, int rmin, int rmax, int geometry)
{
    
    /* Adjusting VF corners correctly: Needs to be correctly implemented for for nprocs > 8 */
    Field(0) = Field(Nx-1) = Field((Ny-1)*Nx) = Field((Ny-1)*Nx + Nx-1) = 1;
    Field((Nz-1)*Nx*Ny) = Field((Nz-1)*Nx*Ny+Nx-1) = Field((Nz-1)*Nx*Ny+(Ny-1)*Nx) = Field((Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1) = 0;
    
    
    Dm->CommunicateMeshHalo(Field);
    pmmc_MeshGradient(Field,Field_x,Field_y,Field_z,Nx,Ny,Nz);
    Dm->CommunicateMeshHalo(Field_x);
    Dm->CommunicateMeshHalo(Field_y);
    Dm->CommunicateMeshHalo(Field_z);
    
    
 
    size_t n;

    Point A,B,C,P,Q,D,E,F;
    
    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    
    int x,y,z;
    
    double ax,ay,az,amag;


    visMyMesh vis_n;

    PhaseVolume.fill(0);

    int count = 0;
    
    
    for (int k=1; k<Nz-1; k++)
    for (int j=1; j<Ny-1; j++)
    for (int i=1; i<Nx-1; i++) {
        n = i + j*Nx + k*Nx*Ny;
        n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
        n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;

        pmmc_ConstructLocalCube(Field,  zero,  zero,  zero,  zero, solid_isovalue, fluid_isovalue,nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
                                local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                                n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                                n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg, i, j, k, Nx, Ny, Nz);
        
        ax = ay = az = amag = 0;
        
        if (n_ws_tris > 0) {
            /* Allocate new mesh and vertex iterator */
            wsMyMesh m_ws;
            wsMyMesh::VertexIterator vi;
            /* Create convex hull mesh object */
            chMyMesh CH_ws;

            /* Add points to the mesh */
            for (int r=0;r<n_ws_tris;r++) {
                A = ws_pts(ws_tris(0,r));
                B = ws_pts(ws_tris(1,r));
                C = ws_pts(ws_tris(2,r));

                double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                double normal_x = Field_x(i,j,k);
                double normal_y = Field_y(i,j,k);
                double normal_z = Field_z(i,j,k);

                if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                    P = A;
                    A = C;
                    C = P;
                }

                // Remap the points
                A.x += 1.0*Dm->iproc()*(Nx-2);
                A.y += 1.0*Dm->jproc()*(Nx-2);
                A.z += 1.0*Dm->kproc()*(Nx-2);
                B.x += 1.0*Dm->iproc()*(Nx-2);
                B.y += 1.0*Dm->jproc()*(Nx-2);
                B.z += 1.0*Dm->kproc()*(Nx-2);
                C.x += 1.0*Dm->iproc()*(Nx-2);
                C.y += 1.0*Dm->jproc()*(Nx-2);
                C.z += 1.0*Dm->kproc()*(Nx-2);
                wsMyMesh::FaceIterator fi = vcg::tri::Allocator<wsMyMesh>::AddFace(m_ws,wsMyMesh::CoordType ( A.x, A.y, A.z),wsMyMesh::CoordType ( B.x, B.y, B.z),wsMyMesh::CoordType ( C.x, C.y, C.z));

                vcg::tri::UpdateNormal<wsMyMesh>::PerFaceNormalized(m_ws);
                vcg::tri::Stat<wsMyMesh>::ComputeMeshArea(m_ws);
                double area = 0.5*DoubleArea(*fi);
                double nx = (*fi).N()[0];
                double ny = (*fi).N()[1];
                double nz = (*fi).N()[2];

                // multiply by area?
                ax += nx*area;
                ay += ny*area;
                az += nz*area;

            }

            //if (UpdateNormalToSolid) {
//                amag = sqrt(ax*ax + ay*ay + az*az);
//                if (amag == 0) amag = 1;
//                ax/=amag;
//                ay/=amag;
//                az/=amag;
          //  printf("yolo");
            GradPhiX(i,j,k) = ax;
            GradPhiY(i,j,k) = ay;
            GradPhiZ(i,j,k) = az;
           // }


            /* Add corner points to surface for volume calculations */
            for (int p = 0; p < 8; p++) {
                int m = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
                if ( Field(m) < 0 ) {
                    z = m/(Nx*Ny); y = (m-Nx*Ny*z)/Nx; x = m-Nx*Ny*z-Nx*y;
                    Q.z = z; Q.z += 1.0*Dm->kproc()*(Nx-2);
                    Q.y = y; Q.y += 1.0*Dm->jproc()*(Nx-2);
                    Q.x = x; Q.x += 1.0*Dm->iproc()*(Nx-2);

                    vi = vcg::tri::Allocator<wsMyMesh>::AddVertices(m_ws,1);
                    vi->P()=wsMyMesh::CoordType ( Q.x, Q.y, Q.z); vi++;
                }
            }

            vcg::tri::Clean<wsMyMesh>::RemoveDuplicateFace(m_ws);
            vcg::tri::Clean<wsMyMesh>::RemoveDuplicateVertex(m_ws,false);

            vcg::tri::ConvexHull<wsMyMesh,chMyMesh>::ComputeConvexHull(m_ws,CH_ws);
            VolumeFraction(i,j,k) =  vcg::tri::Stat<chMyMesh>::ComputeMeshVolume(CH_ws);
        }
    } // i



for (int k=1; k<Nz-1; k++) {
    for (int j=1; j<Ny-1; j++) {
        for (int i=1; i<Nx-1; i++) {
            n = i + j*Nx + k*Nx*Ny;
            if (Field(n) < 0.0 && VolumeFraction(n) == 0 ){
                VolumeFraction(n) = 1.0;
            }
        }
    }
}

Dm->CommunicateMeshHalo(VolumeFraction);
Dm->CommunicateMeshHalo(GradPhiX);
Dm->CommunicateMeshHalo(GradPhiY);
Dm->CommunicateMeshHalo(GradPhiZ);
    
    Dm->CommunicateMeshHalo(VolumeFraction);
//    char LocalRankString[8];
//    char LocalRankFilename[40];
//    sprintf(LocalRankString,"%05d",Dm->rank());
//    sprintf(LocalRankFilename,"%s%s%s","visdata.",LocalRankString,".off");
//    vcg::tri::io::ExporterOFF<visMyMesh>::Save(vis_n,LocalRankFilename, 1);
    
    
//    double DVALUE = 0;
//    printf("Convex Hull result:\n");
//    for (int i=3;i<4;i++){
//        for (int j=1;j<Ny-1;j++){
//            for (int k=1;k<Nz-1;k++){
//                int n=k*(Nx)*(Ny)+j*(Nz)+i;
//                DVALUE = VolumeFraction(i,j,k);
//                printf("%.2f ",DVALUE);
//            }
//            printf("\n");
//        }
//        printf("\n\n");
//    }
    
    
//    for (int k=1; k<Nz-1; k++)
//    for (int j=1; j<Ny-1; j++)
//    for (int i=1; i<Nx-1; i++) {
//        n = i + j*Nx + k*Nx*Ny;
//        VolumeFraction(n) = 1-PhaseVolume(n);
//    }
//
    
//    VolumeFraction = PhaseVolume;
}






void TwoPhase::ComputeLocal()
{
    nwRootMesh m_nw;
    nwRootMesh  m_n;
    
    nsMyMesh m_ns;
    wsMyMesh m_ws;
    
    nwMyMesh m_ni1;
    nwsMyMesh m_nws_2;
    
    Point P,A,B,C;
    Point Q,D,E,F;
    Point R,G,H,I;
    
    
    Dm->CommunicateMeshHalo(Phase);
    pmmc_MeshGradient( Phase, SDn_x, SDn_y, SDn_z,Nx,Ny,Nz);
    
    
    
    
   
    
//    pmmc_MeshGradient( SDn, SDn_x, SDn_y, SDn_z,Nx,Ny,Nz);
//    //...........................................................................
//    // Gradient of the phase indicator field
//    //...........................................................................
    Dm->CommunicateMeshHalo( SDn_x);
    //...........................................................................
    Dm->CommunicateMeshHalo( SDn_y);
    //...........................................................................
    Dm->CommunicateMeshHalo( SDn_z);
    //...........................................................................
    
    
    pmmc_MeshGradient( SDs, SDs_x, SDs_y, SDs_z,Nx,Ny,Nz);
    //...........................................................................
    // Gradient of the phase indicator field
    //...........................................................................
    Dm->CommunicateMeshHalo( SDs_x);
    //...........................................................................
    Dm->CommunicateMeshHalo( SDs_y);
    //...........................................................................
    Dm->CommunicateMeshHalo( SDs_z);
    //...........................................................................

    int amin = Dm->amin;
    int amax = Dm->amax;

    double pqx,pqy,pqz;
    double avgpqx,avgpqy,avgpqz;
    const int N2x = ( Nx-2 ) / subdivide[0];
    const int N2y = ( Ny-2 ) / subdivide[1];
    const int N2z = ( Nz-2 ) / subdivide[2];
    
    DoubleArray nwOwsO_tmp(3), nwOwnO_tmp(3);
    DoubleArray nw_tmp(3);
    DoubleArray wwn_tmp(3);
    DoubleArray vwn_tmp(3), vwns_tmp(3), vws_tmp(3), vns_tmp(3), Gwn_tmp(6), Gws_tmp(6), Gns_tmp(6), Gwns_tmp(6);
    
    Phase(0) = Phase(Nx-1) = Phase((Ny-1)*Nx) = Phase((Ny-1)*Nx + Nx-1) = NAN;
    Phase((Nz-1)*Nx*Ny) = Phase((Nz-1)*Nx*Ny+Nx-1) = Phase((Nz-1)*Nx*Ny+(Ny-1)*Nx) = Phase((Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1) = NAN;
    
    SDs(0) = SDs(Nx-1) = SDs((Ny-1)*Nx) = SDs((Ny-1)*Nx + Nx-1) = NAN;
    SDs((Nz-1)*Nx*Ny) = SDs((Nz-1)*Nx*Ny+Nx-1) = SDs((Nz-1)*Nx*Ny+(Ny-1)*Nx) = SDs((Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1) = NAN;
        
    for (int k=0; k<amax; k++) {
        const int k2 = std::min<int>( ( k - 1) / N2z, subdivide[2] - 1 );
        for (int j=0; (size_t)j<Ny-1; j++){
            const int j2 = std::min<int>( ( j - 1 ) / N2y, subdivide[1] - 1 );
            for (int i=0; (size_t)i<Nx-1; i++) {
                const int i2 = std::min<int>( ( i - 1 ) / N2x, subdivide[0] - 1 );
                
                vwn_tmp.fill( 0 );  vws_tmp.fill( 0 );  vns_tmp.fill( 0 );  vwns_tmp.fill( 0 );
                Gwn_tmp.fill( 0 );  Gws_tmp.fill( 0 );  Gns_tmp.fill( 0 );  Gwns_tmp.fill( 0 );
                wwn_tmp.fill( 0 ); nwOwnO_tmp.fill( 0 ); nwOwsO_tmp.fill( 0 );
                
                /*
                 Fix loop for range of 1 to Nx-1....
                 
                 
                 */
                // Compute volume averages
                int n = i + j*Nx + k*Nx*Ny;
                double vs = VFmask(n);

                double phase_val =  SDn(n);
                double delPhi_val = DelPhi(n);
                char id_val = Dm->id[n];
                
                volume_sub(i2,j2,k2) += 1.0;

                if (id_val == 0) {
                    porosity_sub(i2,j2,k2) += 1.0;
                }
                if (vs < 0.5) {
                    if (phase_val <= 0) {
                        vw_sub(i2,j2,k2,0) += Vel_x(n);
                        vw_sub(i2,j2,k2,1) += Vel_y(n);
                        vw_sub(i2,j2,k2,2) += Vel_z(n);

                        wp_volume_sub(i2,j2,k2) += 1.0;
                        aw_sub(i2,j2,k2) += 1.0;

                        if (delPhi_val < 1e-4) {  
                            pw_sub(i2,j2,k2) += Press(n);  
                            vol_w_sub(i2,j2,k2) += 1; 
                            if (Press(n) > 50000.0){
                                printf("Pressure,i,j,k: %f %i %i %i \n",Press(n),i,j,k);
                            }
                        }


                        
                        if (k==amin) { pw_lhs += Press(n); pw_lhs_count += 1; }
                        if (k==(amax-1)) { pw_rhs += Press(n); pw_rhs_count += 1; }
                        
                    }
                }

                if (vs < 0.5) {
                    if (phase_val > 0) {
                        vn_sub(i2,j2,k2,0) += Vel_x(n);
                        vn_sub(i2,j2,k2,1) += Vel_y(n);
                        vn_sub(i2,j2,k2,2) += Vel_z(n);
                        grav_nw += Fz;
                        nwp_volume_sub(i2,j2,k2) += 1.0;
                        an_sub(i2,j2,k2) += 1.0;

                        if (delPhi_val < 1e-4) { pn_sub(i2,j2,k2) += Press(n);   vol_n_sub(i2,j2,k2) += 1; }
                        
                        if (k==amin) { pn_lhs += Press(n); pn_lhs_count += 1; }
                        if (k==(amax-1)) { pn_rhs += Press(n); pn_rhs_count += 1; }

                    }
                }
                
              

             
                //...........................................................................
                n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
                n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;
                
                //...........................................................................
                pmmc_ConstructLocalCube(SDs,  Phase,  SDn_x,  SDn_y,  SDn_z, solid_isovalue, fluid_isovalue, nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
                                        local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                                        n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                                        n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
                                        i, j, k, Nx, Ny, Nz);
                
                if (n_nw_tris > 0) {
                    for (int r=0;r<n_nw_tris;r++) {
                        A = nw_pts(nw_tris(0,r));
                        B = nw_pts(nw_tris(1,r));
                        C = nw_pts(nw_tris(2,r));

                        double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                        double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                        double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                        double normal_x =  -SDn_x(i,j,k);
                        double normal_y =  -SDn_y(i,j,k);
                        double normal_z =  -SDn_z(i,j,k);

                        if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                            P = A;
                            A = C;
                            C = P;
                        }
                        A.x += 1.0*Dm->iproc()*(Nx-2);
                        A.y += 1.0*Dm->jproc()*(Nx-2);
                        A.z += 1.0*Dm->kproc()*(Nx-2);
                        B.x += 1.0*Dm->iproc()*(Nx-2);
                        B.y += 1.0*Dm->jproc()*(Nx-2);
                        B.z += 1.0*Dm->kproc()*(Nx-2);
                        C.x += 1.0*Dm->iproc()*(Nx-2);
                        C.y += 1.0*Dm->jproc()*(Nx-2);
                        C.z += 1.0*Dm->kproc()*(Nx-2);
                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                        
                        A.y += offset_distance;
                        B.y += offset_distance;
                        C.y += offset_distance;
                        
                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_nw,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                        
//                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                    }
                }
                
                if (n_ns_tris > 0) {
                    for (int r=0;r<n_ns_tris;r++) {
                        A = ns_pts(ns_tris(0,r));
                        B = ns_pts(ns_tris(1,r));
                        C = ns_pts(ns_tris(2,r));
                        
                        double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                        double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                        double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);
                        
                        double normal_x = -SDs_x(i,j,k);
                        double normal_y = -SDs_y(i,j,k);
                        double normal_z = -SDs_z(i,j,k);
                        
                        if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                            P = A;
                            A = C;
                            C = P;
                        }
                        A.x += 1.0*Dm->iproc()*(Nx-2);
                        A.y += 1.0*Dm->jproc()*(Nx-2);
                        A.z += 1.0*Dm->kproc()*(Nx-2);
                        B.x += 1.0*Dm->iproc()*(Nx-2);
                        B.y += 1.0*Dm->jproc()*(Nx-2);
                        B.z += 1.0*Dm->kproc()*(Nx-2);
                        C.x += 1.0*Dm->iproc()*(Nx-2);
                        C.y += 1.0*Dm->jproc()*(Nx-2);
                        C.z += 1.0*Dm->kproc()*(Nx-2);
                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                       // printf("offset_distance=%f",offset_distance);
                        A.y += offset_distance;
                        B.y += offset_distance;
                        C.y += offset_distance;

                        vcg::tri::Allocator<nsMyMesh>::AddFace(m_ns,nsMyMesh::CoordType ( A.x, A.y, A.z),nsMyMesh::CoordType ( B.x, B.y, B.z),nsMyMesh::CoordType ( C.x, C.y, C.z));
                        
//                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                    }
                }
                if (n_ws_tris > 0) {
               //     printf("n_ws_tris=%d\n",n_ws_tris);
                    for (int r=0;r<n_ws_tris;r++) {
                        A = ws_pts(ws_tris(0,r));
                        B = ws_pts(ws_tris(1,r));
                        C = ws_pts(ws_tris(2,r));
                        
                        double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                        double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                        double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);
                        
                        double normal_x = SDs_x(i,j,k);
                        double normal_y = SDs_y(i,j,k);
                        double normal_z = SDs_z(i,j,k);
                        
                        if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                            P = A;
                            A = C;
                            C = P;
                        }
                        A.x += 1.0*Dm->iproc()*(Nx-2);
                        A.y += 1.0*Dm->jproc()*(Nx-2);
                        A.z += 1.0*Dm->kproc()*(Nx-2);
                        B.x += 1.0*Dm->iproc()*(Nx-2);
                        B.y += 1.0*Dm->jproc()*(Nx-2);
                        B.z += 1.0*Dm->kproc()*(Nx-2);
                        C.x += 1.0*Dm->iproc()*(Nx-2);
                        C.y += 1.0*Dm->jproc()*(Nx-2);
                        C.z += 1.0*Dm->kproc()*(Nx-2);
                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                        A.y += offset_distance;
                        B.y += offset_distance;
                        C.y += offset_distance;

                        vcg::tri::Allocator<wsMyMesh>::AddFace(m_ws,wsMyMesh::CoordType ( A.x, A.y, A.z),wsMyMesh::CoordType ( B.x, B.y, B.z),wsMyMesh::CoordType ( C.x, C.y, C.z));
                    }
                }
                if (n_nws_seg > 0) {
                    for (int r = 0; r < n_nws_seg; r++) {

                        nwsMyMesh m_nws;
                        nwsMyMesh m_nws_edges;
                        Q = nws_pts(nws_seg(0,r));
                        R = nws_pts(nws_seg(1,r));


                        nwMyMesh m_nw;

                        if (n_ns_tris > 0 ) {
                            for (int r=0;r<n_ns_tris;r++) {
                                A = ns_pts(ns_tris(0,r));
                                B = ns_pts(ns_tris(1,r));
                                C = ns_pts(ns_tris(2,r));

                                double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                                double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                                double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                                double normal_x = -SDs_x(i,j,k);
                                double normal_y = -SDs_y(i,j,k);
                                double normal_z = -SDs_z(i,j,k);

                                if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                                    P = A;
                                    A = C;
                                    C = P;
                                }

                                if (Q == A) {
                                    if (R == B) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                         vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));

                                    }
                                    if (R == C) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                }
                                if (Q == B) {
                                    if (R == A) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                    if (R == C) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                }
                                if (Q == C) {
                                    if (R == B) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                    if (R == A) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                }
                            }
                        }
                        if (n_nw_tris > 0 ) {
                          //  std::cout << "n_nw_tris for nws=" << n_nw_tris << std::endl;
                            for (int r=0;r<n_nw_tris;r++) {
                                A = nw_pts(nw_tris(0,r));
                                B = nw_pts(nw_tris(1,r));
                                C = nw_pts(nw_tris(2,r));

                                double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                                double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                                double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                                double normal_x =  -SDn_x(i,j,k);
                                double normal_y =  -SDn_y(i,j,k);
                                double normal_z =  -SDn_z(i,j,k);

                                if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                                    P = A;
                                    A = C;
                                    C = P;
                                }
                                if (Q == A) {
                                    if (R == B) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                             
                                        
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));
                                       vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                    if (R == C) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                                     
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));
                                        vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                }
                                if (Q == B) {
                                    if (R == A) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                                  
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));
                                        vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                    if (R == C) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                                 
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));
                                        vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                }
                                if (Q == C) {
                                    if (R == B) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                           
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));
                                        vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                    if (R == A) {
                                        A.x += 1.0*Dm->iproc()*(Nx-2);
                                        A.y += 1.0*Dm->jproc()*(Nx-2);
                                        A.z += 1.0*Dm->kproc()*(Nx-2);
                                        B.x += 1.0*Dm->iproc()*(Nx-2);
                                        B.y += 1.0*Dm->jproc()*(Nx-2);
                                        B.z += 1.0*Dm->kproc()*(Nx-2);
                                        C.x += 1.0*Dm->iproc()*(Nx-2);
                                        C.y += 1.0*Dm->jproc()*(Nx-2);
                                        C.z += 1.0*Dm->kproc()*(Nx-2);
                                        A.x +=0.5; A.y +=0.5; A.z+=0.5;
                                        B.x +=0.5; B.y +=0.5; B.z+=0.5;
                                        C.x +=0.5; C.y +=0.5; C.z+=0.5;
                        
                                        A.y += offset_distance;
                                        B.y += offset_distance;
                                        C.y += offset_distance;
                                        vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));
                                        vcg::tri::Allocator<nwsMyMesh>::AddFace(m_nws,nwsMyMesh::CoordType ( A.x, A.y, A.z),nwsMyMesh::CoordType ( B.x, B.y, B.z),nwsMyMesh::CoordType ( C.x, C.y, C.z));
                                        
//                                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                                    }
                                }
                            }
                        }
                        vcg::tri::Append<nwsMyMesh, nwsMyMesh>::Mesh(m_nws_2, m_nws, false /* selected */, false /* adjFlag */);
                        // Get segment vector
                       
                        nwsComputeQuantities(m_nws,m_nw,Vel_x,Vel_y,Vel_z,pqx,pqy,pqz,avgpqx,avgpqy,avgpqz);
                    }
                }
                
                
               
                
                
                if ( n_local_nws_pts > 0 ) {

                    wwnsdnwn_sub(i2,j2,k2) += pmmc_CommonCurveSpeed(CubeValues, dPdt, vwns_tmp,
                                                                     SDn_x,  SDn_y,  SDn_z,
                                                                    SDs_x,SDs_y,SDs_z,
                                                                    local_nws_pts,i,j,k,n_local_nws_pts);
                    vwns_sub(i2,j2,k2,0) += vwns_tmp(0);
                    vwns_sub(i2,j2,k2,1) += vwns_tmp(1);
                    vwns_sub(i2,j2,k2,2) += vwns_tmp(2);

                    pmmc_CurveCurvature( SDn, SDs,
                                         SDn_x,  SDn_y,  SDn_z,
                                        SDs_x,SDs_y,SDs_z,
                                        KNwns_values, KGwns_values, KNwns, KGwns,
                                        nws_pts, n_nws_pts, i, j, k);
                    KNwns_sub(i2,j2,k2) += KNwns;
                    KGwns_sub(i2,j2,k2) += KGwns;
                    lwns_sub(i2,j2,k2) +=  pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
                    pmmc_CurveOrientation(Gwns_tmp,local_nws_pts,n_local_nws_pts,i,j,k);
                    for (int m=0; m<6; m++)  Gwns_sub(i2,j2,k2,m) += Gwns_tmp(m);
                    
                }
                if (n_ws_tris > 0) {
                  //  std::cout << "n_ws_tris > 0" << std::endl;
                    aws_sub(i2,j2,k2) += pmmc_CubeSurfaceOrientation(Gws_tmp,ws_pts,ws_tris,n_ws_tris);
                    for (int m=0; m<6; m++) { Gws_sub(i2,j2,k2,m) += Gws_tmp(m); }
                    pmmc_InterfaceSpeed(dPdt,  SDn_x,  SDn_y,  SDn_z, CubeValues, ws_pts, ws_tris, NormalVector, InterfaceSpeed, vws_tmp, i, j, k, n_ws_pts, n_ws_tris);
                    vws_sub(i2,j2,k2,0) += vws_tmp(0);
                    vws_sub(i2,j2,k2,1) += vws_tmp(1);
                    vws_sub(i2,j2,k2,2) += vws_tmp(2);
                   
                }
                if (n_ns_tris > 0) {
                  //  std::cout << "n_ns_tris > 0" << std::endl;
                    ans_sub(i2,j2,k2) += pmmc_CubeSurfaceOrientation(Gns_tmp,ns_pts,ns_tris,n_ns_tris);
                    for (int m=0; m<6; m++) {  Gns_sub(i2,j2,k2,m) += Gns_tmp(m); }
                    pmmc_InterfaceSpeed(dPdt,  SDn_x,  SDn_y,  SDn_z,  CubeValues, ns_pts,ns_tris, NormalVector, InterfaceSpeed, vns_tmp, i, j, k, n_ns_pts, n_ns_tris);
                    vns_sub(i2,j2,k2,0) += vns_tmp(0);
                    vns_sub(i2,j2,k2,1) += vns_tmp(1);
                    vns_sub(i2,j2,k2,2) += vns_tmp(2);
                    euler_sub(i2,j2,k2) += geomavg_EulerCharacteristic(ns_pts,ns_tris,n_ns_pts,n_ns_tris,i,j,k);
                }
                
                if (n_nw_tris > 0) {
                    awn_sub(i2,j2,k2) += pmmc_CubeSurfaceOrientation(Gwn_tmp,nw_pts,nw_tris,n_nw_tris);
                    for (int m=0; m<6; m++) {  Gwn_sub(i2,j2,k2,m) += Gwn_tmp(m); }
                    pmmc_InterfaceSpeed(dPdt,  SDn_x,  SDn_y,  SDn_z,  CubeValues, nw_pts,nw_tris, NormalVector, InterfaceSpeed, vwn_tmp, i, j, k, n_nw_pts, n_nw_tris);
                    vwn_sub(i2,j2,k2,0) += vwn_tmp(0);
                    vwn_sub(i2,j2,k2,1) += vwn_tmp(1);
                    vwn_sub(i2,j2,k2,2) += vwn_tmp(2);
                     euler_sub(i2,j2,k2) += geomavg_EulerCharacteristic(nw_pts,nw_tris,n_nw_pts,n_nw_tris,i,j,k);
                    pwn_sub(i2,j2,k2) += pmmc_CubeSurfaceInterpValue(CubeValues,Press,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);
                }
            }
        }
    }
        
    char TimeString[8];
    char OffsetString[8];
    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(OffsetString,"%.1f",offset_distance);
    sprintf(TimeString,"%.05d",time_step);
    sprintf(LocalRankString,"%05d",Dm->rank());
    sprintf(LocalRankFilename,"%s%s%s%s%s%s",OffsetString,"_",TimeString,"_nwsdata.",LocalRankString,".off");
    vcg::tri::io::ExporterOFF<nwsMyMesh>::Save(m_nws_2,LocalRankFilename, 1);
    
    awn_tcenter = 0;
    Jwn_tcenter = 0;
    Kwn_tcenter = 0;
    vwndnn_tcenter = 0;

    ComputeAccurateNWInterfaceCurvatures(m_nw,&Jwn_tcenter,&Kwn_tcenter,&awn_tcenter,&vwndnn_tcenter, ewnwwn, wwnJn, vwndnnIIprime, Vel_x,Vel_y,Vel_z);
    nsComputeQuantities(m_ns,Vel_x,Vel_y,Vel_z);
    wsComputeQuantities(m_ws,Vel_x,Vel_y,Vel_z);
  
    if (pw_lhs_count > 0) pw_lhs /= double(pw_lhs_count); else pw_lhs = 0;
    if (pw_rhs_count > 0) pw_rhs /= double(pw_rhs_count); else pw_rhs = 0;
    if (pn_lhs_count > 0) pn_lhs /= double(pn_lhs_count); else pn_lhs = 0;
    if (pn_rhs_count > 0) pn_rhs /= double(pn_rhs_count); else pn_rhs = 0;
   
    if (isnan(pw_lhs)) pw_lhs = 0.0;
    if (isnan(pw_rhs)) pw_rhs = 0.0;
    if (isnan(pn_lhs)) pn_lhs = 0.0;
    if (isnan(pn_rhs)) pn_rhs = 0.0;

    //std::cout << "awn_tcenter" << awn_tcenter << std::endl;
   
    volume = volume_sub.sum();

    porosity = porosity_sub.sum();
    //volume = volume_sub.sum();
    an = an_sub.sum();
    aw = aw_sub.sum();
    
    awn = awn_sub.sum();
    ans = ans_sub.sum();
    aws = aws_sub.sum();
    
  
    vol_n = vol_n_sub.sum();
    vol_w = vol_w_sub.sum();
    euler = euler_sub.sum();
  
    pn = pn_sub.sum();
    pw = pw_sub.sum();

    pwn = pwn_sub.sum();
    wp_volume = wp_volume_sub.sum();
    nwp_volume = nwp_volume_sub.sum();

    KGwns_sub.sum();
    KNwns_sub.sum();
   

    lwns = lwns_sub.sum();
 
    
    for ( size_t i=0; i<subdivide[0]; i++) {
        for ( size_t j=0; j<subdivide[1]; j++) {
            for ( size_t k=0; k<subdivide[2]; k++) {
                for ( int m=0; m<3; m++) {
                    vwn(m) += vwn_sub(i,j,k,m);
                    vws(m) += vws_sub(i,j,k,m);
                    vns(m) += vns_sub(i,j,k,m);
                    vwns(m) += vwns_sub(i,j,k,m);
                    vw(m) += vw_sub(i,j,k,m);
                    vn(m) += vn_sub(i,j,k,m);
                }
                for ( int m=0; m<6; m++) {
                    Gwn(m) += Gwn_sub(i,j,k,m);
                    Gws(m) += Gws_sub(i,j,k,m);
                    Gns(m) += Gns_sub(i,j,k,m);
                    Gwns(m) += Gwns_sub(i,j,k,m);
                }
            }
        }
    }
    
    //std::cout << "vwn(2)=" << vwn(2) << std::endl;
 
}





void TwoPhase::AssignComponentLabels()
{
    NumberComponents_NWP = ComputeGlobalBlobIDs(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2,Dm->rank_info,SDs, SDn,solid_isovalue,fluid_isovalue,Label_NWP,Dm->Comm);
}



void TwoPhase::ReduceGlobal()
{

    int amin = Dm->amin;
    int amax = Dm->amax;
    nodeVolume = 1.0*(double(amax-amin)*double(Nx-2)*double(Ny-2));
    MPI_Allreduce(&nodeVolume,&Volume,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    double iVol_global = 1.0/double(Volume);

    // Sum local averages to global
    MPI_Barrier(Dm->Comm);
    MPI_Allreduce(&volume,&volume_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&porosity,&porosity_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    
    MPI_Allreduce(&nsEuler,&nsEuler_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&nwEuler,&nwEuler_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&euler,&euler_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    nGlobalEuler = nwEuler_global;
  
    MPI_Allreduce(&Jwn_tcenter,&Jwn_tcenter_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Kwn_tcenter,&Kwn_tcenter_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&awn_tcenter,&awn_tcenter_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    //std::cout << "awn_tcenter_global=" << awn_tcenter_global << std::endl;
    
    /* Volume Reductions */
    MPI_Allreduce(&nwp_volume,&nwp_volume_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&wp_volume,&wp_volume_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&an,&an_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&aw,&aw_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm); 
    
    MPI_Allreduce(&Eta,&Eta_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Eta1,&Eta_global1,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&ans_tcenter,&ans_tcenter_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&aws,&aws_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&euler,&euler_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&vol_w,&vol_w_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vol_n,&vol_n_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&local_nw_phase_volume,&global_nw_phase_volume,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&mass_local,&mass_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&Jns_tcenter,&Jns_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Kns_tcenter,&Kns_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&wwnx_tcenter,&wwnx_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&wwny_tcenter,&wwny_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&wwnz_tcenter,&wwnz_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    if (awn_tcenter_global == 0.0) {
        Jwn_tcenter_global = 0.0;
        Kwn_tcenter_global = 0.0;
        wwnx_global = 0;
        wwny_global = 0;
        wwnz_global = 0;
    }
    else {
        Jwn_tcenter_global /= awn_tcenter_global;
        Kwn_tcenter_global /= awn_tcenter_global;
        wwnx_global /= awn_tcenter_global;
        wwny_global /= awn_tcenter_global;
        wwnz_global /= awn_tcenter_global;
    }
    if (aws_global == 0.0){
        Jws_global = 0.0;
        Kws_global = 0.0;
    }
    if (ans_tcenter_global == 0.0) {
        Jns_global = 0.0;
        Kns_global = 0.0;
    }
    else {
        Jns_global /= ans_tcenter_global;
        Kns_global /= ans_tcenter_global;
    }
    
    /* Common curve curvatures  */
    MPI_Allreduce(&KGwns,&KGwns_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&KNwns,&KNwns_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    efawns_global_count = 0;
    efawns_global = 0.0;
    /* Common Curve Reduction */
    MPI_Allreduce(&efawns,&efawns_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    //    efawns_global /= double(npx * npy * npz);
    MPI_Allreduce(&efawns_count,&efawns_global_count,1,MPI_INT,MPI_SUM,Dm->Comm);
    
    if (efawns_count == 0) { efawns = 0.0; }
    else { efawns /= double(efawns_count);}
    if (efawns_global_count == 0) { efawns_global = 0.0; }
    else { efawns_global /= double(efawns_global_count);  }

   
    //std::cout << "efawns=" << efawns << " efawns_count=" << efawns_count << std::endl;
    
    /* Pressure */
    MPI_Allreduce(&pw,&pw_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&pn,&pn_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&pwn,&pwn_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    if (vol_n_global > 0) {
        pn_global /= vol_n_global;
    }
    else {
        pn_global = 0.0;
    }
    
    if (vol_w_global > 0) {
        pw_global /= vol_w_global;
    }
    else {
        pw_global = 0.0;
    }
//
    if (awn_tcenter_global > 0 ) {
        pwn_global /= awn_tcenter_global;
    }
    else {
        pwn_global= 0.0;
    }
    
    /* Velocity-Related Reductions */
    MPI_Allreduce(&vn(0),&vn_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vw(0),&vw_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&wwn(0),&wwn_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&vwn(0),&vwn_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vws(0),&vws_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vns(0),&vns_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vwns(0),&vwns_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vwndnw,&vwndnw_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&wwndnw,&wwndnw_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&wwnsdnwn,&wwnsdnwn_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    MPI_Allreduce(&nGaussianCurvature,&nGaussianCurvature_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    
    MPI_Allreduce(&nwOwnO(0),&nwOwnO_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    
    if (wp_volume_global > 0.0){
        vw_global(0) /= wp_volume_global;
        vw_global(1) /= wp_volume_global;
        vw_global(2) /= wp_volume_global;
    }
    else {
        vw_global(0) = 0.0;
        vw_global(1) = 0.0;
        vw_global(2) = 0.0;
        
    }
    if (nwp_volume_global > 0.0){
        vn_global(0) /= nwp_volume_global;
        vn_global(1) /= nwp_volume_global;
        vn_global(2) /= nwp_volume_global;
    }
    else {
        vn_global(0) = 0.0;
        vn_global(1) = 0.0;
        vn_global(2) = 0.0;
    }
    
    if (awn_tcenter_global > 0.0) {
        for (int i=0; i<3; i++) vwn_global(i)/= awn_tcenter_global;
    }
    else {
        for (int i=0; i<3; i++) vwn_global(i) = 0;
    }
    if (awn_tcenter_global > 0) {
        wwndnw_global /= awn_tcenter_global;
    }
    else {
        wwndnw_global = 0;
    }
    
    
    /* Orientation Reductions */
    MPI_Allreduce(&Gwn(0),&Gwn_global(0),6,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Gns(0),&Gns_global(0),6,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Gws(0),&Gws_global(0),6,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Gwns(0),&Gwns_global(0),6,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    if (awn_tcenter_global > 0.0) {
        for (int i=0; i<6; i++) Gwn_global(i)/= awn_tcenter_global;
    }
    else {
        for (int i=0; i<6; i++) Gwn_global(i) = 0.0;
    }
    
    if (ans_tcenter_global > 0.0) {
        for (int i=0; i<3; i++) vns_global(i) /= ans_tcenter_global;
        for (int i=0; i<6; i++) Gns_global(i) /= ans_tcenter_global;
        

    }
    else {
        for (int i=0; i<3; i++) vns_global(i) = 0.0;
        for (int i=0; i<6; i++) Gns_global(i) = 0.0;
    }
    if (aws_global > 0.0) {
        for (int i=0; i<3; i++) vws_global(i) /= aws_global;
        for (int i=0; i<6; i++) Gws_global(i) /= aws_global;
        
       
    }
    else {
        for (int i=0; i<3; i++) vws_global(i) = 0.0;
        for (int i=0; i<6; i++) Gws_global(i) = 0.0;
    }
    

    /* Curvature Reductions */
    MPI_Allreduce(&lwns,&lwns_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    if (lwns_global > 0.0){
        // efawns_global /= lwns_global;
        KNwns_global /= lwns_global;
        KGwns_global /= lwns_global;
        for (int i=0; i<3; i++) vwns_global(i) /= lwns_global;
        for (int i=0; i<6; i++) Gwns_global(i) /= lwns_global;
    }
    else {
        efawns_global = 1.0;
        KNwns_global = 0.0;
        KGwns_global = 0.0;
        for (int i=0; i<3; i++) vwns_global(i) = 0.0;
        for (int i=0; i<6; i++) Gwns_global(i) = 0.0;
    }
    
    MPI_Allreduce(&nwOwsO(0),&nwOwsO_global(0),3,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Barrier(Dm->Comm);
    // <v_wn . n_w >_OwnO
    vwndnw_global = vwndnw_global*iVol_global;

    porosity_global *= iVol_global;
    porosity_global = 1.0 - porosity_global;
    eulerW_global *= iVol_global;
    eulerS_global *= iVol_global;
    
    
    //    awn2_global *=iVol_global;
    /* Do these last */
    awn_tcenter_global = awn_tcenter_global*iVol_global;
    ans_tcenter_global = ans_tcenter_global*iVol_global;
    aws_global = aws_global*iVol_global;
    lwns_global = lwns_global*iVol_global;
    an_global *= iVol_global;
    aw_global *= iVol_global;

}

void TwoPhase::NonDimensionalize(double D, double viscosity, double IFT)
{
    NULL_USE( viscosity );
    NULL_USE( IFT );
    awn_global *= D;
    ans_tcenter_global *= D;
    ans_tcenter_global *= D;
    lwns_global *= D*D;
}

void TwoPhase::PrintAll(int timestep)
{
    if ( TIMELOG_global ){

        Ca_global = abs(vn_global(2) *average_mu)/(5.796*alpha*1.0);
        Re_global = abs(vn_global(2) *1.0 * sphere_diameter)/(average_mu);
        
        if (Ca_global < 1E-15) Ca_global = 0.0;
        if (Re_global < 1E-15) Re_global = 0.0;
        
        double sat_w = wp_volume_global/(nwp_volume_global+wp_volume_global);

        fprintf(TIMELOG_global,"%i %-15.15E %-15.15E ",timestep,Volume,porosity_global);
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",sat_w,pw_global,pn_global);
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",awn_tcenter_global,ans_tcenter_global,aws_global);
        fprintf(TIMELOG_global,"%-15.15E %-15.15E ",Jwn_tcenter_global, Kwn_tcenter_global);
        fprintf(TIMELOG_global,"%-15.15E ",lwns_global);
        fprintf(TIMELOG_global,"%.1f ",efawns_global);
        fprintf(TIMELOG_global,"%-15.15E %-15.15E ",KNwns_global, KGwns_global);
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",vw_global(0),vw_global(1),vw_global(2));
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",vn_global(0),vn_global(1),vn_global(2));
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",vwn_global(0),vwn_global(1),vwn_global(2));
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",vwns_global(0),vwns_global(1),vwns_global(2));
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gwn_global(0),Gwn_global(1),Gwn_global(2),Gwn_global(3),Gwn_global(4),Gwn_global(5));
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gws_global(0),Gws_global(1),Gws_global(2),Gws_global(3),Gws_global(4),Gws_global(5));
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gns_global(0),Gns_global(1),Gns_global(2),Gns_global(3),Gns_global(4),Gns_global(5));
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E ", pwn_global, euler_global, aw_global, an_global);
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gwns_global(0),Gwns_global(1),Gwns_global(2),Gwns_global(3),Gwns_global(4),Gwns_global(5));
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                vws_global(0),vws_global(1),vws_global(2),vns_global(0),vns_global(1),vns_global(2));
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E ",
                wwnx_global,wwny_global,wwnz_global);
        fprintf(TIMELOG_global,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E\n", Ca_global, Re_global, Jns_global, Kns_global, 
                pw_tminus_global, pn_tminus_global, aw_tminus_global, an_tminus_global, awn_tminus_global, Jwn_tminus_global, Gwn_tminus_global(0));
        fflush(TIMELOG_global);
        
        
//        if (timestep > 1800) {
            printf("timestep=%d resulting contact angle=%f\n",timestep,efawns_global);
            string filename("/Users/cpf/Desktop/build_LBPM_UNC_042522/Testing/Results.txt");
            fstream file;

            file.open(filename, std::ios_base::app | std::ios_base::in);
            if (file.is_open() && offset_distance != 28 ) {
                file << efawns_global << ", ";
            }
            if (file.is_open() && offset_distance == 28) {
                file << efawns_global << endl;
            }
//        }
    }
    if ( TIMELOG_local ){

        double iVol_nodal;
        if (nodeVolume == 0.0){
            iVol_nodal = 0.0;
        } else{
            iVol_nodal = 1.0/nodeVolume;
        }
        
        //porosity = 0.0;
        if (iVol_nodal > 0.0){
            porosity *= iVol_nodal;
            porosity = 1.0 - porosity;
        }
      
        if (vol_n > 0) {
            pn /= vol_n;
        }
        else {
            pn = 0.0;
        }
        
        if (vol_w > 0) {
            pw /= vol_w;
        }
        else {
            pw = 0.0;
        }
        
        
        if (awn_tcenter > 0.0) {
            pwn /= awn;
        }
        else {
            Jwn = 0;
            Kwn = 0;
            pwn = 0;
        }
        
        
        if (lwns > 0.0) {
            KNwns /= lwns;
            KGwns /= lwns;
            // efawns /= lwns;
        }
        else {
            KNwns = 0;
            KGwns = 0;
            // efawns = 1.0;
        }
        
        
        for (int m = 0; m < 3; m++) {
            if (wp_volume > 0.0) {vw(m) /= wp_volume;   } else { vw(m) = 0; }
            if (nwp_volume > 0.0) {vn(m) /= nwp_volume; } else { vn(m) = 0; }
            if (awn_tcenter > 0.0) {vwn(m) /= awn_tcenter;} else { vwn(m) = 0; }
            if (aws > 0.0) {vws(m) /= aws;} else { vws(m) = 0; }
            if (ans > 0.0) {vns(m) /= ans;} else { vns(m) = 0; }
            if (lwns > 0.0) {vwns(m) /= lwns;} else { vwns(m) = 0; }
        }
        
        
        if (ans > 0.0) { Kns_tcenter /= ans;} else { Kns_tcenter = 0; }
        if (ans > 0.0) { Jns_tcenter /= ans;  } else { Jns_tcenter = 0; }
        
        for (int m = 0; m < 6; m++ ) {
            if (awn_tcenter > 0.0) {Gwn(m) /= awn_tcenter; } else { Gwn(m) = 0; }
            if (ans > 0.0) {Gns(m) /= ans; } else { Gns(m) = 0; }
            if (aws > 0.0) {Gws(m) /= aws; } else { Gws(m) = 0; }
            if (lwns > 0.0){Gwns(m) /= lwns; } else { Gwns(m) = 0; }
        }
        
        //        wwn(0) *= iVol_nodal;
        if (awn_tcenter > 0) { wwn(0) /= awn_tcenter;  }  else{ wwn(0) = 0;}
        //        wwn(1) *= iVol_nodal;
        if (awn_tcenter > 0) { wwn(1) /= awn_tcenter; }  else{ wwn(1) = 0;}
        //        wwn(2) *= iVol_nodal;
        if (awn_tcenter > 0) { wwn(2) /= awn_tcenter; }  else{ wwn(2) = 0;}
        
        // nEuler = nsEuler + nwEuler;
        // nEuler /= (4*3.14159265359);
        
        // nEuler *= iVol_nodal;
        
        //nGaussianCurvature /= (2*PI);
        
        if (iVol_nodal > 0.0){
            aw *= iVol_nodal;
            aw_tminus *= iVol_nodal;

            an *= iVol_nodal;
            an_tminus *= iVol_nodal;

            awn_tcenter *= iVol_nodal;
            awn_tminus *= iVol_nodal;

            aws *= iVol_nodal;

            ans *= iVol_nodal;
            ans_tcenter *= iVol_nodal;

            lwns *= iVol_nodal;
        }

        //nw(0) = nw(1) = nw(2) = 0;
        
        
        double Ca = fabs(vn(2) *average_mu)/(5.796*alpha*1.0);
        double Re = fabs(vn(2) *1.0 * sphere_diameter)/(average_mu);
        
        double sat_w = 0.0;
        if (nwp_volume+wp_volume > 0.0){
            sat_w = 1.0 - nwp_volume/(nwp_volume+wp_volume);
        }
       
        fprintf(TIMELOG_local,"%i %-15.15E %-15.15E ",timestep, nodeVolume, porosity);
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",sat_w,pw,pn);
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",awn_tcenter,ans_tcenter,aws);
        fprintf(TIMELOG_local,"%-15.15E %-15.15E ",Jwn_tcenter, Kwn_tcenter);
        fprintf(TIMELOG_local,"%-15.15E ",lwns);
        fprintf(TIMELOG_local,"%-15.15E ",efawns);
        fprintf(TIMELOG_local,"%-15.15E %-15.15E ",KNwns, KGwns);
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",vw(0),vw(1),vw(2));
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",vn(0),vn(1),vn(2));
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",vwn(0),vwn(1),vwn(2));
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",vwns(0),vwns(1),vwns(2));
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gwn(0),Gwn(1),Gwn(2),Gwn(3),Gwn(4),Gwn(5));
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gws(0),Gws(1),Gws(2),Gws(3),Gws(4),Gws(5));
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gns(0),Gns(1),Gns(2),Gns(3),Gns(4),Gns(5));
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E ", pwn, euler, aw, an);
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                Gwns(0),Gwns(1),Gwns(2),Gwns(3),Gwns(4),Gwns(5));
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E ",
                vws(0),vws(1),vws(2),vns(0),vns(1),vns(2));
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E ",
                wwnx_tcenter,wwny_tcenter,wwnz_tcenter);
        fprintf(TIMELOG_local,"%-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E %-15.15E\n",Ca, Re, Jns_tcenter, Kns_tcenter, pw_rhs, pw_lhs, pn_rhs, pn_lhs,
                aw_tminus,an_tminus,awn_tminus);
        fflush(TIMELOG_local);
    }

}

void TwoPhase::AggregateLabels( const std::string& filename )
{

	int nx = Dm->Nx;
	int ny = Dm->Ny;
	int nz = Dm->Nz;
	
	// assign the ID from the phase indicator field
	for (int k=0; k<nz; k++){
		for (int j=0; j<ny; j++){
			for (int i=0; i<nx; i++){
				size_t n = k*nx*ny+j*nx+i;
				char local_id_val = Dm->id[n]; 
				if (local_id_val > 0){
					double value = SDn(i,j,k);
					if (value > 0.0)	local_id_val = 1;
					else 				local_id_val = 2;
				}
				Dm->id[n] = local_id_val;
			}
		}
	}

	MPI_Barrier(Dm->Comm);

	Dm->AggregateLabels( filename );

}



void TwoPhase::ComputeEwn()
{
    nwRootMesh m_nw;
    nwRootMesh  m_n;
    nwRootMesh m_ns;
    wsMyMesh m_ws;
    nwMyMesh m_ni1;
    nwsMyMesh m_nws_2;
    
    Point P,A,B,C;
    Point Q,D,E,F;
    Point R,G,H,I;

    Dm->CommunicateMeshHalo(SDs);
    
    pmmc_MeshGradient(SDs,SDs_x,SDs_y,SDs_z,Nx,Ny,Nz);
    
    // Gradient of the phase indicator field
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_z);

    
    Dm->CommunicateMeshHalo(SDn);
    
    pmmc_MeshGradient( SDn, SDn_x, SDn_y, SDn_z,Nx,Ny,Nz);
    //...........................................................................
    // Gradient of the phase indicator field
    //...........................................................................
    Dm->CommunicateMeshHalo( SDn_x);
    //...........................................................................
    Dm->CommunicateMeshHalo( SDn_y);
    //...........................................................................
    Dm->CommunicateMeshHalo( SDn_z);
    //...........................................................................

    int amin = Dm->amin;
    int amax = Dm->amax;
    nodeVolume = 1.0*(double(amax-amin)*double(Nx-2)*double(Ny-2));
    MPI_Allreduce(&nodeVolume,&Volume,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    const int N2x = ( Nx-2 ) / subdivide[0];
    const int N2y = ( Ny-2 ) / subdivide[1];
    const int N2z = ( Nz-2 ) / subdivide[2];

    //subNodeVolume = 1.0*(N2x*N2y*N2z);
    double iVol_global = 1.0/double(Volume);

    DoubleArray Gwn_tmp(6);
        
    for (int k=amin; k<amax; k++) {
        const int k2 = std::min<int>( ( k - 1) / N2z, subdivide[2] - 1 );
        for (int j=1; (size_t)j<Ny-1; j++){
            const int j2 = std::min<int>( ( j - 1 ) / N2y, subdivide[1] - 1 );
            for (int i=1; (size_t)i<Nx-1; i++) {
                const int i2 = std::min<int>( ( i - 1 ) / N2x, subdivide[0] - 1 );

                Gwn_tmp.fill( 0 );
                
                // Compute volume averages
                int n = i + j*Nx + k*Nx*Ny;
                double vs = VFmask(n);
                double phase_val =  SDn(n);
                double delPhi_val = DelPhi(n);

                //volume_sub(i2,j2,k2) += 1.0;

                if (vs < 0.5) {
                    if (phase_val <= 0) {

                        //wp_volume_sub(i2,j2,k2) += 1.0;
                        aw_sub(i2,j2,k2) += 1.0;

                        if (delPhi_val < 1e-4) {  
                            pw_sub(i2,j2,k2) += Press(n); 
                            vol_w_sub(i2,j2,k2) += 1;  
                        }

                    }
                }
             
                if (vs < 0.5) {
                    if (phase_val > 0) {

                        //nwp_volume_sub(i2,j2,k2) += 1.0;
                        an_sub(i2,j2,k2) += 1.0;

                        if (delPhi_val < 1e-4) { 
                            pn_sub(i2,j2,k2) += Press(n);   
                            vol_n_sub(i2,j2,k2) += 1; 
                        }

                    }
                }

                //...........................................................................
                n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
                n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;
                
                //...........................................................................
                pmmc_ConstructLocalCube(SDs,  SDn,  SDn_x,  SDn_y,  SDn_z, solid_isovalue, fluid_isovalue, nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
                                        local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                                        n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                                        n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
                                        i, j, k, Nx, Ny, Nz);
                
                if (n_nw_tris > 0) {
                    for (int r=0;r<n_nw_tris;r++) {
                        A = nw_pts(nw_tris(0,r));
                        B = nw_pts(nw_tris(1,r));
                        C = nw_pts(nw_tris(2,r));

                        double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                        double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                        double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                        double normal_x =  -SDn_x(i,j,k);
                        double normal_y =  -SDn_y(i,j,k);
                        double normal_z =  -SDn_z(i,j,k);

                        if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                            P = A;
                            A = C;
                            C = P;
                        }
                        
                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_nw,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                        
                    }
                }
                
                if (n_nws_seg > 0) {
                    for (int r = 0; r < n_nws_seg; r++) {

                        nwsMyMesh m_nws;
                        nwsMyMesh m_nws_edges;
                        Q = nws_pts(nws_seg(0,r));
                        R = nws_pts(nws_seg(1,r));

                        nwMyMesh m_nw;

                        if (n_nw_tris > 0 ) {

                            for (int r=0;r<n_nw_tris;r++) {
                                A = nw_pts(nw_tris(0,r));
                                B = nw_pts(nw_tris(1,r));
                                C = nw_pts(nw_tris(2,r));

                                double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                                double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                                double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                                double normal_x =  -SDn_x(i,j,k);
                                double normal_y =  -SDn_y(i,j,k);
                                double normal_z =  -SDn_z(i,j,k);

                                if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                                    P = A;
                                    A = C;
                                    C = P;
                                }
                                if (Q == A) {
                                    if (R == B) {

                                       vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));

                                    }
                                    if (R == C) {

                                       vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));

                                    }
                                }
                                if (Q == B) {
                                    if (R == A) {

                                      vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));

                                    }
                                    if (R == C) {

                                        vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));

                                    }
                                }
                                if (Q == C) {
                                    if (R == B) {

                                         vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));

                                    }
                                    if (R == A) {

                                       vcg::tri::Allocator<nwMyMesh>::AddFace(m_nw,nwMyMesh::CoordType ( A.x, A.y, A.z),nwMyMesh::CoordType ( B.x, B.y, B.z),nwMyMesh::CoordType ( C.x, C.y, C.z));

                                    }
                                }
                            }
                        }
                    }
                }

                if (n_nw_tris > 0) {
                    pmmc_CubeSurfaceOrientation(Gwn_tmp,nw_pts,nw_tris,n_nw_tris);
                    for (int m=0; m<6; m++) {  Gwn_sub(i2,j2,k2,m) += Gwn_tmp(m); }
                }

            }
        }
    }

    pn_tminus = 0;
    pw_tminus = 0;
    vol_n_tminus = 0;
    vol_w_tminus = 0;

    an_tminus = 0;
    aw_tminus = 0;
    awn_tcenter = 0;
    awn_tminus = 0;

    Jwn_tcenter = 0;
    Jwn_tminus = 0;
    Kwn_tcenter = 0;
    vwndnn_tcenter = 0;

    ComputeAccurateNWInterfaceCurvatures(m_nw,&Jwn_tcenter,&Kwn_tcenter,&awn_tcenter,&vwndnn_tcenter, ewnwwn, wwnJn, vwndnnIIprime, Vel_x,Vel_y,Vel_z);

    awn_tminus = awn_tcenter;
    Jwn_tminus = Jwn_tcenter;
    awn_tcenter = Jwn_tcenter = 0;

    pn_tminus = pn_sub.sum();
    pw_tminus = pw_sub.sum();
    pn_sub.fill(0);
    pw_sub.fill(0);

    vol_n_tminus = vol_n_sub.sum();
    vol_w_tminus = vol_w_sub.sum();
    vol_n_sub.fill(0);
    vol_w_sub.fill(0);

    an_tminus = an_sub.sum();
    aw_tminus = aw_sub.sum();
    an_sub.fill(0);
    aw_sub.fill(0);

    for ( size_t i=0; i<subdivide[0]; i++) {
        for ( size_t j=0; j<subdivide[1]; j++) {
            for ( size_t k=0; k<subdivide[2]; k++) {
                for ( int m=0; m<6; m++) {
                    Gwn_tminus(m) += Gwn_sub(i,j,k,m);
                }
            }
        }
    }
    Gwn_sub.fill(0);

    //Find global sums
    MPI_Allreduce(&pw_tminus,&pw_tminus_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&pn_tminus,&pn_tminus_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    MPI_Allreduce(&vol_w_tminus,&vol_w_tminus_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&vol_n_tminus,&vol_n_tminus_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    MPI_Allreduce(&an_tminus,&an_tminus_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&aw_tminus,&aw_tminus_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&awn_tminus,&awn_tminus_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    MPI_Allreduce(&Jwn_tminus,&Jwn_tminus_global,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
    MPI_Allreduce(&Gwn_tminus(0),&Gwn_tminus_global(0),6,MPI_DOUBLE,MPI_SUM,Dm->Comm);

    // //Compute local averages
    // if (vol_n_tminus == 0.0) {
    //     pn_tminus = 0.0;
    // } else {
    //     pn_tminus /= vol_n_tminus;
    // }
    
    // if (vol_w_tminus == 0.0) {
    //     pw_tminus = 0.0;
    // } else {
    //     pw_tminus /= vol_w_tminus;
    // }

    //Compute global averages
    if (vol_w_tminus_global == 0.0){
        pw_tminus_global = 0.0;
    } else {
        pw_tminus_global /= vol_w_tminus_global;
    }

    if (vol_n_tminus_global == 0.0){
        pn_tminus_global = 0.0;
    } else {
        pn_tminus_global /= vol_n_tminus_global;
    }

    if (awn_tminus_global == 0.0) {
        Jwn_tminus_global = 0.0;
        for (int i=0; i<6; i++) Gwn_tminus_global(i) = 0.0;
    } else {
        Jwn_tminus_global /= awn_tminus_global;
        for (int i=0; i<6; i++) Gwn_tminus_global(i)/= awn_tminus_global;
    }

    awn_tminus_global *= iVol_global;
    an_tminus_global *= iVol_global;
    aw_tminus_global *= iVol_global;

}

void TwoPhase::ComponentAverages()
{

    int i,j,k;
    // int kmin,kmax;
	//int LabelWP,LabelNWP;
    int LabelNWP;
	//double TempLocal;

    nwRootMesh m_nw;
    nwRootMesh  m_n;
    nwRootMesh m_ns;
    wsMyMesh m_ws;
    nwMyMesh m_ni1;
    nwsMyMesh m_nws_2;
    
    Point P,A,B,C;
    Point Q,D,E,F;
    Point R,G,H,I;

    DoubleArray vwn_tmp(3), vns_tmp(3), Gwn_tmp(6), Gns_tmp(6);

    //const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    
    // ComponentAverages_WP.resize(BLOB_AVG_COUNT,NumberComponents_WP);
    ComponentAverages_NWP.resize(BLOB_AVG_COUNT,NumberComponents_NWP);
    
    // ComponentAverages_WP.fill(0.0);
    ComponentAverages_NWP.fill(0.0);
    
    //if (Dm->rank()==0){
    //    printf("Number of wetting phase components is %i \n",NumberComponents_WP);
    //    printf("Number of non-wetting phase components is %i \n",NumberComponents_NWP);
    //}

    Dm->CommunicateMeshHalo(SDn);
    
    pmmc_MeshGradient( SDn, SDn_x, SDn_y, SDn_z,Nx,Ny,Nz);
    //...........................................................................
    // Gradient of the phase indicator field
    //...........................................................................
    Dm->CommunicateMeshHalo( SDn_x);
    //...........................................................................
    Dm->CommunicateMeshHalo( SDn_y);
    //...........................................................................
    Dm->CommunicateMeshHalo( SDn_z);
    //...........................................................................

    int amin = Dm->amin;
    int amax = Dm->amax;

	for (k=amin; k<amax; k++){
		for (j=1; (size_t)j<Ny-1; j++){
			for (i=1; (size_t)i<Nx-1; i++){
				//LabelWP=GetCubeLabel(i,j,k,Label_WP);
				LabelNWP=GetCubeLabel(i,j,k,Label_NWP);
				n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
				n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;
				// Initialize the averaged quantities
				awn = aws = ans = lwns = 0.0;
				vwn(0) = vwn(1) = vwn(2) = 0.0;
				vwns(0) = vwns(1) = vwns(2) = 0.0;
				Gwn(0) = Gwn(1) = Gwn(2) = 0.0;
				Gwn(3) = Gwn(4) = Gwn(5) = 0.0;
				Gws(0) = Gws(1) = Gws(2) = 0.0;
				Gws(3) = Gws(4) = Gws(5) = 0.0;
				Gns(0) = Gns(1) = Gns(2) = 0.0;
				Gns(3) = Gns(4) = Gns(5) = 0.0;
				Jwn = Kwn = efawns = 0.0;
				//trawn=trJwn=0.0;
				euler=0.0;

                int n = i + j*Nx + k*Nx*Ny;
                double vs = VFmask(n);

                double phase_val =  SDn(n);
                double delPhi_val = DelPhi(n);
                //char id_val = Dm->id[n];


                if (vs < 0.5) {
                    if (phase_val > 0) {
                        ComponentAverages_NWP(VX,LabelNWP) += Vel_x(n);
                        ComponentAverages_NWP(VY,LabelNWP) += Vel_y(n);
                        ComponentAverages_NWP(VZ,LabelNWP) += Vel_z(n);

                        ComponentAverages_NWP(VOL,LabelNWP)  += 1.0;

                        if (delPhi_val < 1e-4) {  ComponentAverages_NWP(PRS,LabelNWP) += Press(n);  ComponentAverages_NWP(TRIMVOL,LabelNWP) += 1; }
                        
                    }
                }

                //...........................................................................
                n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
                n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;
                
                //...........................................................................
                pmmc_ConstructLocalCube(SDs,  SDn,  SDn_x,  SDn_y,  SDn_z, solid_isovalue, fluid_isovalue, nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
                                        local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                                        n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                                        n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
                                        i, j, k, Nx, Ny, Nz);

                if (n_nw_tris > 0) {
                    for (int r=0;r<n_nw_tris;r++) {
                        A = nw_pts(nw_tris(0,r));
                        B = nw_pts(nw_tris(1,r));
                        C = nw_pts(nw_tris(2,r));

                        double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                        double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                        double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

                        double normal_x =  -SDn_x(i,j,k);
                        double normal_y =  -SDn_y(i,j,k);
                        double normal_z =  -SDn_z(i,j,k);

                        if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                            P = A;
                            A = C;
                            C = P;
                        }
                        
                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_nw,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                        
//                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                    }
                }

                if (n_ns_tris > 0) {
                    for (int r=0;r<n_ns_tris;r++) {
                        A = ns_pts(ns_tris(0,r));
                        B = ns_pts(ns_tris(1,r));
                        C = ns_pts(ns_tris(2,r));
                        
                        double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                        double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                        double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);
                        
                        double normal_x = -SDs_x(i,j,k);
                        double normal_y = -SDs_y(i,j,k);
                        double normal_z = -SDs_z(i,j,k);
                        
                        if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                            P = A;
                            A = C;
                            C = P;
                        }
                        

                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_ns,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                        
//                        vcg::tri::Allocator<nwRootMesh>::AddFace(m_n,nwRootMesh::CoordType ( A.x, A.y, A.z),nwRootMesh::CoordType ( B.x, B.y, B.z),nwRootMesh::CoordType ( C.x, C.y, C.z));
                    }
                }          

                if (n_ws_tris > 0) {
                    for (int r=0;r<n_ws_tris;r++) {
                        A = ws_pts(ws_tris(0,r));
                        B = ws_pts(ws_tris(1,r));
                        C = ws_pts(ws_tris(2,r));
                        
                        double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
                        double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
                        double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);
                        
                        double normal_x = SDs_x(i,j,k);
                        double normal_y = SDs_y(i,j,k);
                        double normal_z = SDs_z(i,j,k);
                        
                        if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
                            P = A;
                            A = C;
                            C = P;
                        }
                        

                        vcg::tri::Allocator<wsMyMesh>::AddFace(m_ws,wsMyMesh::CoordType ( A.x, A.y, A.z),wsMyMesh::CoordType ( B.x, B.y, B.z),wsMyMesh::CoordType ( C.x, C.y, C.z));
                    }
                } 

                // if (n_ws_tris > 0) {
                //   //  std::cout << "n_ws_tris > 0" << std::endl;
                //     ComponentAverages_NWP(AWS,LabelNWP) += pmmc_CubeSurfaceOrientation(Gws_tmp,ws_pts,ws_tris,n_ws_tris);
                //     //for (int m=0; m<6; m++) { Gws_sub(i2,j2,k2,m) += Gws_tmp(m); }
                //     // pmmc_InterfaceSpeed(dPdt,  SDn_x,  SDn_y,  SDn_z, CubeValues, ws_pts, ws_tris, NormalVector, InterfaceSpeed, vws_tmp, i, j, k, n_ws_pts, n_ws_tris);
                //     // ComponentAverages_NWP(VWSX,LabelNWP) += vws_tmp(0);
                //     // ComponentAverages_NWP(VWSY,LabelNWP) += vws_tmp(1);
                //     // ComponentAverages_NWP(VWSZ,LabelNWP) += vws_tmp(2);
                   
                // }
                if (n_ns_tris > 0) {
                  //  std::cout << "n_ns_tris > 0" << std::endl;
                    ComponentAverages_NWP(ANS,LabelNWP) += pmmc_CubeSurfaceOrientation(Gns_tmp,ns_pts,ns_tris,n_ns_tris);
                    //for (int m=0; m<6; m++) {  Gns_sub(i2,j2,k2,m) += Gns_tmp(m); }
                    pmmc_InterfaceSpeed(dPdt,  SDn_x,  SDn_y,  SDn_z,  CubeValues, ns_pts,ns_tris, NormalVector, InterfaceSpeed, vns_tmp, i, j, k, n_ns_pts, n_ns_tris);
                    ComponentAverages_NWP(VNSX,LabelNWP) += vns_tmp(0);
                    ComponentAverages_NWP(VNSY,LabelNWP) += vns_tmp(1);
                    ComponentAverages_NWP(VNSZ,LabelNWP) += vns_tmp(2);
                    //euler_sub(i2,j2,k2) += geomavg_EulerCharacteristic(ns_pts,ns_tris,n_ns_pts,n_ns_tris,i,j,k);
                }
                
                if (n_nw_tris > 0) {
                   ComponentAverages_NWP(AWN,LabelNWP) += pmmc_CubeSurfaceOrientation(Gwn_tmp,nw_pts,nw_tris,n_nw_tris);
                    //for (int m=0; m<6; m++) {  Gwn_sub(i2,j2,k2,m) += Gwn_tmp(m); }
                    pmmc_InterfaceSpeed(dPdt,  SDn_x,  SDn_y,  SDn_z,  CubeValues, nw_pts,nw_tris, NormalVector, InterfaceSpeed, vwn_tmp, i, j, k, n_nw_pts, n_nw_tris);
                    ComponentAverages_NWP(VWNX,LabelNWP) += vwn_tmp(0);
                    ComponentAverages_NWP(VWNY,LabelNWP) += vwn_tmp(1);
                    ComponentAverages_NWP(VWNZ,LabelNWP) += vwn_tmp(2);
                    //euler_sub(i2,j2,k2) += geomavg_EulerCharacteristic(nw_pts,nw_tris,n_nw_pts,n_nw_tris,i,j,k);
                    //pwn_sub(i2,j2,k2) += pmmc_CubeSurfaceInterpValue(CubeValues,Press,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);
                }
			}
		}
	}

	MPI_Barrier(Dm->Comm);
	/*if (Dm->rank()==0){
		printf("Component averages computed locally -- reducing result... \n");
	}*/
	// Globally reduce the non-wetting phase averages
	RecvBuffer.resize(BLOB_AVG_COUNT,NumberComponents_NWP);

/*	for (int b=0; b<NumberComponents_NWP; b++){
		MPI_Barrier(Dm->Comm);
		MPI_Allreduce(&ComponentAverages_NWP(0,b),&RecvBuffer(0),BLOB_AVG_COUNT,MPI_DOUBLE,MPI_SUM,Dm->Comm);
		for (int idx=0; idx<BLOB_AVG_COUNT; idx++) ComponentAverages_NWP(idx,b)=RecvBuffer(idx);
	}
	*/
	MPI_Barrier(Dm->Comm);
	MPI_Allreduce(ComponentAverages_NWP.data(),RecvBuffer.data(),BLOB_AVG_COUNT*NumberComponents_NWP,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	//	MPI_Reduce(ComponentAverages_NWP.data(),RecvBuffer.data(),BLOB_AVG_COUNT,MPI_DOUBLE,MPI_SUM,0,Dm->Comm);
  
	/*if (Dm->rank()==0){
		printf("rescaling... \n");
	}*/

	for (int b=0; b<NumberComponents_NWP; b++){
		for (int idx=0; idx<BLOB_AVG_COUNT; idx++) ComponentAverages_NWP(idx,b)=RecvBuffer(idx,b);
	}

	for (int b=0; b<NumberComponents_NWP; b++){
        
		if (ComponentAverages_NWP(VOL,b) > 0.0){
			//double Vn,pn,awn,ans,Jwn,Kwn,lwns,cwns,vsq;
            double Vn,pn,awn,ans;

			Vn = ComponentAverages_NWP(VOL,b);
			awn = ComponentAverages_NWP(AWN,b);
			ans = ComponentAverages_NWP(ANS,b);
			vn(0) = ComponentAverages_NWP(VX,b)/Vn;
			vn(1) = ComponentAverages_NWP(VY,b)/Vn;
			vn(2) = ComponentAverages_NWP(VZ,b)/Vn;
            // NULL_USE(ans);

			if (ComponentAverages_NWP(TRIMVOL,b) > 0.0){
				pn = ComponentAverages_NWP(PRS,b)/ComponentAverages_NWP(TRIMVOL,b);
			}
			else pn = 0.0;

			if (awn != 0.0){
				// Jwn = ComponentAverages_NWP(JWN,b)/awn;
				// Kwn = ComponentAverages_NWP(KWN,b)/awn;
				vwn(0) = ComponentAverages_NWP(VWNX,b)/awn;
				vwn(1) = ComponentAverages_NWP(VWNY,b)/awn;
				vwn(2) = ComponentAverages_NWP(VWNZ,b)/awn;
				// Gwn(0) = ComponentAverages_NWP(GWNXX,b)/awn;
				// Gwn(1) = ComponentAverages_NWP(GWNYY,b)/awn;
				// Gwn(2) = ComponentAverages_NWP(GWNZZ,b)/awn;
				// Gwn(3) = ComponentAverages_NWP(GWNXY,b)/awn;
				// Gwn(4) = ComponentAverages_NWP(GWNXZ,b)/awn;
				// Gwn(5) = ComponentAverages_NWP(GWNYZ,b)/awn;
			}
			// else Jwn=Kwn=0.0;

            if (ans != 0.0){
                vns(0) = ComponentAverages_NWP(VNSX,b)/ans;
				vns(1) = ComponentAverages_NWP(VNSY,b)/ans;
				vns(2) = ComponentAverages_NWP(VNSZ,b)/ans;
                // Gns(0) = ComponentAverages_NWP(GXX,b)/ans;
				// Gns(1) = ComponentAverages_NWP(GYY,b)/ans;
				// Gns(2) = ComponentAverages_NWP(GZZ,b)/ans;
				// Gns(3) = ComponentAverages_NWP(GXY,b)/ans;
				// Gns(4) = ComponentAverages_NWP(GXZ,b)/ans;
				// Gns(5) = ComponentAverages_NWP(GYZ,b)/ans;
            }
            // else Gns(0) = Gns(1) = Gns(2) = Gns(3) = Gns(4) = Gns(5) = 0.0;

			// lwns = ComponentAverages_NWP(LWNS,b);
			// if (lwns != 0.0){
			// 	cwns = ComponentAverages_NWP(CWNS,b)/lwns;
			// 	vwns(0) = ComponentAverages_NWP(VWNSX,b)/lwns;
			// 	vwns(1) = ComponentAverages_NWP(VWNSY,b)/lwns;
			// 	vwns(2) = ComponentAverages_NWP(VWNSZ,b)/lwns;
            //     KNwns = ComponentAverages_NWP(KNWNS,b)/lwns;
            //     KGwns = ComponentAverages_NWP(KGWNS,b)/lwns;
			// }
			// else  {
            //     cwns=0.0; 
            //     KNwns = 0.0; 
            //     KGwns = 0.0;
            // }

			ComponentAverages_NWP(PRS,b) = pn;
			ComponentAverages_NWP(VX,b) = vn(0);
			ComponentAverages_NWP(VY,b) = vn(1);
			ComponentAverages_NWP(VZ,b) = vn(2);
			// ComponentAverages_NWP(VSQ,b) = vsq;

			// ComponentAverages_NWP(JWN,b) = Jwn;
			// ComponentAverages_NWP(KWN,b) = Kwn;
			ComponentAverages_NWP(VWNX,b) = vwn(0);
			ComponentAverages_NWP(VWNY,b) = vwn(1);
			ComponentAverages_NWP(VWNZ,b) = vwn(2);

            ComponentAverages_NWP(VNSX,b) = vns(0);
			ComponentAverages_NWP(VNSY,b) = vns(1);
			ComponentAverages_NWP(VNSZ,b) = vns(2);
			
			// ComponentAverages_NWP(GWNXX,b) = Gwn(0);
			// ComponentAverages_NWP(GWNYY,b) = Gwn(1);
			// ComponentAverages_NWP(GWNZZ,b) = Gwn(2);
			// ComponentAverages_NWP(GWNXY,b) = Gwn(3);
			// ComponentAverages_NWP(GWNXZ,b) = Gwn(4);
			// ComponentAverages_NWP(GWNYZ,b) = Gwn(5);

            // ComponentAverages_NWP(GXX,b) = Gns(0);
			// ComponentAverages_NWP(GYY,b) = Gns(1);
			// ComponentAverages_NWP(GZZ,b) = Gns(2);
			// ComponentAverages_NWP(GXY,b) = Gns(3);
			// ComponentAverages_NWP(GXZ,b) = Gns(4);
			// ComponentAverages_NWP(GYZ,b) = Gns(5);

			// ComponentAverages_NWP(CWNS,b) = cwns;
			// ComponentAverages_NWP(VWNSX,b) = vwns(0);
			// ComponentAverages_NWP(VWNSY,b) = vwns(1);
			// ComponentAverages_NWP(VWNSZ,b) = vwns(2);
            // ComponentAverages_NWP(KNWNS,b) = KNwns;
            // ComponentAverages_NWP(KGWNS,b) = KGwns;

			// ComponentAverages_NWP(CMX,b) /= Vn;
			// ComponentAverages_NWP(CMY,b) /= Vn;
			// ComponentAverages_NWP(CMZ,b) /= Vn;

			// ComponentAverages_NWP(EULER,b) /= (2*PI);

		}

	}
}

inline int TwoPhase::GetCubeLabel(int i, int j, int k, IntArray &BlobLabel)
{
    int label;
    label=BlobLabel(i,j,k);
    label=max(label,BlobLabel(i+1,j,k));
    label=max(label,BlobLabel(i,j+1,k));
    label=max(label,BlobLabel(i+1,j+1,k));
    label=max(label,BlobLabel(i,j,k+1));
    label=max(label,BlobLabel(i+1,j,k+1));
    label=max(label,BlobLabel(i,j+1,k+1));
    label=max(label,BlobLabel(i+1,j+1,k+1));
    
    return label;
}

void TwoPhase::SortBlobs()
{
    //printf("Sorting the blobs based on volume \n");
    //printf("-----------------------------------------------\n");
    int TempLabel,a,aa,bb,i,j,k,idx;
    double TempValue;
    //.......................................................................
    // Sort NWP components by volume
    //.......................................................................
    IntArray OldLabel(NumberComponents_NWP);
    for (a=0; a<NumberComponents_NWP; a++)	OldLabel(a) = a;
    // Sort the blob averages based on volume
    for (aa=0; aa<NumberComponents_NWP-1; aa++){
        for ( bb=aa+1; bb<NumberComponents_NWP; bb++){
            if (ComponentAverages_NWP(VOL,aa) < ComponentAverages_NWP(VOL,bb)){
                // Exchange location of blobs aa and bb
                //printf("Switch blob %i with %i \n", OldLabel(aa),OldLabel(bb));
                // switch the label
                TempLabel = OldLabel(bb);
                OldLabel(bb) = OldLabel(aa);
                OldLabel(aa) = TempLabel;
                // switch the averages
                for (idx=0; idx<BLOB_AVG_COUNT; idx++){
                    TempValue = ComponentAverages_NWP(idx,bb);
                    ComponentAverages_NWP(idx,bb) = ComponentAverages_NWP(idx,aa);
                    ComponentAverages_NWP(idx,aa) = TempValue;
                }
            }
        }
    }
    IntArray NewLabel(NumberComponents_NWP);
    for (aa=0; aa<NumberComponents_NWP; aa++){
        // Match the new label for original blob aa
        bb=0;
        while (OldLabel(bb) != aa)	bb++;
        NewLabel(aa) = bb;
    }
    // Re-label the blob ID
    //	printf("Re-labeling the blobs, now indexed by volume \n");
    for (k=0; (size_t)k<Nz; k++){
        for (j=0; (size_t)j<Ny; j++){
            for (i=0; (size_t)i<Nx; i++){
                if (Label_NWP(i,j,k) > -1){
                    TempLabel = NewLabel(Label_NWP(i,j,k));
                    Label_NWP(i,j,k) = TempLabel;
                }
            }
        }
    }
    //.......................................................................
    // Sort WP components by volume
    //.......................................................................
    // OldLabel.resize(NumberComponents_WP);
    // for (a=0; a<NumberComponents_WP; a++)	OldLabel(a) = a;
    // // Sort the blob averages based on volume
    // for (aa=0; aa<NumberComponents_WP-1; aa++){
    //     for ( bb=aa+1; bb<NumberComponents_WP; bb++){
    //         if (ComponentAverages_WP(VOL,aa) < ComponentAverages_WP(VOL,bb)){
    //             // Exchange location of blobs aa and bb
    //             //printf("Switch blob %i with %i \n", OldLabel(aa),OldLabel(bb));
    //             // switch the label
    //             TempLabel = OldLabel(bb);
    //             OldLabel(bb) = OldLabel(aa);
    //             OldLabel(aa) = TempLabel;
    //             // switch the averages
    //             for (idx=0; idx<BLOB_AVG_COUNT; idx++){
    //                 TempValue = ComponentAverages_WP(idx,bb);
    //                 ComponentAverages_WP(idx,bb) = ComponentAverages_WP(idx,aa);
    //                 ComponentAverages_WP(idx,aa) = TempValue;
    //             }
    //         }
    //     }
    // }
    // NewLabel.resize(NumberComponents_WP);
    // for (aa=0; aa<NumberComponents_WP; aa++){
    //     // Match the new label for original blob aa
    //     bb=0;
    //     while (OldLabel(bb) != aa)	bb++;
    //     NewLabel(aa) = bb;
    // }
    // // Re-label the blob ID
    // //	printf("Re-labeling the blobs, now indexed by volume \n");
    // for (k=0; k<Nz; k++){
    //     for (j=0; j<Ny; j++){
    //         for (i=0; i<Nx; i++){
    //             if (Label_WP(i,j,k) > -1){
    //                 TempLabel = NewLabel(Label_WP(i,j,k));
    //                 Label_WP(i,j,k) = TempLabel;
    //             }
    //         }
    //     }
    // }  
}

void TwoPhase::PrintComponents(int timestep)
{
    if (Dm->rank()==0){
        printf("Print %i nonwetting component averages at timestep %i \n",(int)ComponentAverages_NWP.size(1),timestep);
        for (int b=0; b<NumberComponents_NWP; b++){
            //if (ComponentAverages_NWP(TRIMVOL,b) > 0.0){
            fprintf(NWPLOG,"%i ",timestep);
            if ( Label_NWP_map.empty() ) {
                // print index
                fprintf(NWPLOG,"%i ",b);
            } else {
                // print final global id
                fprintf(NWPLOG,"%i ",Label_NWP_map[b]);
            }
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VOL,b));
            //			fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(TRIMVOL,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(PRS,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(AWN,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(ANS,b));
            /*fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(JWN,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(KWN,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(LWNS,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(CWNS,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(KNWNS,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(KGWNS,b));*/
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VWNX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VWNY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VWNZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VNSX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VNSY,b));
            fprintf(NWPLOG,"%-15.15E\n",ComponentAverages_NWP(VNSZ,b));
            /*fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VWNSX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VWNSY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(VWNSZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNXX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNYY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNZZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNXY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNXZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GWNYZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GXX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GYY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GZZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GXY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GXZ,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(GYZ,b));
            //fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(TRAWN,b));
            //fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(TRJWN,b));
            //fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(INTCURV,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(EULER,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(CMX,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(CMY,b));
            fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(CMZ,b));
            fprintf(NWPLOG,"%-15.15E\n",ComponentAverages_NWP(VSQ,b));*/
            //				fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(NVERT,b));
            //			fprintf(NWPLOG,"%-15.15E ",ComponentAverages_NWP(NSIDE,b));
            //		fprintf(NWPLOG,"%-15.15E\n",ComponentAverages_NWP(NFACE,b));
            //}
        }
        fflush(NWPLOG);
        
        // for (int b=0; b<NumberComponents_WP; b++){
        //     if (ComponentAverages_WP(TRIMVOL,b) > 0.0 ){
        //         fprintf(WPLOG,"%i ",timestep);
        //         fprintf(WPLOG,"%i ",b);
        //         fprintf(WPLOG,"%-15.15E\n",ComponentAverages_WP(VOL,b));
                //			fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(TRIMVOL,b));
                /*fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(PRS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(AWN,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(AWNS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(JWN,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(KWN,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(LWNS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(CWNS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(KNWNS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(KGWNS,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VWNX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VWNY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VWNZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VWNSX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VWNSY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(VWNSZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNXX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNYY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNZZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNXY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNXZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GWNYZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GXX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GYY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GZZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GXY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GXZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(GYZ,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(CMX,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(CMY,b));
                fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(CMZ,b));
                fprintf(WPLOG,"%-15.15E\n",ComponentAverages_WP(VSQ,b));*/
                //fprintf(WPLOG,"%-15.15E ",ComponentAverages_WP(TRAWN,b));
                //fprintf(WPLOG,"%-15.15E\n",ComponentAverages_WP(TRJWN,b));
                
        //     }
        // }
        // fflush(WPLOG);
    }
}
