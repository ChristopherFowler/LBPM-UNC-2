// Header file for two-phase averaging class
#ifndef TwoPhase_INC
#define TwoPhase_INC

#include <array>
#include <vector>

#include "analysis/pmmc.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "analysis/analysis.h"

#include "shared_ptr.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"


#include <wrap/io_trimesh/export_off.h>
#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/selection.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/curvature_fitting.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <vcg/complex/algorithms/create/ball_pivoting.h>
#include<vcg/complex/algorithms/inertia.h>
#include <vcg/complex/algorithms/convex_hull.h>
#include <vcg/complex/append.h>

#include <vcg/complex/algorithms/hole.h>




class visMyEdge;
class visMyFace;
class visMyVertex;
struct visMyUsedTypes : public vcg::UsedTypes< vcg::Use<visMyVertex>::AsVertexType, vcg::Use<visMyEdge>::AsEdgeType,vcg::Use<visMyFace>::AsFaceType>{};
class visMyVertex:public vcg::Vertex<visMyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDirf, vcg::vertex::Mark>{};
class visMyFace:public vcg::Face< visMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3f, vcg::face::BitFlags, vcg::face::Mark > {};
class visMyEdge:public vcg::Edge<visMyUsedTypes>{};
class visMyMesh:public vcg::tri::TriMesh< std::vector<visMyVertex>, std::vector<visMyFace> , std::vector<visMyEdge>  > {
public:
    visMyMesh() = default;
};

class nwMyEdge;
class nwMyFace;
class nwMyVertex;
struct nwMyUsedTypes : public vcg::UsedTypes< vcg::Use<nwMyVertex>::AsVertexType, vcg::Use<nwMyEdge>::AsEdgeType,vcg::Use<nwMyFace>::AsFaceType>{};
class nwMyVertex:public vcg::Vertex<nwMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class nwMyFace:public vcg::Face< nwMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class nwMyEdge:public vcg::Edge<nwMyUsedTypes>{};
class nwMyMesh: public vcg::tri::TriMesh< std::vector<nwMyVertex>, std::vector<nwMyFace> , std::vector<nwMyEdge>  > {};

class nwRootEdge;
class nwRootFace;
class nwRootVertex;
struct nwRootUsedTypes : public vcg::UsedTypes< vcg::Use<nwRootVertex>::AsVertexType, vcg::Use<nwRootEdge>::AsEdgeType,vcg::Use<nwRootFace>::AsFaceType>{};
class nwRootVertex:public vcg::Vertex<nwRootUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class nwRootFace:public vcg::Face< nwRootUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class nwRootEdge:public vcg::Edge<nwRootUsedTypes, vcg::edge::EEAdj>{};
class nwRootMesh: public vcg::tri::TriMesh< std::vector<nwRootVertex>, std::vector<nwRootFace> , std::vector<nwRootEdge>  > {};

class nwsMyEdge;
class nwsMyFace;
class nwsMyVertex;
struct nwsMyUsedTypes : public vcg::UsedTypes< vcg::Use<nwsMyVertex>::AsVertexType, vcg::Use<nwsMyEdge>::AsEdgeType,vcg::Use<nwsMyFace>::AsFaceType>{};
class nwsMyVertex:public vcg::Vertex<nwsMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class nwsMyFace:public vcg::Face< nwsMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class nwsMyEdge:public vcg::Edge<nwsMyUsedTypes>{};
class nwsMyMesh:public vcg::tri::TriMesh< std::vector<nwsMyVertex>, std::vector<nwsMyFace> , std::vector<nwsMyEdge>  > {};

class nsMyEdge;
class nsMyFace;
class nsMyVertex;
struct nsMyUsedTypes : public vcg::UsedTypes< vcg::Use<nsMyVertex>::AsVertexType, vcg::Use<nsMyEdge>::AsEdgeType,vcg::Use<nsMyFace>::AsFaceType>{};
class nsMyVertex:public vcg::Vertex<nsMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class nsMyFace:public vcg::Face< nsMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class nsMyEdge:public vcg::Edge<nsMyUsedTypes, vcg::edge::EEAdj>{};
class nsMyMesh: public vcg::tri::TriMesh< std::vector<nsMyVertex>, std::vector<nsMyFace> , std::vector<nsMyEdge>  > {};

class wsMyEdge;
class wsMyFace;
class wsMyVertex;
struct wsMyUsedTypes : public vcg::UsedTypes< vcg::Use<wsMyVertex>::AsVertexType, vcg::Use<wsMyEdge>::AsEdgeType,vcg::Use<wsMyFace>::AsFaceType>{};
class wsMyVertex:public vcg::Vertex<wsMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class wsMyFace:public vcg::Face< wsMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class wsMyEdge:public vcg::Edge<wsMyUsedTypes, vcg::edge::EEAdj>{};
class wsMyMesh: public vcg::tri::TriMesh< std::vector<wsMyVertex>, std::vector<wsMyFace> , std::vector<wsMyEdge>  > {};


/* CONVEX HULL*/
class chMyEdge;
class chMyFace;
class chMyVertex;
struct chMyUsedTypes : public vcg::UsedTypes< vcg::Use<chMyVertex>::AsVertexType, vcg::Use<chMyEdge>::AsEdgeType,vcg::Use<chMyFace>::AsFaceType>{};
class chMyVertex:public vcg::Vertex<chMyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags, vcg::vertex::VFAdj,  vcg::vertex::CurvatureDird, vcg::vertex::Mark>{};
class chMyFace:public vcg::Face< chMyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef,vcg::face::Normal3d, vcg::face::BitFlags, vcg::face::Mark > {};
class chMyEdge:public vcg::Edge<chMyUsedTypes>{};
class chMyMesh:public vcg::tri::TriMesh< std::vector<chMyVertex>, std::vector<chMyFace> , std::vector<chMyEdge>  > {};

class TwoPhase{

public:
    
    int wsGlobal_dS;
    int nsGlobal_dS;
    int nwGlobal_dS;
    
    int wsGlobal_CCs;
    int nsGlobal_CCs;
    int nwGlobal_CCs;
    
    //...........................................................................
    int n_nw_pts,n_ns_pts,n_ws_pts,n_nws_pts,n_local_sol_pts,n_local_nws_pts;
    int n_nw_tris,n_ns_tris,n_ws_tris,n_nws_seg,n_local_sol_tris;
    //...........................................................................
private:
    int nc;
    int kstart,kfinish;

    double fluid_isovalue, solid_isovalue;
    double Volume;
    double nodeVolume;
    double subNodeVolume;
    // initialize lists for vertices for surfaces, common line
    DTMutableList<Point> nw_pts;
    DTMutableList<Point> ns_pts;
    DTMutableList<Point> ws_pts;
    DTMutableList<Point> nws_pts;
    DTMutableList<Point> local_sol_pts;
    DTMutableList<Point> local_nws_pts;
    DTMutableList<Point> tmp;

public:
    // initialize triangle lists for surfaces
    IntArray w_tris;
    IntArray nw_tris;
    IntArray ns_tris;
    IntArray ws_tris;
    IntArray nws_seg;
    IntArray local_sol_tris;

private:
    // Temporary storage arrays
    DoubleArray CubeValues;
    DoubleArray Values;
    DoubleArray DistanceValues;
    DoubleArray KGwns_values;
    DoubleArray KNwns_values;
    DoubleArray InterfaceSpeed;
    DoubleArray NormalVector;

    char *TempID;

    DoubleArray RecvBuffer;

    // CSV / text file where time history of averages is saved
    Array<FILE*> TIMELOG_subdomain;
    FILE *TIMELOG_local;
    FILE *TIMELOG_global;
    FILE *NWPLOG;
    FILE *WPLOG;
    
    FILE *GEOMETRY_local;
    FILE *GEOMETRY_global;

public: // Internal data
    //...........................................................................
    std::shared_ptr <Domain> Dm;
    int NumberComponents_WP,NumberComponents_NWP;

public: // local averages (to each MPI process)
        std::shared_ptr<IO::TriList> wn_mesh;//( new IO::TriList() );
        std::shared_ptr<IO::TriList> ns_mesh;//( new IO::TriList() );
        std::shared_ptr<IO::TriList> ws_mesh;//( new IO::TriList() );

        std::shared_ptr<IO::TriList> w_mesh;//( new IO::TriList() );
    size_t npx,npy,npz;
    double del_IIprimevwndnn;
    double deldel_IIprimevwndnn;
    double IIprimevwndnn;
    double trimdist;                    // pixel distance to trim surface for specified averages
    double porosity, poreVol;
    double volume;
    double awn, ans, aws, lwns;
    double wp_volume, nwp_volume;
    double As;
    double vol_w, vol_n;                // volumes the exclude the interfacial region
    double vol_w_tminus, vol_n_tminus;
    double pn, pw;                    // local phase averaged pressure
    double pn_tminus, pw_tminus;
    double efawns;                      // averaged contact angle
    size_t efawns_count;
    double euler,Kn,Jn,An,eulerW,eulerS;
    double nsEuler,wsEuler,nwEuler;
    double nwEuler_global, nsEuler_global;
    double nEuler,wEuler,sEuler;
    double nGlobalEuler, sGlobalEuler, wGlobalEuler;
    double Jwn,Jwn2;                       // average mean curavture - wn interface
    double Kwn,Kwn2;                         // average Gaussian curavture - wn interface
    double JwnW;                         // average mean curavture - wn interface
    double KwnW;                         // average Gaussian curavture - wn interface
    double KNwns;                       // wns common curve normal curavture
    double KGwns;                       // wns common curve geodesic curavture
    double trawn;                       // trimmed interfacial area
    double trJwn;                       // trimmed interfacial area
    double trRwn;                       // trimmed interfacial area
    double wwndnw;
    double wwnsdnwn;
    double Jwnwwndnw;
    double swn;
    
    double gradpz;
    DoubleArray nwPress;
    double nw_nnnnxz, nw_nnnnyz, nw_nnnnzz;
    double Jwn_old, Kwn_old;
    
    double Jws;
    double Kws;
    double Jns;
    double Kns;
    
    double JwsS;
    double KwsS;
    double JnsS;
    double KnsS;
    
    double JnJnnnz;
    double TwoKnnz;
    
    DoubleArray qDistances;
    
    int NSpheresOnRank;
    
    bool ComputeAccurateNWPhaseVolume;
    
    
    
    double Jwn2_global;
    double Kwn2_global;
    
    double Ca_global, Re_global;
    double alpha,sphere_diameter,average_rho_input, average_mu;
    size_t nwBorderVert, nsBorderVert, wsBorderVert;
    size_t nwBorderFace, nsBorderFace, wsBorderFace;
    size_t nwBorderVert_global, nsBorderVert_global, wsBorderVert_global;
    size_t nwBorderFace_global, nsBorderFace_global, wsBorderFace_global;
  
    double aw, an, Xn, pwn;
    double awn2, awn2_global;
    
    double Jnx,Jny,Jnz;
    double integrated_cpx,integrated_cpy,integrated_cpz;
    double gradz_an_tminus,gradz_an_tplus;
    double d_gradz_an_dt;
    double igl_mean_curvature, igl_Gaussian_curvature;

    DoubleArray vn;
    DoubleArray vw;
    
    DoubleArray time_avg_vw;
    DoubleArray time_avg_vn;
       
//    DoubleArray no_interface_time_avg_vw;
//    DoubleArray no_interface_time_avg_vn;
//
//    DoubleArray no_interface_vw;
//    DoubleArray no_interface_vn;
    
        DoubleArray nw;
    DoubleArray vwn;
    DoubleArray vws;
    DoubleArray vns;
    DoubleArray vwns;
    DoubleArray Gwn;
    DoubleArray Gwn_tminus;
    DoubleArray Gns;
    DoubleArray Gws;
    DoubleArray Gwns;
    
    DoubleArray PhaseVolume_old;
    DoubleArray dendt_;

public: // sub-domain averages
    DoubleArray trimdist_sub;               // pixel distance to trim surface for specified averages
    DoubleArray porosity_sub, poreVol_sub;
    DoubleArray volume_sub;
    DoubleArray awn_sub, ans_sub, aws_sub, lwns_sub, tmp_tmp;
    DoubleArray wp_volume_sub, nwp_volume_sub;
    DoubleArray As_sub;
    DoubleArray vol_w_sub, vol_n_sub;       // volumes the exclude the interfacial region
    DoubleArray pn_sub, pw_sub;           // local phase averaged pressure
    DoubleArray efawns_sub;                 // averaged contact angle
    DoubleArray euler_sub, Kn_sub, Jn_sub, An_sub, eulerW_sub,eulerS_sub;
    
    DoubleArray Jwn_sub;                    // average mean curavture - wn interface
    DoubleArray Kwn_sub;                    // average Gaussian curavture - wn interface
    DoubleArray JwnW_sub;                    // average mean curavture - wn interface
    DoubleArray KwnW_sub;                    // average Gaussian curavture - wn interface
    double ans_tcenter;
    double ans_tcenter_global;
    DoubleArray Jws_sub;
    DoubleArray Kws_sub;
    
    DoubleArray Jns_sub;
    DoubleArray Kns_sub;
    
    DoubleArray JwsS_sub;
    DoubleArray KwsS_sub;
    
    DoubleArray JnsS_sub;
    DoubleArray KnsS_sub;
    double nwndvwnsnnz;
    double nwnwwndnnz;
    
    DoubleArray wwnnn_;
    double deld_wwndnn;
    double wwnnn;
    IntArray ID;
    double HOvwnz;
    DoubleArray NxArray;
    DoubleArray NyArray;
    DoubleArray NzArray;
    DoubleArray HOvwnz_;
    DoubleArray rankSDs;
    DoubleArray nw_nnnnzz_;
    DoubleArray vwnz_;
    
    double ewnwwny;
    double deld_ewnwwn2;
    DoubleArray ewnwwnx_;
    DoubleArray ewnwwny_;
    
    DoubleArray potential;//.resize(Nx,Ny,Nz);  gravity_potential.fill(0);
    DoubleArray gravity;//.resize(Nx,Ny,Nz);  gravity.fill(0);
    DoubleArray potential_nnw;//.resize(Nx,Ny,Nz);   potential_nnw.fill(0);
    double grad_potential;
    double pot_nnw;
    double grav;
    // DoubleArray ewnwwnz;
    double ewnwwnz;
    DoubleArray Jnwwnx;
    DoubleArray Jnwwny;
    DoubleArray Jnwwnz;
    
    DoubleArray vwndnn_;

    
    DoubleArray TotalMass;
    DoubleArray DensityA, DensityB;

    DoubleArray Eta_sub; // Interfacial speed
    double Eta;
    double Eta_global;
    

    DoubleArray SDn_tminus;
    DoubleArray SDn_tplus;
    
    DoubleArray KNwns_sub;                  // wns common curve normal curavture
    DoubleArray KGwns_sub;                  // wns common curve geodesic curavture
    DoubleArray trawn_sub;                  // trimmed interfacial area
    DoubleArray trJwn_sub;                  // trimmed interfacial area
    DoubleArray trRwn_sub;                  // trimmed interfacial area
    DoubleArray wwndnw_sub;
    DoubleArray wwnsdnwn_sub;
    DoubleArray Jwnwwndnw_sub;
    DoubleArray swn_sub;
     DoubleArray vn_sub;
    DoubleArray vws_sub;
    DoubleArray vns_sub;
    DoubleArray vw_sub;
    DoubleArray vwn_sub;
    DoubleArray vwns_sub;
    DoubleArray Gwn_sub;
    DoubleArray Gns_sub;
    DoubleArray Gws_sub;
    DoubleArray Gwns_sub;
    
    DoubleArray Jnvwndnn_;
    DoubleArray Kwn_;
    
    
    double MAXZ, MINZ;
    
    double vwndnw;
    double vwndnw_global;
    DoubleArray vwndnw_sub;
    
    DoubleArray wwn;
    DoubleArray wwn_global;
    DoubleArray wwn_sub;
    
    DoubleArray nwOwnO;
    DoubleArray nwOwnO_global;
    DoubleArray nwOwnO_sub;
    
    DoubleArray nwOwnO2;
    DoubleArray nwOwnO_global2;
    DoubleArray nwOwnO_sub2;
    
    DoubleArray nwOwsO;
    DoubleArray nwOwsO_global;
    DoubleArray nwOwsO_sub;
    
    DoubleArray nwOwsO2;
    DoubleArray nwOwsO_global2;
    DoubleArray nwOwsO_sub2;
    
    DoubleArray nwOnsO;
    DoubleArray nwOnsO_global;
    DoubleArray nwOnsO_sub;
    
    DoubleArray nwOnsO2;
    DoubleArray nwOnsO_global2;
    DoubleArray nwOnsO_sub2;
   
        DoubleArray nw_sub;

    
    DoubleArray IIprimezzvwndnn;
 
    vector<double> MagVel;
    
    DoubleArray time_avg_vw_sub;
    DoubleArray time_avg_vn_sub;
       
    DoubleArray no_interface_time_avg_vw_sub;
    DoubleArray no_interface_time_avg_vn_sub;
       
    DoubleArray no_interface_vw_sub;
    DoubleArray no_interface_vn_sub;
    DoubleArray nw_speed_;
    double nw_speed;
    DoubleArray maggradphi;
    DoubleArray aw_sub, an_sub, Xn_sub, pwn_sub;
    
    DoubleArray Kns_;
    DoubleArray Jns_;
    
    DoubleArray an_sub1, an_sub2, an_sub3;
    DoubleArray Eta_sub1;
    
    double deld_nnz;
    
    double Kns_tcenter;
    double Jns_tcenter;
    
    DoubleArray nnz_;
    
    DoubleArray speed_;

public: // Global averages (all processes)
    
    
    
    
    double LEFT_an;
    double RIGHT_an;

    double dJwndt;
    double Jnvwndnn;
    
    double nnz;
    double nnz_global;
  
    
    
    double Knvwndnn;
    
   
    double gradzgradzIIprimezzvwndnn;
    
    
    double mass_local;
    double mass_global;
    double porosity_global;
    double volume_global;
    double pn_global,pw_global;               // local phase averaged pressure
    double pn_tminus_global, pw_tminus_global;
    double vol_w_global, vol_n_global;          // volumes the exclude the interfacial region
    double vol_w_tminus_global, vol_n_tminus_global;
    double awn_global,ans_global,aws_global;
    double lwns_global;
    double efawns_global;
    size_t efawns_global_count;                    // averaged contact angle
    double euler_global,Kn_global,Jn_global,An_global,eulerW_global,eulerS_global;
    double Jwn_global;                          // average mean curavture - wn interface
    double Kwn_global;                          // average Gaussian curavture - wn interface
    double JwnW_global;                          // average mean curavture - wn interface
    double KwnW_global;                          // average Gaussian curavture - wn interface
    
    double an1, an2, an3;
    double Eta1;
    
    double an_global1, an_global2,an_global3;
    double Eta_global1;
    
    double affinity;
    
    double an_old;
    double an_new;
    double an_global_old;
    
    double an_old3;
    double an_global_old3;
    
    double dendt;
    double dendt_global;
    
    double dendt3;
    double dendt_global3;
    
    double Jws_global;
    double Kws_global;
    double Jns_global;
    double Kns_global;
    
    double JwsS_global;
    double KwsS_global;
    double JnsS_global;
    double KnsS_global;
    
    double func, deriv;
    
    double local_nw_phase_volume;
    double global_nw_phase_volume;
    
    DoubleArray sincos;
    DoubleArray coscos;
    double euler_global_NS;
    
    double Jwn_tcenter_global;
    double Kwn_tcenter_global;
    double awn_tcenter_global;
    double awn_tminus_global;
    double an_tminus_global;
    double aw_tminus_global;
    double Jwn_tminus_global;
    
    double Kwn_tminus, Kwn_tplus, Kwn_tcenter;
    double Jwn_tminus, Jwn_tplus, Jwn_tcenter;
    double awn_tcenter, awn_tplus, awn_tminus;
    double an_tminus, aw_tminus;

    double vwndnn2;
    
    double vwndnn_tcenter;
    double vwndnn_tminus;
    double vwndnn_tplus;
    double wwnz_tcenter;
    double wwnx_tcenter;
    double wwny_tcenter;
    double dawndt;
    
    double pw_lhs,pw_rhs;
    double pn_lhs,pn_rhs;
    
    int pw_lhs_count,pw_rhs_count;
    int pn_lhs_count,pn_rhs_count;
    
    double left_wwnJnx_1;
    double left_wwnJny_1;
    double left_wwnJnz_1;
    
    double left_wwnJnx_2;
    double left_wwnJny_2;
    double left_wwnJnz_2;
    
    double right_wwnJnx_1;
    double right_wwnJny_1;
    double right_wwnJnz_1;
    
    double right_wwnJnx_2;
    double right_wwnJny_2;
    double right_wwnJnz_2;
    
    double left_wwnJnx;
    double left_wwnJny;
    double left_wwnJnz;
    double right_wwnJnx;
    double right_wwnJny;
    double right_wwnJnz;
    
    double right_IIprimezzvwndnn;
    double right_IIprimezzvwndnn_1;
    double right_IIprimezzvwndnn_2;
    
    double center_IIprimezzvwndnn;
    
    double left_IIprimezzvwndnn;
    double left_IIprimezzvwndnn_1;
    double left_IIprimezzvwndnn_2;
   
    double deldelddIIprimezzvwndnn;
    
    double deld_wwnKn;
    
    double left_ewnwwnx_1;
    double left_ewnwwny_1;
    double left_ewnwwnz_1;
    
    double nwndvwns;
     
    double left_ewnwwnx_2;
    double left_ewnwwny_2;
    double left_ewnwwnz_2;
    
    double right_ewnwwnx_1;
    double right_ewnwwny_1;
    double right_ewnwwnz_1;
    
    double right_ewnwwnx_2;
    double right_ewnwwny_2;
    double right_ewnwwnz_2;
    
    double left_ewnwwnx;
    double left_ewnwwny;
    double left_ewnwwnz;
    double right_ewnwwnx;
    double right_ewnwwny;
    double right_ewnwwnz;
    
    double IIprimezzwwnz;
    double deldelddIIprimezznnz;
    
    double gradz_an;
    double deld_wwnJn;
    double deld_ewnwwn;
    
    double vwndnnnnx, vwndnnnny, vwndnnnnz;
    
    double JnJnvwndnn;
    
    double sum_local_Jwn;
    double sum_local_Kwn;
    
    double KNwns_global;                        // wns common curve normal curavture
    double KGwns_global;                        // wns common curve geodesic curavture
    double trawn_global;                        // trimmed interfacial area
    double trJwn_global;                        // trimmed interfacial area
    double trRwn_global;                        // trimmed interfacial area
    double nwp_volume_global;                   // volume for the non-wetting phase
    double wp_volume_global;                    // volume for the wetting phase
    double As_global;
    double wwndnw_global;
    double wwnsdnwn_global;
    double Jwnwwndnw_global;
    double swn_global;
    double dEs,dAwn,dAns;
    double ens;
    double ews;
    
    double wwnx_global;
    double wwny_global;
    double wwnz_global;
    
    double vnsx;
    double vnsy;
    double vnsz;
    
    double nnnn_nszz;
    double nnnn_nsyz;
    double nnnn_nsxz;
    double nwSecondCC;
    double total_volume_global;
   
    double mesh_euler;
    double delpdnn;
    double vnsdnn_tcenter; 
    double deldelddIIprimennz;
    double deldIIprimez;
    double Jnnnz;

        double maxMagVel_global;

    double aw_global, an_global, Xn_global, pwn_global;
                        // Global surface energy (calculated by rank=0)
    DoubleArray vn_global;
    DoubleArray vw_global;
    
    DoubleArray press_tmp;
    
    DoubleArray time_avg_vw_global;
    DoubleArray time_avg_vn_global;
    
    DoubleArray no_interface_time_avg_vw_global;
    DoubleArray no_interface_time_avg_vn_global;
    
    DoubleArray no_interface_vw_global;
    DoubleArray no_interface_vn_global;
    
        DoubleArray nw_global;

    DoubleArray vwn_global;
    DoubleArray vws_global;
    DoubleArray vns_global;
    DoubleArray vwns_global;
    DoubleArray Gwn_global;
    DoubleArray Gwn_tminus_global;
    DoubleArray Gns_global;
    DoubleArray Gws_global;
    DoubleArray Gwns_global;

public: // Internal data
    std::array<size_t,3> subdivide;
    size_t Nx,Ny,Nz;
    IntArray PhaseID;    // Phase ID array (solid=0, non-wetting=1, wetting=2)
    
    BlobIDArray Label_WP;   // Wetting phase label
    BlobIDArray Label_NWP;  // Non-wetting phase label index (0:nblobs-1)
    std::vector<BlobIDType> Label_NWP_map;  // Non-wetting phase label for each index
    DoubleArray SDn;
    DoubleArray TempPhase;
     
    DoubleArray ParticularSDs;
    
    int time_flag;
    
    DoubleArray SDw;
    DoubleArray SDs;
    DoubleArray Phase;
    DoubleArray Press;
    DoubleArray dPdt;
    DoubleArray PhiField;
    DoubleArray BlobID;
    
    DoubleArray SDs_x;        // Gradient of the signed distance
    DoubleArray SDs_y;
    DoubleArray SDs_z;
    DoubleArray SDn_x;        // Gradient of the signed distance
    DoubleArray SDn_y;
    DoubleArray SDn_z;
    DoubleArray SDw_x;              // Gradient of the signed distance
    DoubleArray SDw_y;
    DoubleArray SDw_z;

   DoubleArray normalw_x;
   DoubleArray normalw_y;
   DoubleArray normalw_z;
    
    DoubleArray ewnwwnz_;
    
    IntArray BoolField;

    DoubleArray DelPhi;
    DoubleArray zero;
    
    DoubleArray VFmask, VFboundary, vfmask;
    
    DoubleArray wnCube;
    DoubleArray wsCube;
    DoubleArray nsCube;
    
    DoubleArray InactiveCellFluidCenterOfMass;
    
    IntArray InactiveCellsActiveNeighborInterpolationAbscissa;
    
    DoubleArray ParticularSDs_x;
    DoubleArray ParticularSDs_y;
    DoubleArray ParticularSDs_z;
    
    DoubleArray VolumeFraction;
    
    DoubleArray A;
    DoubleArray Aco;
    DoubleArray D;
    DoubleArray Dco;
    DoubleArray E;
    DoubleArray Qarray;
    
    DoubleArray A2;
    DoubleArray Aco2;
    DoubleArray D2;
    DoubleArray Dco2;
    DoubleArray E2;
    DoubleArray Qarray2;
    
    DoubleArray GradPhiX;
    DoubleArray GradPhiY;
    DoubleArray GradPhiZ;
    DoubleArray CField;
    
    
    DoubleArray PhaseVolume;
    
    DoubleArray Null3;
    DoubleArray Null9;
    
    DoubleArray IntegrationRegionOne;
    
    DoubleArray ewnwwn;
    DoubleArray wwnJn;
    DoubleArray vwndnnIIprime;
    
    DoubleArray IIprimexxvwndnn;
    DoubleArray IIprimexyvwndnn;
    DoubleArray IIprimexzvwndnn;
    DoubleArray IIprimeyxvwndnn;
    DoubleArray IIprimeyyvwndnn;
    DoubleArray IIprimeyzvwndnn;
    DoubleArray IIprimezxvwndnn;
    DoubleArray IIprimezyvwndnn;
    double Fx,Fy,Fz;
    DoubleArray ewn_;
    DoubleArray Jwn_;
    DoubleArray Knvwndnn_;
    
    double nzmagvwn;
    
    double nnnnxx, nnnnxy, nnnnxz;
    double nnnnyx, nnnnyy, nnnnyz;
    double nnnnzx, nnnnzy, nnnnzz;
    double nnsz,vwsz, n_wsz;
    double vwnx, vwny, vwnz;
    double centroid_ewn;
    double centroid_Jwn;
    double centroid_Knvwndnn;
    double grav_nw;
    DoubleArray JwnIIprimeMINUSnnnnxz_;
    DoubleArray JwnIIprimeMINUSnnnnyz_;
    DoubleArray JwnIIprimeMINUSnnnnzz_;
    double deld_JwnIIprimeMINUSnnnnz;
    
    DoubleArray IIprimexxnnz;
    DoubleArray IIprimexynnz;
    DoubleArray IIprimexznnz;
    DoubleArray IIprimeyxnnz;
    DoubleArray IIprimeyynnz;
    DoubleArray IIprimeyznnz;
    DoubleArray IIprimezxnnz;
    DoubleArray IIprimezynnz;
    DoubleArray IIprimezznnz;
    
    DoubleArray Delpvwndnnz;
    
    DoubleArray Knwwnx;
    DoubleArray Knwwny;
    DoubleArray Knwwnz;
    
    DoubleArray IIprimexz_;
    DoubleArray IIprimeyz_;
    DoubleArray IIprimezz_;
 
    DoubleArray PhaseField;
    DoubleArray PhaseField_old;
    
    DoubleArray temp_tplus;
    DoubleArray SDw_tplus;
    DoubleArray SDw_tplusplus;
    DoubleArray SDn_tplusplus;
    DoubleArray Phase_tplus;
    DoubleArray Phase_tminus;
    DoubleArray Vel_x;        // Velocity
    DoubleArray Vel_y;
    DoubleArray Vel_z;
    /* Time-averaged velocities */
    DoubleArray Time_Avg_Vel_x;
    DoubleArray Time_Avg_Vel_y;
    DoubleArray Time_Avg_Vel_z;
    //    Container for averages;
    DoubleArray ComponentAverages_WP;
    DoubleArray ComponentAverages_NWP;
    
    double delddelpvwndnn;
    
    
    double nnz_ns;
    
    
    int numToAvg;

public: // Public interfaces
    TwoPhase(std::shared_ptr <Domain> Dm, std::array<size_t,3> subdivide = { 1, 1, 1 });
    ~TwoPhase();
    void Initialize();
//    void SetupCubes(Domain &Dm);
    
    
    void UpdateMeshValues();
    void UpdateSolid();
    void ComputeDelPhi();
    void ComputeSDw();
    void ColorToSignedDistance(double Beta, DoubleArray &ColorData, DoubleArray &DistData);
    void ComputeLocal();
    void ComputeEwn();
  
    
    double right_gradz, left_gradz;
    double nnz_left, nnz_right, nnz_center;
    double gradzgradz_an;
    double grady_an, gradx_an;
 
    
    void SavePhaseIntoOldPhase(DoubleArray & NewField, DoubleArray & Field );
    
    void Speed();
    

 
    
  
 
  
    void ComputePhiDistance(DoubleArray &Field);
    
    void FindIsolatedLatticeSites();

    double nGaussianCurvature;
    double nGaussianCurvature_global;
  
    
    void ComputeVolumeFraction(DoubleArray& Field, DoubleArray& Field_x, DoubleArray& Field_y, DoubleArray& Field_z, DoubleArray& VolumeFraction, bool UpdateNormalToSolid, int rmin, int rmax, int geometry);
    
    void CommunicateField(int time_flag);
    void FluidCentroidInit(double affinity);
    
    void AssignComponentLabels();
    void ComponentAverages();
    inline void Reduce() { ReduceGlobal(); }
    void WriteSurfaces(int logcount);
    void NonDimensionalize(double D, double viscosity, double IFT);
    double SignedVolumeOfTriangle(Point p1, Point p2, Point p3);
    void PrintAll(int timestep);
    void GeometryOutput(int timestep);
    int GetCubeLabel(int i, int j, int k, IntArray &BlobLabel);
    
    void ComputeEuler(nwRootMesh &  );
    void nsComputeQuantities(nsMyMesh & /* some tri mesh */,DoubleArray &vx, DoubleArray &vy, DoubleArray &vz);
    void wsComputeQuantities(wsMyMesh & /* some tri mesh */,DoubleArray &vx, DoubleArray &vy, DoubleArray &vz);
  
 
    
   

    
    void nwsComputeQuantities(nwsMyMesh & ,nwMyMesh &, DoubleArray & vx, DoubleArray & vy, DoubleArray & vz, double pqx, double pqy, double pqz, double, double, double);

    void ComputeAccurateNWInterfaceCurvatures(nwRootMesh & m_nw, double * avg_mean_curvature,  double * avg_gaussian_curvature, double * avg_area, double * avg_vwndnn, DoubleArray & ewnwwn_, DoubleArray & wwnJn_, DoubleArray & vwndnnIIprime_, DoubleArray &vx, DoubleArray &vy, DoubleArray &vz);
    
  
    
    void SortBlobs();
    void PrintComponents(int timestep);
    void SetParams(double rhoA, double rhoB, double tauA, double tauB, double force_x, double force_y, double force_z, double alpha);
    double Volume_w(){
        return wp_volume_global;
    }
    double Volume_n(){
        return nwp_volume_global;
    }

    void AggregateLabels( const std::string& filename );
    
private: // Private interfaces
    void ReduceGlobal();
};






// bool compute_alpha_shapes(int dim, int numpoints, MeshModel &m, MeshModel &pm, double alpha, bool alphashape);



#endif
