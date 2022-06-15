#include "analysis/TwoPhase.h"

#include <stdint.h>
#include <limits.h>

#include <string>

void TwoPhase::nsComputeQuantities(nsMyMesh & m_ns, DoubleArray & vx, DoubleArray & vy, DoubleArray & vz) {
    
    vcg::tri::Clean<nsMyMesh>::RemoveDuplicateVertex(m_ns,true);
    vcg::tri::Clean<nsMyMesh>::RemoveDegenerateFace(m_ns);
    vcg::tri::Clean<nsMyMesh>::RemoveDuplicateFace(m_ns);
    vcg::tri::UpdateTopology<nsMyMesh>::FaceFace(m_ns);
    vcg::tri::Clean<nsMyMesh>::RemoveDegenerateVertex(m_ns);
    vcg::tri::UpdateTopology<nsMyMesh>::VertexFace(m_ns);
    vcg::tri::Clean<nsMyMesh>::RemoveNonManifoldFace(m_ns);
    std::vector< std::pair<int, typename nsMyMesh::FacePointer> > CCV;
    int TotalCC=vcg::tri::Clean<nsMyMesh>::ConnectedComponents(m_ns, CCV);
//    ws_CCs =  vcg::tri::Clean<nsMyMesh>::RemoveSmallConnectedComponentsSize(m_ns, 0.0);
//    ws_CCs.first = ws_CCs.first - ws_CCs.second;
    std::vector< std::pair<int, typename nsMyMesh::FacePointer> > CCV2;
    TotalCC=vcg::tri::Clean<nsMyMesh>::ConnectedComponents(m_ns, CCV2);
    int totalUnref = vcg::tri::Clean<nsMyMesh>::RemoveUnreferencedVertex(m_ns, true);
    vcg::tri::Clean<nsMyMesh>::RemoveNonManifoldFace(m_ns);
    vcg::tri::Clean<nsMyMesh>::RemoveDegenerateVertex(m_ns);
    vcg::tri::Clean<nsMyMesh>::RemoveDuplicateVertex(m_ns,false);
    vcg::tri::Clean<nsMyMesh>::RemoveUnreferencedVertex(m_ns, true);
    vcg::tri::UpdateFlags<nsMyMesh>::VertexSetB(m_ns);
    vcg::tri::UpdateSelection<nsMyMesh>::VertexFromBorderFlag(m_ns,false  );
    vcg::tri::UpdateTopology<nsMyMesh>::AllocateEdge(m_ns);
    vcg::tri::Allocator<nsMyMesh>::CompactEveryVector(m_ns);
    vcg::tri::UpdateCurvatureFitting<nsMyMesh>::computeCurvature(m_ns);
    vcg::tri::UpdateSelection<nsMyMesh>::VertexClear(m_ns);
    vcg::tri::UpdateSelection<nsMyMesh>::FaceClear(m_ns);
    vcg::tri::Stat<nsMyMesh>::ComputeMeshArea(m_ns);
    vcg::tri::UpdateNormal<nsMyMesh>::PerFaceNormalized(m_ns);
    
    
    double sum_local_area = 0;
    double sum_nnsz = 0;
    double sum_nnnn_nsxx = 0;
    double sum_nnnn_nsyy = 0;
    double sum_nnnn_nszz = 0;
    double sum_nnnn_nsxz = 0;
    double sum_nnnn_nsyz = 0;
    double sum_vnsx = 0;
    double sum_vnsy = 0;
    double sum_vnsz = 0;
    double sum_vnsdnn = 0;
    
    
    int x_rank_min;
    int y_rank_min;
    int z_rank_min;
    int x_rank_max;
    int y_rank_max;
    int z_rank_max;
    
    x_rank_min = 1;
    y_rank_min = 1;
    z_rank_min = 1;
    x_rank_max = Nx-1;
    y_rank_max = Ny-1;
    z_rank_max = Nz-1;
    
    vcg::Point3<double> BottomCorner(x_rank_min,y_rank_min,z_rank_min);
    vcg::Point3<double> TopCorner(x_rank_max,y_rank_max,z_rank_max);
    vcg::Box3<double> RankBox(BottomCorner,TopCorner);
    
    nsMyMesh::CoordType b(0,0,0);
    
    
    
    for(nsMyMesh::FaceIterator fi=m_ns.face.begin(); fi!=m_ns.face.end();++fi) if(!(*fi).IsD()) {
        b = vcg::Barycenter(*fi);
        if(RankBox.IsIn(b) ) {
            (*fi).SetS();
            //  selCnt++;
        }
    }
    
    for(nsMyMesh::VertexIterator vi = m_ns.vert.begin(); vi != m_ns.vert.end(); ++vi) if(!(*vi).IsD()) {
        if(RankBox.IsIn((*vi).P()) ) {
            (*vi).SetS();
        }
    }
    
    TriLinPoly Velx,Vely,Velz;
    Point P;
    double area = 0;
    
    double normx, normy, normz;
    int a;
    int i,j,k;
    nsMyMesh::FaceIterator fi;
    
    ewn_.fill(0);
    Jns_.fill(0);
    Kns_.fill(0);
    
    double avgk1,avgk2,localJns,localKns;
    
    for(a = 0, fi=m_ns.face.begin(); fi!=m_ns.face.end();++fi,++a){
        if(!(*fi).IsD()) {
            if ((*fi).IsS()){
                // VELOCITY MEASUREMENTS
                b = vcg::Barycenter(*fi);
                
                P.x = b[0];
                P.y = b[1];
                P.z = b[2];
                
                i = floor(P.x);
                j = floor(P.y);
                k = floor(P.z);
                
                // AREA MEASUREMENT
                area = 0.5*DoubleArea(*fi);
                ewn_(i,j,k) += area;
                
                
                // NORMAL MEASUREMENT
                normx = (double)(*fi).N()[0];
                normy = (double)(*fi).N()[1];
                normz = (double)(*fi).N()[2];
                
                avgk1 = 0;
                avgk1 = (*(*fi).V(0)).K1();
                avgk1 += (*(*fi).V(1)).K1();
                avgk1 += (*(*fi).V(2)).K1();
                avgk1 *= 0.3333333333333333;
                
                avgk2 = 0;
                avgk2 = (*(*fi).V(0)).K2();
                avgk2 += (*(*fi).V(1)).K2();
                avgk2 += (*(*fi).V(2)).K2();
                avgk2 *= 0.3333333333333333;

                localJns = (avgk1+avgk2);
                localKns = avgk1*avgk2;
                if (isnan(localJns)) localJns = 0;
                if (isnan(localKns)) localKns = 0;
                if (abs(localKns) > 4 ) localKns = 0;
                if (abs(localJns) > 2 ) localJns = 0;
                Jns_(i,j,k) += localJns*area;
                Kns_(i,j,k) += localKns*area;
                
                sum_local_area += area;
            }
        }
    }
    
    ens = sum_local_area;
    
    Dm->CommunicateMeshHalo(Jns_);
    Dm->CommunicateMeshHalo(Kns_);
    
    Kns_tcenter = Jns_tcenter = ans_tcenter = 0;
    int kmin=1; int kmax=Nz-1;
    // if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=126; //5
    // if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-1)) kmax=1;
    // if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-2)) kmax=1; 
    // if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-3)) kmax=Nz-1-135;  //Nz-1-82
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=5; //python says 89, set to +1 from python
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-1)) kmax=Nz-1-5; //python says 93, set to -1 from python (Nz-94)

    int n;
    for (int k=kmin; k<kmax; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++){
                n = i + j*Nx + k*Nx*Ny;
                Jns_tcenter += Jns_(i,j,k);
                Kns_tcenter += Kns_(i,j,k);
                ans_tcenter += ewn_(i,j,k);
            }
        }
    }
    
    //printf("offset_distance=%f",offset_distance);
    char TimeString[8];
    char OffsetString[8];
    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(OffsetString,"%.1f",offset_distance);
    sprintf(TimeString,"%.05d",time_step);
    sprintf(LocalRankString,"%05d",Dm->rank());
    sprintf(LocalRankFilename,"%s%s%s%s%s%s",OffsetString,"_",TimeString,"_nsdata.",LocalRankString,".off");
    vcg::tri::io::ExporterOFF<nsMyMesh>::Save(m_ns,LocalRankFilename, 1);
    
}




void TwoPhase::ComputeAccurateNWInterfaceCurvatures(nwRootMesh & m_root_nw, double * avg_mean_curvature, double * avg_gaussian_curvature, double * avg_area, double * avg_vwndnn, DoubleArray & ewnwwn_, DoubleArray & wwnJn_, DoubleArray & vwndnnIIprime_, DoubleArray &vx, DoubleArray &vy, DoubleArray &vz) {
    
    vcg::tri::Clean<nwRootMesh>::RemoveDuplicateVertex(m_root_nw,true);
    vcg::tri::Clean<nwRootMesh>::RemoveDegenerateFace(m_root_nw);
    vcg::tri::Clean<nwRootMesh>::RemoveDuplicateFace(m_root_nw);
    vcg::tri::UpdateTopology<nwRootMesh>::FaceFace(m_root_nw);
    vcg::tri::Clean<nwRootMesh>::RemoveDegenerateVertex(m_root_nw);
    vcg::tri::UpdateTopology<nwRootMesh>::VertexFace(m_root_nw);
    vcg::tri::Clean<nwRootMesh>::RemoveNonManifoldFace(m_root_nw);
    
    std::vector< std::pair<int, typename nwRootMesh::FacePointer> > CCV;
    int TotalCC=vcg::tri::Clean<nwRootMesh>::ConnectedComponents(m_root_nw, CCV);
    
    std::pair<int,int> ws_CCs;
    double CCsize = double(Nz-2)/3.0;
    ws_CCs =  vcg::tri::Clean<nwRootMesh>::RemoveSmallConnectedComponentsSize(m_root_nw, 0);
    
    ws_CCs.first = ws_CCs.first - ws_CCs.second;
    
    std::vector< std::pair<int, typename nwRootMesh::FacePointer> > CCV2;
    TotalCC=vcg::tri::Clean<nwRootMesh>::ConnectedComponents(m_root_nw, CCV2);
    
    int totalUnref = vcg::tri::Clean<nwRootMesh>::RemoveUnreferencedVertex(m_root_nw, true);
    vcg::tri::Clean<nwRootMesh>::RemoveNonManifoldFace(m_root_nw);
    vcg::tri::Clean<nwRootMesh>::RemoveDegenerateVertex(m_root_nw);
    vcg::tri::Clean<nwRootMesh>::RemoveDuplicateVertex(m_root_nw,false);
    vcg::tri::Clean<nwRootMesh>::RemoveUnreferencedVertex(m_root_nw, true);
    vcg::tri::UpdateFlags<nwRootMesh>::VertexSetB(m_root_nw);
    vcg::tri::UpdateSelection<nwRootMesh>::VertexFromBorderFlag(m_root_nw,false);
    vcg::tri::UpdateTopology<nwRootMesh>::AllocateEdge(m_root_nw);
    vcg::tri::UpdateTopology<nwRootMesh>::FaceFace(m_root_nw);
    vcg::tri::Allocator<nwRootMesh>::CompactEveryVector(m_root_nw);
    vcg::tri::UpdateCurvatureFitting<nwRootMesh>::computeCurvature(m_root_nw);
    vcg::tri::UpdateSelection<nwRootMesh>::VertexClear(m_root_nw);
    vcg::tri::UpdateSelection<nwRootMesh>::FaceClear(m_root_nw);
    vcg::tri::Stat<nwRootMesh>::ComputeMeshArea(m_root_nw);
    vcg::tri::UpdateNormal<nwRootMesh>::PerFaceNormalized(m_root_nw);
    vcg::tri::UpdateNormal<nwRootMesh>::PerVertexNormalized(m_root_nw);
    
    int x_rank_min;
    int y_rank_min;
    int z_rank_min;
    int x_rank_max;
    int y_rank_max;
    int z_rank_max;
    
    x_rank_min = 1;
    y_rank_min = 1;
    z_rank_min = 1;
    x_rank_max = (Nx-1);
    y_rank_max = (Ny-1);
    z_rank_max = (Nz-1);
    
    vcg::Point3<double> BottomCorner(x_rank_min,y_rank_min,z_rank_min);
    vcg::Point3<double> TopCorner(x_rank_max,y_rank_max,z_rank_max);
    vcg::Box3<double> RankBox(BottomCorner,TopCorner);
    
    
    size_t selCnt = 0;
    
    nwMyMesh::CoordType b(0,0,0);
    for(nwRootMesh::FaceIterator fi=m_root_nw.face.begin(); fi!=m_root_nw.face.end();++fi) if(!(*fi).IsD()) {
        b = vcg::Barycenter(*fi);
        if(RankBox.IsIn(b) ) {
            (*fi).SetS();
            selCnt++;
        }
    }
    
    for(nwRootMesh::VertexIterator vi = m_root_nw.vert.begin(); vi != m_root_nw.vert.end(); ++vi) if(!(*vi).IsD()) {
        if(RankBox.IsIn((*vi).P()) ) {
            (*vi).SetS();
        }
    }
    
    
    
    double vwndnn;
    double vnnnnx, vnnnny, vnnnnz;
    double nnnnxx, nnnnxy, nnnnxz;
    double nnnnyx, nnnnyy, nnnnyz;
    double nnnnzx, nnnnzy, nnnnzz;
    
    
    
    double sum_local_Jwn = 0;
    double sum_local_Kwn = 0;
    
    double sum_nnz = 0;
    double sum_local_area = 0;
    double sum_vwndnn = 0;
    double sum_Knvwndnn = 0;
    
    double sum_IIprimezzvnnnnz = 0;
    double sum_wwnz = 0;
    double sum_wwnx = 0;
    double sum_wwny = 0;
    double sum_JnJnvwndnn = 0;
    
    double sum_delpdnn = 0;
    double delpdnnx = 0;
    double delpdnny = 0;
    double delpdnnz = 0;
    double sum_JnJnnnz = 0;
    double sum_TwoKnnz = 0;
    double sum_nnnnxz = 0;
    double sum_nnnnyz = 0;
    //   double sum_nnnnzz = 0;
    double sum_nnnnxx = 0;
    double sum_nnnnyy = 0;
    double sum_vwnz = 0;
    double sum_vwnx = 0;
    double sum_vwny = 0;
    double sum_IIprimevwndnn = 0;
    //   double sum_nw_speed = 0;
    
    sum_local_area=0;
    sum_nnz=0;
    sum_local_Jwn =0;
    sum_local_Kwn =0;
    sum_vwndnn =0;
    sum_Knvwndnn =0;
    sum_JnJnvwndnn =0;
    sum_wwnz =0;
    sum_wwnx =0;
    sum_wwny =0;
    
    
    
    
    Kwn_.fill(0);
    
    wwnnn_.fill(0);
    
    ewn_.fill(0);
    Jwn_.fill(0);
    
    speed_.fill(0);
    
    
    Point P,Q;
    double avgk1, avgk2;
    double localJwn;
    double localKwn;
    
    double normx, normy, normz;
    double thespeed = 0;
    
    vwndnn_tcenter = 0;
    double area = 0;
    int i,j,k,ii,jj,kk;
    
    
    double gp = 0;
    
    TriLinPoly Velx,Vely,Velz,speed;
    
    double sum_Jnnnz = 0;
    double compare_area;
    double dx1,dx2,dx3,cpx,cpy,cpz;
    double icpx, icpy, icpz;
    double fval;
    double vnormx,vnormy,vnormz;
    double qeval;
    double ax,ay,az;
    double dval;
    double v1x,v1y,v1z,v2x,v2y,v2z,v0x,v0y,v0z;
    double sum_nzmagvwn = 0;
    double euler_sum = 0;
    
    
    nwRootMesh::FaceIterator fi;
    int a;
    for(a = 0, fi=m_root_nw.face.begin(); fi!=m_root_nw.face.end();++fi,++a){
        if(!(*fi).IsD()) {
            if ((*fi).IsS()) {
                // VELOCITY MEASUREMENTS
                b = vcg::Barycenter(*fi);
                
                P.x = b[0];
                P.y = b[1];
                P.z = b[2];
                
                i = floor(P.x);
                j = floor(P.y);
                k = floor(P.z);
                
                // AREA MEASUREMENT
                area = 0.5*DoubleArea(*fi);
                
                // CURVATURE MEASUREMENTS
                avgk1 = 0;
                avgk1 = (*(*fi).V(0)).K1();
                avgk1 += (*(*fi).V(1)).K1();
                avgk1 += (*(*fi).V(2)).K1();
                avgk1 *= 0.3333333333333333;
                avgk2 = 0;
                avgk2 = (*(*fi).V(0)).K2();
                avgk2 += (*(*fi).V(1)).K2();
                avgk2 += (*(*fi).V(2)).K2();
                avgk2 *= 0.3333333333333333;
                localJwn = (avgk1+avgk2);
                localKwn = avgk1*avgk2;
                if (isnan(localJwn)) localJwn = 0;
                if (isnan(localKwn)) localKwn = 0;
                if (abs(localKwn) > 4 ) localKwn = 0;
                if (abs(localJwn) > 2 ) localJwn = 0;
                // if (abs(localKwn) > 1.0 ) localKwn = 0;
                // if (abs(localJwn) > 1.0 ) localJwn = 0; 
                // if (abs(localKwn) < 0.001 ) localKwn = 0;
                // if (abs(localJwn) < 0.001 ) localJwn = 0;                
                Jwn_(i,j,k) += localJwn*area;
                Kwn_(i,j,k) += localKwn*area;
                
                
                Velx.assign(vx,i,j,k);
                Vely.assign(vy,i,j,k);
                Velz.assign(vz,i,j,k);
                
                
                vwny = Vely.eval(P);
                vwnz = Velz.eval(P);
                
                
                
                vwndnn = vwnz*normz + vwnx*normx + vwny*normy;
                
                
                vwndnn_(i,j,k) += vwndnn*area;
                
                
                sum_local_area += area;
                
                // Int f dA
                sum_wwnz += vnnnnz*area;
                sum_wwnx += vnnnnx*area;
                sum_wwny += vnnnny*area;
                
                ewn_(i,j,k) += area;
            // std::cout << "area=" << area << std::endl;
                
                
                speed.assign(nw_speed_,i,j,k);
                thespeed = speed.eval(P);
                
                wwnnn_(i,j,k) += (thespeed*normz*normz + thespeed*normy*normy + thespeed*normx*normx)*area;
            }
        }
    }
    
    
    
    
    vwndnn_tcenter = sum_vwndnn;
    vwnz = sum_vwnz;
    vwny = sum_vwny;
    vwnx = sum_vwnx;
    wwnz_tcenter = sum_wwnz;
    wwnx_tcenter = sum_wwnx;
    wwny_tcenter = sum_wwny;
    Dm->CommunicateMeshHalo(vwndnn_);
    Dm->CommunicateMeshHalo(wwnnn_);
    Dm->CommunicateMeshHalo(ewn_);
    Dm->CommunicateMeshHalo(Jwn_);
    Dm->CommunicateMeshHalo(Kwn_);
    
    nw_speed = 0;
    int kmin=1; int kmax=Nz-1;
    // if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=126; //5
    // if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-1)) kmax=1;
    // if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-2)) kmax=1; 
    // if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-3)) kmax=Nz-1-135;  //Nz-1-82
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0) kmin=5; //python says 89, set to +1 from python
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == (Dm->nprocz()-1)) kmax=Nz-1-5; //python says 93, set to -1 from python (Nz-94)
 
    int n;
    for (int k = kmin; k < kmax; k++){
        for (int j=1; j<Ny-1; j++){
            for (int i=1; i<Nx-1; i++) {
                n = i + j*Nx + k*Nx*Ny;
                dendt += dendt_(i,j,k);
                nw_speed += speed_(i,j,k);
                awn_tcenter += ewn_(i,j,k);
                Jwn_tcenter += Jwn_(i,j,k);
                Kwn_tcenter += Kwn_(i,j,k);
                wwnnn += wwnnn_(i,j,k);
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
    sprintf(LocalRankFilename,"%s%s%s%s%s%s",OffsetString,"_",TimeString,"_nwdata.",LocalRankString,".off");
    vcg::tri::io::ExporterOFF<nwRootMesh>::Save(m_root_nw,LocalRankFilename, 1);
    
}


void TwoPhase::wsComputeQuantities(wsMyMesh & m_ws, DoubleArray & vx,DoubleArray & vy,DoubleArray & vz) {
    
    
    char TimeString[8];
    char OffsetString[8];
    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(OffsetString,"%.1f",offset_distance);
    sprintf(TimeString,"%.05d",time_step);
    sprintf(LocalRankString,"%05d",Dm->rank());
    sprintf(LocalRankFilename,"%s%s%s%s%s%s",OffsetString,"_",TimeString,"_wsdata.",LocalRankString,".off");
    vcg::tri::io::ExporterOFF<wsMyMesh>::Save(m_ws,LocalRankFilename, 1);
    
    
}




void TwoPhase::nwsComputeQuantities(nwsMyMesh & m_nws, nwMyMesh & m_nw, DoubleArray & vx, DoubleArray & vy, DoubleArray & vz, double b1, double b2, double b3, double avgpqx, double avgpqy, double avgpqz ) {
    
    if ( m_nw.FN() == 1) {
        
        vcg::tri::UpdateNormal<nwMyMesh>::PerFaceNormalized(m_nw);
        vcg::tri::Clean<nwsMyMesh>::RemoveDuplicateFace(m_nws);
        vcg::tri::Clean<nwsMyMesh>::RemoveDuplicateVertex(m_nws,false);
        vcg::tri::Clean<nwsMyMesh>::RemoveDegenerateFace(m_nws);
        vcg::tri::UpdateTopology<nwsMyMesh>::FaceFace(m_nws);
        vcg::tri::Clean<nwsMyMesh>::RemoveDegenerateVertex(m_nws);
        vcg::tri::UpdateTopology<nwsMyMesh>::VertexFace(m_nws);
        
        double x_rank_min;
        double y_rank_min;
        double z_rank_min;
        double x_rank_max;
        double y_rank_max;
        double z_rank_max;
        
        x_rank_min = 0.5;
        y_rank_min = 0.5 + double(offset_distance);
        z_rank_min = 0.5;
        x_rank_max = (Nx-2)+0.5;
        y_rank_max = (Ny-2)+0.5 + double(offset_distance);
        z_rank_max = (Nz-2)+0.5;
        
        vcg::Point3<double> BottomCorner(x_rank_min,y_rank_min,z_rank_min);
        vcg::Point3<double> TopCorner(x_rank_max,y_rank_max,z_rank_max);
        vcg::Box3<double> RankBox(BottomCorner,TopCorner);
        
        vcg::tri::UpdateSelection<nwsMyMesh>::VertexClear(m_nws);
        vcg::tri::UpdateSelection<nwsMyMesh>::FaceClear(m_nws);
        
        size_t selCnt = 0;
        
        nwsMyMesh::CoordType b(0,0,0);
        for(nwsMyMesh::FaceIterator fi=m_nws.face.begin(); fi!=m_nws.face.end();++fi) if(!(*fi).IsD()) {
            b = vcg::Barycenter(*fi);
            if(RankBox.IsIn(b) ) {
                (*fi).SetS();
                selCnt++;
            }
        }
        
        
        for(nwsMyMesh::VertexIterator vi = m_nws.vert.begin(); vi != m_nws.vert.end(); ++vi) if(!(*vi).IsD()) {
            if(RankBox.IsIn((*vi).P()) ) {
                (*vi).SetS();
            }
        }
        
        double RadAngle = 0;
        double DegAngle = 0;
        size_t DegCount = 0;
        for(nwsMyMesh::FaceIterator fi=m_nws.face.begin(); fi!=m_nws.face.end();++fi) if(!(*fi).IsD()) if ((*fi).IsS())
        {
            for(int z=0; z<(*fi).VN();++z)
            {
                RadAngle = std::abs(vcg::face::DihedralAngleRad(*fi,z));
                
                if (RadAngle > 1E-5) {
                    DegAngle += (RadAngle*(180.0 /3.14159265359));
                    DegCount++;
                }
            }
        }
        
        if (DegCount > 0) {
            efawns_count += DegCount;
            efawns += DegAngle;
        }
        
    }
    
}

void TwoPhase::Speed() {
    
    nw_speed_.fill(0);
    dPdt.fill(0);
    
    
    Dm->CommunicateMeshHalo(SDn_tminus);
    Dm->CommunicateMeshHalo(SDn_tplus);
    double sp;
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++) {
        
        sp = 0.25*(SDn_tplus(i,j,k)) - (SDn_tminus(i,j,k));
        nw_speed_(i,j,k) = sp;
        dPdt(i,j,k) = 0.25*(SDn_tplus(i,j,k) - SDn_tminus(i,j,k));
    }
    
}

