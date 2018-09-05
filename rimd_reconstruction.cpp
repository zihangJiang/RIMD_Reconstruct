#include "rimd_reconstruction.h"
#include <iostream>
#include "util_3drotation_log_exp.h"
#include <queue>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
RIMD_Reconstruction::RIMD_Reconstruction()
{
    ref_mesh_.add_property(rotation_matrixs);
    ref_mesh_.add_property(scaling_matrixs);
    ref_mesh_.add_property(T_matrixs);
    ref_mesh_.add_property(LB_weights);
    ref_mesh_.add_property(log_dRs);
    ref_mesh_.add_property(buffer_for_compute);

//    defor_mesh_.add_property(rotation_matrixs);
//    defor_mesh_.add_property(scaling_matrixs);
//    defor_mesh_.add_property(T_matrixs);
//    defor_mesh_.add_property(LB_weights);
//    defor_mesh_.add_property(log_dRs);
}

void RIMD_Reconstruction::read_ref_mesh_from_file(std::string _filename)
{
    if(!OpenMesh::IO::read_mesh(ref_mesh_, _filename))
    {
        std::cout<<"RIMD_Reconstruction::read_ref_mesh_from_file::read file wrong!!!"<<std::endl;
        return ;
    }
    is_anchor.resize(ref_mesh_.n_vertices(),false);
}

void RIMD_Reconstruction::read_defor_mesh_from_file(std::string _filename)
{
    if(!OpenMesh::IO::read_mesh(defor_mesh_, _filename))
    {
        std::cout<<"RIMD_Reconstruction::read_defor_mesh_from_file::read file wrong!!!"<<std::endl;
        return ;
    }
    TriMesh::FaceIter f_it = defor_mesh_.faces_begin();
    for(;f_it!=defor_mesh_.faces_end();f_it++)
    {
        TriMesh::FaceHandle f_h = *f_it;
        TriMesh::FaceVertexIter fe_it0,fe_it1;
        fe_it0 = ref_mesh_.fv_iter(f_h);
        fe_it1 = defor_mesh_.fv_iter(f_h);
        int v00,v01,v02,v10,v11,v12;
        v00 = (*fe_it0).idx();fe_it0++;
        v01 = (*fe_it0).idx();fe_it0++;
        v02 = (*fe_it0).idx();fe_it0++;

        v10 = (*fe_it1).idx();fe_it1++;
        v11 = (*fe_it1).idx();fe_it1++;
        v12 = (*fe_it1).idx();fe_it1++;

        if(v00 == v10)
        {
            if(v01!=v11||v02!=v12)
            {
                std::cout<<"RIMD defor and ref are not compatible!!!"<<std::endl;
                return;
            }
        }
        else if(v00 == v11)
        {
            if(v01!=v12||v02!=v10)
            {
                std::cout<<"RIMD defor and ref are not compatible!!!"<<std::endl;
                return;
            }
        }
        else if(v00 == v12)
        {
            if(v01!=v10||v02!=v11)
            {
                std::cout<<"RIMD defor and ref are not compatible!!!"<<std::endl;
                return;
            }
        }
        else
        {
            std::cout<<"RIMD defor and ref are not compatible!!!"<<std::endl;
            return;
        }
    }



}

void RIMD_Reconstruction::read_anchor_points_id(const char* _filename)
{
    if(ref_mesh_.n_vertices()==0)
    {
        std::cout<<"RIMD_Reconstruction::read_anchor ref_mesh is empty!!!"<<std::endl;
        return;
    }
    FILE * fin = fopen(_filename, "r");
    if(fin==NULL)
    {
        std::cout<<"RIMD_Reconstruction read anchor points fail!!!"<<std::endl;
        return;
    }
    int idx;
    is_anchor.resize(ref_mesh_.n_vertices(),false);
    while(fscanf(fin, "%d", &idx) != EOF)
    {
        is_anchor[idx]=true;
    }
}

void RIMD_Reconstruction::Preprocess()
{
    compute_ref_LB_weights();
    compute_A_for_globalstep();
    compute_RIMD_of_ref_to_defor();
//    check_RIMD_correct();
}

void RIMD_Reconstruction::Reconstruction()
{
//    compute_RIMD_of_ref_to_defor();
    initial_P_to_ref_mesh();
    initial_Rs_for_reconstruction();
//    Eigen::VectorXd P_temp = P_;
    double last_energy = compute_Reconstruction_energy();
    double energy = last_energy;
    int iter_times=0;
//    std::cout<<"iter:"<<iter_times<<" norm:"<<(P_-P_temp).norm()<<std::endl;
    do{

        compute_b_for_globalstep();
//        std::cout<<b_<<std::endl;
        P_=A_solver_.solve(b_);
        if(A_solver_.info()!=Eigen::Success) {
          // decomposition failed
            std::cout<<"solver compute error:"<<A_solver_.info()<<std::endl;
        }
        iter_times++;
        if(iter_times%100==0)
        {
            energy = compute_Reconstruction_energy();
            std::cout<<"iter:"<<iter_times<<" error:"<<energy<<std::endl;
            last_energy = energy;
        }
        if(iter_times%100==1)
        {
            energy = compute_Reconstruction_energy();
            if(fabs(energy-last_energy)<0.01)
                break;
        }
        if(/*fabs(energy-last_energy)<0.001||*/iter_times>1000)
            break;
//        last_energy = energy;
        compute_local_step();
    }while(1);

}

void RIMD_Reconstruction::InterlateRIMD(double t, std::string _filename)
{
    TriMesh::VertexIter v_it = ref_mesh_.vertices_begin();
    Eigen::Matrix3d iden = Eigen::Matrix3d::Identity();
    for(;v_it!=ref_mesh_.vertices_end();v_it++)
    {
        Eigen::Matrix3d Scaling = ref_mesh_.property(scaling_matrixs, *v_it);
        ref_mesh_.property(scaling_matrixs, *v_it) = (1.0-t)*iden + t*Scaling;
    }

    TriMesh::HalfedgeIter he_it = ref_mesh_.halfedges_begin();
    Eigen::Matrix3d zero;
    zero.setZero();
    for(;he_it!=ref_mesh_.halfedges_end();he_it++)
    {
        Eigen::Matrix3d logR = ref_mesh_.property(log_dRs,*he_it);
        ref_mesh_.property(log_dRs,*he_it) = (1.0-t)*zero + t*logR;
    }
    check_RIMD_correct();

    // Store S and logdR feature separately.
    double S_feature[56994];
    int v_num=0;
    for(TriMesh::VertexIter v_it=ref_mesh_.vertices_begin(); v_it!=ref_mesh_.vertices_end(); v_it++)
    {
        Eigen::Matrix3d Scale = ref_mesh_.property(scaling_matrixs, *v_it);
        //std::cout<<"Scale "<<Scale<<std::endl;
        S_feature[v_num*6+0] = Scale(0,0);
        S_feature[v_num*6+1] = Scale(0,1);
        S_feature[v_num*6+2] = Scale(0,2);
        S_feature[v_num*6+3] = Scale(1,1);
        S_feature[v_num*6+4] = Scale(1,2);
        S_feature[v_num*6+5] = Scale(2,2);
        v_num++;
    }

    double logdR_feature[170946];
    int half_edge_num=0;
    for(TriMesh::HalfedgeIter h_it=ref_mesh_.halfedges_begin();h_it!=ref_mesh_.halfedges_end();h_it++)
    {
        Eigen::Matrix3d logdR = ref_mesh_.property(log_dRs,*h_it);
        //std::cout<<"logdR "<<std::endl<<logdR<<std::endl;
        logdR_feature[half_edge_num*3+0] = logdR(0,1);
        logdR_feature[half_edge_num*3+1] = logdR(0,2);
        logdR_feature[half_edge_num*3+2] = logdR(1,2);
        half_edge_num++;
    }

    // Write out as binary file.
    std::ofstream f(_filename,std::ios::binary);
    if(!f)
    {
        std::cout << "创建文件失败" <<std::endl;
        return;
    }
    f.write((char*)S_feature, 56994*sizeof(double));
    f.write((char*)logdR_feature, 170946*sizeof(double));
    f.close();

    // Test
    double read_out_data[227940];
    std::ifstream r(_filename, std::ios::binary);
    if(!r)
    {
        std::cout << "读取文件失败" <<std::endl;
        return;
    }
    r.read((char*)read_out_data,227940*sizeof(double));
    for(int i = 0; i < 56994; i++){
        if(read_out_data[i] - S_feature[i] > 0.001)
            std::cout<<"Read-Write Error!"<<std::endl;
    }
    for(int i = 56994; i < 227940; i++){
        if(read_out_data[i] - logdR_feature[i-56994] > 0.001)
            std::cout<<"Read-Write Error!"<<std::endl;
    }
    f.close();
}

void RIMD_Reconstruction::LoadRIMD(std::string _filename)
{
    double rimd_feature[227940];
    std::ifstream r(_filename, std::ios::binary);
    if(!r)
    {
        std::cout << "读取文件失败" <<std::endl;
        return;
    }
    r.read((char*)rimd_feature,227940*sizeof(double));
    int i=0;
    for(TriMesh::VertexIter v_it = ref_mesh_.vertices_begin();v_it!=ref_mesh_.vertices_end();v_it++){
        Eigen::Matrix3d Scale;
        Scale(0, 0)=rimd_feature[i+0];
        Scale(0, 1)=rimd_feature[i+1];
        Scale(1, 0)=rimd_feature[i+1];
        Scale(0, 2)=rimd_feature[i+2];
        Scale(2, 0)=rimd_feature[i+2];
        Scale(1, 1)=rimd_feature[i+3];
        Scale(1, 2)=rimd_feature[i+4];
        Scale(2, 1)=rimd_feature[i+4];
        Scale(2, 2)=rimd_feature[i+5];
        ref_mesh_.property(scaling_matrixs, *v_it) = Scale;
        i=i+6;
    }
    for(TriMesh::HalfedgeIter he_it = ref_mesh_.halfedges_begin();he_it!=ref_mesh_.halfedges_end();he_it++){
        Eigen::Matrix3d logdR;
        logdR(0, 0)=0;
        logdR(1, 1)=0;
        logdR(2, 2)=0;
        logdR(0, 1)=rimd_feature[i+0];
        logdR(1, 0)=-rimd_feature[i+0];
        logdR(0, 2)=rimd_feature[i+1];
        logdR(2, 0)=-rimd_feature[i+1];
        logdR(1, 2)=rimd_feature[i+2];
        logdR(2, 1)=-rimd_feature[i+2];
        ref_mesh_.property(log_dRs,*he_it) = logdR;
        i=i+3;
    }
    r.close();
}

void RIMD_Reconstruction::GetReconstructionMesh(TriMesh &mesh)
{
    mesh.clear();
    TriMesh::VertexIter v_it = ref_mesh_.vertices_begin();
    for(;v_it!=ref_mesh_.vertices_end();v_it++)
    {
        mesh.add_vertex(TriMesh::Point(0.0,0.0,0.0));
    }
    for(v_it = ref_mesh_.vertices_begin();v_it!=ref_mesh_.vertices_end();v_it++)
    {
        int id = (*v_it).idx();
        mesh.set_point(*v_it,TriMesh::Point(P_[3*id],P_[3*id+1],P_[3*id+2]));
    }
    TriMesh::FaceIter f_it = ref_mesh_.faces_begin();
    for(;f_it!=ref_mesh_.faces_end();f_it++)
    {
        TriMesh::FaceVertexIter fv_it = ref_mesh_.fv_iter(*f_it);
        TriMesh::VertexHandle v0,v1,v2;
        v0 = *fv_it;fv_it++;
        v1 = *fv_it;fv_it++;
        v2 = *fv_it;fv_it++;
        mesh.add_face(v0,v1,v2);
    }
    mesh.update_normals();
}

void RIMD_Reconstruction::compute_RIMD_of_ref_to_defor()
{
//    compute_ref_LB_weights();
    compute_ref_to_defor_Tmatrixs();
    compute_rotation_scaling_matrixs();
    compute_logdR_matrixs();
//    int test_id[5]={93,98,122,24317,37749};
//    for(int i=0;i<5;i++)
//    {
//        TriMesh::VertexHandle v(test_id[i]);
//        for(int i=0;i<3;i++)
//            for(int j=0;j<3;j++)
//                if(_isnan(ref_mesh_.property(rotation_matrixs,v)(i,j)))
//                {std::cout<<v.idx()<<std::endl;}
//        TriMesh::VertexVertexIter vv_it = ref_mesh_.vv_iter(*vv_it);
//        for(;vv_it.is_valid();vv_it++)
//        {
//            for(int i=0;i<3;i++)
//                for(int j=0;j<3;j++)
//                    if(_isnan(ref_mesh_.property(rotation_matrixs,*vv_it)(i,j)))
//                        std::cout<<(*vv_it).idx()<<std::endl;
//        }
//    }
//    std::cout<<"---------"<<std::endl;
}

void RIMD_Reconstruction::compute_ref_LB_weights()
{
    TriMesh::EdgeIter e_it, e_end(ref_mesh_.edges_end());
    TriMesh::HalfedgeHandle    h0, h1, h2;
    TriMesh::VertexHandle      v0, v1;
    TriMesh::Point             p0, p1, p2, d0, d1;
    TriMesh::Scalar w;
    for (e_it=ref_mesh_.edges_begin(); e_it!=e_end; ++e_it)
    {
        w  = 0.0;
        if(ref_mesh_.is_boundary(*e_it))
        {
            h0 = ref_mesh_.halfedge_handle(e_it.handle(), 0);
            if(ref_mesh_.is_boundary(h0))
                h0 = ref_mesh_.halfedge_handle(e_it.handle(),1);

            v0 = ref_mesh_.to_vertex_handle(h0);
            v1 = ref_mesh_.from_vertex_handle(h0);
            p0 = ref_mesh_.point(v0);
            p1 = ref_mesh_.point(v1);
            h1 = ref_mesh_.next_halfedge_handle(h0);
            p2 = ref_mesh_.point(ref_mesh_.to_vertex_handle(h1));
            d0 = (p0-p2).normalize();
            d1 = (p1-p2).normalize();
            w += 2.0 / tan(acos(std::min(0.99, std::max(-0.99, (d0|d1)))));
            w = std::max(0.0, w);//w小于0的时候还要仔细思考一下怎么处理
            ref_mesh_.property(LB_weights,e_it) = w;
            // Debug
            //ref_mesh_.property(LB_weights,e_it) = 1.0;
            continue;
        }
        h0 = ref_mesh_.halfedge_handle(e_it.handle(), 0);
        v0 = ref_mesh_.to_vertex_handle(h0);
        p0 = ref_mesh_.point(v0);

        h1 = ref_mesh_.halfedge_handle(e_it.handle(), 1);
        v1 = ref_mesh_.to_vertex_handle(h1);
        p1 = ref_mesh_.point(v1);

        h2 = ref_mesh_.next_halfedge_handle(h0);
        p2 = ref_mesh_.point(ref_mesh_.to_vertex_handle(h2));
        d0 = (p0 - p2).normalize();
        d1 = (p1 - p2).normalize();
        w += 1.0 / tan(acos(std::min(0.99, std::max(-0.99, (d0|d1)))));

        h2 = ref_mesh_.next_halfedge_handle(h1);
        p2 = ref_mesh_.point(ref_mesh_.to_vertex_handle(h2));
        d0 = (p0 - p2).normalize();
        d1 = (p1 - p2).normalize();
        w += 1.0 / tan(acos(std::min(0.99, std::max(-0.99, (d0|d1)))));

        w = std::max(0.0, w);//w小于0的时候还要仔细思考一下怎么处理
        ref_mesh_.property(LB_weights,e_it) = w;
        // Debug
        //ref_mesh_.property(LB_weights,e_it) = 1.0;
    }
}

void RIMD_Reconstruction::compute_ref_to_defor_Tmatrixs()
{
    TriMesh::VertexIter v_it, v_to_it;
    // 这里要求ref和defor的网格不仅拓扑一致，顶点的顺序也是要对应好的
    for(v_it=ref_mesh_.vertices_begin(),v_to_it=defor_mesh_.vertices_begin()
        ;v_it!=ref_mesh_.vertices_end()&&v_to_it!=defor_mesh_.vertices_end()
        ;v_it++,v_to_it++)
    {
//        if((*v_it).idx()!=(*v_to_it).idx())
//            std::cout<<"RIMD::compute_ref_to_defor_Tmatrixs different topology!!!"<<std::endl;
        compute_Ti(*v_it,*v_to_it);
//        if((*v_it).idx()<5)
//            std::cout<<ref_mesh_.property(T_matrixs,*v_it)<<std::endl;
    }
}

void RIMD_Reconstruction::compute_rotation_scaling_matrixs()
{
    TriMesh::VertexIter v_it;
    for(v_it=ref_mesh_.vertices_begin(); v_it!=ref_mesh_.vertices_end(); v_it++)
    {
        Eigen::MatrixXd T = ref_mesh_.property(T_matrixs,*v_it);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(T, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::Matrix3d U,V;
        U=svd.matrixU();
        V=svd.matrixV();
        Eigen::Matrix3d S(svd.singularValues().asDiagonal());
        Eigen::Matrix3d Temp=Eigen::Matrix3d::Identity();
        Temp(2,2) = (U*V.transpose()).determinant();
        Eigen::Matrix3d R=U*Temp*V.transpose();
        Eigen::Matrix3d Scale = V*Temp*S*V.transpose();
//        if(R.determinant()<0)
//        {
////            std::cout<<"ref mesh vertex "<<(*v_it).idx()<<" has reflection while deforming to defor mesh "<<std::endl;
//            R=-R;
//            Scale = -Scale;
//        }
        ref_mesh_.property(rotation_matrixs,*v_it)=R;
        ref_mesh_.property(scaling_matrixs,*v_it)=Scale;
//        if((*v_it).idx()<5)
//        {
//            std::cout<<R<<std::endl;
//            std::cout<<Scale<<std::endl;
//        }

    }
}

void RIMD_Reconstruction::compute_logdR_matrixs()
{
    TriMesh::HalfedgeIter h_it;
    for(h_it=ref_mesh_.halfedges_begin();h_it!=ref_mesh_.halfedges_end();h_it++)
    {
        TriMesh::VertexHandle to_vertex,from_vertex;
        to_vertex = ref_mesh_.to_vertex_handle(*h_it);
        from_vertex = ref_mesh_.from_vertex_handle(*h_it);
        Eigen::Matrix3d Ri = ref_mesh_.property(rotation_matrixs,from_vertex);
        Eigen::Matrix3d Rj = ref_mesh_.property(rotation_matrixs,to_vertex);
        Eigen::Matrix3d dR = Ri.transpose()*Rj;
        dR = rotation_log_exp::log(dR);
        ref_mesh_.property(log_dRs,*h_it) = dR;
//        if((*h_it).idx()<5) std::cout<<dR<<std::endl;
    }
}

void RIMD_Reconstruction::compute_Ti(TriMesh::VertexHandle v_it, TriMesh::VertexHandle v_to_it)
{
    if(v_it.idx()!=v_to_it.idx())
        std::cout<<"RIMD_Reconstruction::compute_Ti correspond is wrong!!!"<<std::endl;
    TriMesh::VertexEdgeIter veiter=ref_mesh_.ve_iter(v_it);
    int v_id = v_it.idx();
    TriMesh::Point p0,p1;
    p0 = ref_mesh_.point(v_it);
    p1 = defor_mesh_.point(v_to_it);
    Eigen::Matrix3d L,RI;
    L.setZero();
    RI.setZero();
    double tolerance = 1.0e-6;
    //三角形太小的话，计算的矩阵行列式就很小了，为了鲁棒性，需要对计算数据进行倍增
    TriMesh::HalfedgeHandle h_e=ref_mesh_.halfedge_handle(v_it);
    TriMesh::VertexHandle test_v = ref_mesh_.to_vertex_handle(h_e);
    TriMesh::Point tp0,tp1;
    tp0 = ref_mesh_.point(v_it); tp1 = ref_mesh_.point(test_v);
    double scale=1.0;
    if(((tp0[0]-tp1[0])*(tp0[0]-tp1[0])+(tp0[1]-tp1[1])*(tp0[1]-tp1[1])+(tp0[2]-tp1[2])*(tp0[2]-tp1[2]))<0.01)
        scale = 100;

    for(;veiter.is_valid();veiter++)
    {
        double weight = sqrt(ref_mesh_.property(LB_weights,(*veiter)));
        int to_id;
        TriMesh::VertexHandle to_v=ref_mesh_.to_vertex_handle(ref_mesh_.halfedge_handle(*veiter, 0));
        if(to_v.idx()==v_id)
            to_v = ref_mesh_.from_vertex_handle(ref_mesh_.halfedge_handle(*veiter, 0));
        to_id = to_v.idx();
        Eigen::Vector3d eij0,eij1;
        TriMesh::Point q0,q1;
        q0 = ref_mesh_.point(to_v);
        q1 = defor_mesh_.point(to_v);
        eij0(0) = p0[0]-q0[0];
        eij0(1) = p0[1]-q0[1];
        eij0(2) = p0[2]-q0[2];
        eij0*=weight*scale;

        eij1(0) = p1[0]-q1[0];
        eij1(1) = p1[1]-q1[1];
        eij1(2) = p1[2]-q1[2];
        eij1*=weight*scale;

        L+=eij1*eij0.transpose();
        RI+=eij0*eij0.transpose();
    }    
    Eigen::MatrixXd T;
    if(fabs(RI.determinant())>tolerance)
         T = L*RI.inverse();
    else
    {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(RI, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::Matrix3d U,V;
        U=svd.matrixU();
        V=svd.matrixV();
        Eigen::Matrix3d S_inv(svd.singularValues().asDiagonal());
        for(int i=0;i<3;i++)
        {
            if(S_inv(i,i)>tolerance)
                S_inv(i,i)=1.0/S_inv(i,i);
            else
                S_inv(i,i)=0.0;
        }
        T = L*V*S_inv*U.transpose();
        std::cout<<"T "<<(v_it).idx()<<" is singular, use pseudo inverse to computation!"<<std::endl;
//        if((v_it).idx()==0)
//        {
//            std::cout<<T<<std::endl;
//            std::cout<<L<<std::endl;
//            std::cout<<RI<<std::endl;
//        }
        //        std::cout<<T<<std::endl;
    }
    ref_mesh_.property(T_matrixs,v_it) = T;
//    if(v_it.idx()==122)
//    {
//        std::cout<<RI.determinant()<<std::endl;
//        std::cout<<L<<std::endl;
//        std::cout<<RI<<std::endl;
//        std::cout<<T<<std::endl;
//    }
}

void RIMD_Reconstruction::compute_A_for_globalstep()
{
    A.resize(3*ref_mesh_.n_vertices(),3*ref_mesh_.n_vertices());
    std::cout<<"A size: "<<A.cols()<<" "<<A.rows()<<std::endl;
//    Eigen::MatrixXd A_test(3*ref_mesh_.n_vertices(),3*ref_mesh_.n_vertices());
    std::vector<Eigen::Triplet<double> >    tripletlist;
    TriMesh::VertexIter v_it = ref_mesh_.vertices_begin();
    for(;v_it!=ref_mesh_.vertices_end();v_it++)
    {
        TriMesh::VertexEdgeIter ve_iter = ref_mesh_.ve_iter(*v_it);
        int center_id = (*v_it).idx();
        double center_val[3]={0.0,0.0,0.0};
        //固定边界测试
        if(is_anchor[(*v_it).idx()])
        {
            for(int i=0;i<3;i++)
            {
                tripletlist.push_back(Eigen::Triplet<double>(3*center_id+i,3*center_id+i,1.0));
//                A_test(3*center_id+i,3*center_id+i)=1.0;
            }
            continue;
        }
        for(;ve_iter.is_valid();ve_iter++)
        {
            double w = ref_mesh_.property(LB_weights,*ve_iter);
            TriMesh::VertexHandle to_v = ref_mesh_.to_vertex_handle(ref_mesh_.halfedge_handle(*ve_iter,0));
            if(to_v.idx()==center_id)
                to_v = ref_mesh_.from_vertex_handle(ref_mesh_.halfedge_handle(*ve_iter,0));
            for(int i=0;i<3;i++)
            {
                center_val[i]+=w;
                tripletlist.push_back(Eigen::Triplet<double>(3*center_id+i,3*to_v.idx()+i,-w));
//                A_test(3*center_id+i,3*to_v.idx()+i) = -w;
            }
        }
        for(int i=0;i<3;i++)
        {
            tripletlist.push_back(Eigen::Triplet<double>(3*center_id+i,3*center_id+i,center_val[i]));
//            A_test(3*center_id+i,3*center_id+i) = center_val[i];
        }
    }
    A.setFromTriplets(tripletlist.begin(),tripletlist.end());
    A_solver_.compute(A);
    if(A_solver_.info()!=Eigen::Success) {
      // decomposition failed
        std::cout<<"A decompose solver compute error:"<<A_solver_.info()<<std::endl;
    }

//    std::cout<<A<<std::endl;
    std::cout<<"A determinant:"<<A_solver_.determinant()<<std::endl;
//    Eigen::VectorXd rb(3*ref_mesh_.n_vertices());
//    rb.setOnes();
//    Eigen::VectorXd x = A_solver_.solve(rb);
//    std::cout<<(A*x-rb).norm()<<std::endl;
}

void RIMD_Reconstruction::compute_b_for_globalstep()
{
    b_.resize(3*ref_mesh_.n_vertices());
    compute_internal_for_globalstep();
    TriMesh::VertexIter v_it = ref_mesh_.vertices_begin();
    int num=0;
    for(;v_it!=ref_mesh_.vertices_end();v_it++)
    {
        TriMesh::VertexEdgeIter ve_it = ref_mesh_.ve_iter(*v_it);
        Eigen::Vector3d temp(0.0,0.0,0.0);
        //固定边界测试
        if(is_anchor[(*v_it).idx()])
        {
            TriMesh::Point p = ref_mesh_.point(*v_it);
            b_.block<3,1>(3*(*v_it).idx(),0) = Eigen::Vector3d(p[0],p[1],p[2]);
            continue;
        }
        for(;ve_it.is_valid();ve_it++)
        {
            TriMesh::VertexHandle v0 = ref_mesh_.to_vertex_handle(ref_mesh_.halfedge_handle(*ve_it,0));
            if(v0.idx()==(*v_it).idx())
                v0 = ref_mesh_.from_vertex_handle(ref_mesh_.halfedge_handle(*ve_it,0));
            TriMesh::Point Pj,Pk;
            Pj = ref_mesh_.point(*v_it);
            Pk = ref_mesh_.point(v0);
            Eigen::Vector3d ejk(Pj[0]-Pk[0],Pj[1]-Pk[1],Pj[2]-Pk[2]);
            double ck,cj;
            cj = 1.0/double(ref_mesh_.valence(*v_it));
            ck = 1.0/double(ref_mesh_.valence(v0));
            double cjk=ref_mesh_.property(LB_weights,*ve_it);
            temp+=cjk*(ck*ref_mesh_.property(buffer_for_compute,v0)+cj*ref_mesh_.property(buffer_for_compute,*v_it))*ejk;
        }
        temp*=0.5;
        b_.block<3,1>(3*(*v_it).idx(),0) = temp;
//        if((*v_it).idx()==42)   std::cout<<temp<<std::endl;
    }
}

void RIMD_Reconstruction::compute_local_step()
{
    compute_Q_for_local_step();
    TriMesh::VertexIter v_it = ref_mesh_.vertices_begin();
    for(;v_it!=ref_mesh_.vertices_end();v_it++)
    {
        Eigen::MatrixXd Qi = ref_mesh_.property(buffer_for_compute,(*v_it));
        //参考“实现此算法的一些参考”中的Least-Squares Rigid Motion Using SVD
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(Qi, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::Matrix3d U,V;
        U=svd.matrixU();
        V=svd.matrixV();
        Eigen::Matrix3d T=Eigen::Matrix3d::Identity();
        T(2,2) = (V*U.transpose()).determinant();
        Eigen::Matrix3d R=V*T*U.transpose();
        ref_mesh_.property(rotation_matrixs,(*v_it))=R;
    }
}

void RIMD_Reconstruction::compute_Q_for_local_step()
{
    std::vector<Eigen::Matrix3d> temps;
    TriMesh::VertexIter v_it = ref_mesh_.vertices_begin();
    for(;v_it!=ref_mesh_.vertices_end();v_it++)
    {
        TriMesh::VertexEdgeIter ve_it = ref_mesh_.ve_iter(*v_it);
        Eigen::Matrix3d temp;
        temp.setZero();
        for(;ve_it.is_valid();ve_it++)
        {
            TriMesh::VertexHandle vk = ref_mesh_.to_vertex_handle(ref_mesh_.halfedge_handle(*ve_it,0));
            if(vk.idx()==(*v_it).idx())
                vk = ref_mesh_.from_vertex_handle(ref_mesh_.halfedge_handle(*ve_it,0));
            TriMesh::Point pj,pk;
            pj = ref_mesh_.point(*v_it);
            pk = ref_mesh_.point(vk);
            Eigen::Vector3d ejk(pj[0]-pk[0],pj[1]-pk[1],pj[2]-pk[2]);
            Eigen::Vector3d ejk_ = P_.block<3,1>(3*(*v_it).idx(),0)-P_.block<3,1>(3*vk.idx(),0);
            double cjk = ref_mesh_.property(LB_weights,*ve_it);
            temp+=cjk*ejk*ejk_.transpose();
        }
        temps.push_back(temp);
    }
    v_it = ref_mesh_.vertices_begin();
    for(;v_it!=ref_mesh_.vertices_end();v_it++)
    {
        TriMesh::VertexEdgeIter ve_it = ref_mesh_.ve_iter(*v_it);
        Eigen::Matrix3d temp;
        temp.setZero();
        for(;ve_it.is_valid();ve_it++)
        {
            TriMesh::HalfedgeHandle h_e = ref_mesh_.halfedge_handle(*ve_it,0);
            TriMesh::VertexHandle vj = ref_mesh_.to_vertex_handle(h_e);
            if(vj.idx()==(*v_it).idx())
            {
                vj = ref_mesh_.from_vertex_handle(h_e);
                h_e = ref_mesh_.opposite_halfedge_handle(h_e);
            }
            double cj = 1.0/double(ref_mesh_.valence(vj));
            temp+=cj*rotation_log_exp::exp(ref_mesh_.property(log_dRs,h_e))*ref_mesh_.property(scaling_matrixs,vj)*temps[vj.idx()];
        }
        ref_mesh_.property(buffer_for_compute,(*v_it))=temp;
    }
}

void RIMD_Reconstruction::compute_internal_for_globalstep()
{
    TriMesh::VertexIter v_it = ref_mesh_.vertices_begin();
    for(;v_it!=ref_mesh_.vertices_end();v_it++)
    {
        TriMesh::VertexEdgeIter ve_it = ref_mesh_.ve_iter(*v_it);
        Eigen::Matrix3d temp;
        temp.setZero();
        for(;ve_it.is_valid();ve_it++)
        {
            TriMesh::HalfedgeHandle h_e = ref_mesh_.halfedge_handle(*ve_it,0);
            if(ref_mesh_.from_vertex_handle(h_e).idx()==(*v_it).idx())
                h_e = ref_mesh_.opposite_halfedge_handle(h_e);
            TriMesh::VertexHandle v0,v1;
            v0 = ref_mesh_.from_vertex_handle(h_e);
            v1 = ref_mesh_.to_vertex_handle(h_e);
            temp+=ref_mesh_.property(rotation_matrixs,v0)*rotation_log_exp::exp(ref_mesh_.property(log_dRs,h_e))*ref_mesh_.property(scaling_matrixs,v1);
        }
//        if(_isnan(temp(0,0)))
//            std::cout<<(*v_it).idx()<<std::endl;
        ref_mesh_.property(buffer_for_compute,*v_it)=temp;
    }
}

void RIMD_Reconstruction::initial_P_to_ref_mesh()
{
    TriMesh::VertexIter v_it = ref_mesh_.vertices_begin();
    P_.resize(3*ref_mesh_.n_vertices());
    for(;v_it!=ref_mesh_.vertices_end();v_it++)
    {
        TriMesh::Point p = ref_mesh_.point(*v_it);
        int id = (*v_it).idx();
        for(int i=0;i<3;i++)
            P_[3*id+i]=p[i];
    }
}
//先暂时用边的遍历的方式来初始化，跟用顶点还是稍有不同的
void RIMD_Reconstruction::initial_Rs_for_reconstruction()
{
    std::vector<bool> is_edge_in_queue;
    is_edge_in_queue.resize(ref_mesh_.n_edges(),false);
    int start_id = ref_mesh_.n_vertices()/2;
    TriMesh::VertexHandle v_start(start_id);
    TriMesh::VertexEdgeIter ve_it = ref_mesh_.ve_iter(v_start);
//    Eigen::Matrix3d R_start;
//    R_start<<0,1,0,-1,0,0,0,0,1;
    ref_mesh_.property(rotation_matrixs,v_start)=Eigen::Matrix3d::Identity()/*R_start*/;
    std::queue<TriMesh::EdgeHandle> edges_q;
    std::queue<int> cid_q;
    for(;ve_it.is_valid();ve_it++)
    {
        edges_q.push(*ve_it);
        is_edge_in_queue[(*ve_it).idx()]=true;
        cid_q.push(start_id);
    }
    while(!edges_q.empty())
    {
        TriMesh::EdgeHandle e_h = edges_q.front();edges_q.pop();
        int v_id = cid_q.front(); cid_q.pop();
        TriMesh::HalfedgeHandle h_e = ref_mesh_.halfedge_handle(e_h,0);
        if(ref_mesh_.from_vertex_handle(h_e).idx()!=v_id)
            h_e = ref_mesh_.halfedge_handle(e_h,1);
        TriMesh::VertexHandle vf,vt;
        vf = ref_mesh_.from_vertex_handle(h_e);
        vt = ref_mesh_.to_vertex_handle(h_e);
        ref_mesh_.property(rotation_matrixs,vt)
                = ref_mesh_.property(rotation_matrixs,vf)*rotation_log_exp::exp(ref_mesh_.property(log_dRs,h_e));

        ve_it = ref_mesh_.ve_iter(vt);
        for(;ve_it.is_valid();ve_it++)
        {
            if(!is_edge_in_queue[(*ve_it).idx()])
            {
                edges_q.push(*ve_it);
                is_edge_in_queue[(*ve_it).idx()]=true;
                cid_q.push(vt.idx());
            }
        }
    }
}

double RIMD_Reconstruction::compute_Reconstruction_energy()
{
    TriMesh::VertexIter v_it = ref_mesh_.vertices_begin();
    double sum = 0.0;
    for(;v_it!=ref_mesh_.vertices_end();v_it++)
    {
        TriMesh::VertexEdgeIter vej_it = ref_mesh_.ve_iter(*v_it);
        TriMesh::VertexHandle vi = *v_it;
        double sub1_sum = 0.0;
        for(;vej_it.is_valid();vej_it++)
        {
            TriMesh::HalfedgeHandle h_eij=ref_mesh_.halfedge_handle(*vej_it,0);
            TriMesh::VertexHandle vj = ref_mesh_.to_vertex_handle(h_eij);
            if(vj.idx() == vi.idx())
            {
                vj = ref_mesh_.from_vertex_handle(ref_mesh_.halfedge_handle(*vej_it,0));
                h_eij = ref_mesh_.opposite_halfedge_handle(h_eij);
            }
            double cj = 1.0/double(ref_mesh_.valence(vj));
            TriMesh::VertexEdgeIter vek_it = ref_mesh_.ve_iter(vj);
            double sub2_sum = 0.0;
            for(;vek_it.is_valid();vek_it++)
            {
                TriMesh::VertexHandle vk = ref_mesh_.to_vertex_handle(ref_mesh_.halfedge_handle(*vek_it,0));
                if(vk.idx()==vj.idx())
                    vk = ref_mesh_.from_vertex_handle(ref_mesh_.halfedge_handle(*vek_it,0));
                double cjk = ref_mesh_.property(LB_weights,*vek_it);
                Eigen::Vector3d e1jk;
                e1jk = P_.block<3,1>(3*vj.idx(),0)-P_.block<3,1>(3*vk.idx(),0);
                TriMesh::Point pj,pk;
                pj = ref_mesh_.point(vj);
                pk = ref_mesh_.point(vk);
                Eigen::Vector3d e0jk(pj[0]-pk[0],pj[1]-pk[1],pj[2]-pk[2]);
                Eigen::Vector3d temp = e1jk -
                        ref_mesh_.property(rotation_matrixs,vi)*rotation_log_exp::exp(ref_mesh_.property(log_dRs,h_eij))*
                        ref_mesh_.property(scaling_matrixs,vj)*e0jk;
                sub2_sum += cjk*(temp*temp.transpose())(0,0);
            }
            sub1_sum += cj*sub2_sum;
        }
        sum += sub1_sum;
    }
    return sum;
}

void RIMD_Reconstruction::check_RIMD_correct()
{
    TriMesh::EdgeIter e_it = ref_mesh_.edges_begin();
    for(;e_it!=ref_mesh_.edges_end();e_it++)
    {
        TriMesh::HalfedgeHandle h0,h1;
        h0 = ref_mesh_.halfedge_handle(*e_it,0);
        h1 = ref_mesh_.halfedge_handle(*e_it,1);
        Eigen::Matrix3d logdR0,logdR1;
        logdR0 = ref_mesh_.property(log_dRs,h0);
        logdR1 = ref_mesh_.property(log_dRs,h1);
        double angle0,angle1;
        Eigen::Vector3d axis0,axis1;
//        rotation_log_exp::so3_to_angle_axis(logdR0,angle0,axis0);
//        rotation_log_exp::so3_to_angle_axis(logdR1,angle1,axis1);
//        if((axis0+axis1).norm()>0.00001)
//        {
//            if(fabs(angle0-M_PI)>0.00001||fabs(angle1-M_PI)>0.00001)
//                std::cout<<(*e_it).idx()<<" edge logdR are not opposite!!!"<<std::endl;
//        }
//        if(fabs(angle0-angle1)>0.00001)
//            std::cout<<(*e_it).idx()<<" edge logdR angles are not equal!!!"<<std::endl;
        if((logdR0+logdR1).norm()>0.00001)
            std::cout<<(*e_it).idx()<<" edge logdR are not compatible!!!"<<std::endl;



//        if((*e_it).idx()==77076||(*e_it).idx()==193268)
//        {
//            std::cout<<angle0<<axis0<<std::endl;
//            std::cout<<angle1<<axis1<<std::endl;
//        }
    }
}
