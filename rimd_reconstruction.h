#ifndef RIMD_RECONSTRUCTION_H
#define RIMD_RECONSTRUCTION_H
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "trimesh_types_hh.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
class RIMD_Reconstruction
{
public:
    RIMD_Reconstruction();
    void read_ref_mesh_from_file(const char* _filename);
    void read_defor_mesh_from_file(const char* _filename);
    void read_anchor_points_id(const char* _filename);
    void Preprocess();
    void Reconstruction();
    // (1-t)ref + t defor
    void InterlateRIMD(double t);
    // Load RIMD feature from file
    void LoadRIMD(std::string _filename);
    void GetReconstructionMesh(TriMesh &mesh);
private:
    void compute_RIMD_of_ref_to_defor();
    void compute_ref_LB_weights();
    void compute_ref_to_defor_Tmatrixs();
    void compute_rotation_scaling_matrixs();
    void compute_logdR_matrixs();

    void compute_Ti(TriMesh::VertexHandle v_it,TriMesh::VertexHandle v_to_it);

    void compute_A_for_globalstep();
    void compute_b_for_globalstep();

    void compute_local_step();
    void compute_Q_for_local_step();

    void compute_internal_for_globalstep();

    void initial_P_to_ref_mesh();
    void initial_Rs_for_reconstruction();

    double compute_Reconstruction_energy();

    void check_RIMD_correct();
private:
    //ref and defor must have the same topology
    TriMesh ref_mesh_;
    TriMesh defor_mesh_;
    OpenMesh::VPropHandleT<Eigen::Matrix3d> rotation_matrixs;
    OpenMesh::VPropHandleT<Eigen::Matrix3d> scaling_matrixs;
    OpenMesh::VPropHandleT<Eigen::Matrix3d> T_matrixs;
    OpenMesh::EPropHandleT<double>  LB_weights;
    OpenMesh::HPropHandleT<Eigen::Matrix3d> log_dRs;
    OpenMesh::VPropHandleT<Eigen::Matrix3d> buffer_for_compute;
    Eigen::SparseMatrix<double> A;
    Eigen::SparseLU<Eigen::SparseMatrix<double> > A_solver_;
    Eigen::VectorXd b_;

    Eigen::VectorXd P_;

    std::vector<bool> is_anchor;
};

#endif // RIMD_RECONSTRUCTION_H
