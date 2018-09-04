#include "util_3drotation_log_exp.h"
#include <iostream>
#include "rimd_reconstruction.h"
#include <queue>
int main(int argc, char *argv[])
{
    RIMD_Reconstruction rimd_r;
//    rimd_r.read_ref_mesh_from_file("Debug_ref.obj");
//    rimd_r.read_defor_mesh_from_file("Debug_def.obj");
//    rimd_r.read_anchor_points_id("");
//    rimd_r.read_anchor_points_id("debug_anchors.txt");

//    rimd_r.read_ref_mesh_from_file("face0_triangle_version.obj");
//    rimd_r.read_defor_mesh_from_file("face22_triangle_version.obj");
//    rimd_r.read_anchor_points_id("face_anchors.txt");
//    rimd_r.read_anchor_points_id("");

//    rimd_r.read_ref_mesh_from_file("no_detail.off");
//    rimd_r.read_defor_mesh_from_file("detail.off");
//    rimd_r.read_anchor_points_id("");

    rimd_r.read_ref_mesh_from_file("/home/chern/Project/RIMD/RIMD_Reconstruct/shape_0.obj"); // read ref_mesh
    rimd_r.read_defor_mesh_from_file("/home/chern/Project/RIMD/RIMD_Reconstruct/shape_29.obj");  // read defor_mesh
//    rimd_r.read_anchor_points_id("");  // fix some points
    rimd_r.read_anchor_points_id("/home/chern/Project/RIMD/RIMD_Reconstruct/one_anchor.txt");

    rimd_r.Preprocess();
    //rimd_r.InterlateRIMD(0.9);  // get new RIMD;
    rimd_r.LoadRIMD("/home/chern/Project/RIMD/RIMD_Reconstruct/shape_29.dat");
    rimd_r.Reconstruction();
    TriMesh mesh;
    rimd_r.GetReconstructionMesh(mesh);
    OpenMesh::IO::write_mesh(mesh,"defor_mesh.obj");
    std::cout<<"done"<<std::endl;

//    TriMesh mesh;
//    OpenMesh::IO::read_mesh(mesh,"large_detail.obj");
//    TriMesh::VertexIter v_it = mesh.vertices_begin();
//    FILE *fin = fopen("detail_anchors.txt","w");
//    for(;v_it!=mesh.vertices_end();v_it++)
//    {
//        if(mesh.is_boundary(*v_it))
//        {
//            fprintf(fin,"%d\n",(*v_it).idx());
//        }
//    }
//    fclose(fin);
//    std::cout<<"done"<<std::endl;
//    OpenMesh::IO::write_mesh(mesh,"large_detail.obj");
}
