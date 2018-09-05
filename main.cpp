#include "util_3drotation_log_exp.h"
#include <iostream>
#include "rimd_reconstruction.h"
#include <queue>
#include <fstream>
#include <string>

int main(int argc, char *argv[])
{
    RIMD_Reconstruction rimd_r;


    /*
     * read_ref_mesh_from_file()
     * read_defor_mesh_from_file()
     * read_anchor_points_id()
     * Preprocess() //Computes LB weights, solve svd decompositions, etc.
     * Interlate(double, string) //Generate new RIMD feature(linear combination) and store those features to a binary file.
     * LoadRIMD(string) //Load RIMD feature to compute the reconstructed mesh.
     * Reconstruction() //Compute reconstructed mesh via RIMD feature.
     * GetReconstructionMesh()
     * write_mesh(mesh, string)
     */
    /*
    rimd_r.read_ref_mesh_from_file("/home/chern/Data/Dog_new_expressions/1_new.obj");
    rimd_r.read_defor_mesh_from_file("/home/chern/Data/Dog_new_expressions/12_new.obj");
    rimd_r.read_anchor_points_id("/home/chern/Project/RIMD/RIMD_Reconstruct/one_anchor.txt");
    rimd_r.Preprocess();
    //rimd_r.InterlateRIMD(1.0,"/home/chern/Data/Dog_new_expressions/12_new.dat");  // get new RIMD;
    rimd_r.LoadRIMD("/home/chern/Data/Dog_new_expressions/12_new.dat");
    rimd_r.Reconstruction();
    TriMesh mesh;
    rimd_r.GetReconstructionMesh(mesh);
    OpenMesh::IO::write_mesh(mesh,"/home/chern/Data/Dog_new_expressions/12_reconstruction.obj");
    std::cout<<"done"<<std::endl;
    */

    std::string ref_mesh_name="/home/chern/Data/Dog_new_expressions/1_new.obj";
    std::string def_mesh_head="/home/chern/Data/Dog_new_expressions/";
    std::string def_mesh_tail="_new.obj";
    std::string data_tail=".dat";
    for(int i=1;i<48;i++){
        std::ostringstream oss_def_mesh_name;
        std::ostringstream oss_data_name;

        oss_def_mesh_name<<def_mesh_head<<i<<def_mesh_tail;
        std::string def_mesh_name=oss_def_mesh_name.str();

        oss_data_name<<def_mesh_head<<i<<data_tail;
        std::string data_name=oss_data_name.str();

        rimd_r.read_ref_mesh_from_file(ref_mesh_name);
        rimd_r.read_defor_mesh_from_file(def_mesh_name);
        rimd_r.read_anchor_points_id("/home/chern/Project/RIMD/RIMD_Reconstruct/one_anchor.txt");
        rimd_r.Preprocess();
        rimd_r.InterlateRIMD(1.0, data_name);

    }

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
