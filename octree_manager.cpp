#include <iostream>
#include "octree_manager.h"

Octree_Manager:: Octree_Manager(Global_Data *g){

    gdata = g;

}


void Octree_Manager:: build_octree(){

    int mpiret;
    gdata->mpicomm = sc_MPI_COMM_WORLD;

    mpiret = sc_MPI_Comm_size (gdata->mpicomm, &gdata->mpisize);

    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Comm_rank (gdata->mpicomm, &gdata->mpirank);

    SC_CHECK_MPI (mpiret);
  
    
    sc_init (gdata->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  
    p4est_init (NULL, SC_LP_DEFAULT);

//build octree connectivity
    gdata->conn = p8est_connectivity_new_unitcube ();
//build octree
    gdata->p4est = p8est_new_ext (gdata->mpicomm, gdata->conn, 0,
                            0 , 1,
                            sizeof (octant_data_t), NULL, gdata);

}



void Octree_Manager:: destroy_octree(){


    p4est_destroy (gdata->p4est);
  
    gdata->p4est = NULL;
  
    p4est_connectivity_destroy (gdata->conn);
  
    gdata->conn = NULL;
  
    sc_finalize ();

}


void Octree_Manager:: partition_octree(int allow_for_coarsening,p4est_weight_t weight_fn){

   p4est_partition(gdata->p4est,allow_for_coarsening,weight_fn); 

}

        
void Octree_Manager:: refine_octree(int recursive, p4est_refine_t refine_fn, p4est_init_t init_fn, p4est_replace_t replace_fn){
   
    p4est_refine_ext(gdata->p4est,recursive,-1,refine_fn,init_fn,replace_fn);

}


