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
                            gdata->initlevel , 1,
                            sizeof (octant_data_t), NULL, gdata);

}



void Octree_Manager:: destroy_octree(){


    p4est_destroy (gdata->p4est);
  
    gdata->p4est = NULL;
  
    p4est_connectivity_destroy (gdata->conn);
  
    gdata->conn = NULL;
  
    sc_finalize ();

}


