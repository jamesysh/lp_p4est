#include <iostream>
#include "octree_manager.h"
using namespace std;
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
    gdata->p8est = p8est_new_ext (gdata->mpicomm, gdata->conn, 0,
                            0 , 1,
                            sizeof (octant_data_t), NULL, gdata);

}



void Octree_Manager:: destroy_octree(){


    p8est_destroy (gdata->p8est);
  
    gdata->p8est = NULL;
  
    p8est_connectivity_destroy (gdata->conn);
  
    gdata->conn = NULL;
  
    sc_finalize ();

}


void Octree_Manager:: partition_octree(int allow_for_coarsening,p8est_weight_t weight_fn){

   p8est_partition(gdata->p8est,allow_for_coarsening,weight_fn); 

}

        
void Octree_Manager:: refine_octree(int recursive, p8est_refine_t refine_fn, p8est_init_t init_fn, p8est_replace_t replace_fn){
   
    p8est_refine_ext(gdata->p8est,recursive,-1,refine_fn,init_fn,replace_fn);

}

int Octree_Manager:: adapt_coarsen (p8est_t * p8est, p4est_topidx_t which_tree,
               p8est_quadrant_t * quadrants[])
{
    int i;
    p4est_locidx_t remain = 0, receive = 0;
    octant_data_t *oud;
    Global_Data *g = (Global_Data *)p8est->user_pointer;
    
    if (quadrants[1] == NULL ||
      quadrants[0]->level == g->minlevel) {
    
        oud = (octant_data_t *) quadrants[0]->p.user_data;
        g->ireindex += oud->premain;
        g->irvindex += oud->preceive;
    
        return 0;
  }
   
    remain = receive = 0;
  
    for (i = 0; i < P8EST_CHILDREN; ++i) {
        oud = (octant_data_t *) quadrants[i]->p.user_data;
        remain += oud->premain;
        receive += oud->preceive;
  }
  if ((double) (remain + receive) < .5 * g->elem_particles) {
    /* we will coarsen and adjust ireindex, irvindex in adapt_replace */
    g->qremain = remain;
    g->qreceive = receive;
    return 1;
  }
  else {
    /* we will not coarsen and proceed with next quadrant */
    oud = (octant_data_t *) quadrants[0]->p.user_data;
    g->ireindex += oud->premain;
    g->irvindex += oud->preceive;
    return 0;
  }


}

void Octree_Manager:: adapt_replace (p8est_t * p8est, p4est_topidx_t which_tree,
               int num_outgoing, p8est_quadrant_t * outgoing[],
               int num_incoming, p8est_quadrant_t * incoming[]){

  int                 wx, wy, wz;
  double              lxyz[3], hxyz[3], dxyz[3];
  sc_array_t          iview, *arr;
  p4est_locidx_t      irem, ibeg;
  p8est_quadrant_t  **pchild;
  octant_data_t          *oud;
  Global_Data      *g = (Global_Data *) p8est->user_pointer;
  if (num_outgoing == P8EST_CHILDREN) {
    // we are coarsening 
    
      oud = (octant_data_t *) incoming[0]->p.user_data;
    g->ireindex += (oud->premain = g->qremain);
    g->irvindex += (oud->preceive = g->qreceive);
  }

  else {
    // we are refining 
    // access parent quadrant 
    g->loopquad (which_tree, outgoing[0],lxyz,hxyz,dxyz);

    // recover window onto remaining particles for the new family 
    ibeg = g->ire2;
    irem = g->ireindex - ibeg;
    sc_array_init_view (&iview, g->iremain, ibeg, irem);

    // sort remaining particles into the children 
    pchild = incoming;
    g->split_by_coord ( &iview, g->klh, PA_MODE_REMAIN, 2, lxyz, dxyz);
    
    for (wz = 0; wz < 2; ++wz) {
    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord (g->jlh[wy], g->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (octant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
      }
    }
    }

    // recover window onto received particles for the new family 
    ibeg = g->irv2;
    irem = g->irvindex - ibeg;
    P4EST_ASSERT (irem >= 0);
    sc_array_init_view (&iview, g->ireceive, ibeg, irem);
    P4EST_ASSERT (qod->preceive == irem);

    // sort received particles into the children 
    pchild = incoming;
    g->split_by_coord ( &iview, g->klh, PA_MODE_RECEIVE, 2, lxyz, dxyz);
    for (wz = 0; wz < 2; ++wz) {
    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_RECEIVE, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord ( g->jlh[wy], g->ilh, PA_MODE_RECEIVE, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->ireceive, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (octant_data_t *) (*pchild++)->p.user_data;
        P4EST_ASSERT (oud->u.lpend == -1);
        ibeg += (oud->preceive = (p4est_locidx_t) arr->elem_count);
      }
    }
    }
    P4EST_ASSERT (ibeg == g->irvindex);
    P4EST_ASSERT (pchild == incoming + P4EST_CHILDREN);

  }

}


int Octree_Manager:: adapt_refine (p8est_t * p8est, p4est_topidx_t which_tree,
              p8est_quadrant_t * quadrant)
{
 
  Global_Data      *g = (Global_Data *) p8est->user_pointer;


  
  octant_data_t          *oud = (octant_data_t *) quadrant->p.user_data;

  /* we have set this to -1 in adapt_coarsen */

  if ((double) (oud->premain + oud->preceive) > g->elem_particles) {
    /* we are trying to refine, we will possibly go into the replace function */
    g->ire2 = g->ireindex;
    g->ireindex += oud->premain;
    g->irv2 = g->irvindex;
    g->irvindex += oud->preceive;
    return 1;
  }
  else {
    /* maintain cumulative particle count for next quadrant */
    g->ireindex += oud->premain;
    g->irvindex += oud->preceive;
    return 0;
  }
}

void Octree_Manager:: adapt_octree(){

    p8est_t *p8est = gdata->p8est;
    gdata->ireindex = gdata->irvindex = 0;
    p8est_coarsen_ext (p8est, 0, 1, adapt_coarsen, NULL, adapt_replace);
    
    
    gdata->ireindex = gdata->ire2 = 0;
    gdata->irvindex = gdata->irv2 = 0;
    p8est_refine_ext (p8est, 0, gdata->maxlevel, adapt_refine, NULL, adapt_replace);
    
}





