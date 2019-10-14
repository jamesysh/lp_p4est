#include "pellet_solver.h"


PelletSolver::PelletSolver(Global_Data*g){
    gdata = g;
    sc_array_copy(particle_data_copy, gdata->particle_data);

    }


void PelletSolver::build_quadtree(){
    conn = p4est_connectivity_new_unitsquare();
    p4est_heating = p4est_new_ext(gdata->mpicomm, conn,0,0,1,sizeof(octant_data_t), NULL, this);  
    
    }


void PelletSolver::prerun(){
    
   for (int i = 0; i < 2; ++i) {
    ilh[i] = sc_array_new (sizeof (p4est_locidx_t));
    jlh[i] = sc_array_new (sizeof (p4est_locidx_t));
    klh[i] = NULL;
   
   }
    
    }


void PelletSolver::resetOctantData2d(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est_heating->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qud->premain = qud->preceive = 0;
    }
  }

}
