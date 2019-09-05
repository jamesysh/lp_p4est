#include<iostream>

#include "p4est_to_p8est.h"
#include "mpi.h"
#include "particle_data.h"
#include "initializer.h"
#include "octree_manager.h"
using namespace std;

int main(){

    Initializer initializer = Initializer();
    Initializer *init=&initializer;

    Global_Data global_data = Global_Data(init);

    Global_Data *gdata = &global_data;

    int mpiret;
    mpiret = sc_MPI_Init (NULL, NULL);

    SC_CHECK_MPI (mpiret);


    Octree_Manager octree_manager = Octree_Manager(gdata);
    Octree_Manager *octree = &octree_manager;

    octree->build_octree();


p4est_quadrant_t   *quad;

   p4est_topidx_t      tt;
  p4est_locidx_t      lq;
  p4est_tree_t       *tree;
  double len; 
  for (tt = gdata->p4est->first_local_tree; tt <= gdata->p4est->last_local_tree;
         ++tt) {
      tree = p4est_tree_array_index (gdata->p4est->trees, tt);
      for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
 
        quad = p4est_quadrant_array_index (&tree->quadrants, lq);
       p4est_qcoord_t      qh;

        qh = P4EST_QUADRANT_LEN (quad->level);
         len = (double)qh/P4EST_ROOT_LEN*gdata->domain_len;
        //   cout<<qh<<endl;
      } 
    }
  

  octree->destroy_octree();
    mpiret = sc_MPI_Finalize ();
    return 0;
}
