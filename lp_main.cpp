#include<iostream>
#include "initializer.h"
#include "mpi.h"
#include "particle_data.h"
#include "octree_manager.h"
#include "lp_solver.h"
using namespace std;

static  int
refine_init (p8est_t * p8est, p4est_topidx_t which_tree,
           p8est_quadrant_t * quadrant)
{
    Global_Data *g = (Global_Data *) p8est->user_pointer;
    int initlevel = g->initlevel;
    if(quadrant->level >= initlevel)

        return 0;
    else 
        return 1;
}



int main(){

    Initializer *init = new Initializer();

    Global_Data *gdata = new Global_Data(init);


    int mpiret;
    mpiret = sc_MPI_Init (NULL, NULL);

    SC_CHECK_MPI (mpiret);

    Octree_Manager *octree = new Octree_Manager(gdata);

    octree->build_octree();

 
    octree->refine_octree(1,refine_init,NULL,NULL);  //initial refinement of octree


    octree->partition_octree(1,NULL);
 
    gdata->initFluidParticles();

    
    gdata->ireindex = gdata->irvindex = 0;
    p8est_coarsen_ext (gdata->p8est, 0, 1, octree->adapt_coarsen, NULL, octree->adapt_replace);
    LPSolver * lpsolver = new LPSolver(gdata);

  //  gdata->writeVTKFiles();
    
   // lpsolver->moveParticlesByG(lpsolver->dt);
    
    gdata->cleanUpArrays(); 
    
    octree->destroy_octree();
   
    
    delete gdata;
    delete octree;
    delete lpsolver;
    delete init;

    mpiret = sc_MPI_Finalize ();
    return 0;
}
