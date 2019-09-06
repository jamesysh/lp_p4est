#include<iostream>
#include "initializer.h"
#include "mpi.h"
#include "particle_data.h"
#include "octree_manager.h"
#include "lp_solver.h"
using namespace std;

static  int
refine_init (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
    if(quadrant->level>4)

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


    octree->partition_octree(0,NULL);
 
    gdata->initFluidParticles();

    LPSolver * lpsolver = new LPSolver(gdata);

    gdata->writeVTKFiles();
    
    lpsolver->moveParticlesByG(lpsolver->dt);


    
    gdata->cleanUpArrays(); 
    
    octree->destroy_octree();
   
    
    delete gdata;
    delete octree;
    delete init;

    mpiret = sc_MPI_Finalize ();
    return 0;
}
