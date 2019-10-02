#include<iostream>
#include "initializer.h"
#include "mpi.h"
#include "particle_data.h"
#include "octree_manager.h"
#include "lp_solver.h"
#include "particle_viewer.h"
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


static  int
refine_init2d (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
    Global_Data *g = (Global_Data *) p4est->user_pointer;
    int initlevel = g->initlevel;
    if(quadrant->level >= initlevel)

        return 0;
    else 
        return 1;
}



int main(){
    double t1 = MPI_Wtime();
    Initializer *init = new Initializer();

    Global_Data *gdata = new Global_Data(init);

    ParticleViewer * viewer = new ParticleViewer(gdata,"output_",4);
    int mpiret;
    mpiret = sc_MPI_Init (NULL, NULL);

    SC_CHECK_MPI (mpiret);

    Octree_Manager *octree = new Octree_Manager(gdata);

    octree->build_octree();

    if(gdata->dimension == 3) 
        octree->refine_octree(1,refine_init,NULL,NULL);  //initial refinement of octree
    else if(gdata->dimension == 2)
        octree->refine_octree2d(1,refine_init2d,NULL,NULL);  //initial refinement of octree

    octree->partition_octree(1);
    gdata->prerun(); 
    gdata->initFluidParticles_distributed();
  //  gdata->initFluidParticles_hexagonal();
    if(gdata->dimension == 3)
       gdata->resetOctantData(); 
    else if(gdata->dimension == 2)
        gdata->resetOctantData2d();
    LPSolver * lpsolver = new LPSolver(gdata,octree,viewer);
    if(gdata->dimension == 3 )
        lpsolver->solve_3d();
    else if(gdata->dimension == 2)
        lpsolver->solve_2d();
    
    
    gdata->cleanUpArrays(); 
    
    octree->destroy_octree();
   
    
    delete gdata;
    delete octree;
    delete lpsolver;
    delete init;

    mpiret = sc_MPI_Finalize ();

    double t2 = MPI_Wtime();
    printf( "Elapsed time is %f\n", t2 - t1 ); 
    return 0;
}
