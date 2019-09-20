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



int main(){
    Initializer *init = new Initializer();

    Global_Data *gdata = new Global_Data(init);

    ParticleViewer * viewer = new ParticleViewer(gdata,"output_",4);
    int mpiret;
    mpiret = sc_MPI_Init (NULL, NULL);

    SC_CHECK_MPI (mpiret);

    Octree_Manager *octree = new Octree_Manager(gdata);

    octree->build_octree();

 
    octree->refine_octree(1,refine_init,NULL,NULL);  //initial refinement of octree


    octree->partition_octree(1,NULL);
    gdata->prerun(); 
    gdata->initFluidParticles();
    gdata->boundary->generateBoundaryParticle(gdata,gdata->eos,gdata->initlocalspacing);
    
    octree->adapt_octree(); 
    
    
    sc_array_destroy(gdata->irecumu);
    sc_array_destroy(gdata->irvcumu);
    gdata->resetOctantData(); 
   sc_array_destroy(gdata->ireceive);
    sc_array_destroy(gdata->iremain);
    LPSolver * lpsolver = new LPSolver(gdata);

     
    double tstart = 0;
    double tend = 0.01;
    double nextwritetime = 0;
    while(tstart<tend)
    {
    tstart += lpsolver->dt;
    
    
    
    lpsolver->moveParticlesByG(lpsolver->dt);
    
    gdata->boundary->UpdateInflowBoundary(gdata,gdata->eos,lpsolver->dt,gdata->initlocalspacing);
    
    gdata->presearch();
        
    gdata->packParticles();
    
    if(gdata->gpnum == 0)
        
    {
      sc_array_destroy_null (&gdata->recevs);
      sc_hash_destroy_null (&gdata->psend);
      sc_array_destroy_null(&gdata->iremain);

      gdata->psend = NULL;
      sc_mempool_destroy (gdata->psmem);
      gdata->psmem = NULL;
        break;
    }
    gdata->communicateParticles();
    gdata->postsearch();
    octree->adapt_octree(); 

    octree->balance_octree(NULL,octree->balance_replace);
    gdata->regroupParticles(); 
 
    
    gdata->createViewForOctant();
    

    gdata->searchNeighbourOctant();
    

    gdata->searchNeighbourParticle();
    
    gdata->searchUpwindNeighbourParticle(); 
    gdata->generateGhostParticle();
   // gdata->testquad(); 
   // lpsolver->solve_upwind(0);
    if(tstart  >= nextwritetime)
    
    {
        nextwritetime += 0.001;    
        viewer->writeResult(tstart);
        //viewer->writeGhost(tstart);
    }
   
    gdata->cleanForTimeStep();
   
    gdata->partitionParticles();
    
    }
    
    
    gdata->cleanUpArrays(); 
    
    octree->destroy_octree();
   
    
    delete gdata;
    delete octree;
    delete lpsolver;
    delete init;

    mpiret = sc_MPI_Finalize ();
    return 0;
}
