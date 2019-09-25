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
    gdata->initFluidParticles();
    //gdata->boundary->generateBoundaryParticle(gdata,gdata->eos,gdata->initlocalspacing);
    if(gdata->dimension == 3 ) 
        octree->adapt_octree(); 
    else if(gdata->dimension == 2)
        octree->adapt_octree2d();
    
    sc_array_destroy(gdata->irecumu);
    sc_array_destroy(gdata->irvcumu);
    if(gdata->dimension == 3)
       gdata->resetOctantData(); 
    else if(gdata->dimension == 2)
        gdata->resetOctantData2d();
   sc_array_destroy(gdata->ireceive);
    sc_array_destroy(gdata->iremain);
    LPSolver * lpsolver = new LPSolver(gdata);
    
     
    viewer->writeResult(0);
    double tstart = 0;
    double tend = 0.0005;
    double nextwritetime = 0;
    while(tstart<tend)
    {
    tstart += lpsolver->cfldt;
    
    
    
    
    //gdata->boundary->UpdateInflowBoundary(gdata,gdata->eos,lpsolver->dt,gdata->initlocalspacing);
    if(gdata->dimension == 3) 
        gdata->presearch();
    else if(gdata->dimension == 2)
        gdata->presearch2d();
        
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
    if(gdata->dimension == 3)
        gdata->postsearch();
    else if(gdata->dimension == 2)
        gdata->postsearch2d();

    if(gdata->dimension == 3 ) 
        octree->adapt_octree(); 
    else if(gdata->dimension == 2)
        octree->adapt_octree2d();
    
    if(gdata->dimension == 3)
        octree->balance_octree(NULL,octree->balance_replace);
    else if(gdata->dimension == 2)
        octree->balance_octree2d(NULL,octree->balance_replace2d);
    
    if(gdata->dimension == 3)
        gdata->regroupParticles(); 
    else if(gdata->dimension == 2) 
        gdata->regroupParticles2d(); 
    
    if(gdata->dimension == 3)
        gdata->partitionParticles();
    else if(gdata->dimension == 2)
        gdata->partitionParticles2d();
    
    if(gdata->dimension == 3)
        gdata->createViewForOctant();
    else if(gdata->dimension == 2)
        gdata->createViewForOctant2d();
    
    
    if(gdata->dimension == 3)
        gdata->searchNeighbourOctant();
    else if(gdata->dimension == 2)
        gdata->searchNeighbourOctant2d();


    if(gdata->dimension == 3)
        gdata->searchNeighbourParticle();
    else if(gdata->dimension == 2)
        gdata->searchNeighbourParticle2d();
    if(gdata->dimension == 3)
        gdata->searchUpwindNeighbourParticle(); 
    else if(gdata->dimension == 2)
        gdata->searchUpwindNeighbourParticle2d(); 
    
    if(gdata->dimension == 3)
        gdata->generateGhostParticle();
    else if(gdata->dimension == 2)
        gdata->generateGhostParticle2d();

        // gdata->testquad();
    lpsolver->computeCFLCondition();
    
    for(int phase = 0;phase< lpsolver->totalphase;phase++){
        lpsolver->solve_upwind(phase);
        MPI_Barrier(gdata->mpicomm);
        if(gdata->dimension == 3)
            gdata->updateViewForOctant(phase);
        else if(gdata->dimension == 2)
            gdata->updateViewForOctant2d(phase);
        MPI_Barrier(gdata->mpicomm); 
    }
    gdata->updateParticleStates();
   
    lpsolver->updateLocalSpacing();
     
    lpsolver->moveParticle();
    if(tstart  >= nextwritetime)
    
    {
        nextwritetime += lpsolver->cfldt;    
        viewer->writeResult(tstart);
    //    viewer->writeGhost(tstart);
    }
   
    MPI_Barrier(gdata->mpicomm); 
    
    if(gdata->dimension == 3)
        gdata->cleanForTimeStep();
    else if(gdata->dimension == 2)
        gdata->cleanForTimeStep2d();
    
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
