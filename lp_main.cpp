#include<iostream>
#include "tool_fn.h"
#include "initializer.h"
#include "mpi.h"
#include "particle_data.h"
#include "octree_manager.h"
using namespace std;

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

   
    //gdata->writeVTKFiles();



    
    gdata->cleanUpArrays(); 
    
    octree->destroy_octree();
   
    
    delete gdata;
    delete octree;
    delete init;

    mpiret = sc_MPI_Finalize ();
    return 0;
}
