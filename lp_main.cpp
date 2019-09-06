#include<iostream>

#include "initializer.h"
#include "mpi.h"
#include "particle_data.h"
#include "octree_manager.h"
#include "tool_fn.h"
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

 
    octree->refine_octree(1,refine_init,NULL,NULL);  //initial refinement of octree


    octree->partition_octree(0,NULL);
 
    gdata->initFluidParticles();

    octree->destroy_octree();
    mpiret = sc_MPI_Finalize ();
    return 0;
}
