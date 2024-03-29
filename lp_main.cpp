#include<iostream>

#include <sys/stat.h> // mkdir()
#include <cstring> // strerror()
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



int main(int argc, const char* argv[]){
    
	string inputfileName;
	string datafileName; // restart datafile
	string outputfileNameFluid;		
	string debugfileName;
    bool ifDebug = false;
    
    int mpiret;
    mpiret = sc_MPI_Init (NULL, NULL);

    SC_CHECK_MPI (mpiret);
    int mpirank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	for(int i=1; i<argc; i++) { // argv[0] is name of the executable 
		if(!strcmp(argv[i], "-i")) { // input; strcmp returns 0 if two c_str are the same
			if(i+1 >= argc || argv[i+1][0]=='-') { // no inputfile name following -i option
				cerr<<"ERROR: Needs to provide inputfile name!"<<endl;
				exit(1);
			}
			inputfileName = argv[i+1];
			cout<<"input file name: "<<inputfileName<<endl;
			if(i+2 < argc && argv[i+2][0]!='-') { // there is second input file (restart data file)
				datafileName = argv[i+2];
				cout<<"resatrt data file name: "<<datafileName<<endl;
			}
		}
		else if(!strcmp(argv[i], "-o")) { // output
			if(i+1 >= argc || argv[i+1][0]=='-') { // no outputfile name following -o option
				cerr<<"ERROR: Needs to provide outputfile name!"<<endl;
				exit(1);
			}
            if(mpirank == 0){
			if(mkdir(argv[i+1], 0777) == -1) { //mkdir takes const char*
				cerr<<"ERROR:"<<strerror(errno)<<endl;
				exit(1);
			}
            }
			outputfileNameFluid = string(argv[i+1]);
			debugfileName = string(argv[i+1]);
			cout<<"output file name: "<<outputfileNameFluid<<endl;
		}
		else if(!strcmp(argv[i], "-d")) { // debug
			ifDebug = true;
			if(i+1 >= argc || argv[i+1][0]=='-') // no debugfile name following -d option is fine, use default name
			{
                debugfileName = debugfileName + "/" + "debug";
                
                cout<<"debug file name: "<<debugfileName<<endl;

                continue;
            }
			debugfileName = debugfileName + "/" + argv[i+1];
			cout<<"debug file name: "<<debugfileName<<endl;
		}	
	}
	if(inputfileName.empty()) {
		cerr<<"ERROR: Needs to provide -i option!"<<endl;
		exit(1);
	}
	if( outputfileNameFluid.empty()) {
		cerr<<"ERROR: Needs to provide -o option!"<<endl;
		exit(1);
	}
    
    
    
    
    
    double t1 = MPI_Wtime();
    Initializer *init;

	if(datafileName.empty())
		init = new Initializer(inputfileName, ifDebug, debugfileName);
	else // restart to do
		init = new Initializer(inputfileName, datafileName, ifDebug, debugfileName);
    
    Global_Data *gdata = new Global_Data(init);

    ParticleViewer * viewer = new ParticleViewer(gdata,outputfileNameFluid,4);

    Octree_Manager *octree = new Octree_Manager(gdata);

    octree->build_octree();

    if(gdata->dimension == 3) 
        octree->refine_octree(1,refine_init,NULL,NULL);  //initial refinement of octree
    else if(gdata->dimension == 2)
        octree->refine_octree2d(1,refine_init2d,NULL,NULL);  //initial refinement of octree

    octree->partition_octree(1);
    gdata->prerun(); 
   // gdata->initFluidParticles_distributed();
    gdata->initFluidParticles_hexagonal();
    if(gdata->dimension == 3)
       gdata->resetOctantData(); 
    else if(gdata->dimension == 2)
        gdata->resetOctantData2d();
    LPSolver * lpsolver = new LPSolver(init, gdata,octree,viewer);
    if(gdata->dimension == 3 )
        lpsolver->solve_3d();
    else if(gdata->dimension == 2)
        lpsolver->solve_2d();
    
    
    gdata->cleanUpArrays(); 
  

    double t2 = MPI_Wtime();
   
    P4EST_GLOBAL_ESSENTIALF ("Elapsed time is %f.\n", t2-t1);
    octree->destroy_octree();
   
    
    delete gdata;
    delete octree;
    delete lpsolver;
    delete init;

    mpiret = sc_MPI_Finalize ();

   
    return 0;
}
