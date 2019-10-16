#ifndef __PELLET_SOLVER_H__
#define __PELLET_SOLVER_H__
#include "particle_data.h"


typedef struct quadrant_data
{   
    
    p4est_locidx_t      lpend;
//    p4est_locidx_t ghostneighbourid[300];
  /** counts of local particles remaining on this quadrant and recieved ones */
  p4est_locidx_t      premain, preceive;
}
quadrant_data_t;

class PelletSolver{
    
    public:
    
        PelletSolver(Global_Data *gdata);
        ~PelletSolver();
        
        void build_quadtree();
        void prerun();
        void resetQuadrantData();
        void presearch2d(); 
        void packParticles();
        
        size_t elem_particle_box;
        p4est_t *p4est_heating;
        p4est_connectivity_t *conn;
        sc_array_t *particle_data_copy; //used for pellet problem;
        Global_Data *gdata;
    private:
        
        
        
        
        sc_array_t *ilh[2],*jlh[2],*klh[2];
    };






#endif
