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
        void postsearch2d();
        void regroupParticles2d();
        
        void communicateParticles();
  
        void partitionParticles2d();
        void adaptQuadtree();
        void destoryQuadtree();
  
        void swapXAndZCoordinate();
        void computeDensityIntegral();
        void split_by_coord ( sc_array_t * in,
                sc_array_t * out[2], pa_mode_t mode, int component,
                const double lxyz[3], const double dxyz[3]);
        
        int countNumberinRange(sc_array_t *view, int n, double x, double y); 
        
        
        
        int elem_particle_box;
        int elem_particle_cell;
        sc_array_t *prebuf;
        p4est_t *p4est_heating;
        p4est_connectivity_t *conn;
        sc_array_t *particle_data_copy; //used for pellet problem;
        Global_Data *gdata;
    
        sc_array_t *ilh[2],*jlh[2],*klh[2];
        
        
        
        
    };






#endif
