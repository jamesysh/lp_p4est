#ifndef __PELLET_SOLVER_H__
#define __PELLET_SOLVER_H__
#include "particle_data.h"

class PelletSolver{
    
    public:
    
        PelletSolver(Global_Data *gdata);
        ~PelletSolver();
        
        void build_quadtree();
        void prerun();
        void resetOctantData2d();
        void presearch2d(); 
        p4est_t *p4est_heating;
        p4est_connectivity_t *conn;
        sc_array_t *particle_data_copy; //used for pellet problem;
        Global_Data *gdata;
    private:
        
        
        
        
        sc_array_t *ilh[2],*jlh[2],*klh[2];
    };






#endif
