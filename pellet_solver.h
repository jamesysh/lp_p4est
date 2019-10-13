#ifndef __PELLET_SOLVER_H__
#define __PELLET_SOLVER_H__
#include "particle_data.h"

class PelletSolver{
    
    public:
    
        PelletSolver(Global_Data *gdata);
        ~PelletSolver();
    
    private:
        Global_Data *gdata;
        sc_array_t *particle_data_copy;
        
        
        
        
        p4est_t *p4est_heating;
        p4est_connectivity_t *conn;
    
    };






#endif
