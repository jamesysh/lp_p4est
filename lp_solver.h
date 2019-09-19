#ifndef __LP_SOLVER_H__
#define __LP_SOLVER_H__
#include "initializer.h"
#include "particle_data.h"

class LPSolver {

public:
	
    LPSolver(Global_Data *g);
    ~LPSolver(){}
    
    void moveParticlesByG(double dt);
    void solve_upwind(int phase);

    double dt = 0.0001;
    Global_Data * gdata; 
     
    void setInAndOutPointer(pdata_t *pad, double *inpressure, double *outpressure, double *involume, double *outvolume,
        double* invelocity, double *outvelocity, double *insoundspeed, double *outsoundspeed, int dir, int phase);

    

};


#endif
