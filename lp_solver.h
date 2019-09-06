#ifndef __LP_SOLVER_H__
#define __LP_SOLVER_H__
#include "initializer.h"
#include "particle_data.h"

class LPSolver {

public:
	
    LPSolver(Global_Data *g);
    ~LPSolver(){}
    
    void moveParticlesByG(double dt);

    double dt;
    Global_Data * gdata; 

};


#endif
