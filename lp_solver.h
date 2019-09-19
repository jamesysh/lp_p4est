#ifndef __LP_SOLVER_H__
#define __LP_SOLVER_H__
#include "initializer.h"
#include "particle_data.h"
#include <vector>
class LPSolver {

public:
	
    LPSolver(Global_Data *g);
    
    
    ~LPSolver(){}
    
    void moveParticlesByG(double dt);
    void solve_upwind(int phase);
     
    void setInAndOutPointer(pdata_t *pad, double **inpressure, double **outpressure, double **involume, double **outvolume,
        double** invelocity, double **outvelocity, double **insoundspeed, double **outsoundspeed, int dir, int phase);

    void setNeighbourListPointer(pdata_t *pad, sc_array_t** neilist0, sc_array_t **neilist1, int dir );
    









    double dt = 0.0001;
    Global_Data * gdata; 

	std::vector<std::vector<int> > m_vDirSplitTable; 
    int splitorder;
};


#endif
