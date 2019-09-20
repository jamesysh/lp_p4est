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
    
    void computeSpatialDer(int dir,pdata_t *pad, sc_array_t *neighbourlist, const double* inpressure, const double *invelocity,
        double *vel_d, double *vel_dd, double *p_d, double *p_dd);


    void computeA3D(pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, size_t numcol, double* distance);







    double dt = 0.0001;
    Global_Data * gdata; 

	std::vector<std::vector<int> > m_vDirSplitTable; 
    int splitorder;
};


#endif
