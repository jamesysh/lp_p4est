#include "iostream"
#include "lp_solver.h"
LPSolver::LPSolver(Global_Data *g){
    gdata = g;

}

void LPSolver::moveParticlesByG(double dt){
    double g = 9.8;
    pdata_t *pad;
    p4est_locidx_t li, lpnum = gdata->lpnum;
     
    pad = (pdata_t *) gdata->sc_array_index_begin(gdata->particle_data);
    for(li = 0; li<lpnum; li++){
       pad->oldv[0] = pad->v[0];
       pad->oldv[1] = pad->v[1];
       pad->oldv[2] = pad->v[2];
       
       pad->v[1] -= g*dt;
   
//move
       pad->xyz[0] += 0.5*dt*(pad->oldv[0]+pad->v[0]);
       pad->xyz[1] += 4;// 0.5*dt*(pad->oldv[1]+pad->v[1]);
       pad->xyz[2] += 0.5*dt*(pad->oldv[2]+pad->v[2]);
  
       pad ++;
    }
}
