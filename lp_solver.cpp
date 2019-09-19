#include "iostream"
#include "lp_solver.h"
LPSolver::LPSolver(Global_Data *g){
    gdata = g;

}

void LPSolver::moveParticlesByG(double dt){
    double g = 9.8;
    pdata_t *pad;
    size_t li, lpnum = gdata->particle_data->elem_count;
     
    for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
       double x = pad->xyz[0];
        
       double y = pad->xyz[1];
       double z = pad->xyz[2];
       double r = sqrt(x*x+y*y+z*z);

       pad->oldv[0] = pad->v[0];
       pad->oldv[1] = pad->v[1];
       pad->oldv[2] = pad->v[2];
        
       pad->v[0] = 200*x;
       pad->v[1] = 200*y;
       pad->v[2] = 200*z;


   
//move
       pad->xyz[0] +=  0.5*dt*(pad->oldv[0]+pad->v[0]);
       pad->xyz[1] += 0.5*dt*(pad->oldv[1]+pad->v[1]);
       pad->xyz[2] += 0.5*dt*(pad->oldv[2]+pad->v[2]);
  
    }
}
