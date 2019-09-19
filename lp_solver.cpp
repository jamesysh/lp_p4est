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


//dir: 0 x, 1 y, 2 z
void LPSolver::solve_upwind(int phase){
    sc_array_t *neighbourlist1;
    sc_array_t *neighbourlist2;
    double *inpressure, *outpressure;
    double *involume, *outvolume;
    double *invelocity, *outvelocity;
    double *insoundspeed, *outsoundspeed;
    
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;

  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  for (tt = gdata->p8est->first_local_tree; tt <= gdata->p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (gdata->p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      sc_array_destroy(qud->localneighbourid);
      sc_array_destroy(qud->ghostneighbourid);
    
      lpend = qud->lpend;
      for(int i=offset;i<lpend;i++){
         pad = (pdata_t *)sc_array_index(gdata->particle_data,i);
         if(pad->ifboundary)
             continue;
         
      }
       offset = lpend;  
    
    }
  }

}

void LPSolver::setInAndOutPointer(pdata_t *pad, double *inpressure, double *outpressure, double *involume, double *outvolume,
        double* invelocity, double *outvelocity, double *insoundspeed, double *outsoundspeed, int dir, int phase){

     if(phase == 0){
        inpressure = &pad->pressure;
        outpressure = &pad->pressureT1;
        involume = &pad->volume;
        outvolume = &pad->volumeT1;
        insoundspeed = &pad->soundspeed;
        outsoundspeed = &pad->soundspeedT1;
     }


}
