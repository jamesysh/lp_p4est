#include <iostream>
#include "lp_solver.h"

using namespace std;
LPSolver::LPSolver(Global_Data *g){
    gdata = g;
    splitorder = 0;
    m_vDirSplitTable = vector<vector<int> >
    ({{0,1,2},
      {0,2,1},
      {1,0,2},
      {1,2,0},
      {2,0,1},
      {2,1,0}});
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
    
    const int dir = m_vDirSplitTable[splitorder][phase];
    double realdt = dt/3;
    sc_array_t *neighbourlist0;
    sc_array_t *neighbourlist1;
    double *inpressure, *outpressure;
    double *involume, *outvolume;
    double *invelocity, *outvelocity;
    double *insoundspeed, *outsoundspeed;
    double vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1; // output
    size_t numrow, numcol;
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
    
      lpend = qud->lpend;
      for(int i=offset;i<lpend;i++){
         pad = (pdata_t *)sc_array_index(gdata->particle_data,i);
         if(pad->ifboundary)
             continue;
         setInAndOutPointer(pad, &inpressure, &outpressure, &involume, &outvolume,
                  &invelocity, &outvelocity, &insoundspeed, &outsoundspeed, dir, phase);
      
         setNeighbourListPointer(pad, &neighbourlist0, &neighbourlist1,dir);
         if(*insoundspeed == 0 || *involume == 0){
             pad->schemeorder = 0;
            *outvolume = *involume;
            *outpressure = *inpressure;
            *outvelocity = *invelocity;
            *outsoundspeed = *insoundspeed;
            printf("Detect a particle which has 0 volume or 0 soundspeed!!!.\n"); 
         }
            
         pad->schemeorder = 1;
         assert(neighbourlist0->elem_count >= gdata->numrow1st);
         assert(neighbourlist1->elem_count >= gdata->numrow1st);
      }
       offset = lpend;  
    
    }
  }

} //to do

void LPSolver::setInAndOutPointer(pdata_t *pad, double **inpressure, double **outpressure, double **involume, double **outvolume,
        double** invelocity, double **outvelocity, double **insoundspeed, double **outsoundspeed, int dir, int phase){

     if(phase == 0){
        *inpressure = &pad->pressure;
        *outpressure = &pad->pressureT1;
        *involume = &pad->volume;
        *outvolume = &pad->volumeT1;
        *insoundspeed = &pad->soundspeed;
        *outsoundspeed = &pad->soundspeedT1;
     }

     if(phase == 1){
        *inpressure = &pad->pressureT1;
        *outpressure = &pad->pressureT2;
        *involume = &pad->volumeT1;
        *outvolume = &pad->volumeT2;
        *insoundspeed = &pad->soundspeedT1;
        *outsoundspeed = &pad->soundspeedT2;
     }

     if(phase == 2){
        *inpressure = &pad->pressureT2;
        *outpressure = &pad->pressureT1;
        *involume = &pad->volumeT2;
        *outvolume = &pad->volumeT1;
        *insoundspeed = &pad->soundspeedT2;
        *outsoundspeed = &pad->soundspeedT1;
     }
    
     if(dir == 0){
        *invelocity = &pad->v[0];
        *outvelocity = &pad->oldv[0];
     
     }
     if(dir = 1){
     
        *invelocity = &pad->v[1];
        *outvelocity = &pad->oldv[1];
     }

     if(dir = 2){
     
        *invelocity = &pad->v[2];
        *outvelocity = &pad->oldv[2];
     }

}

void LPSolver::setNeighbourListPointer(pdata_t *pad, sc_array_t** neilist0, sc_array_t **neilist1,int dir){
    if(dir == 0){
        *neilist0 = pad->neighbourfrontparticle;
        *neilist0 = pad->neighbourbackparticle;
    }
    else if(dir == 1){
    
        *neilist0 = pad->neighbourrightparticle;
        *neilist0 = pad->neighbourleftparticle;
    }
    else if(dir == 2){
    
        *neilist0 = pad->neighbourupparticle;
        *neilist0 = pad->neighbourdownparticle;
    }
}



void LPSolver::computeSpatialDer(int dir,pdata_t *pad, sc_array_t *neighbourlist, const double* inpressure, const double *invelocity,
        double *vel_d, double *vel_dd, double *p_d, double *p_dd) {
    size_t numrow = gdata->numrow1st;
    size_t numcol = 3;
    double distance;
    while( true ){
        if(numrow > neighbourlist->elem_count){
            *vel_d = 0;
            *vel_dd = 0;
            *p_d = 0;
            *p_dd = 0;
            pad->schemeorder -= 1./6;
            break;
        }
        double A[numrow*numcol];


    
    }

}// to do



void LPSolver::computeA3D(double *A, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, double distance){
    double x, y, z, x0, y0, z0;
    x = pad->xyz[0];
    y = pad->xyz[1];
    z = pad->xyz[2];

}







