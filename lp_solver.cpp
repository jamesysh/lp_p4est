#include <iostream>
#include "lp_solver.h"

using namespace std;
LPSolver::LPSolver(Global_Data *g){
    gdata = g;
    splitorder = 4;
    invalidpressure = 0;
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
    double realdt = cfldt/3;
    sc_array_t *neighbourlist0;
    sc_array_t *neighbourlist1;
    double *inpressure, *outpressure;
    double *involume, *outvolume;
    double *invelocity, *outvelocity;
    double *insoundspeed, *outsoundspeed;
    double vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1; // output
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
         computeSpatialDer(dir, pad, neighbourlist0, inpressure, invelocity,
        &vel_d_0, &vel_dd_0, &p_d_0, &p_dd_0);
     
        computeSpatialDer(dir, pad, neighbourlist1, inpressure, invelocity,
        &vel_d_1, &vel_dd_1, &p_d_1, &p_dd_1);
        /* 
         if(p_d_1 > 0){

             printf("%f %f %f\n",pad->xyz[0],pad->xyz[1],pad->xyz[2]);
             double x = pad->xyz[0];
             double y = pad->xyz[1];
             double z = pad->xyz[2];
             cout<<sqrt(x*x+y*y+z*z)<<endl;
         }
         
         if(p_d_0 > 0){
             printf("%f %f %f\n",pad->xyz[0],pad->xyz[1],pad->xyz[2]);
             double x = pad->xyz[0];
             double y = pad->xyz[1];
             double z = pad->xyz[2];
             cout<<sqrt(x*x+y*y+z*z)<<endl;
         }*/
        /* 
if(p_d_0>0){
    printf("%f %f %f %f\n",p_d_0,p_d_1,vel_d_0,vel_d_1);
 
    printf("%f %f %f\n",pad->xyz[0],pad->xyz[1],pad->xyz[2]);
}
*/
    timeIntegration( realdt,
	 0,  *involume, *invelocity, *inpressure, *insoundspeed, 
	vel_d_0, vel_dd_0, p_d_0, p_dd_0,
    vel_d_1, vel_dd_1, p_d_1, p_dd_1,
    outvolume, outvelocity, outpressure);
        
         if(*outpressure < invalidpressure)
            pad->schemeorder = 0;
        
        if(std::isnan(*outvolume) || std::isnan(*outvelocity) || std::isnan(*outpressure)) 
            pad->schemeorder = 0;

        if(pad->schemeorder == 0)
        {
            *outvolume = *involume;
            *outpressure = *inpressure;
            *outvelocity = *invelocity;
            *outsoundspeed = *insoundspeed;
        }
        else{
            *outsoundspeed = gdata->eos->getSoundSpeed(*outpressure, 1./(*outvolume));
        
        }

      
      
      
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
     if(dir == 1){
     
        *invelocity = &pad->v[1];
        *outvelocity = &pad->oldv[1];
     }

     if(dir == 2){
     
        *invelocity = &pad->v[2];
        *outvelocity = &pad->oldv[2];
     }

}

void LPSolver::setNeighbourListPointer(pdata_t *pad, sc_array_t** neilist0, sc_array_t **neilist1,int dir){
    if(dir == 0){
        *neilist0 = pad->neighbourfrontparticle;
        *neilist1 = pad->neighbourbackparticle;
    }
    else if(dir == 1){
    
        *neilist0 = pad->neighbourrightparticle;
        *neilist1 = pad->neighbourleftparticle;
    }
    else if(dir == 2){
    
        *neilist0 = pad->neighbourupparticle;
        *neilist1 = pad->neighbourdownparticle;
    }
}



void LPSolver::computeSpatialDer(int dir,pdata_t *pad, sc_array_t *neighbourlist, const double* inpressure, const double *invelocity,
        double *vel_d, double *vel_dd, double *p_d, double *p_dd) {
    size_t numrow = gdata->numrow1st;
    size_t numcol = 3;
    double distance;
    int info;
    int offset = 3; // 2 for 2 dimension
    neighbour_info_t *ninfo;
    while( true ){
        if(numrow > neighbourlist->elem_count){
            *vel_d = 0;
            *vel_dd = 0;
            *p_d = 0;
            *p_dd = 0;
            break;
        }
        ninfo = (neighbour_info_t *)sc_array_index(neighbourlist,0);
        distance = ninfo->distance;
        double A[numrow*numcol];
        computeA3D(A, pad, neighbourlist,numrow,distance);
        double B[numrow];
    
        computeB(B, pad, neighbourlist, numrow, inpressure ,PRESSURE, dir);
    
		QRSolver qrSolver(numrow,numcol,A);
		
		double result[numcol];
	    info = qrSolver.solve(result,B);
        if(info != 0){
            numrow ++;
        } 
        else{
            
            *p_d = result[dir]/distance; //dir=0 (x), dir=1(y), dir=2(z)
            *p_dd = pad->schemeorder == 2?  result[dir+offset]/distance/distance:0;
           /* 
                cout<<"newone"<<*p_d*distance<<endl;
            if(*p_d > 0){
                cout<<"newone"<<*p_d<<endl;
                
                for(int i=0;i<numrow;i++){
                cout<<B[i]<<endl;
                }
            
            }
            */
             computeB(B, pad, neighbourlist, numrow, invelocity ,VELOCITY, dir);
        
			 qrSolver.solve(result,B);
        
			*vel_d = result[dir]/distance; //dir=0 (x), dir=1(y), dir=2(z)
        
            *vel_dd = pad->schemeorder == 2? result[dir+offset]/distance/distance:0;
        
			if(std::isnan(*p_d) || std::isnan(*p_dd) || std::isnan(*vel_d) || std::isnan(*vel_dd) ||
			   std::isinf(*p_d) || std::isinf(*p_dd) || std::isinf(*vel_d) || std::isinf(*vel_dd)) {
                
                numrow ++;
                
            } 
            
            else
            {
                break;
            }
        }
    }

}// to do



void LPSolver::computeA3D(double *A, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, double distance){
    double x, y, z, x0, y0, z0;
    double h,k,l;
    pdata_copy_t * padnei;
    x = pad->xyz[0];
    y = pad->xyz[1];
    z = pad->xyz[2];
   
    if(pad->schemeorder == 1){
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist,i);
            x0 = padnei->xyz[0];
            y0 = padnei->xyz[1];
            z0 = padnei->xyz[2];
            h = (x0-x)/distance;
            k = (y0-y)/distance;
            l = (z0-z)/distance;
//           printf("%f %f %f \n",h,k,l);
            //  cout<<k<<endl;
          //  if(k>0)
          //      cout<<"warning"<<endl;
            A[i]            = h;
            A[i + 1*numrow] = k;
            A[i + 2*numrow] = l;
   
        }
    
    }
    if(pad->schemeorder == 2){
        for(size_t i=0; i<numrow; i++){
        
            gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist,i);
            x0 = padnei->xyz[0];
            y0 = padnei->xyz[1];
            z0 = padnei->xyz[2];
            h = (x0-x)/distance;
            k = (y0-y)/distance;
            l = (z0-z)/distance;
        
			A[i]            = h;
			A[i + 1*numrow] = k;
			A[i + 2*numrow] = l;
			A[i + 3*numrow] = 0.5*h*h;
			A[i + 4*numrow] = 0.5*k*k;
			A[i + 5*numrow] = 0.5*l*l;
			A[i + 6*numrow] = h*k;
			A[i + 7*numrow] = h*l;
			A[i + 8*numrow] = k*l;
        } 
    
    
    }


}

void LPSolver::computeB(double *B, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, const double* indata, indata_t datatype, int dir){
    pdata_copy_t *padnei;
    if(datatype == PRESSURE){
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist,i);
            B[i] = padnei->pressure-*indata;  
//            cout<<B[i]<<endl;
        }   
    }

    else if(datatype == VELOCITY){
    
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle(pad,&padnei,neighbourlist,i);
            B[i] = padnei->v[dir] - *indata;  
        }   
    
    }
}


void LPSolver::timeIntegration(
	double realDt,
	double gravity, double inVolume, double inVelocity, double inPressure, double inSoundSpeed, 
	double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0,
	double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1,
	double* outVolume, double* outVelocity, double* outPressure){

    double multiplier1st = 3.; // 2 for 2 dimension;
    double multiplier2nd = 3./4;


	double K = inSoundSpeed*inSoundSpeed/inVolume/inVolume; 

	double Pt1st = -0.5*inVolume*K*(vel_d_0+vel_d_1) + 0.5*inVolume*sqrt(K)*(p_d_0-p_d_1);
//	printf("%f %f %f %f\n",p_d_0,p_d_1,vel_d_0,vel_d_1);
    double Pt2nd = -inVolume*inVolume*pow(K,1.5)*(vel_dd_0-vel_dd_1) + inVolume*inVolume*K*(p_dd_0+p_dd_1);
	double Pt = multiplier1st*Pt1st + multiplier2nd*Pt2nd;
	// Vt
	double Vt = -Pt/K;

	// VELt
	double VELt1st = 0.5*inVolume*sqrt(K)*(vel_d_0-vel_d_1) - 0.5*inVolume*(p_d_0+p_d_1);
	double VELt2nd = inVolume*inVolume*K*(vel_dd_0+vel_dd_1) - inVolume*inVolume*sqrt(K)*(p_dd_0-p_dd_1);
	double VELt = multiplier1st*VELt1st + multiplier2nd*VELt2nd;

	(*outVolume)   = inVolume   + realDt*Vt;
	(*outPressure) = inPressure + realDt*Pt;
	(*outVelocity) = inVelocity + realDt*(VELt+gravity);	

	if(std::isnan(*outVolume) || std::isinf(*outVolume) || 
	   std::isnan(*outPressure) || std::isinf(*outPressure) ||
	   std::isnan(*outVelocity) || std::isinf(*outVelocity)) {
	    printf("%.16g,%.16g,%.16g\n",*outVolume,*outPressure,*outVelocity);
        assert(false);   
	}

}

void LPSolver:: computeCFLCondition(){

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;

  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  double lmindt = 100, gmindt = 0;
  double dt;
  double speed,sc;
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
         sc = pad->soundspeed;
         speed = sqrt(pad->v[0]*pad->v[0]+pad->v[1]*pad->v[1]+pad->v[2]*pad->v[2]);
         dt = pad->localspacing/max(sc,speed);
         if(dt < lmindt || lmindt == 100){
            lmindt = dt;
         }
      }
       offset = lpend;  
    
    }
  }

    int mpiret = MPI_Allreduce(&lmindt,&gmindt,1,MPI_DOUBLE,MPI_MIN,gdata->mpicomm);

    SC_CHECK_MPI (mpiret);

    cfldt = gmindt;

    P4EST_GLOBAL_ESSENTIALF ("MINCFL timestep is %f. \n", cfldt);

}

