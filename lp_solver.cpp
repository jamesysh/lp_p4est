#include <iostream>
#include "lp_solver.h"

using namespace std;
LPSolver::LPSolver(Global_Data *g, Octree_Manager *o, ParticleViewer *v){
    gdata = g;
    octree = o;
    viewer = v;
    splitorder = 0;
    cflcoefficient = 0.5;
    invalidpressure = 0;
    if(gdata->dimension == 3){
        m_vDirSplitTable = vector<vector<int> >
        ({{0,1,2},
          {0,2,1},
          {1,0,2},
          {1,2,0},
          {2,0,1},
          {2,1,0}});
        totalphase = 3;
    }
    else if(gdata->dimension == 2){
        m_vDirSplitTable = vector<vector<int> >
        ({{0,1},
          {1,0}});
        totalphase = 2;
    }
}

void LPSolver::moveParticle(){
    pdata_t *pad;
    double dt = cfldt;
    size_t li, lpnum = gdata->particle_data->elem_count;
     
    for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
       pad->xyz[0] += 0.5*dt*(pad->oldv[0]+pad->v[0]);
       pad->xyz[1] += 0.5*dt*(pad->oldv[1]+pad->v[1]);
       if(gdata->dimension == 3)
           pad->xyz[2] += 0.5*dt*(pad->oldv[2]+pad->v[2]);
  
    }
}

void LPSolver::updateLocalSpacing(){

    pdata_t *pad;
    size_t li, lpnum = gdata->particle_data->elem_count;
     
    for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
       if(gdata->dimension == 3) 
           pad->localspacing *= cbrt(pad->volume/pad->volumeT1); 
       else if(gdata->dimension == 2)
           pad->localspacing *= sqrt(pad->volume/pad->volumeT2); 
    }

}

//dir: 0 x, 1 y, 2 z
void LPSolver::solve_upwind(int phase){
    
    const int dir = m_vDirSplitTable[splitorder][phase];
    double realdt = gdata->dimension == 3? cfldt/3 : cfldt/2;
    sc_array_t *neighbourlist0;
    sc_array_t *neighbourlist1;
    double *inpressure, *outpressure;
    double *involume, *outvolume;
    double *invelocity, *outvelocity;
    double *insoundspeed, *outsoundspeed;
    double vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1; // output
  


  pdata_t * pad;

  size_t li, lpnum = gdata->particle_data->elem_count;
     
  for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
         bool redo = false;
         
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
        
         if(*outpressure < invalidpressure || 1./(*outvolume)<0)
         {
             redo  = true;
         }
        if(std::isnan(*outvolume) || std::isnan(*outvelocity) || std::isnan(*outpressure)) 
        {
                redo = true;}


        if(redo == true){

            pad->redocount++;
            
            if(pad->redocount >4){
                pad->schemeorder = 0;
            }
            else{
                li--;
                cout<<"redo upwind for a particle"<<endl;
            }

            
        }
        else{
        
        
            *outsoundspeed = gdata->eos->getSoundSpeed(*outpressure, 1./(*outvolume));
            pad->redocount = 0;        
        }

        if(pad->schemeorder == 0)
        {
            pad->redocount = 0;
            *outvolume = *involume;
            *outpressure = *inpressure;
            *outvelocity = *invelocity;
            *outsoundspeed = *insoundspeed;
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
    size_t numrow = pad->redocount + gdata->dimension == 3? gdata->numrow1st:gdata->numrow1st2d;
    size_t numcol = gdata->dimension == 3? 3:2;
    double distance;
    int info;
    int offset = gdata->dimension == 3? 3:2; // 2 for 2 dimension
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
        
        if(gdata->dimension == 3)
            computeA3D(A, pad, neighbourlist,numrow,distance);
        else if(gdata->dimension == 2)
            computeA2D(A, pad, neighbourlist,numrow,distance);
        
        double B[numrow];
    
        if(gdata->dimension == 3)
            computeB3d(B, pad, neighbourlist, numrow, inpressure ,PRESSURE, dir);
        else if(gdata->dimension == 2)
            computeB2d(B, pad, neighbourlist, numrow, inpressure ,PRESSURE, dir);
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
             
        if(gdata->dimension == 3)
            computeB3d(B, pad, neighbourlist, numrow, invelocity ,VELOCITY, dir);
        else if(gdata->dimension == 2)
            computeB2d(B, pad, neighbourlist, numrow, invelocity ,VELOCITY, dir);
        
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


void LPSolver::computeA2D(double *A, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, double distance){
    double x, y,  x0, y0;
    double h,k;
    pdata_copy_t * padnei;
    x = pad->xyz[0];
    y = pad->xyz[1];
   
    if(pad->schemeorder == 1){

        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle2d(pad,&padnei,neighbourlist,i);
            x0 = padnei->xyz[0];
            y0 = padnei->xyz[1];
            h = (x0-x)/distance;
            k = (y0-y)/distance;
//           printf("%f %f %f \n",h,k,l);
            //  cout<<k<<endl;
          //  if(k>0)
          //      cout<<"warning"<<endl;
            A[i]            = h;
            A[i + 1*numrow] = k;
   
        }
    
    }
    if(pad->schemeorder == 2){
        for(size_t i=0; i<numrow; i++){
        
            gdata->fetchNeighbourParticle2d(pad,&padnei,neighbourlist,i);
            x0 = padnei->xyz[0];
            y0 = padnei->xyz[1];
            h = (x0-x)/distance;
            k = (y0-y)/distance;
        
            A[i]            = h;
			A[i + 1*numrow] = k;
			A[i + 2*numrow]	= 0.5*h*h;
			A[i + 3*numrow]	= 0.5*k*k;
			A[i + 4*numrow]	= h*k;
	
        } 
    
    
    }


}


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

void LPSolver::computeB3d(double *B, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, const double* indata, indata_t datatype, int dir){
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
void LPSolver::computeB2d(double *B, pdata_t *pad, sc_array_t *neighbourlist, size_t numrow, const double* indata, indata_t datatype, int dir){
    pdata_copy_t *padnei;
    
    
    if(datatype == PRESSURE){
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle2d(pad,&padnei,neighbourlist,i);
            B[i] = padnei->pressure-*indata;  
//            cout<<B[i]<<endl;
        }   
    }

    else if(datatype == VELOCITY){
    
        for(size_t i=0; i<numrow; i++){
            gdata->fetchNeighbourParticle2d(pad,&padnei,neighbourlist,i);
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

    double multiplier1st = gdata->dimension == 3? 3.:2.; // 2 for 2 dimension;
    double multiplier2nd = gdata->dimension == 3? 3./4*cfldt : 1./2*cfldt;


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



    pdata_t *pad;
    size_t li, lpnum = gdata->particle_data->elem_count;
     
    double lmindt = 100, gmindt = 0;
    double dt;
    double speed,sc;
    for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           continue;
         sc = pad->soundspeed;
         speed = sqrt(pad->v[0]*pad->v[0]+pad->v[1]*pad->v[1]+pad->v[2]*pad->v[2]);
         dt = pad->localspacing/max(sc,speed);
         if(dt < lmindt || lmindt == 100){
            lmindt = dt;
        }
    }
    int mpiret = MPI_Allreduce(&lmindt,&gmindt,1,MPI_DOUBLE,MPI_MIN,gdata->mpicomm);

    SC_CHECK_MPI (mpiret);

    cfldt = gmindt * cflcoefficient;

    P4EST_GLOBAL_ESSENTIALF ("MINCFL timestep is %f. \n", cfldt);

}

void LPSolver:: solve_2d(){


    computeLocalBoundaryAndFluidNum();
    viewer->writeResult(0);
    double tstart = 0;
    double tend = 1;
    double nextwritetime = 0;
    while(tstart<tend)
    {
    
        tstart += cfldt;
    
    
    
        gdata->boundary->generateBoundaryParticle(gdata,gdata->eos,gdata->initlocalspacing);
    
        
    //gdata->boundary->UpdateInflowBoundary(gdata,gdata->eos,lpsolver->dt,gdata->initlocalspacing);
        
        gdata->presearch2d();
        gdata->packParticles();
    
        if(gdata->gpnum == 0)
            
        {
          sc_array_destroy_null (&gdata->recevs);
          sc_hash_destroy_null (&gdata->psend);
          sc_array_destroy_null(&gdata->iremain);

          gdata->psend = NULL;
          sc_mempool_destroy (gdata->psmem);
          gdata->psmem = NULL;
            break;
        }
        gdata->communicateParticles();
        gdata->postsearch2d();

        octree->adapt_octree2d();
    
        octree->balance_octree2d(NULL,octree->balance_replace2d);
    
        gdata->regroupParticles2d(); 
    
        gdata->partitionParticles2d();
    
        gdata->createViewForOctant2d();
    
    
        gdata->searchNeighbourOctant2d();

        gdata->searchNeighbourParticle2d();
        gdata->searchUpwindNeighbourParticle2d(); 
    
        gdata->generateGhostParticle2d();
        
    //gdata->testquad2d();
        computeCFLCondition();
        P4EST_GLOBAL_ESSENTIALF ("Current Time: %f .\n", tstart);
        splitorder = (int)rand()%2;
        MPI_Bcast(&splitorder,1,MPI_INT,0,gdata->mpicomm);
        for(int phase = 0;phase < totalphase;phase++){
        solve_upwind(phase);
        MPI_Barrier(gdata->mpicomm);
        gdata->updateViewForOctant2d(phase);
        MPI_Barrier(gdata->mpicomm); 
    }
    gdata->updateParticleStates();
   
    updateLocalSpacing();
   
     moveParticle();
    if(tstart  >= nextwritetime)
    
    {

        computeLocalBoundaryAndFluidNum();
        nextwritetime += 0.01;    
        viewer->writeResult(tstart);
//        viewer->writeGhost(tstart);
    }
   
    
        gdata->cleanForTimeStep2d();
    
        gdata->switchFlagDelete();
    }


}


void LPSolver::solve_3d(){


    computeLocalBoundaryAndFluidNum();
    viewer->writeResult(0);
    double tstart = 0;
    double tend = 0.0005;
    double nextwritetime = 0;
    while(tstart<tend)
    {
    
    
    
    
    //gdata->boundary->UpdateInflowBoundary(gdata,gdata->eos,lpsolver->dt,gdata->initlocalspacing);
    gdata->presearch();
        
    gdata->packParticles();
    
    if(gdata->gpnum == 0)
        
    {
      sc_array_destroy_null (&gdata->recevs);
      sc_hash_destroy_null (&gdata->psend);
      sc_array_destroy_null(&gdata->iremain);

      gdata->psend = NULL;
      sc_mempool_destroy (gdata->psmem);
      gdata->psmem = NULL;
        break;
    }
    gdata->communicateParticles();
        gdata->postsearch();

        octree->adapt_octree(); 
    
        octree->balance_octree(NULL,octree->balance_replace);
    
        gdata->regroupParticles(); 
    
        gdata->partitionParticles();
    
        gdata->createViewForOctant();
    
        gdata->searchNeighbourOctant();


        gdata->searchNeighbourParticle();
        gdata->searchUpwindNeighbourParticle(); 
    
        gdata->generateGhostParticle();

    
        computeCFLCondition();
    
        tstart += cfldt;
    for(int phase = 0;phase< totalphase;phase++){
        solve_upwind(phase);
        MPI_Barrier(gdata->mpicomm);
        gdata->updateViewForOctant(phase);
        MPI_Barrier(gdata->mpicomm); 
    }
    gdata->updateParticleStates();
   
    updateLocalSpacing();
     
    moveParticle();
    if(tstart  >= nextwritetime)
    
    {

        computeLocalBoundaryAndFluidNum();
        nextwritetime += cfldt;    
        viewer->writeResult(tstart);
//        viewer->writeGhost(tstart);
    }
   
    MPI_Barrier(gdata->mpicomm); 
    
        gdata->cleanForTimeStep();

        gdata->switchFlagDelete();
    }



}


void LPSolver:: computeLocalBoundaryAndFluidNum(){

    pdata_t *pad;
    size_t li, lpnum = gdata->particle_data->elem_count;
    size_t boundarycount = 0, fluidcount = 0; 

    for(li = 0; li<lpnum; li++){
       pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
       if(pad->ifboundary)
           boundarycount ++;
       else
           fluidcount ++;
    }
    gdata->lfluidnum = fluidcount;
    gdata->lboundarynum = boundarycount;
}

