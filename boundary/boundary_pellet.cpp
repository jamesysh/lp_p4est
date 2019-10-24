#include "boundary_pellet.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include "sc.h"
#include "particle_data.h"
#include "pellet_solver.h"
using namespace std;

PelletInflowBoundary::PelletInflowBoundary():Pinflow(30),Uinflow(0),Vinflow(100){}

void PelletInflowBoundary::generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx, double dt){
    computeMassFlowRate(g,dx);
    
    static double mass_fix = dx*dx*dx/Vinflow/sqrt(2);
    
    double pr = 0.2;
    
    double Ts = 700;
    double gamma = 1.67;
    double R = 83.1446/20.18;
    
    //double pv  = Vinflow*massflowrate/4/M_PI/pr/pr;
    
    double pv = sqrt(gamma*R*Ts)/10;
    if(massflowrate != 0){
        Vinflow = 4*M_PI*pr*pr*pv/massflowrate;
        Pinflow = R*Ts/Vinflow;
    }
    
    
    pelletvelocity = pv;
    double xcen = 0;
    double ycen = 0;
    double zcen = 0;
    int n = 4.0*3.1416*pr*pr*pr/dx/dx/dx*sqrt(2.0)/5;

     int numberofNewFluid = massflowrate*dt/mass_fix;
     double actualdx = sqrt(4*M_PI*pr*pr/numberofNewFluid)/2.5;

    if(dx<actualdx && actualdx<0.1)
        dx = actualdx;
    
    computeRadialDerivative(g,dx);
    
    double newpir = pr * 4/5;
    g->gpnum += n;
    pdata_t*  pad;
    
     
    double x, y, z, d_x, d_y, d_z, dr;
    for(int i=0;i<n;i++)
    if(g->mpirank == 0){	
    {

        pad = (pdata_t*)sc_array_push_count(g->particle_data,1);
        double tx=1,ty=1,tz=1,tr;
		
        while(tx*tx+ty*ty+tz*tz>1)
		
        {
				tx=2*(double)rand()/(double)RAND_MAX-1;
				ty=2*(double)rand()/(double)RAND_MAX-1;
				tz=2*(double)rand()/(double)RAND_MAX-1;
		
    }
        tr=sqrt(tx*tx+ty*ty+tz*tz);
        tx=tx/tr,ty=ty/tr,tz=tz/tr;
        tr=newpir+(double)rand()/(double)RAND_MAX*(pr-newpir);
        x = tx*tr;
        y = ty*tr;
        z = tz*tr;
        d_x = tx*tr;
        d_y = ty*tr;
        d_z = tz*tr;
        dr = sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
	    if(dr>0.2)
            cout<<"warning"<<endl;
        /*	vx[inflowEndIndex] = m_voldv[pi]*d_x/dr;
           	vy[inflowEndIndex] = m_voldv[pi]*d_y/dr;
           	vz[inflowEndIndex] = m_voldv[pi]*d_z/dr;
		    */
    
   	    double newpv = pv + (dr-pr)*ux;    
   	       
        pad->v[0] = newpv*d_x/dr;
        pad->v[1] = newpv*d_y/dr;
        pad->v[2] = newpv*d_z/dr;
        pad->xyz[0] = x;
        pad->xyz[1] = y;
        pad->xyz[2] = z;

        pad->volume = Vinflow;
        pad->pressure = Pinflow + (dr-pr)*px;
        pad->localspacing = dx;
        pad->mass = mass_fix;
        pad->ifboundary = true;
        pad->flagdelete = !g->flagdelete;
        pad->soundspeed = m_pEOS->getSoundSpeed(pad->pressure,1./pad->volume);
        }

    }


    
     P4EST_GLOBAL_ESSENTIALF("Generate %d new fluid particles.\n", numberofNewFluid);
     g->gpnum += numberofNewFluid;

     int rnd = rand();

     double dis = pv*dt;
     double offset = 2./numberofNewFluid;
     double increment = M_PI*(3-sqrt(5));
     if(g->mpirank == 0){
         for(int i=0;i<numberofNewFluid;i++){
               double r_random = ((double)rand()/(double)(RAND_MAX))*dis; 
               double  y_tmp =  ((i*offset-1)+offset/2);
               pad = (pdata_t*)sc_array_push_count(g->particle_data,1);
               y = y_tmp * (pr+r_random) + ycen;
               
               pad->xyz[1] = y;
               double r = sqrt(1-y_tmp*y_tmp);
               double phi = ((i+rnd)%numberofNewFluid) * increment;

               x = cos(phi)*r*(pr+r_random) + xcen;
               
               pad->xyz[0] = x;
               z = sin(phi)*r*(pr+r_random) + zcen;
               
               pad->xyz[2] = z;

               double d_x = x;
               double d_y = y;
               double d_z = z;
               double dr = sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
               
   	           double newpv = pv + (dr-pr)*ux;    
               pad->v[0] = newpv*d_x/dr;
               pad->v[1] = newpv*d_y/dr;
               pad->v[2] = newpv*d_z/dr;
               
               pad->volume = Vinflow;
               pad->pressure = Pinflow + (dr-pr)*px;
               pad->localspacing = dx;
               pad->mass = mass_fix;
               pad->ifboundary = false; 
               pad->soundspeed = m_pEOS->getSoundSpeed(pad->pressure,1./pad->volume);
		       pad->redocount = 0;
           }
           
           }


}
void PelletInflowBoundary::UpdateInflowBoundary(Global_Data* g, EOS* m_pEOS, double dt, double dx){
    
    double massflowrate = 0.6;

    double mass_fix = dx*dx*dx/Vinflow/sqrt(2);

    int numberofNewFluid = massflowrate*dt/mass_fix/10;
    g->gpnum += numberofNewFluid;
    pdata_t*  pad;

     int rnd = rand();

     double pr = 0.2;
     double pv  = Vinflow*massflowrate/4/M_PI/pr/pr;
     double x,y,z;
     double dis = dx/5;
     double offset = 2./numberofNewFluid;
     double increment = M_PI*(3-sqrt(5));
     for(int i=0;i<numberofNewFluid;i++){
           if(g->mpirank == i%g->mpisize){
           double r_random = ((double)rand()/(double)(RAND_MAX))*dis; 
           double  y_tmp =  ((i*offset-1)+offset/2);
           pad = (pdata_t*)sc_array_push_count(g->particle_data,1);
           y = y_tmp * (pr+r_random) + 0;
           
           pad->xyz[1] = y;
           double r = sqrt(1-y_tmp*y_tmp);
           double phi = ((i+rnd)%numberofNewFluid) * increment;

           x = cos(phi)*r*(pr+r_random) + 0;
           
           pad->xyz[0] = x;
           z = sin(phi)*r*(pr+r_random) + 0;
           
           pad->xyz[2] = z;

           double d_x = x;
           double d_y = y;
           double d_z = z;
           double dr = sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
           pad->v[0] = pv*d_x/dr;
           pad->v[1] = pv*d_y/dr;
           pad->v[2] = pv*d_z/dr;
           pad->volume = Vinflow;
           pad->pressure = Pinflow;
           pad->localspacing = dx;
           pad->mass = mass_fix;
           pad->ifboundary = false; 
           pad->soundspeed = m_pEOS->getSoundSpeed(pad->pressure,1./pad->volume);
           }
           
           }
     

}


void PelletInflowBoundary::computeMassFlowRate(Global_Data *g,double dx){
   static bool first=true;
   if(first){
       massflowrate = 0;
       first = false;
       return;
       }
   pdata_t *pad;
   size_t li, lpnum = g->particle_data->elem_count;
   double x, y, z, dr;
   double pr = 0.2;
   double xcen=0;
   double ycen=0;
   double zcen=0;
   int counter_nei = 0;
   double qsum = 0;
   int counter_g = 0;
   double qsum_g = 0;
   PelletSolver *p = g->pellet_solver;
   double sublimationenergy = p->sublimationenergy;
   for(li = 0; li<lpnum; li++){
       pad = (pdata_t *)sc_array_index(g->particle_data,li);
       if(pad->ifboundary)
           continue;
       x = pad->xyz[0];
       y = pad->xyz[1];
       z = pad->xyz[2];
       dr = sqrt((x-xcen)*(x-xcen)+(y-ycen)*(y-ycen)+(z-zcen)*(z-zcen))-pr;
       if(dr<dx){
            qsum += pad->qplusminus;
            counter_nei ++; 
           } 
       }
       
       MPI_Barrier(g->mpicomm);
       
       MPI_Allreduce (&counter_nei, &counter_g, 1, MPI_INT,
                                MPI_SUM, g->mpicomm);
       if(counter_g == 0){
           assert(false);
           }
       MPI_Allreduce(&qsum, &qsum_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       
       massflowrate = qsum_g/counter_g*4*M_PI*pr*pr/sublimationenergy*2/M_PI; 

       P4EST_GLOBAL_ESSENTIALF("The massflowrate is %.16g.\n", massflowrate);
}

void PelletInflowBoundary::computeRadialDerivative(Global_Data *g,double dx){

   pdata_t *pad;
   size_t li, lpnum = g->particle_data->elem_count;
   double x, y, z, dr;
   double pr = 0.2;
   double xcen=0;
   double ycen=0;
   double zcen=0;
   int counter_nei = 0;
   double psum = 0;
   double usum = 0;
   int counter_g = 0;
   double psum_g = 0;
   double usum_g = 0;
   double vx, vy, vz;
   for(li = 0; li<lpnum; li++){
       pad = (pdata_t *)sc_array_index(g->particle_data,li);
       if(pad->ifboundary)
           continue;
       x = pad->xyz[0];
       y = pad->xyz[1];
       z = pad->xyz[2];
       dr = sqrt((x-xcen)*(x-xcen)+(y-ycen)*(y-ycen)+(z-zcen)*(z-zcen))-pr;
       if(dr>0*dx && dr < 1*dx){
            vx = pad->v[0];
            vy = pad->v[1];
            vz = pad->v[2];
            psum += pad->pressure;
            usum += (vx*(x-xcen)+vy*(y-ycen)+vz*(z-zcen))/(dr+pr); 
            counter_nei ++; 
           } 
       }
       
       MPI_Barrier(g->mpicomm);
       
       MPI_Allreduce (&counter_nei, &counter_g, 1, MPI_INT,
                                MPI_SUM, g->mpicomm);
       if(counter_g == 0){
           ux = 0;
           px = 0;
           return;
           }

       avg_dis = .5*dx;
       MPI_Allreduce(&psum, &psum_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       
       px = (psum_g/counter_g-Pinflow)/avg_dis;

       MPI_Allreduce(&usum, &usum_g,1,MPI_DOUBLE, MPI_SUM, g->mpicomm);
       
       ux = (usum_g/counter_g-pelletvelocity)/avg_dis;
        
}

