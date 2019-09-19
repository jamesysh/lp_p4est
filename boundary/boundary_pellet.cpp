#include "boundary_pellet.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include "sc.h"
#include "particle_data.h"
using namespace std;

PelletInflowBoundary::PelletInflowBoundary():Pinflow(30),Uinflow(0),Vinflow(100){}

void PelletInflowBoundary::generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx){
    
    double pr = 0.2;
    size_t n = 4.0*3.1416*pr*pr*pr/dx/dx/dx*sqrt(2.0)/5;

    double newpir = pr * 4/5;

    g->gpnum += n;
    pdata_t*  pad;
    
     
    double massflowrate = 0.6;
    double pv  = Vinflow*massflowrate/4/M_PI/pr/pr;
    double x, y, z, d_x, d_y, d_z, dr;
    for(int i=0;i<n;i++)
    if(g->mpirank == i%g->mpisize){	
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
    
   	       
        pad->v[0] = pv*d_x/dr;
        pad->v[1] = pv*d_y/dr;
        pad->v[2] = pv*d_z/dr;
        pad->xyz[0] = x;
        pad->xyz[1] = y;
        pad->xyz[2] = z;

        pad->volume = Vinflow;
        pad->pressure = Pinflow;
        pad->localspacing = dx;
        pad->mass = 19950523;
        pad->ifboundary = true;

        pad->soundspeed = m_pEOS->getSoundSpeed(pad->pressure,1./pad->volume);
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

