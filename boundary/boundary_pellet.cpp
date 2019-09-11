#include "boundary_pellet.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include "sc.h"
#include "particle_data.h"
using namespace std;

PelletInflowBoundary::PelletInflowBoundary():Pinflow(30),Uinflow(0),Vinflow(100){}


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
     if(g->mpirank == 0){
     double offset = 2./numberofNewFluid;
     double increment = M_PI*(3-sqrt(5));
     for(int i=0;i<numberofNewFluid;i++){
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
           
           pad->soundspeed = m_pEOS->getSoundSpeed(pad->pressure,1./pad->volume);
	
           
           }
     }

}

