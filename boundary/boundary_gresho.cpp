#include "boundary_gresho.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include "sc.h"
#include "particle_data.h"
using namespace std;


Gresho2DSolidBoundary::Gresho2DSolidBoundary():radius(1.), thickness(0.3) {
	bo = radius-thickness;
}

int Gresho2DSolidBoundary:: operator()(double x, double y, double z, double pressure, double vx, double vy, double vz, 
	double& xb, double& yb, double& zb, 
	double& pressureb, double& vxb, double& vyb, double& vzb){
	
	double dist = sqrt(x*x+y*y);

	if(dist < bo) return 0; // inside
	
	if(dist > radius) return 0; // outside	
	
	if(dist==0) return 0; // origin (center of circle)

	double factor = (2.*radius-dist)/dist;
//	double normal_vx = x/dist;
//	double normal_vy = y/dist;

	xb = factor*x;
	yb = factor*y;
//	pressureb.push_back(pressure);
	pressureb = 3+4*log(2);
	vxb = 0;
	vyb = 0;
	return 1;	

}


void Gresho2DSolidBoundary::generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx){
    
    size_t n, i,totalcount = 0 ,count = 0;
    double x, y;
    double vx, vy;
    double xb,yb,zb,vxb,vyb,vzb,pressureb;
    int ifin;
    pdata_t *pad, *padb;
    
    n = g->particle_data->elem_count;
        
    for(i = 0; i<n; i++){
        pad = (pdata_t *) sc_array_index(g->particle_data,i);
        if(pad->ifboundary)
            continue;
        x = pad->xyz[0];
        y = pad->xyz[1];
        vx = pad->v[0];
        vy = pad->v[1];
       
        ifin = operator()(x,y,0,pad->pressure,vx,vy,0,xb,yb,zb,pressureb,vxb,vyb,vzb);
        if(ifin){
            padb = (pdata_t *)sc_array_push(g->particle_data);
            padb->xyz[0] = xb;
            padb->xyz[1] = yb;
            padb->v[0] = vxb;
            padb->v[1] = vyb;
            padb->pressure = pressureb;
            padb->ifboundary = true;
            padb->flagdelete = !g->flagdelete;
            count ++;        
        }
    }


    sc_MPI_Allreduce (&count, &totalcount, 1, MPI_INT, sc_MPI_SUM, g->mpicomm);
    g->gpnum += totalcount;

    P4EST_GLOBAL_ESSENTIALF ("Created %lld solid boundary particles \n",   (long long) totalcount);

}


Yee2DSolidBoundary::Yee2DSolidBoundary():radius(5.), thickness(3.0) {
        bo = radius-thickness;
}

int Yee2DSolidBoundary:: operator()(double x, double y, double z, double pressure, double vx, double vy, double vz, 
	double& xb, double& yb, double& zb, 
	double& pressureb, double& vxb, double& vyb, double& vzb){
        double dist = sqrt(x*x+y*y);

        if(dist < bo) return 0; // inside

        if(dist > radius) return 0; // outside  

        if(dist==0) return 0; // origin (center of circle)

        double factor = (2.*radius-dist)/dist;

        xb = factor*x;
        yb = factor*y;

	    double r=2.*radius-dist;
	
        pressureb = (1.0);
        vxb = (2.5/M_PI*exp(0.5-0.5*r*r)*(-factor*y));
        vyb = (2.5/M_PI*exp(0.5-0.5*r*r)*(factor*x));
        return 1;

}




void Yee2DSolidBoundary::generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx){
    
    size_t n, i,totalcount = 0 ,count = 0;
    double x, y;
    double vx, vy;
    double xb,yb,zb,vxb,vyb,vzb,pressureb;
    int ifin;
    pdata_t *pad, *padb;
    n = g->particle_data->elem_count;
        
    for(i = 0; i<n; i++){
        pad = (pdata_t *) sc_array_index(g->particle_data,i);
        if(pad->ifboundary)
            continue;
        x = pad->xyz[0];
        y = pad->xyz[1];
        vx = pad->v[0];
        vy = pad->v[1];
       
        ifin = operator()(x,y,0,pad->pressure,vx,vy,0,xb,yb,zb,pressureb,vxb,vyb,vzb);
        if(ifin){
            padb = (pdata_t *)sc_array_push(g->particle_data);
            padb->xyz[0] = xb;
            padb->xyz[1] = yb;
            padb->v[0] = vxb;
            padb->v[1] = vyb;
            padb->pressure = pressureb;
            padb->ifboundary = true;
            padb->flagdelete = !g->flagdelete;
            count ++;        
        }
    }


    sc_MPI_Allreduce (&count, &totalcount, 1, MPI_INT, sc_MPI_SUM, g->mpicomm);
    g->gpnum += totalcount;

    P4EST_GLOBAL_ESSENTIALF ("Created %lld solid boundary particles \n",   (long long) totalcount);

}














