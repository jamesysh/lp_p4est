#include "boundary_gresho.h"
#include <iostream>
#include <cmath>
#include <cassert>
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
