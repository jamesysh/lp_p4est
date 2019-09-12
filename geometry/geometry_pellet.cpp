#include "geometry_pellet.h"
#include <iostream>
#include <time.h>
using namespace std;
PelletLayer::PelletLayer(){
    
    xcen = 0;
    ycen = 0;
    zcen = 0;
    innerradius = 0.2;
    outerradius = 0.24;

}


bool PelletLayer::operator()(double x, double y, double z) const{
	double r2=(x-xcen)*(x-xcen)+(y-ycen)*(y-ycen)+(z-zcen)*(z-zcen);
	return (innerradius*innerradius<r2 && r2<outerradius*outerradius);
}


void PelletLayer::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax){
	xmin = xcen-outerradius;
	xmax = xcen+outerradius;
	ymin = ycen-outerradius;
	ymax = ycen+outerradius;
	zmin = zcen-outerradius;
	zmax = zcen+outerradius;
}

	
MultiPelletLayer::MultiPelletLayer(){
    srand(time(NULL)); 
    for(int i=0;i<100;i++){
        xcen[i] = 50*(2*(rand()/(double)RAND_MAX)-1); 
        ycen[i] = 50*(2*(rand()/(double)RAND_MAX)-1); 
        zcen[i] = 50*(2*(rand()/(double)RAND_MAX)-1); 
    }


    innerradius = 0.2;
    outerradius = 0.24;
}


bool MultiPelletLayer::operator()(double x, double y, double z) const{
   
	for(int i=0;i<100;i++){
        	double r2=(x-xcen[i])*(x-xcen[i])+(y-ycen[i])*(y-ycen[i])+(z-zcen[i])*(z-zcen[i]);
        	if(innerradius*innerradius<r2 && r2<outerradius*outerradius)
            {
                return true;
            }
	}
	return false;
    
}

void MultiPelletLayer::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax){
	xmin = -51;
	xmax = 51;
	ymin = -51;
	ymax = 51;
	zmin = -51;
	zmax = 51;
}

