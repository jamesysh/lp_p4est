#include "geometry_pellet.h"

PelletLayer::PelletLayer(){
    
    xcen = 0;
    ycen = 0;
    zcen = 0;
    innerradius = 0.;
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

