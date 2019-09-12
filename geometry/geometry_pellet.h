#ifndef __GEOMETRY_PELLET_H
#define __GEOMETRY_PELLET_H

#include "geometry.h"
#include <math.h>

class PelletLayer: public Geometry {
public:
	PelletLayer();
	virtual ~PelletLayer() {}
	virtual bool operator()(double x, double y, double z) const;
  	virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
	double xcen;
	double ycen;
	double zcen;
	double innerradius;
	double outerradius;
};


class MultiPelletLayer: public Geometry {
public:
	double number = 10;
    MultiPelletLayer();
	virtual ~MultiPelletLayer() {}
	virtual bool operator()(double x, double y, double z) const;
  	virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
    double xcen[100];
	double ycen[100];
	double zcen[100];
	double innerradius;
	double outerradius;
};



#endif
