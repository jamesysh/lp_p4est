#include <iostream>

#include "sc.h"
#include "mpi.h"
#include "particle_data.h"
#include "geometry_pellet.h"
using namespace std;

Global_Data:: Global_Data(Initializer* init){

    initlevel = init->initlevel;
    maxlevel = init->maxlevel;
    elem_particles = init->elem_particles;
    geometry = GeometryFactory::instance().createGeometry("pelletlayer"); 
    geometry->getBoundingBox(bb[0],bb[1],bb[2],bb[3],bb[4],bb[5]);
}


Global_Data:: ~Global_Data(){}
