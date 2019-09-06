#ifndef __INITIALIZER_H__
#define __INITIALIZER_H__

#include "p4est_to_p8est.h"

#ifndef P4_TO_P8
#include <p4est.h>
#else
#include <p8est.h>
#endif 

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_build.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#else
#include <p4est_to_p8est.h>
#include <p8est_bits.h>
#include <p8est_build.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#endif

#include "geometry.h"
#include "state.h"

#include <sc.h>

class Initializer {

    public:
    
        Initializer();
        ~Initializer(){};

        
        int initlevel ; //init level of octree
        int maxlevel ;
        int elem_particles; //max number of particles per octant
        double cfl_coefficient;
        double endt;
        double initlocalspacing;
        double initperturbation;
};

#endif // __INITIALIZER_H__
