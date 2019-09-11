#ifndef __INITIALIZER_H__
#define __INITIALIZER_H__

#include "eos.h"
#include <p4est.h>
#include <p8est.h>

#include <p4est_bits.h>
#include <p4est_build.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p8est_bits.h>
#include <p8est_build.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>

#include "geometry.h"
#include "state.h"
#include <sc.h>
#include <cassert>
class Initializer {

    public:
    
        Initializer();
        ~Initializer(){};

        
        int initlevel ; //init level of octree
        int maxlevel ;
        int minlevel;
        int elem_particles; //max number of particles per octant
        double cfl_coefficient;
        double endt;
        double initlocalspacing;
        double initperturbation;
        int eoschoice;
        int pelletmaterial;
        int gamma;

};

#endif // __INITIALIZER_H__
