#ifndef __OCTREE_MANAGER__
#define __OCTREE_MANAGER__

#include "p4est_to_p8est.h"

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_build.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#else

#include <p8est_bits.h>
#include <p8est_build.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#endif

#include "sc.h"

#include "particle_data.h"
class Octree_Manager{
    
    public:
        Octree_Manager(Global_Data *g);
        ~Octree_Manager(){};

        Global_Data *gdata;
    
        void build_octree();

        void destroy_octree();



};



#endif