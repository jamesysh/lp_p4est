#ifndef __OCTREE_MANAGER__
#define __OCTREE_MANAGER__


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

#include "sc.h"

#include "particle_data.h"
class Octree_Manager{
    
    public:
        Octree_Manager(Global_Data *g);
        ~Octree_Manager(){};

        Global_Data *gdata;
    
        void build_octree();

        void destroy_octree();
        
        void partition_octree(int allow_for_coarsening,p8est_weight_t weight_fn);
        
        void refine_octree(int recursive, p8est_refine_t refine_fn, p8est_init_t init_fn, p8est_replace_t replace_fn);
      
        void adapt_octree();
        static int adapt_coarsen (p8est_t * p8est, p4est_topidx_t which_tree,
               p8est_quadrant_t * quadrants[]);


        static int adapt_refine (p8est_t * p8est, p4est_topidx_t which_tree,
              p8est_quadrant_t * quadrant);
        static void adapt_replace (p8est_t * p8est, p4est_topidx_t which_tree,
               int num_outgoing, p8est_quadrant_t * outgoing[],
               int num_incoming, p8est_quadrant_t * incoming[]);

        static void balance_replace (p8est_t * p8est, p4est_topidx_t which_tree,
               int num_outgoing, p8est_quadrant_t * outgoing[],
               int num_incoming, p8est_quadrant_t * incoming[]);

};



#endif
