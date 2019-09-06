#include <iostream>

#include "sc.h"
#include "mpi.h"
#include "particle_data.h"
#include "geometry_pellet.h"
using namespace std;

static bool ifPointInsideBox(double x, double y, double z, double bb[6]) {
    
    if(x>bb[0] && x<bb[1] && y>bb[2] && y<bb[3] 

#ifdef P4_TO_P8
            && z>bb[4] && z<bb[5]
#endif
            )
        return true;
    else
        return false;
}

static void adjustCoordByDomain(double dl, double xyz[3]){
    
    for(int i=0;i<3;i++){
        xyz[i] *= dl;
        xyz[i] -= dl/2;
    }
    return;
}

static bool ifOctantInsectBox(double lxyz[3],double bb[6],double l) //l:lenth of octant, lxyz:coord of octant
{
    int i,j,k;
    bool ifin;
    for(i = 0; i<2; i++){
        for( j = 0;j<2;j++){
            for( k = 0; k<2; k++){
                ifin = ifPointInsideBox(lxyz[0]+i*l,lxyz[1]+j*l,lxyz[2]+k*l,bb);
                if(ifin)
                    return true;
            }
        }
    
    }
    return false;
}

static void createParticlesInOctant(p4est_iter_volume_info_t * info, void *user_data){
    double l; //actuall length of a octant
    bool iffill; //if fill octant with fluid particles;
    Global_Data* g = (Global_Data*) user_data;
    
    p4est_topidx_t      tt = info->treeid;  /**< the tree containing \a quad */
    p4est_quadrant_t   *quad = info->quad;
    double domain_len = g->domain_len;
    double* bb = g->bb;
    double bb_quad[6];
    octant_data_t *qud = (octant_data_t *)quad->p.user_data;
    p4est_qcoord_t qh;
    
     
    qh = P4EST_QUADRANT_LEN (quad->level);
    l = qh/(double)P4EST_ROOT_LEN*domain_len;
    p4est_qcoord_to_vertex (g->conn, tt, quad->x, quad->y,
#ifdef P4_TO_P8
                          quad->z,
#endif
                          g->lxyz);
   
    adjustCoordByDomain(domain_len,g->lxyz);
    
    iffill = ifOctantInsectBox(g->lxyz,bb, l);
    if(!iffill)
        return;
   //TO DO FILL WITH PARTICLE DATA 

};

Global_Data:: Global_Data(Initializer* init){

    initlevel = init->initlevel;
    maxlevel = init->maxlevel;
    initlocalspacing = init->initlocalspacing;
    elem_particles = init->elem_particles;
    geometry = GeometryFactory::instance().createGeometry("pelletlayer"); 
    geometry->getBoundingBox(bb[0],bb[1],bb[2],bb[3],bb[4],bb[5]);
}


Global_Data:: ~Global_Data(){}


void Global_Data::initFluidParticles(){

  p4est_iterate(p4est,NULL,(void *)this,createParticlesInOctant,NULL,NULL,NULL); 

}
