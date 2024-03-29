#include <iostream>
#include <time.h>
#include <fstream>
#include "sc.h"
#include "mpi.h"
#include "particle_data.h"
#include "geometry_pellet.h"
#include "sc_notify.h"
#include "hexagonal_packing.h"
#include <math.h>
#include <cassert>
using namespace std;


Global_Data:: Global_Data(Initializer* init){
    
    lpnum = 0;
    gpnum = 0;
    gplost = 0; 
    flagrefine = 1;
    domain_len = init->getDomainlength();
    dimension = init->getDimension();
    initlevel = init->getInitLevel();
    timesearchingradius = init->getTimeSearchRadius();
    maxlevel = init->getMaxLevel();
    minlevel = init->getMinLevel(); 
    initlocalspacing = init->getInitParticleSpacing();
    numrow1st = init->getNumRow1stOrder();
    numrow2nd = init->getNumRow2ndOrder();
    numcol1st = init->getNumCol1stOrder();
    numcol2nd = init->getNumCol2ndOrder(); 
    initperturbation = init->getInitialPerturbation();
    elem_particle = init->getElemParticle();
    string temp;
    temp = init->getGeometryName();
    geometry = GeometryFactory::instance().createGeometry(temp); 
    geometry->getBoundingBox(bb[0],bb[1],bb[2],bb[3],bb[4],bb[5]);
    temp = init->getStateName();
    state = StateFactory::instance().createState(temp);
    boundarynumber = init->getBoundaryNumber();
    vector<string> boundarynames = init->getBoundaryNames();
    for(size_t i=0;i<boundarynumber;i++){

    	m_vBoundary.push_back(BoundaryFactory::instance().createBoundary(boundarynames[i]));

    }
    eoschoice = init->getEOSChoice();
    gamma = init->getGamma();
    pinf = init->getPinf();
    einf = init->getEinf();
    setEOS();    

    invalidpressure = init->getInvalidPressure();
    invaliddensity = init->getInvalidDensity();
    iffreeboundary = init->getIfFreeBoundary();
    flagdelete = true;
    
    pelletnumber = init->getPelletDistribution();

}


Global_Data:: ~Global_Data(){

}
static void testcornerside2d( p4est_iter_corner_info_t * info, void *user_data){
    
    bool treeboundary = info->tree_boundary;
    octant_data_t *ghost_data = (octant_data_t *)user_data;
    sc_array_t         *sides = &(info->sides);
    size_t sidescount = sides->elem_count;
    p4est_iter_corner_side_t *sidedest, *sidesrc; 
    p4est_quadrant_t *qdest, *qsrc;
    octant_data_t *ouddest, *oudsrc;
    p4est_locidx_t quadid;
    p4est_locidx_t *neighbourid;
    for(size_t i=0;i<sidescount;i++){
       sidedest = p4est_iter_cside_array_index_int(sides,i); 
       if(sidedest->is_ghost)
           continue;
       qdest = sidedest->quad;
       ouddest = (octant_data_t *)qdest->p.user_data;
       
        if(treeboundary){
            ouddest->flagboundary = true;
            }
       for(size_t j=0;j<sidescount;j++){
       
           sidesrc = p4est_iter_cside_array_index_int(sides,j); 
           if(sidesrc->is_ghost){
                oudsrc = &ghost_data[sidesrc->quadid]; 
                if(oudsrc->poctant == 0){
                    ouddest->flagboundary = true;
                    continue;
                }
                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->ghostneighbourid,1);
                *neighbourid = sidesrc->quadid;

           }

           else{
           
                qsrc = sidesrc->quad;
                oudsrc = (octant_data_t *)qsrc->p.user_data;
                if(oudsrc->poctant == 0){
                    ouddest->flagboundary = true;
                    continue;
                }
                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->localneighbourid,1);
                *neighbourid = sidesrc->quadid;
           
           }
       
       
       }

    }

}
static void testcornerside( p8est_iter_corner_info_t * info, void *user_data){
    bool treeboundary = info->tree_boundary;
    
    octant_data_t *ghost_data = (octant_data_t *)user_data;
    sc_array_t         *sides = &(info->sides);
    size_t sidescount = sides->elem_count;
    p8est_iter_corner_side_t *sidedest, *sidesrc; 
    p8est_quadrant_t *qdest, *qsrc;
    octant_data_t *ouddest, *oudsrc;
    p4est_locidx_t quadid;
    p4est_locidx_t *neighbourid;
    for(size_t i=0;i<sidescount;i++){
    
       sidedest = p8est_iter_cside_array_index_int(sides,i); 
       if(sidedest->is_ghost)
           continue;
       qdest = sidedest->quad;
       ouddest = (octant_data_t *)qdest->p.user_data;
       
        if(treeboundary){
            ouddest->flagboundary = true;
            }
       for(size_t j=0;j<sidescount;j++){
       
           sidesrc = p8est_iter_cside_array_index_int(sides,j); 
           if(sidesrc->is_ghost){
                oudsrc = &ghost_data[sidesrc->quadid]; 
                if(oudsrc->poctant == 0){
                    ouddest->flagboundary = true;

                    continue;
                }
                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->ghostneighbourid,1);
                *neighbourid = sidesrc->quadid;

           }

           else{
           
                qsrc = sidesrc->quad;
                oudsrc = (octant_data_t *)qsrc->p.user_data;
                if(oudsrc->poctant == 0){
                    ouddest->flagboundary = true;
                    
                    continue;
                }
                
                neighbourid = (p4est_locidx_t *)sc_array_push_count(ouddest->localneighbourid,1);
                *neighbourid = sidesrc->quadid;
           
           }
       
       
       }

    }

}

static void initNeighbourArray2d(p4est_iter_volume_info_t *info, void*user_data){

    p4est_quadrant_t   *q = info->quad;
    octant_data_t       *oud = (octant_data_t *) q->p.user_data;
    oud->localneighbourid = sc_array_new(sizeof(p4est_locidx_t));
    oud->ghostneighbourid = sc_array_new(sizeof(p4est_locidx_t));

}
static void initNeighbourArray(p8est_iter_volume_info_t *info, void*user_data){

    p8est_quadrant_t   *q = info->quad;
    octant_data_t       *oud = (octant_data_t *) q->p.user_data;
    oud->localneighbourid = sc_array_new(sizeof(p4est_locidx_t));
    oud->ghostneighbourid = sc_array_new(sizeof(p4est_locidx_t));

}

static int
slocal_quad2d (p4est_t * p4est, p4est_topidx_t which_tree,
             p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
             void *point)
{
  Global_Data      *g = (Global_Data *) p4est->user_pointer;

  /* compute coordinate range of this quadrant */
  g->loopquad2d ( which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

  /* always return 1 to search particles individually */
  return 1;
}

static int
slocal_quad (p8est_t * p4est, p4est_topidx_t which_tree,
             p8est_quadrant_t * quadrant, p4est_locidx_t local_num,
             void *point)
{
  Global_Data      *g = (Global_Data *) p4est->user_pointer;

  /* compute coordinate range of this quadrant */
  g->loopquad ( which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

  /* always return 1 to search particles individually */
  return 1;
}


static int
slocal_point2d (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
              void *point)
{
  int                 i;
  char               *cf;
  size_t              zp;
  Global_Data      *g = (Global_Data *) p4est->user_pointer;
  octant_data_t          *qud;
  double             *x;
  pdata_t          *pad = (pdata_t *) point;

  /* access location of particle to be searched */
  x = pad->xyz;

  /* due to roundoff we call this even for a local leaf */
  for (i = 0; i < P4EST_DIM; ++i) {
    if (!(g->lxyz[i] <= x[i] && x[i] <= g->hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }

  if (local_num >= 0) {
    /* quadrant is a local leaf */
    /* first local match counts (due to roundoff there may be multiple) */
    zp = sc_array_position (g->prebuf, point);
    cf = (char *) sc_array_index (g->cfound, zp);
    if (!*cf) {
      /* make sure this particle is not found twice */
      *cf = 1;

      /* count this particle in its target quadrant */
      *(p4est_locidx_t *) sc_array_push (g->ireceive) = (p4est_locidx_t) zp;
      qud = (octant_data_t *) quadrant->p.user_data;
      ++qud->preceive;
    }

    /* return value will have no effect */
    return 0;
  }

  /* the leaf for this particle has not yet been found */
  return 1;
}
static int
slocal_point (p8est_t * p4est, p4est_topidx_t which_tree,
              p8est_quadrant_t * quadrant, p4est_locidx_t local_num,
              void *point)
{
  int                 i;
  char               *cf;
  size_t              zp;
  Global_Data      *g = (Global_Data *) p4est->user_pointer;
  octant_data_t          *qud;
  double             *x;
  pdata_t          *pad = (pdata_t *) point;

  /* access location of particle to be searched */
  x = pad->xyz;

  /* due to roundoff we call this even for a local leaf */
  for (i = 0; i < P8EST_DIM; ++i) {
    if (!(g->lxyz[i] <= x[i] && x[i] <= g->hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }

  if (local_num >= 0) {
    /* quadrant is a local leaf */
    /* first local match counts (due to roundoff there may be multiple) */
    zp = sc_array_position (g->prebuf, point);
    cf = (char *) sc_array_index (g->cfound, zp);
    if (!*cf) {
      /* make sure this particle is not found twice */
      *cf = 1;

      /* count this particle in its target quadrant */
      *(p4est_locidx_t *) sc_array_push (g->ireceive) = (p4est_locidx_t) zp;
      qud = (octant_data_t *) quadrant->p.user_data;
      ++qud->preceive;
    }

    /* return value will have no effect */
    return 0;
  }

  /* the leaf for this particle has not yet been found */
  return 1;
}
static int
psearch_quad (p8est_t * p4est, p4est_topidx_t which_tree,
              p8est_quadrant_t * quadrant, int pfirst, int plast,
              p4est_locidx_t local_num, void *point)
{
  Global_Data      *g = (Global_Data *) p4est->user_pointer;

  /* compute coordinate range of this quadrant */
  g->loopquad ( which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

  /* always return 1 to search particles individually */
  return 1;
}

static int
psearch_quad2d (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant, int pfirst, int plast,
              p4est_locidx_t local_num, void *point)
{
  Global_Data      *g = (Global_Data *) p4est->user_pointer;

  /* compute coordinate range of this quadrant */
  g->loopquad2d ( which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

  /* always return 1 to search particles individually */
  return 1;
}

static int
psearch_point2d (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * quadrant, int pfirst, int plast,
               p4est_locidx_t local_num, void *point)
{
  int                 i;
  int                *pfn;
  size_t              zp;
  double             *x;
  Global_Data      *g = (Global_Data *) p4est->user_pointer;
  octant_data_t          *qud;
  pdata_t          *pad = (pdata_t *) point;
  /* access location of particle to be searched */
  if(pad->ifboundary){
    if(pad->flagdelete == g->flagdelete)
    { 
        return 0;
    }
  }
  x = pad->xyz;

  /* due to roundoff we call this even for a local leaf */
  for (i = 0; i < P4EST_DIM; ++i) {
    if (!(g->lxyz[i] <= x[i] && x[i] <= g->hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }

  /* convention for entries of pfound:
     -1              particle has not yet been found
     [0 .. mpisize)  particle found on that rank, me or other
   */

  /* find process/quadrant for this particle */
  if (local_num >= 0) {
    /* quadrant is a local leaf */
    zp = sc_array_position (g->particle_data, point);
    pfn = (int *) sc_array_index (g->pfound, zp);
    /* first local match counts (due to roundoff there may be multiple) */
    if (*pfn != g->mpirank) {
      /* particle was either yet unfound, or found on another process */
      /* bump counter of particles in this local quadrant */

      *pfn = g->mpirank;
      *(p4est_locidx_t *) sc_array_push (g->iremain) = (p4est_locidx_t) zp;
      qud = (octant_data_t *) quadrant->p.user_data;
      ++qud->premain;
    }
    /* return value will have no effect, but we must return */
    return 0;
  }
  if (pfirst == plast) {
    if (pfirst == g->mpirank) {
      /* continue recursion for local branch quadrant */
      return 1;
    }
    
    /* found particle on a remote process */
    zp = sc_array_position (g->particle_data, point);
    pfn = (int *) sc_array_index (g->pfound, zp);
    /* only count match if it has not been found locally or on lower rank */
    if (*pfn < 0 || (*pfn != g->mpirank && pfirst < *pfn)) {
        *pfn = pfirst;
    }

    /* return value will have no effect, but we must return */
    return 0;
  }

  /* the process for this particle has not yet been found */
  return 1;
}


static int
psearch_point (p8est_t * p4est, p4est_topidx_t which_tree,
               p8est_quadrant_t * quadrant, int pfirst, int plast,
               p4est_locidx_t local_num, void *point)
{
  int                 i;
  int                *pfn;
  size_t              zp;
  double             *x;
  Global_Data      *g = (Global_Data *) p4est->user_pointer;
  octant_data_t          *qud;
  pdata_t          *pad = (pdata_t *) point;

  if(pad->ifboundary){
    if(pad->flagdelete == g->flagdelete)
    { 
        return 0;
    }
  }
  /* access location of particle to be searched */
  x = pad->xyz;

  /* due to roundoff we call this even for a local leaf */
  for (i = 0; i < P8EST_DIM; ++i) {
    if (!(g->lxyz[i] <= x[i] && x[i] <= g->hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }

  /* convention for entries of pfound:
     -1              particle has not yet been found
     [0 .. mpisize)  particle found on that rank, me or other
   */

  /* find process/quadrant for this particle */
  if (local_num >= 0) {
    /* quadrant is a local leaf */
    zp = sc_array_position (g->particle_data, point);
    pfn = (int *) sc_array_index (g->pfound, zp);
    /* first local match counts (due to roundoff there may be multiple) */
    if (*pfn != g->mpirank) {
      /* particle was either yet unfound, or found on another process */
      /* bump counter of particles in this local quadrant */

      *pfn = g->mpirank;
      *(p4est_locidx_t *) sc_array_push (g->iremain) = (p4est_locidx_t) zp;
      qud = (octant_data_t *) quadrant->p.user_data;
      ++qud->premain;
    }
    /* return value will have no effect, but we must return */
    return 0;
  }
  if (pfirst == plast) {
    if (pfirst == g->mpirank) {
      /* continue recursion for local branch quadrant */
      return 1;
    }
    
    /* found particle on a remote process */
    zp = sc_array_position (g->particle_data, point);
    pfn = (int *) sc_array_index (g->pfound, zp);
    /* only count match if it has not been found locally or on lower rank */
    if (*pfn < 0 || (*pfn != g->mpirank && pfirst < *pfn)) {
        *pfn = pfirst;
    }

    /* return value will have no effect, but we must return */
    return 0;
  }

  /* the process for this particle has not yet been found */
  return 1;
}

static int
part_weight2d (p4est_t * p4est,
             p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  p4est_locidx_t      ilem_particles;
  Global_Data      *g = (Global_Data *) p4est->user_pointer;
  octant_data_t          *qud = (octant_data_t *) quadrant->p.user_data;


  ilem_particles = qud->lpend - g->prevlp;

  g->prevlp = qud->lpend;
  *(int *) sc_array_index (g->src_fixed, g->qcount++) =
    (int) (ilem_particles * sizeof (pdata_t));
  return 1+ qud->fluidnum;
}
static int
part_weight (p8est_t * p4est,
             p4est_topidx_t which_tree, p8est_quadrant_t * quadrant)
{
  p4est_locidx_t      ilem_particles;
  Global_Data      *g = (Global_Data *) p4est->user_pointer;
  octant_data_t          *qud = (octant_data_t *) quadrant->p.user_data;


  ilem_particles = qud->lpend - g->prevlp;

  g->prevlp = qud->lpend;
  *(int *) sc_array_index (g->src_fixed, g->qcount++) =
    (int) (ilem_particles * sizeof (pdata_t));
  return 1 + qud->fluidnum;
  
}
static int
comm_prank_compare (const void *v1, const void *v2)
{
  return sc_int_compare (&((const comm_prank_t *) v1)->rank,
                         &((const comm_prank_t *) v2)->rank);
}

static unsigned
psend_hash (const void *v, const void *u)
{
  const comm_psend_t *ps = (const comm_psend_t *) v;


  return ps->rank;
}

static int
psend_equal (const void *v1, const void *v2, const void *u)
{
  const comm_psend_t *ps1 = (const comm_psend_t *) v1;
  const comm_psend_t *ps2 = (const comm_psend_t *) v2;


  return ps1->rank == ps2->rank;
}

static bool ifPointInsideBox(double x, double y, double z, double bb[6]) {
    
    if(x>bb[0] && x<bb[1] && y>bb[2] && y<bb[3] 

            && z>bb[4] && z<bb[5]
            )
        return true;
    else
        return false;
}

static bool ifPointInsideBox2d(double x, double y, double bb[6]) {
    
    if(x>bb[0] && x<bb[1] && y>bb[2] && y<bb[3] 

            )
        return true;
    else
        return false;
}
void Global_Data:: adjustCoordByDomain( double xyz[3]){
    double dl = domain_len;
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


static bool ifOctantInsectBox2d(double lxyz[3],double bb[6],double l) //l:lenth of octant, lxyz:coord of octant
{
    int i,j;
    bool ifin;
    for(i = 0; i<2; i++){
        for( j = 0;j<2;j++){
                ifin = ifPointInsideBox2d(lxyz[0]+i*l,lxyz[1]+j*l,bb);
                if(ifin)
                    return true;
        }
    
    }
    return false;
}
static void createParticlesInOctant2d(p4est_iter_volume_info_t * info, void *user_data){
    double l; //actuall length of a octant
    bool iffill; //if fill octant with fluid particles;
    int nump;
    Global_Data* g = (Global_Data*) user_data;
     
    Geometry * geom = g->geometry;
    State* state = g->state;
    EOS* eos = g->eos;
    double x,y;
    int i,j;
    double ls = g->initlocalspacing;
    double initr = g->initperturbation;
    p4est_locidx_t *lpnum = &g->lpnum; 
    p4est_locidx_t oldnum = *lpnum;
    p4est_topidx_t      tt = info->treeid;  /**< the tree containing \a quad */
    p4est_quadrant_t   *quad = info->quad;
    double domain_len = g->domain_len;
    double* bb = g->bb;
    pdata_t *pd;
    octant_data_t *oud = (octant_data_t *)quad->p.user_data;
    p4est_qcoord_t qh;

    oud->premain = oud->preceive =  oud->poctant = 0;
    qh = P4EST_QUADRANT_LEN (quad->level);
    l = qh/(double)P4EST_ROOT_LEN*domain_len;
    p4est_qcoord_to_vertex (g->conn2d, tt, quad->x, quad->y,
                          g->lxyz);
   
    g->adjustCoordByDomain(g->lxyz);
    
    iffill = ifOctantInsectBox2d(g->lxyz,bb, l);
    if(!iffill){
        oud->lpend = -1;
        return;
    }
    //TO DO FILL WITH PARTICLE DATA 
    nump =(int)(l/ls);
    
    for(i = 0;i<nump;i++){
        for(j=0;j<nump;j++){
                x = g->lxyz[0] + (i+rand()/(double)RAND_MAX*initr)*ls;
                y = g->lxyz[1] + (j+rand()/(double)RAND_MAX*initr)*ls;
            if(geom->operator()(x,y,0)){
       
                pd = (pdata_t *) sc_array_push_count (g->particle_data,1);
                pd->xyz[0] = x;
                pd->xyz[1] = y;
                pd->xyz[2] = 0;
                state->velocity(x,y,0,pd->v[0],pd->v[1],pd->v[2]);                
                pd->volume = 1./state->density(x,y,0);
                pd->pressure = state->pressure(x,y,0);
                pd->localspacing = ls;
                pd->mass = ls*ls/pd->volume/2*sqrt(3); 
                pd->soundspeed = eos->getSoundSpeed(pd->pressure,1./pd->volume);
                pd->ifboundary = false;
                pd->redocount = 0;
                (*lpnum) ++;
                }
        }
    }
    //   pd = (pdata_t *) sc_array_push_count(g->particle_data,nump); 
   oud->lpend = *(lpnum); 
   oud->premain = *(lpnum) - oldnum; 
   oud->poctant = oud->premain;
};


static void createParticlesInOctant(p8est_iter_volume_info_t * info, void *user_data){
    double l; //actuall length of a octant
    bool iffill; //if fill octant with fluid particles;
    int nump;
    Global_Data* g = (Global_Data*) user_data;
     
    Geometry * geom = g->geometry;
    State* state = g->state;
    EOS* eos = g->eos;
    double x,y,z;
    int i,j,k;
    double ls = g->initlocalspacing *sqrt(3)/2;
    double initr = g->initperturbation;
    p4est_locidx_t *lpnum = &g->lpnum; 
    p4est_locidx_t oldnum = *lpnum;
    p4est_topidx_t      tt = info->treeid;  /**< the tree containing \a quad */
    p8est_quadrant_t   *quad = info->quad;
    double domain_len = g->domain_len;
    double* bb = g->bb;
    pdata_t *pd;
    octant_data_t *oud = (octant_data_t *)quad->p.user_data;
    p4est_qcoord_t qh;
   


    oud->premain = oud->preceive =  oud->poctant = 0;
    qh = P8EST_QUADRANT_LEN (quad->level);
    l = qh/(double)P8EST_ROOT_LEN*domain_len;
    p8est_qcoord_to_vertex (g->conn, tt, quad->x, quad->y,
                          quad->z,
                          g->lxyz);
   
    g->adjustCoordByDomain(g->lxyz);
    
    iffill = ifOctantInsectBox(g->lxyz,bb, l);
    if(!iffill){
        oud->lpend = -1;
        return;
    }
    //TO DO FILL WITH PARTICLE DATA 
    nump =(int)(l/ls);
    ls = (double)l/nump; //corrected locapspacing;
    for(i = 0;i<nump;i++){
        for(j=0;j<nump;j++){
            for(k=0;k<nump;k++){
                x = g->lxyz[0] + (i+rand()/(double)RAND_MAX*initr)*ls;
                y = g->lxyz[1] + (j+rand()/(double)RAND_MAX*initr)*ls;
                z = g->lxyz[2] + (k+rand()/(double)RAND_MAX*initr)*ls;
            if(geom->operator()(x,y,z)){
       
                pd = (pdata_t *) sc_array_push_count (g->particle_data,1);
                pd->xyz[0] = x;
                pd->xyz[1] = y;
                pd->xyz[2] = z;
                state->velocity(x,y,z,pd->v[0],pd->v[1],pd->v[2]);                
                pd->volume = 1./state->density(x,y,z);
                pd->pressure = state->pressure(x,y,z);
                pd->localspacing = ls;
                pd->mass = ls*ls*ls/pd->volume/sqrt(2); 
                pd->soundspeed = eos->getSoundSpeed(pd->pressure,1./pd->volume);
                pd->ifboundary = false;
                pd->redocount = 0;
                (*lpnum) ++;
                }
            }
        }
    }
    //   pd = (pdata_t *) sc_array_push_count(g->particle_data,nump); 
   oud->lpend = *(lpnum); 
   oud->premain = *(lpnum) - oldnum; 
   oud->poctant = oud->premain;
};

void Global_Data::initFluidParticles_hexagonal(){

    pdata_t *pd;
	double h_r = 0.5*initlocalspacing;
    Geometry *geom  = geometry;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    lpnum = 0;
 //   p4est_locidx_t *remainid;   
    srand(1);
    if(mpirank == 0){
    if(dimension == 2){
        xmin = bb[0];
        xmax = bb[1];
        ymin = bb[2];
        ymax = bb[3];
	    HexagonalPacking2D hex2D(xmin,xmax,ymin,ymax,h_r);
    
	    size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
        hex2D.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);
        hex2D.getInitialPerturbation(initperturbation);
    
        for(size_t j=m0; j<=m1; j++) { 
            if((j+1)%2 != 0) { // odd-numbered rows 
                for(size_t k=n0_odd; k<=n1_odd; k++) { 
                    double x = hex2D.computeX(0,k);
                    double y = hex2D.computeY(j);
                    if(!geom->operator()(x,y,0)) continue;
                    
                    pd = (pdata_t *) sc_array_push_count (particle_data,1);
                    pd->xyz[0] = x;
                    pd->xyz[1] = y;
                    pd->xyz[2] = 0;
                    state->velocity(x,y,0,pd->v[0],pd->v[1],pd->v[2]);                
                    pd->volume = 1./state->density(x,y,0);
                    pd->pressure = state->pressure(x,y,0);
                    pd->localspacing = initlocalspacing;
                    pd->mass = sqrt(3)*2*h_r*h_r/pd->volume; 
                    pd->soundspeed = eos->getSoundSpeed(pd->pressure,1./pd->volume);
                    pd->ifboundary = false;
                    pd->redocount = 0;
                    
                  //  remainid = (p4est_locidx_t *) sc_array_push_count(iremain,1);
                  //  *remainid = lpnum;
                    lpnum ++;
                }
            }
            else{ // even-numbered rows
                for(size_t k=n0_even; k<=n1_even; k++) {
                    double x = hex2D.computeX(1,k);
//                                                if(x-0.6*h_r>0) x=x-0.6*h_r;
//                                                double multi=1.0;
//                                                if(x>0) multi=8.0;
//                                                x*=multi;
                    double y = hex2D.computeY(j);
                    if(!geom->operator()(x,y,0)) continue; 	
//                                                if((!objs[p]->operator()(8*x,y,0))&&x>0)        continue;
                
                    pd = (pdata_t *) sc_array_push_count (particle_data,1);
                    pd->xyz[0] = x;
                    pd->xyz[1] = y;
                    pd->xyz[2] = 0;
                    state->velocity(x,y,0,pd->v[0],pd->v[1],pd->v[2]);                
                    pd->volume = 1./state->density(x,y,0);
                    pd->pressure = state->pressure(x,y,0);
                    pd->localspacing = initlocalspacing;
                    pd->mass = sqrt(3)*2*h_r*h_r/pd->volume; 
                    pd->soundspeed = eos->getSoundSpeed(pd->pressure,1./pd->volume);
                    pd->ifboundary = false;
                    pd->redocount = 0;
                    
                 //   remainid = (p4est_locidx_t *) sc_array_push_count(iremain,1);
                  //  *remainid = lpnum;
                    lpnum ++;
                    }
                }
            }
         }
   //TO do 3d 
        else if(dimension == 3){
            xmin = bb[0];
            xmax = bb[1];
            ymin = bb[2];
            ymax = bb[3];
            zmin = bb[4];
            zmax = bb[5];
			HexagonalPacking3D hex3D(xmin, xmax, ymin, ymax, zmin, zmax, h_r);
			hex3D.getInitialPerturbation(initperturbation);
			size_t l0,l1;
			size_t m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
			size_t nn0_odd, nn1_odd, nn0_even, nn1_even; 
			hex3D.getParameters(l0, l1, m0_odd, m1_odd, m0_even, m1_even, 
								n0_odd, n1_odd, n0_even, n1_even, 
								nn0_odd, nn1_odd, nn0_even, nn1_even);	
        
			for(size_t i=l0; i<=l1; i++) { 
				if((i+1)%2 != 0) { //odd-numbered layers
					for(size_t j=m0_odd; j<=m1_odd; j++) { 
						if((j+1)%2 != 0) { //odd-numbered rows 
							for(size_t k=n0_odd; k<=n1_odd; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);
//                                                if(x-0.3*h_r>-0.0001) x=x-0.3*h_r;
								if(!geom->operator()(x,y,z)) continue;	
								
                                    pd = (pdata_t *) sc_array_push_count (particle_data,1);
                                    pd->xyz[0] = x;
                                    pd->xyz[1] = y;
                                    pd->xyz[2] = z;
                                    state->velocity(x,y,z,pd->v[0],pd->v[1],pd->v[2]);                
                                    pd->volume = 1./state->density(x,y,z);
                                    pd->pressure = state->pressure(x,y,z);
                                    pd->localspacing = initlocalspacing;
                                    pd->mass = sqrt(2)*4*h_r*h_r*h_r/pd->volume; 
                                    pd->soundspeed = eos->getSoundSpeed(pd->pressure,1./pd->volume);
                                    pd->ifboundary = false;
                                    pd->redocount = 0;
                                    lpnum ++;
							}
						} 
						else{ //even-numbered rows
							for(size_t k=n0_even; k<=n1_even; k++) {
								double x = hex3D.computeX(1,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);
//                                                if(x-0.3*h_r>-0.0001) x=x-0.3*h_r;
								if(!geom->operator()(x,y,z)) continue;	
                                    pd = (pdata_t *) sc_array_push_count (particle_data,1);
                                    pd->xyz[0] = x;
                                    pd->xyz[1] = y;
                                    pd->xyz[2] = z;
                                    state->velocity(x,y,z,pd->v[0],pd->v[1],pd->v[2]);                
                                    pd->volume = 1./state->density(x,y,z);
                                    pd->pressure = state->pressure(x,y,z);
                                    pd->localspacing = initlocalspacing;
                                    pd->mass = sqrt(2)*4*h_r*h_r*h_r/pd->volume; 
                                    pd->soundspeed = eos->getSoundSpeed(pd->pressure,1./pd->volume);
                                    pd->ifboundary = false;
                                    pd->redocount = 0;
                                    lpnum ++;
							}
						}
					}
						
				} 
				else { //even-numbered layers
					for(size_t j=m0_even; j<=m1_even; j++) { 
						if((j+1)%2 != 0) { //odd-numbered rows
							for(size_t k=nn0_odd; k<=nn1_odd; k++) { 
								double x = hex3D.computeX(1,k);
								double y = hex3D.computeY(1,j);
								double z = hex3D.computeZ(i);
//                                                if(x-0.3*h_r>-0.0001) x=x-0.3*h_r;
								if(!geom->operator()(x,y,z)) continue;
								
                                    pd = (pdata_t *) sc_array_push_count (particle_data,1);
                                    pd->xyz[0] = x;
                                    pd->xyz[1] = y;
                                    pd->xyz[2] = z;
                                    state->velocity(x,y,z,pd->v[0],pd->v[1],pd->v[2]);                
                                    pd->volume = 1./state->density(x,y,z);
                                    pd->pressure = state->pressure(x,y,z);
                                    pd->localspacing = initlocalspacing;
                                    pd->mass = sqrt(2)*4*h_r*h_r*h_r/pd->volume; 
                                    pd->soundspeed = eos->getSoundSpeed(pd->pressure,1./pd->volume);
                                    pd->ifboundary = false;
                                    pd->redocount = 0;
                                    lpnum ++;
							}
						} 
						else { //even-numbered rows
							for(size_t k=nn0_even; k<=nn1_even; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(1,j);
								double z = hex3D.computeZ(i);
//                                                if(x-0.3*h_r>-0.0001) x=x-0.3*h_r;
								if(!geom->operator()(x,y,z)) continue; 	
								
                                    pd = (pdata_t *) sc_array_push_count (particle_data,1);
                                    pd->xyz[0] = x;
                                    pd->xyz[1] = y;
                                    pd->xyz[2] = z;
                                    state->velocity(x,y,z,pd->v[0],pd->v[1],pd->v[2]);                
                                    pd->volume = 1./state->density(x,y,z);
                                    pd->pressure = state->pressure(x,y,z);
                                    pd->localspacing = initlocalspacing;
                                    pd->mass = sqrt(2)*4*h_r*h_r*h_r/pd->volume; 
                                    pd->soundspeed = eos->getSoundSpeed(pd->pressure,1./pd->volume);
                                    pd->ifboundary = false;
                                    pd->redocount = 0;
                                    lpnum ++;
							}
						}
					}	
				}    
			}
        
        
        
        }
    

    gpnum = lpnum;

        }

    MPI_Bcast(&gpnum,1,MPI_INT,0,mpicomm);

    P4EST_GLOBAL_ESSENTIALF ("Created %lld fluid particles \n",   (long long) gpnum);
    }

void Global_Data::initFluidParticles_distributed(){
   
   int mpiret;
   
   p4est_locidx_t      li;
   pdata_t          *pad;
   p4est_gloidx_t  gpoffset;
   srand(time(NULL));   

   
   if(dimension == 3) 
       p8est_iterate(p8est,NULL,(void *)this,createParticlesInOctant,NULL,NULL,NULL); 
   else if(dimension == 2)
       p4est_iterate(p4est,NULL,(void *)this,createParticlesInOctant2d,NULL,NULL); 
   
    p4est_gloidx_t gnum = 0,lnum = (p4est_gloidx_t)lpnum; 
    

    mpiret = sc_MPI_Allreduce (&lnum, &gnum, 1, P4EST_MPI_GLOIDX, sc_MPI_SUM, mpicomm);
    SC_CHECK_MPI (mpiret);
    gpnum = gnum;
    P4EST_GLOBAL_ESSENTIALF ("Created %lld fluid particles \n",   (long long) gpnum);


//create global index for particles

    mpiret = sc_MPI_Exscan (&lnum, &gpoffset, 1, P4EST_MPI_GLOIDX,
                          sc_MPI_SUM, mpicomm);
    SC_CHECK_MPI (mpiret);

    if(mpirank == 0)
          gpoffset = 0;
    pad = (pdata_t *) sc_array_index_begin (particle_data);
    for(li = 0; li<lpnum; li++){
    
        (pad++)->id = gpoffset + li;
    }

}
void Global_Data:: cleanUpArrays(){


   sc_array_destroy_null(&particle_data);

//   sc_array_destroy_null(&iremain);
   
 //  sc_array_destroy_null(&ireceive);
//sc_array_destroy_null(&pfound);

  for (int i = 0; i < 2; ++i) {
    sc_array_destroy_null (&ilh[i]);
    sc_array_destroy_null (&jlh[i]);
    if(dimension == 3)
        sc_array_destroy_null (&klh[i]);
  }
}


void Global_Data:: writeVTKFiles(){
    static bool FIRST = true;
    static int timestep = 0;
    size_t lpnum = particle_data->elem_count;
    size_t li;
    pdata_t *pad;
    string filename = "output_" +to_string(mpirank)+"_"+to_string(timestep)+".vtk";
    
	FILE *outfile;
	outfile = fopen(filename.c_str(), "w");
	if(outfile==nullptr) {
		printf("Unable to open file: %s\n",filename.c_str()); 
		return;
	}
	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.16g\n",0.0);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	
	fprintf(outfile,"POINTS %ld double\n",(long int)lpnum);
    
    pad = (pdata_t *)sc_array_index_begin(particle_data);
    for(li = 0; li<lpnum; li++){
    fprintf(outfile,"%.16g %.16g %.16g\n",pad->xyz[0],pad->xyz[1],pad->xyz[2]);
    pad++ ;
    }
  
    fprintf(outfile,"POINT_DATA %ld\n",(long int)lpnum);

	fprintf(outfile,"VECTORS Velocity double\n");
    
    pad = (pdata_t *)sc_array_index_begin(particle_data);
    for(li = 0; li<lpnum; li++){
    fprintf(outfile,"%.16g %.16g %.16g\n",pad->v[0],pad->v[1],pad->v[2]);
    pad++ ;
    }

	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
    
    pad = (pdata_t *)sc_array_index_begin(particle_data);
    for(li = 0; li<lpnum; li++){
    fprintf(outfile,"%.16g\n",pad->pressure);
    pad++ ;
    }
    
    fclose(outfile);
    if(mpirank == 0){
    FILE *visitfile;
    string fname = "output_data.visit";
    visitfile = fopen(fname.c_str(),"a");
   
	if(visitfile==nullptr) {
		printf("Unable to open file: %s\n",fname.c_str()); 
		return;
	}
    if(FIRST){
        FIRST = false;
        fprintf(visitfile,"!NBLOCKS %d\n",mpisize);
    
        }
        for(int i=0;i<mpisize;i++){
            string name ="output_" + to_string(i) + "_" + to_string(timestep) +".vtk"+"\n";
            fprintf(visitfile,name.c_str());
            
            }
    
    fclose(visitfile);
    }   
    timestep ++;
}



void Global_Data::sc_array_paste (sc_array_t * dest, sc_array_t * src)
{
  P4EST_ASSERT (dest->elem_size == src->elem_size);
  P4EST_ASSERT (dest->elem_count == src->elem_count);

  memcpy (dest->array, src->array, src->elem_count * src->elem_size);
}

void * Global_Data::sc_array_index_begin (sc_array_t * arr)
{
  P4EST_ASSERT (arr != NULL);

  if (arr->elem_count == 0) {
    return NULL;
  }

  P4EST_ASSERT (arr->array != NULL);
  return (void *) arr->array;
}

void Global_Data::setEOS(){

	if(eoschoice == 1) // Polytropic gas EOS
    {   
         eos = new PolytropicGasEOS(gamma);
	     std::vector<double> eos_parameters;
	     eos->getParameters(eos_parameters);
    }

    else {
		cout<<"The choice of EOS does not exist!!! Please correct the input file."<<endl;
		assert(false);
	}

}


void Global_Data::loopquad2d (p4est_topidx_t tt, p4est_quadrant_t * quad,double lxyz[3], double hxyz[3], double dxyz[3]){

    
  int                 i;
  p4est_qcoord_t      qh;
  qh = P4EST_QUADRANT_LEN (quad->level);
  p4est_qcoord_to_vertex (conn2d, tt, quad->x, quad->y,    lxyz);
  p4est_qcoord_to_vertex (conn2d, tt, quad->x + qh, quad->y + qh,  hxyz);
  adjustCoordByDomain(lxyz);
  adjustCoordByDomain(hxyz);
  for (i = 0; i < 2; ++i) {

    dxyz[i] = hxyz[i] - lxyz[i];
   
  }

}


void Global_Data::loopquad (p4est_topidx_t tt, p8est_quadrant_t * quad,double lxyz[3], double hxyz[3], double dxyz[3]){

    
  int                 i;
  p4est_qcoord_t      qh;
  qh = P8EST_QUADRANT_LEN (quad->level);
  p8est_qcoord_to_vertex (conn, tt, quad->x, quad->y,  quad->z,   lxyz);
  p8est_qcoord_to_vertex (conn, tt, quad->x + qh, quad->y + qh,quad->z + qh,  hxyz);
  adjustCoordByDomain(lxyz);
  adjustCoordByDomain(hxyz);
  for (i = 0; i < 3; ++i) {

    dxyz[i] = hxyz[i] - lxyz[i];
   
  }

}


void Global_Data::split_by_coord ( sc_array_t * in,
                sc_array_t * out[2], pa_mode_t mode, int component,
                const double lxyz[3], const double dxyz[3])
{
  p4est_locidx_t      ppos;
  const double       *x;
  size_t              zz, znum;
  pdata_t          *pad;

  sc_array_truncate (out[0]);
  sc_array_truncate (out[1]);

  znum = in->elem_count;
  for (zz = 0; zz < znum; ++zz) {
    ppos = *(p4est_locidx_t *) sc_array_index (in, zz);
    if (mode == PA_MODE_REMAIN) {
      pad = (pdata_t *) sc_array_index (particle_data, ppos);
      x = pad->xyz;
    }
    else if (mode == PA_MODE_RECEIVE) {
      pad = (pdata_t *) sc_array_index (prebuf, ppos);
      x = pad->xyz;
    }
    else {
      P4EST_ASSERT (mode == PA_MODE_LOCATE);
      pad = (pdata_t *) sc_array_index (particle_data, ppos);
      x = pad->xyz;
    }
    if (x[component] <= lxyz[component] + .5 * dxyz[component]) {
      *(p4est_locidx_t *) sc_array_push (out[0]) = ppos;
    }
    else {
      *(p4est_locidx_t *) sc_array_push (out[1]) = ppos;
    }
  }
}
void Global_Data::prerun(){

   


   particle_data = sc_array_new(sizeof( pdata_t ));
   for (int i = 0; i < 2; ++i) {
    ilh[i] = sc_array_new (sizeof (p4est_locidx_t));
    jlh[i] = sc_array_new (sizeof (p4est_locidx_t));
    if(dimension == 3) 
        klh[i] = sc_array_new (sizeof (p4est_locidx_t));
    else 
        (klh[i] = NULL);
   }
}

void Global_Data::presearch(){
        
  pfound = sc_array_new_count (sizeof (int), particle_data->elem_count);
  
  iremain = sc_array_new (sizeof (p4est_locidx_t));
  sc_array_memset (pfound, -1);

  p8est_search_all (p8est, 0, psearch_quad, psearch_point, particle_data);



}


void Global_Data::presearch2d(){
        
  pfound = sc_array_new_count (sizeof (int), particle_data->elem_count);
  
  iremain = sc_array_new (sizeof (p4est_locidx_t));
  sc_array_memset (pfound, -1);

  p4est_search_all (p4est, 0, psearch_quad2d, psearch_point2d, particle_data);



}


void Global_Data::packParticles(){
  int                 mpiret;
  int                 retval;
  int                *pfn;
  size_t              zz, numz;
  void              **hfound;
  p4est_locidx_t      lremain, lsend, llost;
  p4est_gloidx_t      loclrs[4], glolrs[4];
  comm_psend_t       *cps, *there;
  comm_prank_t       *trank;
  pdata_t          *pad;
  psmem = sc_mempool_new (sizeof (comm_psend_t));
  numz = pfound->elem_count;

  psend = sc_hash_new (psend_hash, psend_equal, NULL, NULL);
  recevs = sc_array_new (sizeof (comm_prank_t));

  lremain = lsend = llost = 0;
  cps = (comm_psend_t *) sc_mempool_alloc (psmem);
  cps->rank = -1;

  for (zz = 0; zz < numz; ++zz) {
    pfn = (int *) sc_array_index (pfound, zz);

    /* ignore those that leave the domain or stay local */
    if (*pfn < 0) {
      assert(*pfn == -1);
      ++llost;
      continue;
    }
    if (*pfn == mpirank) {
      ++lremain;
      continue;
    }
    cps->rank = *pfn;
    retval = sc_hash_insert_unique (psend, cps, &hfound);
  
    there = *((comm_psend_t **) hfound);
  
    if (!retval) {
      /* message for this rank exists already */
      assert (there->message.elem_size == sizeof(pdata_t));
      assert (there->message.elem_count > 0);
    }
  
    else {
      /* message is added for this rank */
      assert (there == cps);
      trank = (comm_prank_t *) sc_array_push (recevs);
      trank->rank = there->rank;
      trank->psend = there;
      sc_array_init (&there->message, sizeof(pdata_t));
      cps = (comm_psend_t *) sc_mempool_alloc (psmem);
      cps->rank = -1;
    }
  
    pad = (pdata_t *) sc_array_push (&there->message);
    memcpy (pad, sc_array_index (particle_data, zz), sizeof (pdata_t));
  
    ++lsend;
  
  } 
  
  sc_mempool_free (psmem, cps);
  sc_array_sort (recevs, comm_prank_compare);

  loclrs[0] = (p4est_gloidx_t) lremain;
  loclrs[1] = (p4est_gloidx_t) lsend;
  loclrs[2] = (p4est_gloidx_t) llost;
  loclrs[3] = (p4est_gloidx_t) recevs->elem_count;
  mpiret = sc_MPI_Allreduce (loclrs, glolrs, 4, P4EST_MPI_GLOIDX,
                                sc_MPI_SUM, mpicomm);
  SC_CHECK_MPI (mpiret);

  P4EST_GLOBAL_ESSENTIALF
    ("In current timestep from %lld remain %lld sent %lld lost %lld number of particles, each processer receives data from average of %f peers\n",
      (long long) gpnum, (long long) glolrs[0],
     (long long) glolrs[1], (long long) glolrs[2],
     glolrs[3] / (double)mpisize);
  assert (glolrs[0] + glolrs[1] + glolrs[2] == gpnum);
  gplost += glolrs[2];
  gpnum -= glolrs[2];


  sc_array_destroy_null (&pfound);

}



void Global_Data::communicateParticles(){

  int                 mpiret;
  int                 i;
  int                 num_receivers;
  int                 num_senders;
  int                 count, cucount;
  int                 msglen;
  sc_MPI_Request     *reqs;
  sc_array_t         *notif, *payl;
  sc_array_t         *arr;
  comm_psend_t       *cps;
  comm_prank_t       *trank;

  num_receivers = (int) recevs->elem_count;

  notif = sc_array_new_count (sizeof (int), num_receivers);
  payl = sc_array_new_count (sizeof (int), num_receivers);   //payload


  for (i = 0; i < num_receivers; ++i) {

    trank = (comm_prank_t *) sc_array_index_int (recevs, i);

    *(int *) sc_array_index_int (notif, i) = trank->rank;

    cps = trank->psend;
  
    assert(cps->rank == trank->rank);

    arr = &cps->message;
 
    *(int *) sc_array_index_int (payl, i) = (int) arr->elem_count;
 
  
  }


  sc_notify_ext (notif, NULL, payl, NULL, mpicomm);

  assert (payl->elem_count == notif->elem_count);

  num_senders = (int) notif->elem_count;

  cucount = 0;
  for (i = 0; i < num_senders; ++i) {
    cucount += *(int *) sc_array_index_int (payl, i);
  }
  prebuf = sc_array_new_count (sizeof(pdata_t), cucount);

  /* post non-blocking receive */
  recv_req = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);

  cucount = 0;

  for (i = 0; i < num_senders; ++i) {
    count = *(int *) sc_array_index_int (payl, i);
    msglen = count * (int) sizeof(pdata_t);
    mpiret = sc_MPI_Irecv
      (sc_array_index (prebuf, cucount), msglen, sc_MPI_BYTE,
       *(int *) sc_array_index_int (notif, i), COMM_TAG_PART, mpicomm,
       (sc_MPI_Request *) sc_array_index_int (recv_req, i));
    SC_CHECK_MPI (mpiret);
    cucount += count;
  }
  assert(cucount == (int) prebuf->elem_count);

  sc_array_destroy_null (&notif);
  sc_array_destroy_null (&payl);

    send_req = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
    for (i = 0; i < num_receivers; ++i) {
      trank = (comm_prank_t *) sc_array_index_int (recevs, i);
      cps = trank->psend;
      arr = &cps->message;
      msglen = (int) (arr->elem_count * arr->elem_size);
      mpiret = sc_MPI_Isend
        (arr->array, msglen, sc_MPI_BYTE, cps->rank, COMM_TAG_PART,
         mpicomm, (sc_MPI_Request *) sc_array_index_int (send_req, i));
      SC_CHECK_MPI (mpiret);
    }

  reqs = (sc_MPI_Request *) sc_array_index_begin (recv_req);
  mpiret = sc_MPI_Waitall (num_senders, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&recv_req);

  num_receivers = (int) recevs->elem_count;
  reqs = (sc_MPI_Request *) sc_array_index_begin (send_req),
    mpiret = sc_MPI_Waitall (num_receivers, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&send_req);

  for (i = 0; i < num_receivers; ++i) {
    trank = (comm_prank_t *) sc_array_index_int (recevs, i);
    cps = trank->psend;
    sc_array_reset (&cps->message);
  }
  sc_array_destroy_null (&recevs);
  sc_hash_destroy (psend);

  psend = NULL;
  sc_mempool_destroy (psmem);
  psmem = NULL;
}

void Global_Data::postsearch(){
  ireceive = sc_array_new (sizeof (p4est_locidx_t));
  cfound = sc_array_new_count (sizeof (char), prebuf->elem_count);

  sc_array_memset (cfound, 0);

  p8est_search_local (p8est, 0, slocal_quad, slocal_point, prebuf);
  
  sc_array_destroy_null (&cfound);

/*
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      if(qud->preceive>0 )
          printf("%d %d %d \n",qud->premain,qud->preceive,mpirank);
        
    }
  }
*/
}

void Global_Data::postsearch2d(){
  ireceive = sc_array_new (sizeof (p4est_locidx_t));
  cfound = sc_array_new_count (sizeof (char), prebuf->elem_count);

  sc_array_memset (cfound, 0);

  p4est_search_local (p4est, 0, slocal_quad2d, slocal_point2d, prebuf);
  
  sc_array_destroy_null (&cfound);

/*
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      if(qud->preceive>0 )
          printf("%d %d %d \n",qud->premain,qud->preceive,mpirank);
        
    }
  }
*/
}

void Global_Data::regroupParticles2d(){

  sc_array_t         *newpa;
  p4est_topidx_t      tt;
  p4est_locidx_t      newnum;
  p4est_locidx_t      ppos;
  p4est_locidx_t      lq, prev;
  p4est_locidx_t      qboth, li;
  p4est_locidx_t     *premain, *preceive;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  pdata_t          *pad;

  newnum =
    (p4est_locidx_t) (iremain->elem_count + ireceive->elem_count);
  lpnum = newnum;

  premain = (p4est_locidx_t *) sc_array_index_begin (iremain);
  preceive = (p4est_locidx_t *) sc_array_index_begin (ireceive);
  newpa = sc_array_new_count (sizeof (pdata_t), newnum);
  pad = (pdata_t *) sc_array_index_begin (newpa);
  prev = 0;
  
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    
      tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qboth = qud->premain + qud->preceive;
      
      qud->fluidnum = 0;
      if (qboth == 0) {
        qud->lpend = prev;
        qud->premain = qud->preceive = 0;
        qud->poctant = 0;
        continue;
      }
      prev += qboth;
      for (li = 0; li < qud->premain; ++li) {
        ppos = *premain++;
        memcpy (pad, sc_array_index (particle_data, ppos), sizeof (pdata_t));
        if(!pad->ifboundary)
            qud->fluidnum ++;
        pad ++;
      }
      for (li = 0; li < qud->preceive; ++li) {
        ppos = *preceive++;
        memcpy (pad, sc_array_index (prebuf, ppos), sizeof (pdata_t));
        if(!pad->ifboundary)
            qud->fluidnum ++;
        pad ++;
      }
      qud->lpend = prev;
      qud->poctant = qboth;
      qud->premain = qud->preceive = 0;
    }
    
  } 
   
  sc_array_destroy_null (&iremain);

  sc_array_destroy_null (&prebuf);
  sc_array_destroy_null (&ireceive);
  sc_array_destroy (particle_data);

  particle_data = newpa;
}

void Global_Data::regroupParticles(){

  sc_array_t         *newpa;
  p4est_topidx_t      tt;
  p4est_locidx_t      newnum;
  p4est_locidx_t      ppos;
  p4est_locidx_t      lq, prev;
  p4est_locidx_t      qboth, li;
  p4est_locidx_t     *premain, *preceive;
  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  pdata_t          *pad;

  newnum =
    (p4est_locidx_t) (iremain->elem_count + ireceive->elem_count);
  lpnum = newnum;

  premain = (p4est_locidx_t *) sc_array_index_begin (iremain);
  preceive = (p4est_locidx_t *) sc_array_index_begin (ireceive);
  newpa = sc_array_new_count (sizeof (pdata_t), newnum);
  pad = (pdata_t *) sc_array_index_begin (newpa);
  prev = 0;
  
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    
      tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qboth = qud->premain + qud->preceive;
      qud->fluidnum = 0;
      if (qboth == 0) {
        qud->lpend = prev;
        qud->premain = qud->preceive = 0;
        qud->poctant = 0;
        continue;
      }
      prev += qboth;
      for (li = 0; li < qud->premain; ++li) {
        ppos = *premain++;
        memcpy (pad, sc_array_index (particle_data, ppos), sizeof (pdata_t));
        if(!pad->ifboundary)
            qud->fluidnum ++;
        pad ++;
      }
      for (li = 0; li < qud->preceive; ++li) {
        ppos = *preceive++;
        memcpy (pad, sc_array_index (prebuf, ppos), sizeof (pdata_t));
       
        if(!pad->ifboundary)
            qud->fluidnum ++;
        pad ++;
      }
      qud->lpend = prev;
      qud->poctant = qboth;
      qud->premain = qud->preceive = 0;
    }
    
  } 
   
  sc_array_destroy_null (&iremain);

  sc_array_destroy_null (&prebuf);
  sc_array_destroy_null (&ireceive);
  sc_array_destroy (particle_data);

  particle_data = newpa;
}

void Global_Data:: partitionParticles2d(){


  sc_array_t         *dest_data;
  p4est_topidx_t      tt;
  p4est_locidx_t      ldatasize, lcount;
  p4est_locidx_t      dest_quads, src_quads;
  p4est_locidx_t      dest_parts;
  p4est_locidx_t      lquad, lq;
  p4est_locidx_t      lpnum;
  p4est_gloidx_t      gshipped;
  p4est_gloidx_t     *src_gfq;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
   
  if(mpisize == 1)
  
      return;

  src_gfq = P4EST_ALLOC (p4est_gloidx_t, mpisize + 1);

  memcpy (src_gfq, p4est->global_first_quadrant,
          (mpisize + 1) * sizeof (p4est_gloidx_t));


  src_quads = p4est->local_num_quadrants;

  assert(src_quads == src_gfq[mpirank+1]-src_gfq[mpirank]);

  src_fixed = sc_array_new_count (sizeof (int), src_quads);

  qcount = 0;
  prevlp = 0;

  gshipped = p4est_partition_ext (p4est, 1, part_weight2d);
  dest_quads = p4est->local_num_quadrants;

  if (gshipped == 0) {
    sc_array_destroy_null (&src_fixed);
    P4EST_FREE (src_gfq);
    return;
  }

  dest_fixed = sc_array_new_count (sizeof (int), dest_quads);

  p4est_transfer_fixed (p4est->global_first_quadrant, src_gfq,
                        mpicomm, COMM_TAG_FIXED,
                        (int *) dest_fixed->array,
                        (const int *) src_fixed->array, sizeof (int));

  ldatasize = (p4est_locidx_t) sizeof (pdata_t);

  dest_parts = 0;

  for (lq = 0; lq < dest_quads; ++lq) {
    dest_parts += *(int *) sc_array_index (dest_fixed, lq);
  }
  assert(dest_parts % ldatasize == 0); 
  dest_parts /= ldatasize;
  dest_data = sc_array_new_count (sizeof (pdata_t), dest_parts);
  p4est_transfer_custom (p4est->global_first_quadrant, src_gfq,
                         mpicomm, COMM_TAG_CUSTOM,
                         (pdata_t *) dest_data->array,
                         (const int *) dest_fixed->array,
                         (const pdata_t *) particle_data->array,
                         (const int *) src_fixed->array);

  sc_array_destroy_null (&src_fixed);

  P4EST_FREE (src_gfq);
  sc_array_destroy (particle_data);
  particle_data = dest_data;
  lpnum = 0;
  lquad = 0;
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      /* access quadrant */
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;

      /* back out particle count in quadrant from data size */
      lcount = *(int *) sc_array_index (dest_fixed, lquad);
      assert (lcount % ldatasize == 0);
      lcount /= ldatasize;
      lpnum += lcount;
      qud->lpend = lpnum;
      ++lquad;
    }
  }
  sc_array_destroy_null (&dest_fixed);


}

void Global_Data:: partitionParticles(){


  sc_array_t         *dest_data;
  p4est_topidx_t      tt;
  p4est_locidx_t      ldatasize, lcount;
  p4est_locidx_t      dest_quads, src_quads;
  p4est_locidx_t      dest_parts;
  p4est_locidx_t      lquad, lq;
  p4est_locidx_t      lpnum;
  p4est_gloidx_t      gshipped;
  p4est_gloidx_t     *src_gfq;
  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
   
  if(mpisize == 1)
  
      return;

  src_gfq = P4EST_ALLOC (p4est_gloidx_t, mpisize + 1);

  memcpy (src_gfq, p8est->global_first_quadrant,
          (mpisize + 1) * sizeof (p4est_gloidx_t));


  src_quads = p8est->local_num_quadrants;

  assert(src_quads == src_gfq[mpirank+1]-src_gfq[mpirank]);

  src_fixed = sc_array_new_count (sizeof (int), src_quads);

  qcount = 0;
  prevlp = 0;

  gshipped = p8est_partition_ext (p8est, 1, part_weight);
  dest_quads = p8est->local_num_quadrants;

  if (gshipped == 0) {
    sc_array_destroy_null (&src_fixed);
    P4EST_FREE (src_gfq);
    return;
  }

  dest_fixed = sc_array_new_count (sizeof (int), dest_quads);

  p8est_transfer_fixed (p8est->global_first_quadrant, src_gfq,
                        mpicomm, COMM_TAG_FIXED,
                        (int *) dest_fixed->array,
                        (const int *) src_fixed->array, sizeof (int));

  ldatasize = (p4est_locidx_t) sizeof (pdata_t);

  dest_parts = 0;

  for (lq = 0; lq < dest_quads; ++lq) {
    dest_parts += *(int *) sc_array_index (dest_fixed, lq);
  }
  assert(dest_parts % ldatasize == 0); 
  dest_parts /= ldatasize;
  dest_data = sc_array_new_count (sizeof (pdata_t), dest_parts);
  p8est_transfer_custom (p8est->global_first_quadrant, src_gfq,
                         mpicomm, COMM_TAG_CUSTOM,
                         (pdata_t *) dest_data->array,
                         (const int *) dest_fixed->array,
                         (const pdata_t *) particle_data->array,
                         (const int *) src_fixed->array);

  sc_array_destroy_null (&src_fixed);

  P4EST_FREE (src_gfq);
  sc_array_destroy (particle_data);
  particle_data = dest_data;
  lpnum = 0;
  lquad = 0;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      /* access quadrant */
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;

      /* back out particle count in quadrant from data size */
      lcount = *(int *) sc_array_index (dest_fixed, lquad);
      assert (lcount % ldatasize == 0);
      lcount /= ldatasize;
      lpnum += lcount;
      qud->lpend = lpnum;
      ++lquad;
    }
  }
  sc_array_destroy_null (&dest_fixed);


}

void Global_Data::testquad2d(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad, *quad2;
  octant_data_t          *qud,*qud2;
  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  pdata_copy_t *pad2;
  neighbour_info_t * ninfo;
/*  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      if(qud->poctant && lq>2000){
          qud->flagboundary = 2000;
        size_t cc = qud->localneighbourid->elem_count;
        for(size_t i=0;i<cc;i++){
            p4est_locidx_t* qid = (p4est_locidx_t *)sc_array_index(qud->localneighbourid,i);
            if(*qid == lq)
                continue;
            quad2 = p8est_quadrant_array_index(&tree->quadrants,*qid);
            
            qud2 = (octant_data_t *) quad2->p.user_data;
            qud2->flagboundary = 1000;
        }
        break;
      }

    }
  }
*/
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
    lpend = qud->lpend;
    if(lq>=0){
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        
        if(pad->ifboundary)
            continue;
        
        printf("dest%f %f %f\n",pad->xyz[0],pad->xyz[1],0);
        // cout<<pad->ifhasghostneighbour<<" ifhasghost "<<endl; 
         // pad->flagboundary = 100;
        size_t cc = pad->neighbourleftparticle->elem_count;
        for(size_t i=0;i<cc;i++){
              
           fetchNeighbourParticle2d(pad, &pad2, pad->neighbourleftparticle, i);

           ninfo = (neighbour_info_t *)sc_array_index(pad->neighbourleftparticle,i); 
            /*
            int qid = ninfo->quadid;
            int pid = ninfo->parid;
            quad2 = p8est_quadrant_array_index(&tree->quadrants,qid);
            qud2= (octant_data_t *)quad2->p.user_data;
            pad2 = &qud2->localparticle[pid];*/
           // pad2->flagboundary = 10000;
            if(ninfo->ifghost) 
            printf("%f %f %f %d\n",pad2->xyz[0],pad2->xyz[1],0,9999);
            else if(ninfo->ifremote)
            printf("%f %f %f %d\n",pad2->xyz[0],pad2->xyz[1],0,2);
            else
            printf("%f %f %f %d\n",pad2->xyz[0],pad2->xyz[1],0,3);
        
        }
      //  pad->flagboundary = (double)qud->flagboundary;
      }
    }
       offset = lpend;  
    
    }
  }

}

void Global_Data::testquad(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad, *quad2;
  octant_data_t          *qud,*qud2;
  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  pdata_copy_t *pad2;
  neighbour_info_t * ninfo;
/*  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      if(qud->poctant && lq>2000){
          qud->flagboundary = 2000;
        size_t cc = qud->localneighbourid->elem_count;
        for(size_t i=0;i<cc;i++){
            p4est_locidx_t* qid = (p4est_locidx_t *)sc_array_index(qud->localneighbourid,i);
            if(*qid == lq)
                continue;
            quad2 = p8est_quadrant_array_index(&tree->quadrants,*qid);
            
            qud2 = (octant_data_t *) quad2->p.user_data;
            qud2->flagboundary = 1000;
        }
        break;
      }

    }
  }
*/
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
    lpend = qud->lpend;
    if(lq>=0){
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        if(pad->ifboundary)
            continue;
        if(pad->ifhasghostneighbour == 0) 
            continue;
        
        printf("dest%f %f %f\n",pad->xyz[0],pad->xyz[1],pad->xyz[2]);
        
        // cout<<pad->ifhasghostneighbour<<" ifhasghost "<<endl; 
         // pad->flagboundary = 100;
        
        size_t cc = pad->neighbourleftparticle->elem_count;
        for(size_t i=0;i<cc;i++){
           fetchNeighbourParticle(pad, &pad2, pad->neighbourleftparticle, i);

           ninfo = (neighbour_info_t *)sc_array_index(pad->neighbourleftparticle,i); 
            /*
            int qid = ninfo->quadid;
            int pid = ninfo->parid;
            quad2 = p8est_quadrant_array_index(&tree->quadrants,qid);
            qud2= (octant_data_t *)quad2->p.user_data;
            pad2 = &qud2->localparticle[pid];*/
           // pad2->flagboundary = 10000;
            if(ninfo->ifghost) 
            printf("%f %f %f %d\n",pad2->xyz[0],pad2->xyz[1],pad2->xyz[2],9999);
            else if(ninfo->ifremote)
            printf("%f %f %f %d\n",pad2->xyz[0],pad2->xyz[1],pad2->xyz[2],2);
            else
            printf("%f %f %f %d\n",pad2->xyz[0],pad2->xyz[1],pad2->xyz[2],3);
        
        }
      //  pad->flagboundary = (double)qud->flagboundary;
      }
    }
       offset = lpend;  
    
    }
  }

}


void Global_Data::resetOctantData2d(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qud->premain = qud->preceive = 0;
    }
  }

}

void Global_Data::resetOctantData(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qud->premain = qud->preceive = 0;
    }
  }

}


void Global_Data::createViewForOctant2d(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;
    p4est_locidx_t      offset = 0;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  pdata_t *pads;
  pdata_copy_t *padd;
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qud->flagboundary = false;  
      qud->poctant = qud->lpend - offset;
        if(qud->poctant >= 400){
            printf("This octant has more than 400 particles, please enlarge size of localparticle!!\n");
            assert(false);
        }
//      qud->particle_data_view = sc_array_new_count(sizeof(pdata_t),(size_t)qud->poctant);
//      padd = (pdata_t *) sc_array_index_begin(qud->particle_data_view);
      
      pads = (pdata_t *) sc_array_index(particle_data,(size_t)offset);
      for(size_t i=0;i<(size_t) qud->poctant;i++){
        padd = &qud->localparticle[i];
        copyParticle( padd, pads);
        pads++;
      
      } 
          offset = qud->lpend;
    
    }
  }

}

void Global_Data::setFlagBoundaryForParticle(){
    
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;

  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      lpend = qud->lpend;
      for(int i=offset;i<lpend;i++){
         pad = (pdata_t *)sc_array_index(particle_data,i);
         if(pad->ifboundary)
             continue;
         pad->flagboundary = qud->flagboundary;
            }
       offset = lpend;  
    
        }
    }
    
}


void Global_Data::createViewForOctant(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;
    p4est_locidx_t      offset = 0;
  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  pdata_t *pads;
  pdata_copy_t *padd;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qud->flagboundary = false;  
      qud->poctant = qud->lpend - offset;
        if(qud->poctant >= 400){
            printf("This octant has more than 400 particles, please enlarge size of localparticle!!\n");
            assert(false);
        }
//      qud->particle_data_view = sc_array_new_count(sizeof(pdata_t),(size_t)qud->poctant);
//      padd = (pdata_t *) sc_array_index_begin(qud->particle_data_view);
      
      pads = (pdata_t *) sc_array_index(particle_data,(size_t)offset);
      for(size_t i=0;i<(size_t) qud->poctant;i++){
        padd = &qud->localparticle[i];
        copyParticle( padd, pads);
        pads++;
      
      } 
          offset = qud->lpend;
    
    }
  }

}


void Global_Data::cleanForTimeStep2d(){

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;

  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      sc_array_destroy(qud->localneighbourid);
      sc_array_destroy(qud->ghostneighbourid);
    
      lpend = qud->lpend;
      for(int i=offset;i<lpend;i++){
         pad = (pdata_t *)sc_array_index(particle_data,i);
         if(pad->ifboundary)
             continue;
         sc_array_destroy(pad->neighbourparticle);
         sc_array_destroy(pad->neighbourrightparticle);
         sc_array_destroy(pad->neighbourleftparticle);
         sc_array_destroy(pad->neighbourfrontparticle);
         sc_array_destroy(pad->neighbourbackparticle);
         if(iffreeboundary)
            sc_array_destroy(pad->ghostneighbour); 
         }
       offset = lpend;  
    
    }
  }

    p4est_ghost_destroy (ghost2d);
    P4EST_FREE (ghost_data);
    ghost2d = NULL;
    ghost_data = NULL;
}


void Global_Data::cleanForTimeStep(){

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;

  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      sc_array_destroy(qud->localneighbourid);
      sc_array_destroy(qud->ghostneighbourid);
    
      lpend = qud->lpend;
      for(int i=offset;i<lpend;i++){
         pad = (pdata_t *)sc_array_index(particle_data,i);
         if(pad->ifboundary)
             continue;
         sc_array_destroy(pad->neighbourparticle);

 
         sc_array_destroy(pad->neighbourupparticle);
         sc_array_destroy(pad->neighbourdownparticle);
         sc_array_destroy(pad->neighbourrightparticle);
         sc_array_destroy(pad->neighbourleftparticle);
         sc_array_destroy(pad->neighbourfrontparticle);
         sc_array_destroy(pad->neighbourbackparticle);
         if(iffreeboundary)
            sc_array_destroy(pad->ghostneighbour); 
         }
       offset = lpend;  
    
    }
  }

    p8est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;
}
static void copyNeighbourInfo(neighbour_info_t * d, neighbour_info_t *s){
    d->ifremote = s->ifremote;
    d->distance = s->distance;
    d->ifghost = s->ifghost;
    d->theta = s->theta;
    d->phi = s->phi;
    d->sigma = s->sigma;
    d->parid = s->parid;
    d->quadid = s->quadid;

}

void Global_Data::copyParticle(pdata_copy_t *d, pdata_t *s){
    d->xyz[0] = s->xyz[0];
    d->xyz[1] = s->xyz[1];
    d->xyz[2] = s->xyz[2];
    d->v[0] = s->v[0];
    d->v[1] = s->v[1];
    d->v[2] = s->v[2];
    d->pressure = s->pressure;
    d->soundspeed = s->soundspeed;
    d->temperature = s->temperature;
    d->volume = s->volume;
    d->mass = s->mass;
    d->localspacing = s->localspacing;
    d->ifboundary = s->ifboundary;
    d->flagboundary = s->flagboundary;
}


static int
compareoctant (const void *p1, const void *p2)
{
  p4est_locidx_t                 i1 = *(p4est_locidx_t *) p1;
  p4est_locidx_t                 i2 = *(p4est_locidx_t *) p2;

  return i1 - i2;
}
static int compareneighbour_info(const void *p1,const void *p2)
{
    neighbour_info_t         *  i1 = (neighbour_info_t *)p1;

    neighbour_info_t         *  i2 = (neighbour_info_t *)p2;
    double d1 = i1->distance;
    double d2 = i2->distance;
    if(d1 >d2)
        return 1;
    else
        return 0;
}

void Global_Data::searchNeighbourOctant2d(){

  ghost2d = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  ghost_data = P4EST_ALLOC (octant_data_t, ghost2d->ghosts.elem_count);
  p4est_ghost_exchange_data (p4est, ghost2d, ghost_data);
  p4est_iterate(p4est,ghost2d,(void*)ghost_data,initNeighbourArray2d,NULL,testcornerside2d);

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;

  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      sc_array_sort(qud->localneighbourid,compareoctant);
      sc_array_sort(qud->ghostneighbourid,compareoctant);
      sc_array_uniq(qud->localneighbourid,compareoctant);
      sc_array_uniq(qud->ghostneighbourid,compareoctant);
   
    }
  
  }

  
}

void Global_Data::searchNeighbourOctant(){

  ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FULL);
  ghost_data = P4EST_ALLOC (octant_data_t, ghost->ghosts.elem_count);
  p8est_ghost_exchange_data (p8est, ghost, ghost_data);
  p8est_iterate(p8est,ghost,(void*)ghost_data,initNeighbourArray,NULL,NULL,testcornerside);

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;

  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      sc_array_sort(qud->localneighbourid,compareoctant);
      sc_array_sort(qud->ghostneighbourid,compareoctant);
      sc_array_uniq(qud->localneighbourid,compareoctant);
      sc_array_uniq(qud->ghostneighbourid,compareoctant);
   
    }
  
  }

  
}

void Global_Data:: searchNeighbourParticle2d(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad, *quadnei;
  octant_data_t          *qud,*qudnei;
  p4est_locidx_t   offset = 0,lpend;
  p4est_locidx_t   *localneiid, *ghostneiid;
  pdata_t * pad;
  pdata_copy_t *padnei;
  size_t localsize, ghostsize;
  double *position ;
  double x,y,x0,y0,dx,dy;
  double radius, dissquared;
  neighbour_info_t * nei_info;
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      lpend = qud->lpend;
      localsize = qud->localneighbourid->elem_count;
      ghostsize = qud->ghostneighbourid->elem_count;
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        if(pad->ifboundary)
            continue;
        pad->neighbourparticle = sc_array_new(sizeof(neighbour_info_t));
        radius = pad->localspacing * timesearchingradius; 
        position = pad->xyz;
        x = position[0];
        y = position[1];
        
       // printf("end%f %f %f\n",x,y,z);
        for(size_t i=0; i<localsize; i++){ //iterate through local neighbour octants
            localneiid = (p4est_locidx_t *)sc_array_index(qud->localneighbourid,i);
            quadnei = p4est_quadrant_array_index(&tree->quadrants,*localneiid);
            qudnei = (octant_data_t *) quadnei->p.user_data;
            size_t nump = qudnei->poctant;
            for(size_t pid = 0;pid<nump; pid++){
                padnei = &qudnei->localparticle[pid];
                x0 = padnei->xyz[0]; 
                y0 = padnei->xyz[1]; 
                dx = x-x0;
                dy = y-y0;
                dissquared = dx*dx + dy*dy;
                if(dissquared == 0){
                    continue;  } // the neighbour is itslef
                if(dissquared <= radius*radius){
                    nei_info = (neighbour_info_t *) sc_array_push(pad->neighbourparticle);
                    nei_info->ifremote = false;
                    nei_info->ifghost = false;
                    nei_info->quadid = *localneiid;
                    nei_info->parid = pid;
                    nei_info->distance = sqrt(dissquared);

                    nei_info->phi = acos(dy/sqrt(dissquared));
                    nei_info->sigma = acos(dx/sqrt(dissquared));
                    //printf("%f %f %f\n",x0,y0,z0);
                }

            
            }  
        }
        
        for(size_t i=0; i<ghostsize; i++){ //iterate through local neighbour octants
            ghostneiid = (p4est_locidx_t *)sc_array_index(qud->ghostneighbourid,i);
          //  quadnei = p8est_quadrant_array_index(&tree->quadrants,*localneiid);
            qudnei = &ghost_data[*ghostneiid];
            size_t nump = qudnei->poctant;
            for(size_t pid = 0;pid<nump; pid++){
                padnei = &qudnei->localparticle[pid];
                x0 = padnei->xyz[0]; 
                y0 = padnei->xyz[1]; 
                dx = x-x0;
                dy = y-y0;
                dissquared = dx*dx + dy*dy;
                if(dissquared == 0)
                    continue;   // the neighbour is itslef
                if(dissquared <= radius*radius){
                    nei_info = (neighbour_info_t *) sc_array_push(pad->neighbourparticle);
                    nei_info->ifremote = true;
                    nei_info->ifghost =false;
                    nei_info->quadid = *ghostneiid;
                    nei_info->parid = pid;
                    nei_info->distance = sqrt(dissquared);
                    nei_info->phi = acos(dy/sqrt(dissquared));
                    nei_info->sigma = acos(dx/sqrt(dissquared));
                    // printf("%f %f %f\n",x0,y0,z0);
                }

           
            }  
        }
        sc_array_sort(pad->neighbourparticle,compareneighbour_info); 
    }
       offset = lpend;  
    }
  }

}


void Global_Data:: searchNeighbourParticle(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad, *quadnei;
  octant_data_t          *qud,*qudnei;
  p4est_locidx_t   offset = 0,lpend;
  p4est_locidx_t   *localneiid, *ghostneiid;
  pdata_t * pad;
  pdata_copy_t *padnei;
  size_t localsize, ghostsize;
  double *position ;
  double x,y,z,x0,y0,z0,dx,dy,dz;
  double radius, dissquared;
  neighbour_info_t * nei_info;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      lpend = qud->lpend;
      localsize = qud->localneighbourid->elem_count;
      ghostsize = qud->ghostneighbourid->elem_count;
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        if(pad->ifboundary)
            continue;
        pad->neighbourparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->ifhasghostneighbour = false;
        radius = pad->localspacing * timesearchingradius; 
        position = pad->xyz;
        x = position[0];
        y = position[1];
        z = position[2];
        
       // printf("end%f %f %f\n",x,y,z);
        for(size_t i=0; i<localsize; i++){ //iterate through local neighbour octants
            localneiid = (p4est_locidx_t *)sc_array_index(qud->localneighbourid,i);
            quadnei = p8est_quadrant_array_index(&tree->quadrants,*localneiid);
            qudnei = (octant_data_t *) quadnei->p.user_data;
            size_t nump = qudnei->poctant;
            for(size_t pid = 0;pid<nump; pid++){
                padnei = &qudnei->localparticle[pid];
                x0 = padnei->xyz[0]; 
                y0 = padnei->xyz[1]; 
                z0 = padnei->xyz[2]; 
                dx = x-x0;
                dy = y-y0;
                dz = z-z0;
                dissquared = dx*dx + dy*dy + dz*dz;
                if(dissquared == 0){
                    continue;  } // the neighbour is itslef
                if(dissquared <= radius*radius){
                    pad->ifhasghostneighbour = pad->ifhasghostneighbour || padnei->ifboundary;
                    nei_info = (neighbour_info_t *) sc_array_push(pad->neighbourparticle);
                    nei_info->ifremote = false;
                    nei_info->ifghost = false;
                    nei_info->quadid = *localneiid;
                    nei_info->parid = pid;
                    nei_info->distance = sqrt(dissquared);

                    nei_info->phi = acos(dy/sqrt(dissquared));
                    nei_info->theta = acos(dz/sqrt(dissquared));
                    nei_info->sigma = acos(dx/sqrt(dissquared));
                    //printf("%f %f %f\n",x0,y0,z0);
                }

            
            }  
        }
        
        for(size_t i=0; i<ghostsize; i++){ //iterate through ghost neighbour octants
            ghostneiid = (p4est_locidx_t *)sc_array_index(qud->ghostneighbourid,i);
          //  quadnei = p8est_quadrant_array_index(&tree->quadrants,*localneiid);
            qudnei = &ghost_data[*ghostneiid];
            size_t nump = qudnei->poctant;
            for(size_t pid = 0;pid<nump; pid++){
                padnei = &qudnei->localparticle[pid];
                x0 = padnei->xyz[0]; 
                y0 = padnei->xyz[1]; 
                z0 = padnei->xyz[2]; 
                dx = x-x0;
                dy = y-y0;
                dz = z-z0;
                dissquared = dx*dx + dy*dy + dz*dz;
                if(dissquared == 0)
                    continue;   // the neighbour is itslef
                if(dissquared <= radius*radius){
                    
                    pad->ifhasghostneighbour = pad->ifhasghostneighbour || padnei->ifboundary;
                    nei_info = (neighbour_info_t *) sc_array_push(pad->neighbourparticle);
                    nei_info->ifremote = true;
                    nei_info->ifghost =false;
                    nei_info->quadid = *ghostneiid;
                    nei_info->parid = pid;
                    nei_info->distance = sqrt(dissquared);
                    nei_info->phi = acos(dy/sqrt(dissquared));
                    nei_info->theta = acos(dz/sqrt(dissquared));
                    nei_info->sigma = acos(dx/sqrt(dissquared));
                    // printf("%f %f %f\n",x0,y0,z0);
                }

           
            }  
        }
        sc_array_sort(pad->neighbourparticle,compareneighbour_info); 
    }
       offset = lpend;  
    }
  }

}

void Global_Data::searchUpwindNeighbourParticle2d(){

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  neighbour_info_t * nei_info, *nei_info2;
  double phi, sigma;
  size_t numnei, neiid;
  double anglemin = 1.3734;
  double anglemax = M_PI-anglemin;
  pdata_copy_t *padcopy;

  sc_array_t *lf = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *lb = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *rf = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *rb = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *fl = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *fr = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *bl = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *br = sc_array_new(sizeof(neighbour_info_t));
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      lpend = qud->lpend;
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        if(pad->ifboundary)
            continue;
        pad->neighbourrightparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourleftparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourfrontparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourbackparticle = sc_array_new(sizeof(neighbour_info_t));
        
        numnei = pad->neighbourparticle->elem_count;
        for(neiid = 0; neiid<numnei; neiid++){
            
            nei_info = (neighbour_info_t *)sc_array_index(pad->neighbourparticle,neiid); 
            phi = nei_info->phi;
            sigma = nei_info->sigma;
            fetchParticle2d(pad, &padcopy,nei_info);
       //y- 
            if(phi >= 0 && phi <= anglemin){
                if(padcopy->xyz[0] > pad->xyz[0])  {
                    nei_info2 = (neighbour_info_t *)sc_array_push(lf);
                    copyNeighbourInfo(nei_info2,nei_info);
                }
                else
                {
                
                    nei_info2 = (neighbour_info_t *)sc_array_push(lb);
                    copyNeighbourInfo(nei_info2,nei_info);
                }
            }
           //y+
            if(phi >= anglemax && phi <=M_PI){
                if(padcopy->xyz[0] > pad->xyz[0])  {
                    nei_info2 = (neighbour_info_t *)sc_array_push(rf);
                    copyNeighbourInfo(nei_info2,nei_info);
                }
                else
                {
                
                    nei_info2 = (neighbour_info_t *)sc_array_push(rb);
                    copyNeighbourInfo(nei_info2,nei_info);
                }
            }
        //x-
            if(sigma >= 0 && sigma <= anglemin){
                if(padcopy->xyz[1] > pad->xyz[1])  {
                    nei_info2 = (neighbour_info_t *)sc_array_push(br);
                    copyNeighbourInfo(nei_info2,nei_info);
                }
                else
                {
                
                    nei_info2 = (neighbour_info_t *)sc_array_push(bl);
                    copyNeighbourInfo(nei_info2,nei_info);
                }
            }
         //x+   
            if(sigma >= anglemax && sigma <=M_PI){
                if(padcopy->xyz[1] > pad->xyz[1])  {
                    nei_info2 = (neighbour_info_t *)sc_array_push(fr);
                    copyNeighbourInfo(nei_info2,nei_info);
                }
                else
                {
                
                    nei_info2 = (neighbour_info_t *)sc_array_push(fl);
                    copyNeighbourInfo(nei_info2,nei_info);
                }
            }
        }
    
        setUpwindNeighbourList2d(lf, lb, pad->neighbourleftparticle);
        setUpwindNeighbourList2d(rf, rb, pad->neighbourrightparticle);
        setUpwindNeighbourList2d(fl, fr, pad->neighbourfrontparticle);
        setUpwindNeighbourList2d(bl, br, pad->neighbourbackparticle);
        sc_array_reset(lf); 
        sc_array_reset(lb); 
        sc_array_reset(rf); 
        sc_array_reset(rb); 
        sc_array_reset(bl); 
        sc_array_reset(br); 
        sc_array_reset(fl); 
        sc_array_reset(fr); 
    }
     
    offset = lpend; 
    }
  
  }
    sc_array_destroy(lf);
    sc_array_destroy(lb);
    sc_array_destroy(rf);
    sc_array_destroy(rb);
    sc_array_destroy(fl);
    sc_array_destroy(fr);
    sc_array_destroy(bl);
    sc_array_destroy(br);
}
/*
void Global_Data::searchUpwindNeighbourParticle2d(){

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  neighbour_info_t * nei_info, *nei_info2;
  double phi, sigma;
  size_t numnei, neiid;

  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      lpend = qud->lpend;
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        if(pad->ifboundary)
            continue;
        pad->neighbourrightparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourleftparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourfrontparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourbackparticle = sc_array_new(sizeof(neighbour_info_t));
        
        numnei = pad->neighbourparticle->elem_count;
        for(neiid = 0; neiid<numnei; neiid++){
            
            nei_info = (neighbour_info_t *)sc_array_index(pad->neighbourparticle,neiid); 
            phi = nei_info->phi;
            sigma = nei_info->sigma;
        
            if(phi >= 0 && phi <= 0.96){
                nei_info2 = (neighbour_info_t *)sc_array_push(pad->neighbourleftparticle);
                copyNeighbourInfo(nei_info2,nei_info);
            }
            if(phi >= 2.19 && phi <=M_PI){
                nei_info2 = (neighbour_info_t *)sc_array_push(pad->neighbourrightparticle);
                copyNeighbourInfo(nei_info2,nei_info);
            }
        
            if(sigma >= 0 && sigma <= 0.96){
                nei_info2 = (neighbour_info_t *)sc_array_push(pad->neighbourbackparticle);
                copyNeighbourInfo(nei_info2,nei_info);
            }
            if(sigma >= 2.19 && sigma <=M_PI){
                nei_info2 = (neighbour_info_t *)sc_array_push(pad->neighbourfrontparticle);
                copyNeighbourInfo(nei_info2,nei_info);
            }
        }
    }
  
    offset = lpend; 
    }
  
  }
}
*/

void Global_Data::searchUpwindNeighbourParticle(){

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  neighbour_info_t * nei_info, *nei_info2;
  double theta, phi, sigma;
  size_t numnei, neiid;
  double angle_min = 1.374;
  double angle_max = M_PI-angle_min;
  double anglemin, anglemax;
  double anglefactor;
  pdata_copy_t *padcopy;

  sc_array_t *fru = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *frd = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *flu = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *fld = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *bru = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *brd = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *blu = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *bld = sc_array_new(sizeof(neighbour_info_t));

  sc_array_t *rfu = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *rfd = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *lfu = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *lfd = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *rbu = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *rbd = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *lbu = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *lbd = sc_array_new(sizeof(neighbour_info_t));
  
  sc_array_t *ufr = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *dfr = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *ufl = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *dfl = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *ubr = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *dbr = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *ubl = sc_array_new(sizeof(neighbour_info_t));
  sc_array_t *dbl = sc_array_new(sizeof(neighbour_info_t));
  
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      lpend = qud->lpend;
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        if(pad->ifboundary)
            continue;
        pad->neighbourupparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourdownparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourrightparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourleftparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourfrontparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourbackparticle = sc_array_new(sizeof(neighbour_info_t));
        
        numnei = pad->neighbourparticle->elem_count;
        for(neiid = 0; neiid<numnei; neiid++){
            
            nei_info = (neighbour_info_t *)sc_array_index(pad->neighbourparticle,neiid); 
            anglefactor = max(nei_info->distance/pad->localspacing/3.8*2.,1.);
            anglemin = angle_min/anglefactor;
            anglemax = M_PI-anglemin;
            theta = nei_info->theta;
            phi = nei_info->phi;
            sigma = nei_info->sigma;
            fetchParticle(pad, &padcopy,nei_info);
            
            if(theta >= 0 && theta <=anglemin){
                if(padcopy->xyz[0] <= pad->xyz[0]){
                    if(padcopy->xyz[1] <= pad->xyz[1]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(dbl);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(dbr);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                }
                else{
                    if(padcopy->xyz[1] <= pad->xyz[1]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(dfl);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(dfr);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                } 
            }

            if(theta >= anglemax && theta <=M_PI){
                if(padcopy->xyz[0] <= pad->xyz[0]){
                    if(padcopy->xyz[1] <= pad->xyz[1]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(ubl);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(ubr);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                }
                else{
                    if(padcopy->xyz[1] <= pad->xyz[1]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(ufl);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(ufr);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                } 
            }
        
        
            if(phi >= 0 && phi <= anglemin){
                if(padcopy->xyz[0] <= pad->xyz[0]){
                    if(padcopy->xyz[2] <= pad->xyz[2]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(lbd);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(lbu);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                }
                else{
                    if(padcopy->xyz[2] <= pad->xyz[2]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(lfd);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(lfu);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                } 
            }
            if(phi >= anglemax && phi <=M_PI){
                if(padcopy->xyz[0] <= pad->xyz[0]){
                    if(padcopy->xyz[2] <= pad->xyz[2]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(rbd);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(rbu);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                }
                else{
                    if(padcopy->xyz[2] <= pad->xyz[2]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(rfd);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(rfu);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                } 
            }
        
            if(sigma >= 0 && sigma <= anglemin){
                if(padcopy->xyz[2] <= pad->xyz[2]){
                    if(padcopy->xyz[1] <= pad->xyz[1]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(bld);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(brd);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                }
                else{
                    if(padcopy->xyz[1] <= pad->xyz[1]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(blu);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(bru);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                } 
            }
            if(sigma >= anglemax && sigma <=M_PI){
                if(padcopy->xyz[2] <= pad->xyz[2]){
                    if(padcopy->xyz[1] <= pad->xyz[1]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(fld);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(frd);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                }
                else{
                    if(padcopy->xyz[1] <= pad->xyz[1]){
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(flu);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                    else{
                    
                        nei_info2 = (neighbour_info_t *)sc_array_push(fru);
                        copyNeighbourInfo(nei_info2,nei_info);
                    }
                } 
            }
        }
        
        setUpwindNeighbourList(fru,frd,flu,fld,pad->neighbourfrontparticle);
        setUpwindNeighbourList(bru,brd,blu,bld,pad->neighbourbackparticle);
        setUpwindNeighbourList(rfu,rfd,rbu,rbd,pad->neighbourrightparticle);
        setUpwindNeighbourList(lfu,lfd,lbu,lbd,pad->neighbourleftparticle);
        setUpwindNeighbourList(ufr,ubr,ufl,ubl,pad->neighbourupparticle);
        setUpwindNeighbourList(dfr,dbr,dfl,dbl,pad->neighbourdownparticle);
    
        sc_array_reset(fru);
        sc_array_reset(frd);
        sc_array_reset(flu);
        sc_array_reset(fld);
        sc_array_reset(bru);
        sc_array_reset(brd);
        sc_array_reset(blu);
        sc_array_reset(bld);
        sc_array_reset(rfu);
        sc_array_reset(rfd);
        sc_array_reset(rbu);
        sc_array_reset(rbd);
        sc_array_reset(lfu);
        sc_array_reset(lfd);
        sc_array_reset(lbu);
        sc_array_reset(lbd);
        sc_array_reset(ufr);
        sc_array_reset(ubr);
        sc_array_reset(ufl);
        sc_array_reset(ubl);
        sc_array_reset(dfr);
        sc_array_reset(dbr);
        sc_array_reset(dfl);
        sc_array_reset(dbl);
    }
  
    offset = lpend; 
    }
  
  }

        sc_array_destroy(fru);
        sc_array_destroy(frd);
        sc_array_destroy(flu);
        sc_array_destroy(fld);
        sc_array_destroy(bru);
        sc_array_destroy(brd);
        sc_array_destroy(blu);
        sc_array_destroy(bld);
        sc_array_destroy(rfu);
        sc_array_destroy(rfd);
        sc_array_destroy(rbu);
        sc_array_destroy(rbd);
        sc_array_destroy(lfu);
        sc_array_destroy(lfd);
        sc_array_destroy(lbu);
        sc_array_destroy(lbd);
        sc_array_destroy(ufr);
        sc_array_destroy(ubr);
        sc_array_destroy(ufl);
        sc_array_destroy(ubl);
        sc_array_destroy(dfr);
        sc_array_destroy(dbr);
        sc_array_destroy(dfl);
        sc_array_destroy(dbl);
}

/*
void Global_Data::searchUpwindNeighbourParticle(){

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  neighbour_info_t * nei_info, *nei_info2;
  double theta, phi, sigma;
  size_t numnei, neiid;
  double anglemin = 1.3734;
  double anglemax = M_PI-anglemin;
  
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      lpend = qud->lpend;
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        if(pad->ifboundary)
            continue;
        pad->neighbourupparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourdownparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourrightparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourleftparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourfrontparticle = sc_array_new(sizeof(neighbour_info_t));
        pad->neighbourbackparticle = sc_array_new(sizeof(neighbour_info_t));
        
        numnei = pad->neighbourparticle->elem_count;
        for(neiid = 0; neiid<numnei; neiid++){
            
            nei_info = (neighbour_info_t *)sc_array_index(pad->neighbourparticle,neiid); 
            theta = nei_info->theta;
            phi = nei_info->phi;
            sigma = nei_info->sigma;
            if(theta >= 0 && theta <=anglemin){
                nei_info2 = (neighbour_info_t *)sc_array_push(pad->neighbourdownparticle);
                copyNeighbourInfo(nei_info2,nei_info);
            }

            if(theta >= anglemax && theta <=M_PI){
                nei_info2 = (neighbour_info_t *)sc_array_push(pad->neighbourupparticle);
                copyNeighbourInfo(nei_info2,nei_info);
            }
        
        
            if(phi >= 0 && phi <= anglemin){
                nei_info2 = (neighbour_info_t *)sc_array_push(pad->neighbourleftparticle);
                copyNeighbourInfo(nei_info2,nei_info);
            }
            if(phi >= anglemax && phi <=M_PI){
                nei_info2 = (neighbour_info_t *)sc_array_push(pad->neighbourrightparticle);
                copyNeighbourInfo(nei_info2,nei_info);
            }
        
            if(sigma >= 0 && sigma <= anglemin){
                nei_info2 = (neighbour_info_t *)sc_array_push(pad->neighbourbackparticle);
                copyNeighbourInfo(nei_info2,nei_info);
            }
            if(sigma >= anglemax && sigma <=M_PI){
                nei_info2 = (neighbour_info_t *)sc_array_push(pad->neighbourfrontparticle);
                copyNeighbourInfo(nei_info2,nei_info);
            }
        }
    }
  
    offset = lpend; 
    }
  
  }
}
*/

void Global_Data::generateGhostParticle2d(){


    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  size_t count;
    
  lghostnum = 0;
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      lpend = qud->lpend;
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        if(pad->ifboundary)
            continue;
        pad->ghostneighbour =  sc_array_new(sizeof(pdata_copy_t));
        pad->ifhasghostneighbour = false;


        count = pad->neighbourrightparticle->elem_count;
        if(count<numrow1st)
            fillArrayWithGhostParticle2d(pad->neighbourrightparticle,pad,numrow1st-count,3);

        count = pad->neighbourleftparticle->elem_count;
        if(count<numrow1st)
            fillArrayWithGhostParticle2d(pad->neighbourleftparticle,pad,numrow1st-count,4);

        count = pad->neighbourfrontparticle->elem_count;
        if(count<numrow1st)
            fillArrayWithGhostParticle2d(pad->neighbourfrontparticle,pad,numrow1st-count,5);

        count = pad->neighbourbackparticle->elem_count;
        if(count<numrow1st)
            fillArrayWithGhostParticle2d(pad->neighbourbackparticle,pad,numrow1st-count,6);

        }
  
    
    offset = lpend; 
    
    
    }
    
  
  
  }
    p4est_gloidx_t lgn = (p4est_gloidx_t) lghostnum;
    gghostnum = 0;
    int mpiret = sc_MPI_Allreduce (&lgn, &gghostnum, 1, P4EST_MPI_GLOIDX, sc_MPI_SUM, mpicomm);
    SC_CHECK_MPI (mpiret);

    P4EST_GLOBAL_ESSENTIALF ("Created %lld Ghost Particles \n",   (long long) gghostnum);

}

void Global_Data::generateGhostParticle(){


    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  size_t count;
    
  lghostnum = 0;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      lpend = qud->lpend;
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        if(pad->ifboundary)
            continue;
        pad->ghostneighbour =  sc_array_new(sizeof(pdata_copy_t));
//        pad->ifhasghostneighbour = false;

        count = pad->neighbourupparticle->elem_count;
        if(count<numrow1st)
            fillArrayWithGhostParticle(pad->neighbourupparticle,pad,numrow1st-count,1);

        count = pad->neighbourdownparticle->elem_count;
        if(count<numrow1st)
            fillArrayWithGhostParticle(pad->neighbourdownparticle,pad,numrow1st-count,2);

        count = pad->neighbourrightparticle->elem_count;
        if(count<numrow1st)
            fillArrayWithGhostParticle(pad->neighbourrightparticle,pad,numrow1st-count,3);

        count = pad->neighbourleftparticle->elem_count;
        if(count<numrow1st)
            fillArrayWithGhostParticle(pad->neighbourleftparticle,pad,numrow1st-count,4);

        count = pad->neighbourfrontparticle->elem_count;
        if(count<numrow1st)
            fillArrayWithGhostParticle(pad->neighbourfrontparticle,pad,numrow1st-count,5);

        count = pad->neighbourbackparticle->elem_count;
        if(count<numrow1st)
            fillArrayWithGhostParticle(pad->neighbourbackparticle,pad,numrow1st-count,6);

        }
  
    
    offset = lpend; 
    
    
    }
    
  
  
  }
    p4est_gloidx_t lgn = (p4est_gloidx_t) lghostnum;
    gghostnum = 0;
    int mpiret = sc_MPI_Allreduce (&lgn, &gghostnum, 1, P4EST_MPI_GLOIDX, sc_MPI_SUM, mpicomm);
    SC_CHECK_MPI (mpiret);

    P4EST_GLOBAL_ESSENTIALF ("Created %lld Ghost Particles \n",   (long long) gghostnum);

}

//dir: 1 up 2 down 3 right 4 left 5 front 6 back

void Global_Data::fillArrayWithGhostParticle2d(sc_array_t * neighbourlist, pdata_t * pad, int count, int dir){
      double radius = (timesearchingradius-1) * pad->localspacing;
      double r;
      double dx, dy;
      pdata_copy_t * ghostnei;
      neighbour_info_t *ghostnei_info;
      size_t ghostid;
      pad->ifhasghostneighbour = true;
      
      lghostnum += count;
      for(int i = 0;i<count;i++){
        r = rand()/(double)RAND_MAX * radius/4 + pad->localspacing;
        if( dir == 3 || dir == 5){

            if(dir == 3){
                dy = -r;
                dx = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                ghostnei_info = (neighbour_info_t *)sc_array_push(neighbourlist);
                ghostid = pad->ghostneighbour->elem_count;
                ghostnei = (pdata_copy_t *) sc_array_push(pad->ghostneighbour);
                ghostnei_info->parid = ghostid;
                ghostnei_info->ifremote = false;
                ghostnei_info->ifghost = true;
                ghostnei_info->distance = sqrt(dx*dx+dy*dy);
                addGhostParticle(ghostnei,pad,dx,dy,0);
            }
            if(dir == 5){
                dx = -r;
                dy = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                ghostnei_info = (neighbour_info_t *)sc_array_push(neighbourlist);
                ghostid = pad->ghostneighbour->elem_count;
                ghostnei = (pdata_copy_t *) sc_array_push(pad->ghostneighbour);
                ghostnei_info->parid = ghostid;
                ghostnei_info->ifremote = false;
                ghostnei_info->ifghost = true;
                ghostnei_info->distance = sqrt(dx*dx+dy*dy);
                addGhostParticle(ghostnei,pad,dx,dy,0);
            }
        
        
        } 

        if( dir == 4 || dir == 6){

            if(dir == 4){
                dy = r;
                dx = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                ghostnei_info = (neighbour_info_t *)sc_array_push(neighbourlist);
                ghostid = pad->ghostneighbour->elem_count;
                ghostnei = (pdata_copy_t *) sc_array_push(pad->ghostneighbour);
                ghostnei_info->parid = ghostid;
                ghostnei_info->ifremote = false;
                ghostnei_info->ifghost = true;
                ghostnei_info->distance = sqrt(dx*dx+dy*dy);
                addGhostParticle(ghostnei,pad,dx,dy,0);
            }
            if(dir == 6){
                dx = r;
                
                dy = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                
                ghostnei_info = (neighbour_info_t *)sc_array_push(neighbourlist);
                ghostid = pad->ghostneighbour->elem_count;
                ghostnei = (pdata_copy_t *) sc_array_push(pad->ghostneighbour);
                ghostnei_info->parid = ghostid;
                ghostnei_info->ifremote = false;
                ghostnei_info->ifghost = true;
                ghostnei_info->distance = sqrt(dx*dx+dy*dy);
                addGhostParticle(ghostnei,pad,dx,dy,0);
            }
        
        
        } 
      
      }



}



void Global_Data::fillArrayWithGhostParticle(sc_array_t * neighbourlist, pdata_t * pad, int count, int dir){
      double anglemin = 0.14;
      double anglemax = 3.;
      double radius = (timesearchingradius-1) * pad->localspacing;
      double r;
      double angle;
      double dx, dy, dz;
      double r2;
      pdata_copy_t * ghostnei;
      neighbour_info_t *ghostnei_info;
      size_t ghostid;
      pad->ifhasghostneighbour = true;
      
      lghostnum += count;
      for(int i = 0;i<count;i++){
        r = rand()/(double)RAND_MAX * radius/4 + pad->localspacing;
        if(dir == 1 || dir == 3 || dir == 5){
            angle = anglemax + rand()/(double)RAND_MAX *(M_PI-anglemax);
            if(dir == 1){
                dz = -r;//cos(angle)*r;
                r2 = r*r - dz*dz;
                dx = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                dy = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                ghostnei_info = (neighbour_info_t *)sc_array_push(neighbourlist);
                ghostid = pad->ghostneighbour->elem_count;
                ghostnei = (pdata_copy_t *) sc_array_push(pad->ghostneighbour);
                ghostnei_info->parid = ghostid;
                ghostnei_info->ifremote = false;
                ghostnei_info->ifghost = true;
                ghostnei_info->distance = sqrt(dx*dx+dy*dy+dz*dz);
                addGhostParticle(ghostnei,pad,dx,dy,dz);
            }

            if(dir == 3){
                dy = -r;
                r2 = r*r - dy*dy;
                dz = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                dx = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                ghostnei_info = (neighbour_info_t *)sc_array_push(neighbourlist);
                ghostid = pad->ghostneighbour->elem_count;
                ghostnei = (pdata_copy_t *) sc_array_push(pad->ghostneighbour);
                ghostnei_info->parid = ghostid;
                ghostnei_info->ifremote = false;
                ghostnei_info->ifghost = true;
                ghostnei_info->distance = sqrt(dx*dx+dy*dy+dz*dz);
                addGhostParticle(ghostnei,pad,dx,dy,dz);
            }
            if(dir == 5){
                dx = -r;
                r2 = r*r - dx*dx;
                dz = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                dy = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                ghostnei_info = (neighbour_info_t *)sc_array_push(neighbourlist);
                ghostid = pad->ghostneighbour->elem_count;
                ghostnei = (pdata_copy_t *) sc_array_push(pad->ghostneighbour);
                ghostnei_info->parid = ghostid;
                ghostnei_info->ifremote = false;
                ghostnei_info->ifghost = true;
                ghostnei_info->distance = sqrt(dx*dx+dy*dy+dz*dz);
                addGhostParticle(ghostnei,pad,dx,dy,dz);
            }
        
        
        } 

        if(dir == 2 || dir == 4 || dir == 6){
            angle =  rand()/(double)RAND_MAX *anglemin;
            if(dir == 2){
                dz = r;
                r2 = r*r - dz*dz;
                dx = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                dy = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                ghostnei_info = (neighbour_info_t *)sc_array_push(neighbourlist);
                ghostid = pad->ghostneighbour->elem_count;
                ghostnei = (pdata_copy_t *) sc_array_push(pad->ghostneighbour);
                ghostnei_info->parid = ghostid;
                ghostnei_info->ifremote = false;
                ghostnei_info->ifghost = true;
                ghostnei_info->distance = sqrt(dx*dx+dy*dy+dz*dz);
                addGhostParticle(ghostnei,pad,dx,dy,dz);
            }

            if(dir == 4){
                dy = r;
                r2 = r*r - dy*dy;
                dz = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                dx = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                ghostnei_info = (neighbour_info_t *)sc_array_push(neighbourlist);
                ghostid = pad->ghostneighbour->elem_count;
                ghostnei = (pdata_copy_t *) sc_array_push(pad->ghostneighbour);
                ghostnei_info->parid = ghostid;
                ghostnei_info->ifremote = false;
                ghostnei_info->ifghost = true;
                ghostnei_info->distance = sqrt(dx*dx+dy*dy+dz*dz);
                addGhostParticle(ghostnei,pad,dx,dy,dz);
            }
            if(dir == 6){
                dx = r;
                r2 = r*r - dx*dx;
                
                dz = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                dy = rand()/(double)RAND_MAX *(2*(rand()%2)-1)* r/4;
                
                ghostnei_info = (neighbour_info_t *)sc_array_push(neighbourlist);
                ghostid = pad->ghostneighbour->elem_count;
                ghostnei = (pdata_copy_t *) sc_array_push(pad->ghostneighbour);
                ghostnei_info->parid = ghostid;
                ghostnei_info->ifremote = false;
                ghostnei_info->ifghost = true;
                ghostnei_info->distance = sqrt(dx*dx+dy*dy+dz*dz);
                addGhostParticle(ghostnei,pad,dx,dy,dz);
            }
        
        
        } 
      
      }



}


void Global_Data::fetchParticle(pdata_t* pad, pdata_copy_t **padnei, neighbour_info_t *neiinfo){
  

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  size_t parid;
  size_t quadid;

  tree = p8est_tree_array_index (p8est->trees, 0);

  if(neiinfo->ifghost){         //ghost neighbour
      parid = neiinfo->parid;
      *padnei = (pdata_copy_t*) sc_array_index(pad->ghostneighbour, parid);
      
  
  }

  else if(neiinfo -> ifremote){     //remote neighbour
      quadid = neiinfo->quadid;
      parid = neiinfo->parid;
      qud = &ghost_data[quadid];
      *padnei = &qud->localparticle[parid];     
    
  }
  
  else{              //local neighbour
      quadid = neiinfo->quadid;
      parid = neiinfo->parid;
      quad = p8est_quadrant_array_index(&tree->quadrants,quadid);
      qud = (octant_data_t *) quad->p.user_data;
      *padnei = &qud->localparticle[parid];
  }


} 

void Global_Data::fetchParticle2d(pdata_t* pad, pdata_copy_t **padnei, neighbour_info_t *neiinfo){
  

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  size_t parid;
  size_t quadid;

  tree = p4est_tree_array_index (p4est->trees, 0);

  if(neiinfo->ifghost){         //ghost neighbour
      parid = neiinfo->parid;
      *padnei = (pdata_copy_t*) sc_array_index(pad->ghostneighbour, parid);
      
  
  }

  else if(neiinfo -> ifremote){     //remote neighbour
      quadid = neiinfo->quadid;
      parid = neiinfo->parid;
      qud = &ghost_data[quadid];
      *padnei = &qud->localparticle[parid];     
    
  }
  
  else{              //local neighbour
      quadid = neiinfo->quadid;
      parid = neiinfo->parid;
      quad = p4est_quadrant_array_index(&tree->quadrants,quadid);
      qud = (octant_data_t *) quad->p.user_data;
      *padnei = &qud->localparticle[parid];
  }


} 

void Global_Data::fetchNeighbourParticle2d(pdata_t* pad, pdata_copy_t **padnei ,sc_array_t *neighbourlist, size_t index){

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  neighbour_info_t *neiinfo;
  size_t parid;
  size_t quadid;

  tree = p4est_tree_array_index (p4est->trees, 0);
  neiinfo = (neighbour_info_t *) sc_array_index(neighbourlist, index);

  if(neiinfo->ifghost){         //ghost neighbour
      parid = neiinfo->parid;
      *padnei = (pdata_copy_t*) sc_array_index(pad->ghostneighbour, parid);
      
  
  }

  else if(neiinfo -> ifremote){     //remote neighbour
      quadid = neiinfo->quadid;
      parid = neiinfo->parid;
      qud = &ghost_data[quadid];
      *padnei = &qud->localparticle[parid];     
    
  }
  
  else{              //local neighbour
      quadid = neiinfo->quadid;
      parid = neiinfo->parid;
      quad = p4est_quadrant_array_index(&tree->quadrants,quadid);
      qud = (octant_data_t *) quad->p.user_data;
      *padnei = &qud->localparticle[parid];
  }


} 


void Global_Data::fetchNeighbourParticle(pdata_t* pad, pdata_copy_t **padnei ,sc_array_t *neighbourlist, size_t index){

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  neighbour_info_t *neiinfo;
  size_t parid;
  size_t quadid;

  tree = p8est_tree_array_index (p8est->trees, 0);
  neiinfo = (neighbour_info_t *) sc_array_index(neighbourlist, index);

  if(neiinfo->ifghost){         //ghost neighbour
      parid = neiinfo->parid;
      *padnei = (pdata_copy_t*) sc_array_index(pad->ghostneighbour, parid);
      
  
  }

  else if(neiinfo -> ifremote){     //remote neighbour
      quadid = neiinfo->quadid;
      parid = neiinfo->parid;
      qud = &ghost_data[quadid];
      *padnei = &qud->localparticle[parid];     
    
  }
  
  else{              //local neighbour
      quadid = neiinfo->quadid;
      parid = neiinfo->parid;
      quad = p8est_quadrant_array_index(&tree->quadrants,quadid);
      qud = (octant_data_t *) quad->p.user_data;
      *padnei = &qud->localparticle[parid];
  }


} 



void Global_Data::addGhostParticle(pdata_copy_t * ghostnei, pdata_t *pad, double dx, double dy, double dz){
                double volume;
                double pressure;
                if(pad->flagboundary){
                    volume = 1.0e6;
                    pressure = 0.;
                    }
                else{
                    volume = 1.0e6;//pad->volume;
                    pressure = 0.;//pad->pressure;
                    }
                ghostnei->xyz[0] = pad->xyz[0] - dx; 
                ghostnei->xyz[1] = pad->xyz[1] - dy; 
                ghostnei->xyz[2] = pad->xyz[2] - dz; 
                ghostnei->v[0] = pad->v[0];
                ghostnei->v[1] = pad->v[1];
                ghostnei->v[2] = pad->v[2];
                ghostnei->pressure = pressure; 
                ghostnei->volume = volume; 
                ghostnei->soundspeed = pad->soundspeed; 
                ghostnei->mass = 0; 
                ghostnei->localspacing = pad->localspacing; 

}


void Global_Data::updateViewForOctant2d(int phase){

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;
    p4est_locidx_t      offset = 0;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  pdata_t *pads;
  pdata_copy_t *padd;
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qud->poctant = qud->lpend - offset;
      
      pads = (pdata_t *) sc_array_index(particle_data,(size_t)offset);
      for(size_t i=0;i<(size_t) qud->poctant;i++){
        
        if(pads->ifboundary){
            pads++;
            continue;
        }  
        padd = &qud->localparticle[i];
        if(phase == 0){
            padd->pressure = pads->pressureT1;
            padd->soundspeed = pads->soundspeedT1;
            padd->volume = pads->volumeT1;
        }
        
        if(phase == 1){
            padd->pressure = pads->pressure;
            padd->soundspeed = pads->soundspeed;
            padd->volume = pads->volume;
        }
        pads++;
      
      } 
          offset = qud->lpend;
    
    }
  }

    P4EST_FREE (ghost_data);
    
   // ghost_data = NULL;
    
    ghost_data = P4EST_ALLOC (octant_data_t, ghost2d->ghosts.elem_count);
    p4est_ghost_exchange_data (p4est, ghost2d, ghost_data);
}
void Global_Data::updateViewForOctant(int phase){

    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;
    p4est_locidx_t      offset = 0;
  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  pdata_t *pads;
  pdata_copy_t *padd;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qud->poctant = qud->lpend - offset;
      
      pads = (pdata_t *) sc_array_index(particle_data,(size_t)offset);
      for(size_t i=0;i<(size_t) qud->poctant;i++){
        
        if(pads->ifboundary){
            pads++;
            continue;
        }  
        padd = &qud->localparticle[i];
        if(phase == 0){
            padd->pressure = pads->pressureT1;
            padd->soundspeed = pads->soundspeedT1;
            padd->volume = pads->volumeT1;
        }
        else if(phase == 1 ){
        
            padd->pressure = pads->pressureT2;
            padd->soundspeed = pads->soundspeedT2;
            padd->volume = pads->volumeT2;
        }
        
        else if(phase == 2 ){
        
            padd->pressure = pads->pressure;
            padd->soundspeed = pads->soundspeed;
            padd->volume = pads->volume;
        }
        pads++;
      
      } 
          offset = qud->lpend;
    
    }
  }

    P4EST_FREE (ghost_data);
    
   // ghost_data = NULL;
    
    ghost_data = P4EST_ALLOC (octant_data_t, ghost->ghosts.elem_count);
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
}

void Global_Data::updateParticleStates(){


    pdata_t * pad;
    size_t li, lpnum = particle_data->elem_count;
     
    for(li = 0; li<lpnum; li++){
    
         pad = (pdata_t *)sc_array_index(particle_data,li);
        if(pad->ifboundary)
         {
             continue;
         }
        if(dimension == 3){
             swap(pad->pressure, pad->pressureT1);
             swap(pad->soundspeed,pad->soundspeedT1);
             swap(pad->volume, pad->volumeT1);
        }
        else if(dimension == 2){
             swap(pad->pressure, pad->pressureT2);
             swap(pad->soundspeed,pad->soundspeedT2);
             swap(pad->volume, pad->volumeT2);
        }
         swap(pad->v,pad->oldv);
         
         }

}

void Global_Data::setUpwindNeighbourList(sc_array_t *list0, sc_array_t *list1, sc_array_t *list2, sc_array_t *list3, sc_array_t *neidest){
    size_t n0, n1, n2, n3;
    size_t i0 = 0, i1 = 0, i2 = 0, i3 = 0;
    neighbour_info_t *nei_infos, *nei_infod;
    n0 = list0->elem_count;
    n1 = list1->elem_count;
    n2 = list2->elem_count;
    n3 = list3->elem_count;
    while(true){
       if(i0<n0){
            nei_infos = (neighbour_info_t *)sc_array_index(list0,i0);
            nei_infod = (neighbour_info_t *)sc_array_push(neidest);
            copyNeighbourInfo(nei_infod,nei_infos);
            i0 ++;
       } 
    
       if(i1<n1){
            nei_infos = (neighbour_info_t *)sc_array_index(list1,i1);
            nei_infod = (neighbour_info_t *)sc_array_push(neidest);
            copyNeighbourInfo(nei_infod,nei_infos);
            i1 ++;
       } 
       if(i2<n2){
            nei_infos = (neighbour_info_t *)sc_array_index(list2,i2);
            nei_infod = (neighbour_info_t *)sc_array_push(neidest);
            copyNeighbourInfo(nei_infod,nei_infos);
            i2 ++;
       } 
       if(i3<n3){
            nei_infos = (neighbour_info_t *)sc_array_index(list3,i3);
            nei_infod = (neighbour_info_t *)sc_array_push(neidest);
            copyNeighbourInfo(nei_infod,nei_infos);
            i3 ++;
       } 
       if(i0 == n0 && i1 == n1 && i2 == n2 && i3 == n3){
            
           assert(neidest->elem_count == n0+n1+n2+n3);
           break;

       } 
    }

}

void Global_Data::setUpwindNeighbourList2d(sc_array_t* nei0, sc_array_t *nei1, sc_array_t *neidest){
    
    sc_array_t *neishort;
    sc_array_t *neilong;
    neighbour_info_t *nei_info0, *nei_info1, *nei_infod;
    size_t count;
    if(nei0->elem_count > nei1->elem_count){
        neishort = nei1;
        neilong = nei0;
    }
    else{
    
        neishort = nei0;
        neilong = nei1;
    }

    count = neishort->elem_count;

    for(size_t i=0; i<count; i++){
        nei_info0 = (neighbour_info_t *)sc_array_index(neishort,i);
        nei_info1 = (neighbour_info_t *)sc_array_index(neilong,i);
        if(nei_info0->distance < nei_info1->distance){
        nei_infod = (neighbour_info_t *)sc_array_push(neidest);
        copyNeighbourInfo(nei_infod,nei_info0);
        nei_infod = (neighbour_info_t *)sc_array_push(neidest);
        copyNeighbourInfo(nei_infod,nei_info1);
        }
        else{
        nei_infod = (neighbour_info_t *)sc_array_push(neidest);
        copyNeighbourInfo(nei_infod,nei_info1);
        nei_infod = (neighbour_info_t *)sc_array_push(neidest);
        copyNeighbourInfo(nei_infod,nei_info0);
        }

    }
    
    for(size_t i=count;i<neilong->elem_count;i++){
    
        nei_info1 = (neighbour_info_t *)sc_array_index(neilong,i);
        nei_infod = (neighbour_info_t *)sc_array_push(neidest);
        copyNeighbourInfo(nei_infod,nei_info1);
    
    }

}

void Global_Data::switchFlagDelete(){
    
    flagdelete = !flagdelete;

}

void Global_Data::reorderNeighbourList2d(){
    
    sc_array_t *neilist;
    size_t i, j, count = particle_data->elem_count;
    size_t neisize;
    pdata_t *pad;
    pdata_copy_t *pad_copy;
    neighbour_info_t *ninfo;
    const double *xyz, *xyz0;
    double x,y,x0,y0;
    float *distance;
    float penalty[4];
    float penalty_weight = 10000;
    int region;
    for(i = 0;i<count;i++){
        pad = (pdata_t *) sc_array_index(particle_data,i);
        if(pad->ifboundary) 
            continue;

        xyz0 = pad->xyz;
        x0 = xyz0[0];
        y0 = xyz0[1];
        neilist = pad->neighbourparticle;
        
        neisize = neilist->elem_count;
       
        for(int index=0;index<4;index++){
            penalty[index] = 1;
            }

        for(j=0;j<neisize;j++){
            ninfo = (neighbour_info_t *) sc_array_index(neilist,j);
            fetchParticle2d(pad,&pad_copy,ninfo); 
            xyz = pad_copy->xyz;
            distance = &ninfo->distance;
            x = xyz[0]; 
            y = xyz[1]; 
            if(x0>=x && y0>=y)
                region = 0;
            else if(x0<x && y0>=y)
                region = 1;
            else if(x0>=x && y0<y)
                region = 2;
            else if(x0<x && y0<y )
                region = 3;
            
            *distance = *distance * penalty[region];
            penalty[region] *= penalty_weight;
            } 
    
        sc_array_sort(neilist,compareneighbour_info); 
    }

}

void Global_Data::reorderNeighbourList(){
    
    sc_array_t *neilist;
    size_t i, j, count = particle_data->elem_count;
    size_t neisize;
    pdata_t *pad;
    pdata_copy_t *pad_copy;
    neighbour_info_t *ninfo;
    const double *xyz, *xyz0;
    double x,y,z,x0,y0,z0, dx, dy, dz;
    float *distance;
    float penalty[6];
    float penalty_weight = 100;
    int region;
    for(i = 0;i<count;i++){
        pad = (pdata_t *) sc_array_index(particle_data,i);
        if(pad->ifboundary) 
            continue;

        xyz0 = pad->xyz;
        x0 = xyz0[0];
        y0 = xyz0[1];
        z0 = xyz0[2];
        neilist = pad->neighbourparticle;
        
        neisize = neilist->elem_count;
       
        for(int index=0;index<8;index++){
            penalty[index] = 1;
            }

        for(j=0;j<neisize;j++){
            ninfo = (neighbour_info_t *) sc_array_index(neilist,j);
            fetchParticle(pad,&pad_copy,ninfo); 
            xyz = pad_copy->xyz;
            distance = &ninfo->distance;
            x = xyz[0]; 
            y = xyz[1]; 
            z = xyz[2];
            dx = x-x0;
            dy = y-y0;
            dz = z-z0;

			if ((fabs(dx)>=fabs(dy))&&(fabs(dx)>=fabs(dz)))
			{
				if(dx>0)
					region=0;
				else
					region=1;
			}
            else if ((fabs(dy)>=fabs(dx))&&(fabs(dy)>=fabs(dz)))
            {
                    if(dy>0)
                            region=2;
                    else
                            region=3;
            }
            else
            {
                    if(dz>0)
                            region=4;
                    else
                            region=5;
            }
           
           /* if(x0>=x && y0>=y && z0>=z)
                region = 0;
            else if(x0<x && y0>=y && z0>=z)
                region = 1;
            else if(x0>=x && y0<y && z0>=z)
                region = 2;
            else if(x0>=x && y0>=y && z0<z)
                region = 3;
            else if(x0<x && y0<y && z0>=z)
                region = 4;
            else if(x0<x && y0>=y && z0<z)
                region = 5;
            else if(x0>=x && y0<y && z0<z)
                region = 6;
            else if(x0<x && y0<y && z0<z)
                region = 7;
            */
            
            *distance = (*distance) * penalty[region];
            penalty[region] *= penalty_weight;
            } 
    
        sc_array_sort(neilist,compareneighbour_info); 
    }

}

void Global_Data:: setParticleIDAndRank(){
    pdata_t *pad;
    size_t li, lpnum = particle_data->elem_count;
    for(li = 0; li<lpnum; li++){
        pad = (pdata_t *) sc_array_index(particle_data, li);
        pad->mpirank = mpirank;
        pad->id = li; 
        }
    
    }

