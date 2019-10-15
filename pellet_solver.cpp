#include "pellet_solver.h"
#include <iostream>
using namespace std;
PelletSolver::PelletSolver(Global_Data*g){
    gdata = g;
    particle_data_copy = sc_array_new(sizeof(pdata_t));
    sc_array_copy(particle_data_copy, gdata->particle_data);
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
static void adjustCoordByDomain( double xyz[3], double dl){
    for(int i=0;i<3;i++){
        xyz[i] *= dl;
        xyz[i] -= dl/2;
    }
    return;
}

static void loopquad2d (PelletSolver* p, p4est_topidx_t tt, p4est_quadrant_t * quad,double lxyz[3], double hxyz[3], double dxyz[3]){

    
  int                 i;
  p4est_qcoord_t      qh;
  qh = P4EST_QUADRANT_LEN (quad->level);
  p4est_qcoord_to_vertex (p->conn, tt, quad->x, quad->y,    lxyz);
  p4est_qcoord_to_vertex (p->conn, tt, quad->x + qh, quad->y + qh,  hxyz);
  adjustCoordByDomain(lxyz,p->gdata->domain_len);
  adjustCoordByDomain(hxyz,p->gdata->domain_len);
  for (i = 0; i < 2; ++i) {

    dxyz[i] = hxyz[i] - lxyz[i];
   
  }

}

static int
psearch_quad2d (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant, int pfirst, int plast,
              p4est_locidx_t local_num, void *point)
{
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data *g = p->gdata;
  /* compute coordinate range of this quadrant */
  loopquad2d (p, which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

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
  PelletSolver *p = (PelletSolver *)p4est->user_pointer;
  Global_Data      *g = p->gdata;
  quadrant_data_t          *qud;
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

  /* convention for entries of pfound:
     -1              particle has not yet been found
     [0 .. mpisize)  particle found on that rank, me or other
   */

  /* find process/quadrant for this particle */
  if (local_num >= 0) {
    /* quadrant is a local leaf */
    zp = sc_array_position (p->particle_data_copy, point);
    pfn = (int *) sc_array_index (g->pfound, zp);
    /* first local match counts (due to roundoff there may be multiple) */
    if (*pfn != g->mpirank) {
      /* particle was either yet unfound, or found on another process */
      /* bump counter of particles in this local quadrant */

      *pfn = g->mpirank;
      *(p4est_locidx_t *) sc_array_push (g->iremain) = (p4est_locidx_t) zp;
      qud = (quadrant_data_t *) quadrant->p.user_data;
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
    zp = sc_array_position (p->particle_data_copy, point);
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

static int adapt_refine2d (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant)
{
  
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data *g = p->gdata;

  
  
  quadrant_data_t          *oud = (quadrant_data_t *) quadrant->p.user_data;
  
  
  
        


  
  
  /* we have set this to -1 in adapt_coarsen */

  if ((double) (oud->premain + oud->preceive) > p->elem_particle_box) {
    /* we are trying to refine, we will possibly go into the replace function */
    g->ire2 = g->ireindex;
    g->ireindex += oud->premain;
    g->irv2 = g->irvindex;
    g->irvindex += oud->preceive;
    
    return 1;
  }
  else {

    /* maintain cumulative particle count for next quadrant */
    g->ireindex += oud->premain;
    g->irvindex += oud->preceive;
    p4est_locidx_t *irecumu = (p4est_locidx_t *)sc_array_push(g->irecumu);
    *irecumu = g->ireindex;
    p4est_locidx_t *irvcumu = (p4est_locidx_t *)sc_array_push(g->irvcumu);
    *irvcumu = g->irvindex;
    
    return 0;
  }

}

static void adapt_replace2d (p4est_t * p4est, p4est_topidx_t which_tree,
               int num_outgoing, p4est_quadrant_t * outgoing[],
               int num_incoming, p4est_quadrant_t * incoming[]){
    
  int                 wx, wy, wz;
  double              lxyz[3], hxyz[3], dxyz[3];
  sc_array_t          iview, *arr;
  p4est_locidx_t      irem, ibeg;
  p4est_quadrant_t  **pchild;
  quadrant_data_t          *oud;
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data *g = p->gdata;
  
  if (num_outgoing == P4EST_CHILDREN) {
  
  // we are coarsening 
      oud = (quadrant_data_t *) incoming[0]->p.user_data;
    g->ireindex2 += (oud->premain = g->qremain);
    g->irvindex2 += (oud->preceive = g->qreceive);
  }

  else {
    // we are refining 
    // access parent quadrant 
    
    loopquad2d (p,which_tree, outgoing[0],lxyz,hxyz,dxyz);

    // recover window onto remaining particles for the new family 
    ibeg = g->ire2;
    irem = g->ireindex - ibeg;
    
    sc_array_init_view (&iview, g->iremain, ibeg, irem);

    // sort remaining particles into the children 
    pchild = incoming;
    g->klh[0] = &iview;
    wz = 0;

    g->split_by_coord ( g->klh[wz], g->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord (g->jlh[wy], g->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (quadrant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
      
      }
    }

    g->klh[0] = NULL;

  }

}

void PelletSolver::build_quadtree(){
    p4est_locidx_t *id;
    conn = p4est_connectivity_new_unitsquare();
    p4est_heating = p4est_new_ext(gdata->mpicomm, conn,0,0,1,sizeof(quadrant_data_t), NULL, this);  
    
    resetQuadrantData();
    gdata->ireindex = gdata->ire2 = 0;
    gdata->irvindex = gdata->irv2 = 0;
    gdata->iremain = sc_array_new_count(sizeof(p4est_locidx_t),particle_data_copy->elem_count);
    for(size_t i=0; i<particle_data_copy->elem_count; i++){
        id = (p4est_locidx_t*)sc_array_index(gdata->iremain,i);
        *id = (p4est_locidx_t )i;
        } 
    p4est_refine_ext (p4est_heating, 0, 50 , adapt_refine2d, NULL, adapt_replace2d);
    sc_array_destroy(gdata->iremain);
    
    }


void PelletSolver::prerun(){
    
   for (int i = 0; i < 2; ++i) {
    ilh[i] = sc_array_new (sizeof (p4est_locidx_t));
    jlh[i] = sc_array_new (sizeof (p4est_locidx_t));
    klh[i] = NULL;
   
   }
    
    }


void PelletSolver::resetQuadrantData(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  quadrant_data_t          *qud;
  for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est_heating->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (quadrant_data_t *) quad->p.user_data;
      qud->premain = particle_data_copy->elem_count;
      qud->preceive = 0;
    }
  }

}

void PelletSolver::presearch2d(){
        
  gdata->pfound = sc_array_new_count (sizeof (int), particle_data_copy->elem_count);
  
  gdata->iremain = sc_array_new (sizeof (p4est_locidx_t));
  sc_array_memset (gdata->pfound, -1);

  p4est_search_all (p4est_heating, 0, psearch_quad2d, psearch_point2d, particle_data_copy);



}
/*
void PelletSolver::packParticles(){
  int                 mpiret;
  int                 retval;
  int                *pfn;
  size_t              zz, numz;
  void              **hfound;
  p4est_locidx_t      lremain, lsend, llost;
  comm_psend_t       *cps, *there;
  comm_prank_t       *trank;
  pdata_t          *pad;
  
  sc_array_t *pfound = gdata->pfound; 
  sc_hash_t  *psend = gdata->psend;   
  sc_mempool_t *psmem = gdata->psmem;    
  
  
  
  
  psmem = sc_mempool_new (sizeof (comm_psend_t));
  numz = pfound->elem_count;

  psend = sc_hash_new (psend_hash, psend_equal, NULL, NULL);
  recevs = sc_array_new (sizeof (comm_prank_t));
  lremain = lsend = llost = 0;
  cps = (comm_psend_t *) sc_mempool_alloc (psmem);
  cps->rank = -1;

  for (zz = 0; zz < numz; ++zz) {
    pfn = (int *) sc_array_index (pfound, zz);

    if (*pfn < 0) {
      assert(*pfn == -1);
      ++llost;
      continue;
    }
    if (*pfn == gdata->mpirank) {
      ++lremain;
      continue;
    }
    cps->rank = *pfn;
    retval = sc_hash_insert_unique (psend, cps, &hfound);
  
    there = *((comm_psend_t **) hfound);
  
    if (!retval) {
      assert (there->message.elem_size == sizeof(pdata_t));
      assert (there->message.elem_count > 0);
    }
  
    else {
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


  sc_array_destroy_null (&pfound);

}
*/
