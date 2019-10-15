#include "pellet_solver.h"


PelletSolver::PelletSolver(Global_Data*g){
    gdata = g;
    sc_array_copy(particle_data_copy, gdata->particle_data);

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
  octant_data_t          *qud;
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
void PelletSolver::build_quadtree(){
    conn = p4est_connectivity_new_unitsquare();
    p4est_heating = p4est_new_ext(gdata->mpicomm, conn,0,0,1,sizeof(octant_data_t), NULL, this);  
    
    }


void PelletSolver::prerun(){
    
   for (int i = 0; i < 2; ++i) {
    ilh[i] = sc_array_new (sizeof (p4est_locidx_t));
    jlh[i] = sc_array_new (sizeof (p4est_locidx_t));
    klh[i] = NULL;
   
   }
    
    }


void PelletSolver::resetOctantData2d(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  octant_data_t          *qud;
  for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est_heating->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qud->premain = qud->preceive = 0;
    }
  }

}

void PelletSolver::presearch2d(){
        
  gdata->pfound = sc_array_new_count (sizeof (int), particle_data_copy->elem_count);
  
  gdata->iremain = sc_array_new (sizeof (p4est_locidx_t));
  sc_array_memset (gdata->pfound, -1);

  p4est_search_all (p4est_heating, 0, psearch_quad2d, psearch_point2d, particle_data_copy);



}


