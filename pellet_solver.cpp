#include "pellet_solver.h"
#include <algorithm>
#include "sc_notify.h"
#include <iostream>
using namespace std;

#define c_light 2.99792458e7
float     Bessel_Kn( int n,float  x);
float     Bessel_I0(float  x);
float     Bessel_K1(float  x);
float     Bessel_I1(float  x);


PelletSolver::PelletSolver(Initializer *init,Global_Data*g){
    gdata = g;
    elem_particle_box = init->getQuadtreeResolution();
    elem_particle_cell = init->getBinarytreeResolution();
    magneticfield = init->getMagneticField();
    gdata->pellet_solver = this;
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

  if ((double) (oud->premain + oud->preceive ) > p->elem_particle_box) {
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
    
    return 0;
  }

}

void PelletSolver::split_by_coord ( sc_array_t * in,
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
      pad = (pdata_t *) sc_array_index (particle_data_copy, ppos);
      x = pad->xyz;
    }
    else if (mode == PA_MODE_RECEIVE) {
      pad = (pdata_t *) sc_array_index (prebuf, ppos);
      x = pad->xyz;
    }
    else {
      P4EST_ASSERT (mode == PA_MODE_LOCATE);
      pad = (pdata_t *) sc_array_index (particle_data_copy, ppos);
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
static void  adapt_replace2d (p4est_t * p4est, p4est_topidx_t which_tree,
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
    p->klh[0] = &iview;
    wz = 0;

    p->split_by_coord ( p->klh[wz], p->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      p->split_by_coord (p->jlh[wy], p->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = p->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (quadrant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
      
      }
    }

    // recover window onto received particles for the new family 
    ibeg = g->irv2;
    irem = g->irvindex - ibeg;
    sc_array_init_view (&iview, g->ireceive, ibeg, irem);

    // sort received particles into the children 
    pchild = incoming;
    wz = 0;
    p->split_by_coord ( p->klh[wz], p->jlh, PA_MODE_RECEIVE, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      p->split_by_coord ( p->jlh[wy], p->ilh, PA_MODE_RECEIVE, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = p->ilh[wx];
        sc_array_init_view (&iview, g->ireceive, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (quadrant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->preceive = (p4est_locidx_t) arr->elem_count);
      }
    }
    p->klh[0] = NULL;
    assert (ibeg == g->irvindex);

  }

}
/*
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
    p->klh[0] = &iview;
    wz = 0;

    g->split_by_coord ( p->klh[wz], p->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord (p->jlh[wy], p->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        // we have a set of particles for child 4 * wz + 2 * wy + wx 
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        g->sc_array_paste (&iview, arr);
        oud = (quadrant_data_t *) (*pchild++)->p.user_data;
        ibeg += (oud->premain = (p4est_locidx_t) arr->elem_count);
      
      }
    }

    p->klh[0] = NULL;

  }

}
*/
void PelletSolver::build_quadtree(){

    particle_data_copy = sc_array_new(sizeof(pdata_t));
    sc_array_copy(particle_data_copy,gdata->particle_data);
    swapXAndZCoordinate();
    conn = p4est_connectivity_new_unitsquare();
    p4est_heating = p4est_new_ext(gdata->mpicomm, conn,1,0,1,sizeof(quadrant_data_t), NULL, this);  
    resetQuadrantData();
    
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
      qud->premain = 0;
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

void PelletSolver::packParticles(){
  int                 retval;
  int                *pfn;
  size_t              zz, numz;
  void              **hfound;
  p4est_locidx_t      lremain, lsend, llost;
  comm_psend_t       *cps, *there;
  comm_prank_t       *trank;
  pdata_t          *pad;
  
  
  
  gdata->recevs = sc_array_new (sizeof (comm_prank_t));
  
  
  gdata->psmem = sc_mempool_new (sizeof (comm_psend_t));
  numz = gdata->pfound->elem_count;

  gdata->psend = sc_hash_new (psend_hash, psend_equal, NULL, NULL);
  lremain = lsend = llost = 0;
  cps = (comm_psend_t *) sc_mempool_alloc (gdata->psmem);
  cps->rank = -1;

  for (zz = 0; zz < numz; ++zz) {
    pfn = (int *) sc_array_index (gdata->pfound, zz);

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
    retval = sc_hash_insert_unique (gdata->psend, cps, &hfound);
  
    there = *((comm_psend_t **) hfound);
  
    if (!retval) {
      assert (there->message.elem_size == sizeof(pdata_t));
      assert (there->message.elem_count > 0);
    }
  
    else {
      assert (there == cps);
      trank = (comm_prank_t *) sc_array_push (gdata->recevs);

      trank->rank = there->rank;
      trank->psend = there;
      sc_array_init (&there->message, sizeof(pdata_t));
      cps = (comm_psend_t *) sc_mempool_alloc (gdata->psmem);
      cps->rank = -1;
    }
  
    pad = (pdata_t *) sc_array_push (&there->message);
    memcpy (pad, sc_array_index (particle_data_copy, zz), sizeof (pdata_t));
  
    ++lsend;
  
  } 
  
  sc_mempool_free (gdata->psmem, cps);
  sc_array_sort (gdata->recevs, comm_prank_compare);

  sc_array_destroy_null (&gdata->pfound);

}

void PelletSolver::communicateParticles(){

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
   
  num_receivers = gdata->recevs->elem_count;
 
  
  notif = sc_array_new_count (sizeof (int), num_receivers);
  payl = sc_array_new_count (sizeof (int), num_receivers);   //payload


  for (i = 0; i < num_receivers; ++i) {

    trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);

    *(int *) sc_array_index_int (notif, i) = trank->rank;

    cps = trank->psend;
  
    assert(cps->rank == trank->rank);

    arr = &cps->message;
 
    *(int *) sc_array_index_int (payl, i) = (int) arr->elem_count;
 
  
  }


  sc_notify_ext (notif, NULL, payl, NULL, gdata->mpicomm);

  assert (payl->elem_count == notif->elem_count);

  num_senders = (int) notif->elem_count;

  cucount = 0;
  for (i = 0; i < num_senders; ++i) {
    cucount += *(int *) sc_array_index_int (payl, i);
  }
  prebuf = sc_array_new_count (sizeof(pdata_t), cucount);

  /* post non-blocking receive */
  gdata->recv_req = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);

  cucount = 0;

  for (i = 0; i < num_senders; ++i) {
    count = *(int *) sc_array_index_int (payl, i);
    msglen = count * (int) sizeof(pdata_t);
    mpiret = sc_MPI_Irecv
      (sc_array_index (prebuf, cucount), msglen, sc_MPI_BYTE,
       *(int *) sc_array_index_int (notif, i), 999, gdata->mpicomm,
       (sc_MPI_Request *) sc_array_index_int (gdata->recv_req, i));
    SC_CHECK_MPI (mpiret);
    cucount += count;
  }
  
  assert(cucount == (int) prebuf->elem_count);

  sc_array_destroy_null (&notif);
  sc_array_destroy_null (&payl);

    gdata->send_req = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
    for (i = 0; i < num_receivers; ++i) {
      trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);
      cps = trank->psend;
      arr = &cps->message;
      msglen = (int) (arr->elem_count * arr->elem_size);
      mpiret = sc_MPI_Isend
        (arr->array, msglen, sc_MPI_BYTE, cps->rank, 999,
         gdata->mpicomm, (sc_MPI_Request *) sc_array_index_int (gdata->send_req, i));
      SC_CHECK_MPI (mpiret);
    }

  reqs = (sc_MPI_Request *) sc_array_index (gdata->recv_req,0);
  
  mpiret = sc_MPI_Waitall (num_senders, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&gdata->recv_req);

  num_receivers = (int) gdata->recevs->elem_count;
  reqs = (sc_MPI_Request *) sc_array_index (gdata->send_req,0),
  mpiret = sc_MPI_Waitall (num_receivers, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&gdata->send_req);

  for (i = 0; i < num_receivers; ++i) {
    trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);
    cps = trank->psend;
    sc_array_reset (&cps->message);
  }
  
  sc_array_destroy_null (&gdata->recevs);
  sc_hash_destroy (gdata->psend);

  gdata->psend = NULL;
  sc_mempool_destroy (gdata->psmem);
  gdata->psmem = NULL;
}

static int
slocal_quad2d (p4est_t * p4est, p4est_topidx_t which_tree,
             p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
             void *point)
{
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data      *g = p->gdata;

  /* compute coordinate range of this quadrant */
  loopquad2d (p, which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

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
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data      *g = p->gdata;
  quadrant_data_t          *qud;
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
    zp = sc_array_position (p->prebuf, point);
    cf = (char *) sc_array_index (g->cfound, zp);
    if (!*cf) {
      /* make sure this particle is not found twice */
      *cf = 1;

      /* count this particle in its target quadrant */
      *(p4est_locidx_t *) sc_array_push (g->ireceive) = (p4est_locidx_t) zp;
      qud = (quadrant_data_t *) quadrant->p.user_data;
      ++qud->preceive;
    }

    /* return value will have no effect */
    return 0;
  }

  /* the leaf for this particle has not yet been found */
  return 1;
}

void PelletSolver::postsearch2d(){
  gdata->ireceive = sc_array_new (sizeof (p4est_locidx_t));
  gdata->cfound = sc_array_new_count (sizeof (char), prebuf->elem_count);
    
  sc_array_memset (gdata->cfound, 0);

  p4est_search_local (p4est_heating, 0, slocal_quad2d, slocal_point2d, prebuf);
  
  sc_array_destroy_null (&gdata->cfound);

}

void PelletSolver::regroupParticles2d(){

  sc_array_t         *newpa;
  p4est_topidx_t      tt;
  p4est_locidx_t      newnum;
  p4est_locidx_t      ppos;
  p4est_locidx_t      lq, prev;
  p4est_locidx_t      qboth, li;
  p4est_locidx_t     *premain, *preceive;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  quadrant_data_t          *qud;
  pdata_t          *pad;

  newnum =
    (p4est_locidx_t) (gdata->iremain->elem_count + gdata->ireceive->elem_count);

  premain = (p4est_locidx_t *) sc_array_index (gdata->iremain,0);
  preceive = (p4est_locidx_t *) sc_array_index (gdata->ireceive,0);
  newpa = sc_array_new_count (sizeof (pdata_t), newnum);
  pad = (pdata_t *) sc_array_index (newpa,0);
  prev = 0;
  
  for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
    
      tree = p4est_tree_array_index (p4est_heating->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (quadrant_data_t *) quad->p.user_data;
      qboth = qud->premain + qud->preceive;
      
      if (qboth == 0) {
        qud->lpend = prev;
        qud->premain = qud->preceive = 0;
        continue;
      }
      prev += qboth;
      for (li = 0; li < qud->premain; ++li) {
        ppos = *premain++;
        memcpy (pad, sc_array_index (particle_data_copy, ppos), sizeof (pdata_t));
        pad ++;
      }
      for (li = 0; li < qud->preceive; ++li) {
        ppos = *preceive++;
        memcpy (pad, sc_array_index (prebuf, ppos), sizeof (pdata_t));
        pad ++;
      }
      qud->lpend = prev;
      qud->premain = qud->preceive = 0;
    }
    
  } 
   
  sc_array_destroy_null (&gdata->iremain);

  sc_array_destroy_null (&prebuf);
  sc_array_destroy_null (&gdata->ireceive);
  sc_array_destroy (particle_data_copy);

  particle_data_copy = newpa;
}

void PelletSolver::adaptQuadtree(){
    
    int oldquad = (int)p4est_heating->global_num_quadrants;
    while(true){
    
        gdata->ireindex = gdata->ire2 = 0;
        gdata->irvindex = gdata->irv2 = 0;
        p4est_refine_ext (p4est_heating, 0, 50, adapt_refine2d, NULL, adapt_replace2d);
        if((int)p4est_heating->global_num_quadrants == oldquad)
            break;
        else
            oldquad = p4est_heating->global_num_quadrants;
    }
    
    }

static int
part_weight2d (p4est_t * p4est,
             p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  p4est_locidx_t      ilem_particles;
  PelletSolver      *p = (PelletSolver *) p4est->user_pointer;
  Global_Data      *g = p->gdata;
  quadrant_data_t          *qud = (quadrant_data_t *) quadrant->p.user_data;


  ilem_particles = qud->lpend - g->prevlp;

  g->prevlp = qud->lpend;
  *(int *) sc_array_index (g->src_fixed, g->qcount++) =
    (int) (ilem_particles * sizeof (pdata_t));
  return 1+ ilem_particles;
}
void PelletSolver:: partitionParticles2d(){


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
  quadrant_data_t          *qud;
   
  if(gdata->mpisize == 1)
  
      return;

  src_gfq = P4EST_ALLOC (p4est_gloidx_t, gdata->mpisize + 1);

  memcpy (src_gfq, p4est_heating->global_first_quadrant,
          (gdata->mpisize + 1) * sizeof (p4est_gloidx_t));


  src_quads = p4est_heating->local_num_quadrants;

  assert(src_quads == src_gfq[gdata->mpirank+1]-src_gfq[gdata->mpirank]);

  gdata->src_fixed = sc_array_new_count (sizeof (int), src_quads);

  gdata->qcount = 0;
  gdata->prevlp = 0;

  gshipped = p4est_partition_ext (p4est_heating, 1, part_weight2d);
  dest_quads = p4est_heating->local_num_quadrants;

  if (gshipped == 0) {
    sc_array_destroy_null (&gdata->src_fixed);
    P4EST_FREE (src_gfq);
    return;
  }

  gdata->dest_fixed = sc_array_new_count (sizeof (int), dest_quads);

  p4est_transfer_fixed (p4est_heating->global_first_quadrant, src_gfq,
                        gdata->mpicomm, COMM_TAG_FIXED,
                        (int *) gdata->dest_fixed->array,
                        (const int *) gdata->src_fixed->array, sizeof (int));

  ldatasize = (p4est_locidx_t) sizeof (pdata_t);

  dest_parts = 0;

  for (lq = 0; lq < dest_quads; ++lq) {
    dest_parts += *(int *) sc_array_index (gdata->dest_fixed, lq);
  }
  assert(dest_parts % ldatasize == 0); 
  dest_parts /= ldatasize;
  dest_data = sc_array_new_count (sizeof (pdata_t), dest_parts);
  p4est_transfer_custom (p4est_heating->global_first_quadrant, src_gfq,
                         gdata->mpicomm, COMM_TAG_CUSTOM,
                         (pdata_t *) dest_data->array,
                         (const int *) gdata->dest_fixed->array,
                         (const pdata_t *) particle_data_copy->array,
                         (const int *) gdata->src_fixed->array);

  sc_array_destroy_null (&gdata->src_fixed);

  P4EST_FREE (src_gfq);
  sc_array_destroy (particle_data_copy);
  particle_data_copy = dest_data;
  lpnum = 0;
  lquad = 0;
  for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est_heating->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      /* access quadrant */
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (quadrant_data_t *) quad->p.user_data;

      /* back out particle count in quadrant from data size */
      lcount = *(int *) sc_array_index (gdata->dest_fixed, lquad);
      assert (lcount % ldatasize == 0);
      lcount /= ldatasize;
      lpnum += lcount;
      qud->lpend = lpnum;
      ++lquad;
    }
  }
  sc_array_destroy_null (&gdata->dest_fixed);


}
void PelletSolver::destoryQuadtree(){
    
        p4est_destroy (p4est_heating);
      
        p4est_heating = NULL;
      
        p4est_connectivity_destroy (conn);
      
        conn = NULL;
        
        sc_array_destroy_null(&prebuf_integral);  
        sc_array_destroy_null(&particle_data_copy); 
        for(int i=0;i<2;i++){ 
            sc_array_destroy_null (&ilh[i]);
            sc_array_destroy_null (&jlh[i]);
        }
    
    }


static int
compareXCoordinate (const void *p1, const void *p2)
{
  int t = 0;
  pdata_t                * i1 = (pdata_t *) p1;
  pdata_t                * i2 = (pdata_t *) p2;
  
  if((i1->xyz[2] - i2->xyz[2])>0)
      t = 1;
  return t;
}
void PelletSolver::computeDensityIntegral(){

    p4est_topidx_t      tt;
    p4est_locidx_t      lq, lpend, offset = 0;
    p4est_tree_t       *tree;
    p4est_quadrant_t   *quad;
    quadrant_data_t          *qud;
    pdata_t          *pad;
    double dy, dz;
    double lxyz[3],hxyz[3],dxyz[3];
    double x_min, x_max;
    int i_node, i_nodeend, np_node, node_id;
    double x_lower, x_upper, x_center;
    int division;
    sc_array_t view;  
  
    vector<int> cell_ID;
    vector<double> x_left;
    vector<double> x_right;
    vector<double> x_dividing;
    for (tt = p4est_heating->first_local_tree; tt <= p4est_heating->last_local_tree; ++tt) {
        tree = p4est_tree_array_index (p4est_heating->trees, tt);
        for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
            quad = p4est_quadrant_array_index (&tree->quadrants, lq);
            qud = (quadrant_data_t *) quad->p.user_data;
            lpend = qud->lpend;
            //binary tree construction;            
            if(lpend == offset)
                continue;
            sc_array_init_view(&view,particle_data_copy,offset,lpend-offset);
            sc_array_sort(&view,compareXCoordinate); 
            loopquad2d(this,tt,quad,lxyz,hxyz,dxyz);
            dz = dxyz[0];
            dy = dxyz[1];
            pad = (pdata_t*) sc_array_index_int(&view,0);
            x_min = pad->xyz[2];
            pad = (pdata_t*) sc_array_index_int(&view,lpend-offset-1);
            x_max = pad->xyz[2];
            
            cell_ID.push_back(0);

            x_left.push_back(x_min);
            x_right.push_back(x_max);

            i_node = 0;
            i_nodeend = 1;
            np_node = lpend - offset;
            while(i_node < i_nodeend){
                node_id = cell_ID[i_node];
                x_lower = x_left[node_id];
                x_upper = x_right[node_id];
                np_node = countNumberinRange(&view,lpend-offset,x_lower,x_upper);
                if(np_node <= elem_particle_cell){ 
                    i_node ++;}
                else{
                    x_center = (x_lower+x_upper)/2;
                    cell_ID.push_back(i_nodeend);
                    cell_ID.push_back(i_nodeend+1);
                    x_left.push_back(x_lower);
                    x_left.push_back(x_center);
                    x_right.push_back(x_center);
                    x_right.push_back(x_upper);
                    x_dividing.push_back(x_center);
                    i_node++;
                    i_nodeend += 2;
                    
                    }
                
                
                
                }
            
            x_dividing.push_back(x_max);
            sort(x_dividing.begin(),x_dividing.end());
            division = x_dividing.size();
            vector<double> left_sum(division,0);
            vector<double> right_sum(division,0);
            vector<double> local_sum(division,0) ;
            
            for(int id_particle=0,id_cell=0;id_particle<lpend-offset;)
                {   
                    pad = (pdata_t *) sc_array_index_int(&view,id_particle);
                    if(pad->xyz[2] <= x_dividing[id_cell]){
                            
                        local_sum[id_cell] += pad->mass;
                        if(pad->mass == 0)
                            assert(false);
                        ++id_particle;
                        
                    }
                    else
                    
                        {    
                             ++id_cell;
                          }            
                }
               
             for(int k = 0; k<division; k++){
                for(int p = 0; p<k+1; p++){
                    left_sum[k] += local_sum[p];
                    
                    }
                    left_sum[k] -= local_sum[k]/2;
                    left_sum[k] = left_sum[k]/dy/dz; 
              
                for(int p = k; p<division; p++){
                    right_sum[k] += local_sum[p];
                    }
                    right_sum[k] -= local_sum[k]/2;
                    right_sum[k] = right_sum[k]/dy/dz; 
                
                }
                
                
            for(int id_particle = 0, id_cell = 0;id_particle<lpend - offset; ){
                
                  pad = (pdata_t *) sc_array_index_int(&view,id_particle);
                  if(pad->xyz[2] <= x_dividing[id_cell])   
                        {
                         pad->leftintegral = left_sum[id_cell];
                         pad->rightintegral = right_sum[id_cell];
                         id_particle++;
                         }
                   else
                         {  
                         id_cell++;
                        }
            }
                
            cell_ID.clear();
            x_left.clear();
            x_right.clear();
            x_dividing.clear();
            offset = lpend;
                
             }   
        
        
        }
    
    
    
    }

void PelletSolver::swapXAndZCoordinate(){
    pdata_t *pad;    
    size_t li, lpnum = particle_data_copy->elem_count;
    
    for(li=0; li<lpnum; li++){
        pad = (pdata_t *) sc_array_index(particle_data_copy,li);
        swap(pad->xyz[0],pad->xyz[2]);
        }
    
    } 



int PelletSolver::countNumberinRange(sc_array_t *view, int n, double x, double y) 
    {   
        pdata_t* pad;
        int l = 0, h = n - 1; 
        while (l <= h) { 
            int mid = (l + h) / 2;
            pad = (pdata_t* )sc_array_index_int(view,mid);
            if (pad->xyz[2] >= x) 
                h = mid - 1; 
            else
               l = mid + 1; 
                       } 
        int lower = l;
       l = 0; 
       h = n-1;
         while (l <= h)
         { 
            int mid = (l + h) / 2;
            pad = (pdata_t* )sc_array_index_int(view,mid);
            if (pad->xyz[2] <= y) 
                l = mid + 1; 
            else  
                h = mid - 1;       
                }
        int upper = h;
        return upper-lower+1;
    } 

void PelletSolver::packParticles_phase2(){
    
  int                 retval;
  int                pfn;
  size_t              zz, numz;
  void              **hfound;
  p4est_locidx_t      lremain, lsend;
  comm_psend_t       *cps, *there;
  comm_prank_t       *trank;
  pdata_t          *pad;
  integral_t * integral; 
  
  
  gdata->recevs = sc_array_new (sizeof (comm_prank_t));
  
  gdata->psmem = sc_mempool_new (sizeof (comm_psend_t));
  gdata->psend = sc_hash_new (psend_hash, psend_equal, NULL, NULL);
  numz = particle_data_copy->elem_count;  
    
  lremain = lsend = 0;
  cps = (comm_psend_t *) sc_mempool_alloc (gdata->psmem);
  cps->rank = -1;
  for(zz = 0; zz<numz;zz++){
      pad = (pdata_t *) sc_array_index(particle_data_copy, zz);
      pfn = pad->mpirank;
      if(pfn == gdata->mpirank){
          ++lremain;
          continue;
          }
    
       cps->rank = pfn;

       retval = sc_hash_insert_unique (gdata->psend, cps, &hfound);
            
       there = *((comm_psend_t **) hfound);
      
        if (!retval) {
          assert (there->message.elem_size == sizeof(integral_t));
          assert (there->message.elem_count > 0);
            }
         else{
             assert(there == cps);
             trank = (comm_prank_t *) sc_array_push (gdata->recevs);
             trank->rank = there->rank;
             trank->psend = there;
             sc_array_init (&there->message, sizeof(integral_t));
             cps = (comm_psend_t *) sc_mempool_alloc (gdata->psmem);
             cps->rank = -1;

             }

          integral = (integral_t *) sc_array_push(&there->message);
          integral->leftintegral = pad->leftintegral;
          integral->rightintegral = pad->rightintegral;
          integral->id = pad->id;
          
          ++lsend;
      
      } 
      assert(numz == (size_t)(lsend+lremain));
      sc_mempool_free (gdata->psmem, cps);
      sc_array_sort (gdata->recevs, comm_prank_compare);
    }




void PelletSolver::communicateParticles_phase2(){
    
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
    
  num_receivers = gdata->recevs->elem_count;
  notif = sc_array_new_count (sizeof (int), num_receivers);
  payl = sc_array_new_count (sizeof (int), num_receivers);   //payload
    
  for (i = 0; i < num_receivers; ++i) {

    trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);

    *(int *) sc_array_index_int (notif, i) = trank->rank;

    cps = trank->psend;
  
    assert(cps->rank == trank->rank);

    arr = &cps->message;
 
    *(int *) sc_array_index_int (payl, i) = (int) arr->elem_count;
 
  
    }
    

  sc_notify_ext (notif, NULL, payl, NULL, gdata->mpicomm);
    
  assert (payl->elem_count == notif->elem_count);
    
  num_senders = (int) notif->elem_count;

  cucount = 0;
  for (i = 0; i < num_senders; ++i) {
    cucount += *(int *) sc_array_index_int (payl, i);
  }
  prebuf_integral = sc_array_new_count (sizeof(integral_t), cucount);
    
  gdata->recv_req = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);

  cucount = 0;

  for (i = 0; i < num_senders; ++i) {
    count = *(int *) sc_array_index_int (payl, i);
    msglen = count * (int) sizeof(integral_t);
    mpiret = sc_MPI_Irecv
      (sc_array_index (prebuf_integral, cucount), msglen, sc_MPI_BYTE,
       *(int *) sc_array_index_int (notif, i), 888, gdata->mpicomm,
       (sc_MPI_Request *) sc_array_index_int (gdata->recv_req, i));
    SC_CHECK_MPI (mpiret);
    cucount += count;
  }
    
  assert(cucount == (int) prebuf_integral->elem_count);

  sc_array_destroy_null (&notif);
  sc_array_destroy_null (&payl);
    
    gdata->send_req = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
    for (i = 0; i < num_receivers; ++i) {
      trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);
      cps = trank->psend;
      arr = &cps->message;
      msglen = (int) (arr->elem_count * arr->elem_size);
      mpiret = sc_MPI_Isend
        (arr->array, msglen, sc_MPI_BYTE, cps->rank, 888,
         gdata->mpicomm, (sc_MPI_Request *) sc_array_index_int (gdata->send_req, i));
      SC_CHECK_MPI (mpiret);
    }

  reqs = (sc_MPI_Request *) sc_array_index (gdata->recv_req,0);
  
  mpiret = sc_MPI_Waitall (num_senders, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&gdata->recv_req);

  num_receivers = (int) gdata->recevs->elem_count;
  reqs = (sc_MPI_Request *) sc_array_index (gdata->send_req,0),
  mpiret = sc_MPI_Waitall (num_receivers, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&gdata->send_req);
    
  for (i = 0; i < num_receivers; ++i) {
    trank = (comm_prank_t *) sc_array_index_int (gdata->recevs, i);
    cps = trank->psend;
    sc_array_reset (&cps->message);
  }
    
  sc_array_destroy_null (&gdata->recevs);
  sc_hash_destroy (gdata->psend);

  gdata->psend = NULL;
  sc_mempool_destroy (gdata->psmem);
  gdata->psmem = NULL;
    
    }

void PelletSolver::writeIntegralValue(){
    
    size_t li, num;
    pdata_t *pad, *pad_copy;
    int id;
    integral_t *integral;
    num = prebuf_integral->elem_count;
    size_t counter = 0;
    for(li = 0; li<num ; li++){
        integral = (integral_t *) sc_array_index(prebuf_integral, li);
        id = integral->id;
        pad = (pdata_t *) sc_array_index_int(gdata->particle_data,id);
        pad->leftintegral = integral->leftintegral;
        pad->rightintegral = integral->rightintegral;
        counter ++;
        }
   
    num = particle_data_copy->elem_count;
    for(li = 0; li<num ;li++){
        pad_copy = (pdata_t *) sc_array_index(particle_data_copy,li);
        if(pad_copy->mpirank != gdata->mpirank)
            continue;
        id = pad_copy->id;
        pad = (pdata_t *) sc_array_index_int(gdata->particle_data,id);
        pad->leftintegral = pad_copy->leftintegral;
        pad->rightintegral = pad_copy->rightintegral;
        counter ++;
        }
        
    assert(counter == gdata->particle_data->elem_count);
    
    }

void PelletSolver::computeHeatDeposition( double dt){
    
    pdata_t *pad;
    
    size_t li, lnump = gdata->particle_data->elem_count;
    double e = heatK*(2.99792458e7)/100;
    double lnLambda = log(2*teinf/I*sqrt(exp(1)/2));
    double tauleft; 
    double tauright;
    double tauinf;
    double taueff;
    double uleft; 
    double uright;
    double qinf;
    double guleft;
    double guright; 
    double nt;
    
    for(li = 0; li<lnump ; li++){
        pad = (pdata_t *) sc_array_index(gdata->particle_data,li);
        if(pad->ifboundary)
            continue;
        tauleft = pad->leftintegral/mass*Z;
        tauright = pad->rightintegral/mass*Z;
        tauinf = heatK*heatK*teinf*teinf/(8.0*3.1416*e*e*e*e*lnLambda);
		taueff = tauinf*sqrt(2./(1.+Z));//tauinf/(0.625+0.55*sqrt(one_plus_zstar)); //tauinf*sqrt(2./(1.+ZNe));//tauinf/(0.625+0.55*sqrt(one_plus_Zstar));
        uleft = tauleft/taueff;
        uright = tauright/taueff;                               
        qinf=sqrt(2.0/3.1416/masse)*neinf*pow(heatK*teinf,1.5);
        guleft = sqrt(uleft)*Bessel_K1(sqrt(uleft))/4;
        guright = sqrt(uright)*Bessel_K1(sqrt(uright))/4;
        nt=1.0/pad->volume/mass;

        pad->deltaq = qinf*nt/tauinf*(guleft+guright);
        pad->qplusminus = qinf*0.5*(uleft*Bessel_Kn(2,sqrt(uleft))+uright*Bessel_Kn(2,sqrt(uright)));
       // if(pad->qplusminus == 0)
         //   cout<<pad->leftintegral<<" "<<pad->rightintegral<<endl;
        
        }

}
float           Bessel_I0(
	float           x)
{
        float   p1 = 1.0;
        float   p2 = 3.5156229;
        float   p3 = 3.0899424;
        float   p4 = 1.2067492;
        float   p5 = 0.2659732;
        float   p6 = 0.360768e-1;
        float   p7 = 0.45813e-2;
	
        float   q1 = 0.39894228;
        float   q2 = 0.1328592e-1;
        float   q3 = 0.225319e-2;
        float   q4 = -0.157565e-2;
        float   q5 = 0.916281e-2;
        float   q6 = -0.2057706e-1;
        float   q7 = 0.2635537e-1;
        float   q8 = -0.1647633e-1;
        float   q9 = 0.392377e-2;
	
	float   ax, y, value;
	
	if (fabs(x) < 3.75)
	  {
	    y = (x/3.75)*(x/3.75);//sqr
	    value = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
	  }
	else
	  {
	    ax = fabs(x);
	    y = 3.75/ax;

	    value = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
	  }

	return value;
}


/* Bessel_I1 returns the modifies Bessel function I1(x) of positive real x  */

float           Bessel_I1(
	float           x)
{
        float   p1 = 0.5;
        float   p2 = 0.87890594;
        float   p3 = 0.51498869;
        float   p4 = 0.15084934;
        float   p5 = 0.2658733e-1;
        float   p6 = 0.301532e-2;
        float   p7 = 0.32411e-3;
	
        float   q1 = 0.39894228;
        float   q2 = -0.3988024e-1;
        float   q3 = -0.362018e-2;
        float   q4 = 0.163801e-2;
        float   q5 = -0.1031555e-1;
        float   q6 = 0.2282967e-1;
        float   q7 = -0.2895312e-1;
        float   q8 = 0.1787654e-1;
        float   q9 = -0.420059e-2;
	
	float   ax, y, value;
	
	if (fabs(x) < 3.75)
	  {
	    y = (x/3.75)*(x/3.75);//sqr
	    value = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
	  }
	else
	  {
	    ax = fabs(x);
	    y = 3.75/ax;

	    value = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
	    if (x < 0)
	      value *= -1.0;
	  }
	return value;
}

/* Bessel_K0 returns the modifies Bessel function K0(x) of positive real x  */

float           Bessel_K0(
	float           x)
{
        float   p1 = -0.57721566;
	float   p2 = 0.4227842;
	float   p3 = 0.23069756;
	float   p4 = 0.348859e-1;
	float   p5 = 0.262698e-2;
	float   p6 = 0.1075e-3;
	float   p7 = 0.74e-5;

	float   q1 = 1.25331414;
	float   q2 = -0.7832358e-1;
	float   q3 = 0.2189568e-1;
	float   q4 = -0.1062446e-1;
	float   q5 = 0.587872e-2;
	float   q6 = -0.25154e-2;
	float   q7 = 0.53208e-3;

	float   y, value;

	if (x <= 2.0)
	  {
	    y = x*x/4.0;
	    value = (-log(x/2.0)*Bessel_I0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
	  }
	else
	  {
	    y = 2.0/x;
	    value = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
	  }
	return value;
}

/* Bessel_K1 returns the modifies Bessel function K1(x) of positive real x  */

float           Bessel_K1(
	float           x)
{
        float   p1 = 1.0;
	float   p2 = 0.15443144;
	float   p3 = -0.67278579;
	float   p4 = -0.18156897;
	float   p5 = -0.01919402;
	float   p6 = -0.110404e-2;
	float   p7 = -0.4686e-4;

	float   q1 = 1.25331414;
	float   q2 = 0.23498619;
	float   q3 = -0.3655620e-1;
	float   q4 = 0.1504268e-1;
	float   q5 = -0.780353e-2;
	float   q6 = 0.325614e-2;
	float   q7 = -0.68245e-3;

	float   y, value;

	if (x <= 2.0)
	  {
	    y = x*x/4.0;
	    value = (log(x/2.0)*Bessel_I1(x))+(1.0/x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
	  }
	else
	  {
	    y = 2.0/x;
	    value = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
	  }
	return value;
}
				       		       

/* Bessel_Kn returns the modifies Bessel function Kn(x) of positive real x for n >= 2 */

float           Bessel_Kn(
        int             n,
	float           x)
{
        int    j;
        float  bk, bkm, bkp, tox;

	if (n < 2)
	  {
	    printf("Error in Bessel_Kn(), the order n < 2: n = %d\n",n);
	    assert(false);
	    return 0;
	  }

	tox = 2.0/x;
	bkm = Bessel_K0(x);
	bk = Bessel_K1(x);

	for (j = 1; j < n; j++)
	  {
	    bkp = bkm + j*tox*bk;
	    bkm = bk;
	    bk = bkp;
	  }

	return bk;
}
//Pellet heat deposition calculation






















