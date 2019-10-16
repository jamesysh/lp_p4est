#include "pellet_solver.h"
#include "sc_notify.h"
#include <iostream>
using namespace std;
PelletSolver::PelletSolver(Global_Data*g){
    gdata = g;
    elem_particle_box = 150;
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
    p4est_locidx_t *irecumu = (p4est_locidx_t *)sc_array_push(g->irecumu);
    *irecumu = g->ireindex;
    p4est_locidx_t *irvcumu = (p4est_locidx_t *)sc_array_push(g->irvcumu);
    *irvcumu = g->irvindex;
    
    return 0;
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

    g->split_by_coord ( p->klh[wz], p->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord (p->jlh[wy], p->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
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
    g->split_by_coord ( p->klh[wz], p->jlh, PA_MODE_RECEIVE, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      g->split_by_coord ( p->jlh[wy], p->ilh, PA_MODE_RECEIVE, 0, lxyz, dxyz);
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
       *(int *) sc_array_index_int (notif, i), COMM_TAG_PART, gdata->mpicomm,
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
        (arr->array, msglen, sc_MPI_BYTE, cps->rank, COMM_TAG_PART,
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
    zp = sc_array_position (g->prebuf, point);
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

void PelletSolver::destoryQuadtree(){
    
        p4est_destroy (p4est_heating);
      
        p4est_heating = NULL;
      
        p4est_connectivity_destroy (conn);
      
        conn = NULL;
        sc_array_destroy_null(&particle_data_copy); 
        for(int i=0;i<2;i++){ 
            sc_array_destroy_null (&ilh[i]);
            sc_array_destroy_null (&jlh[i]);
        }
    
    }
