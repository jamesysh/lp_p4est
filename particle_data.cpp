#include <iostream>
#include <time.h>
#include "initializer.h"
#include <fstream>
#include "sc.h"
#include "mpi.h"
#include "particle_data.h"
#include "geometry_pellet.h"
using namespace std;


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

    if(*pfn != g->mpirank && pfirst < *pfn)
        printf("case happened\n");
    /* return value will have no effect, but we must return */
    return 0;
  }

  /* the process for this particle has not yet been found */
  return 1;
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
    double ls = g->initlocalspacing;
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
   
    p4est_locidx_t      *remainid;


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
    
    for(i = 0;i<nump;i++){
        for(j=0;j<nump;j++){
            for(k=0;k<nump;k++){
                x = g->lxyz[0] + (i+rand()/(double)RAND_MAX*initr)*ls;
                y = g->lxyz[1] + (j+rand()/(double)RAND_MAX*initr)*ls;
                z = g->lxyz[2] + (k+rand()/(double)RAND_MAX*initr)*ls;
            if(geom->operator()(x,y,z)){
       
                pd = (pdata_t *) sc_array_push_count (g->particle_data,1);
                remainid = (p4est_locidx_t *) sc_array_push_count(g->iremain,1);
                *remainid = *lpnum;
                pd->xyz[0] = x;
                pd->xyz[1] = y;
                pd->xyz[2] = z;
                state->velocity(x,y,z,pd->v[0],pd->v[1],pd->v[2]);                
                pd->volume = 1./state->density(x,y,z);
                pd->pressure = state->pressure(x,y,z);
                pd->localspacing = ls;
                pd->mass = ls*ls*ls/pd->volume/sqrt(2); 
                pd->soundspeed = eos->getSoundSpeed(pd->pressure,1./pd->volume);
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

Global_Data:: Global_Data(Initializer* init){
    
    lpnum = 0;
    gpnum = 0;
    gplost = 0; 
    flagrefine = 1;
    initlevel = init->initlevel;
    maxlevel = init->maxlevel;
    minlevel = init->minlevel; 
    initlocalspacing = init->initlocalspacing;
    initperturbation = init->initperturbation;
    elem_particles = init->elem_particles;
    geometry = GeometryFactory::instance().createGeometry("pelletlayer"); 
    geometry->getBoundingBox(bb[0],bb[1],bb[2],bb[3],bb[4],bb[5]);
    state = StateFactory::instance().createState("pelletstate");
    eoschoice = init->eoschoice;
    setEOS();    


}


Global_Data:: ~Global_Data(){

}


void Global_Data::initFluidParticles(){
   
   int mpiret;
   
   p4est_locidx_t      li;
   pdata_t          *pad;
   p4est_gloidx_t  gpoffset;
   srand(time(NULL));   

   
   
   p8est_iterate(p8est,NULL,(void *)this,createParticlesInOctant,NULL,NULL,NULL); 
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
    sc_array_destroy_null (&klh[i]);
  }
}


void Global_Data:: writeVTKFiles(){
    static bool FIRST = true;
    static int timestep = 0;
    p4est_locidx_t li;
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
         eos = new PolytropicGasEOS(gamma,pelletmaterial);
	     std::vector<double> eos_parameters;
	     eos->getParameters(eos_parameters);
    }

    else {
		cout<<"The choice of EOS does not exist!!! Please correct the input file."<<endl;
		assert(false);
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

   ireceive = sc_array_new(sizeof(p4est_locidx_t));    

   particle_data = sc_array_new(sizeof( pdata_t ));
   iremain = sc_array_new(sizeof(p4est_locidx_t));
   for (int i = 0; i < 2; ++i) {
    ilh[i] = sc_array_new (sizeof (p4est_locidx_t));
    jlh[i] = sc_array_new (sizeof (p4est_locidx_t));
    klh[i] = sc_array_new (sizeof (p4est_locidx_t));
  }
}

void Global_Data::presearch(){
        
  pfound = sc_array_new_count (sizeof (int), particle_data->elem_count);
  
  iremain = sc_array_new (sizeof (p4est_locidx_t));
  sc_array_memset (pfound, -1);

  p8est_search_all (p8est, 0, psearch_quad, psearch_point, particle_data);

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
    ("In current timestep from %lld remain %lld sent %lld lost %lld number of particles, %.3g of processers are receiving particels\n",
      (long long) gpnum, (long long) glolrs[0],
     (long long) glolrs[1], (long long) glolrs[2],
     glolrs[3] / (double)mpisize);
  assert (glolrs[0] + glolrs[1] + glolrs[2] == gpnum);
  gplost += glolrs[2];
  gpnum -= glolrs[2];


  sc_array_destroy_null (&pfound);
}
