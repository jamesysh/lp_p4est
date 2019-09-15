#include <iostream>
#include <time.h>
#include <fstream>
#include "sc.h"
#include "mpi.h"
#include "particle_data.h"
#include "geometry_pellet.h"
#include "sc_notify.h"
using namespace std;

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


inline bool is_node_intersect_search_region(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z, const double& search_x, const double& search_y, const double& search_z, const double& radius) {
  //The following calculation has been fashioned such that:
  //if final value of squared_dmin == 0 then the point lies inside the node. 
  //if final value of squared_dmin !=0  then the point lies outside the node, AND 
  //             tells the SQUARE of the minimum distance of the point to the points on the node boundary(surface).


  double squared_dmin = 0.0;
  double temp; //Used as a temporary variable to store some intermediate results.
 
  //Process the x cooridinates
  if( search_x < min_x ) { 
    temp = search_x - min_x; 
    squared_dmin += temp*temp;
  } else if( search_x > max_x )	{ 
    temp         = search_x - max_x ; 
    squared_dmin += temp*temp                  ;
  }   


  //Process the Y-coorindtaes
  if( search_y < min_y ) { 
    temp = search_y - min_y ; 
    squared_dmin += temp*temp;
  } else if( search_y > max_y ) { 
    temp = search_y - max_y ; 
    squared_dmin += temp*temp;
  }   

  //Process the Z-coorindtaes
  if( search_z < min_z ) 
  { 
    temp          = search_z - min_z; 
    squared_dmin += temp*temp;
  } else if( search_z > max_z ) { 
    temp          = search_z - max_z; 
    squared_dmin += temp*temp;
  }   

  if (squared_dmin <= radius*radius) return true;
  else return false;

}

static int neighboursearch_point (p8est_t * p4est, p4est_topidx_t which_tree,
               p8est_quadrant_t * quadrant, int pfirst, int plast,
               p4est_locidx_t local_num, void *point){

    pdata_t *pad = (pdata_t*) point;
    double *x = pad->xyz;
  
    Global_Data      *g = (Global_Data *) p4est->user_pointer;
    double timeradius = g->timesearchingradius;
    double ls  = pad->localspacing;
    if(!is_node_intersect_search_region(g->lxyz[0],g->lxyz[1], g->lxyz[2], g->hxyz[0],g->hxyz[1], g->hxyz[2], x[0],x[1], x[2], ls*timeradius) ){
        return 0;
    }
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

    /* return value will have no effect, but we must return */
    return 0;
  }

  /* the process for this particle has not yet been found */
  return 1;
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
  return 1 + ilem_particles;
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
    timesearchingradius = init->timesearchingradius;
    maxlevel = init->maxlevel;
    minlevel = init->minlevel; 
    initlocalspacing = init->initlocalspacing;
    initperturbation = init->initperturbation;
    elem_particles = init->elem_particles;
    geometry = GeometryFactory::instance().createGeometry("pelletlayer"); 
    geometry->getBoundingBox(bb[0],bb[1],bb[2],bb[3],bb[4],bb[5]);
    state = StateFactory::instance().createState("pelletstate");
    boundary = BoundaryFactory::instance().createBoundary("inflowboundary");
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
      if (qboth == 0) {
        qud->lpend = prev;
        qud->premain = qud->preceive = 0;
        qud->poctant = 0;
        continue;
      }
      prev += qboth;
      for (li = 0; li < qud->premain; ++li) {
        ppos = *premain++;
        memcpy (pad++, sc_array_index (particle_data, ppos), sizeof (pdata_t));
      }
      for (li = 0; li < qud->preceive; ++li) {
        ppos = *preceive++;
        memcpy (pad++, sc_array_index (prebuf, ppos), sizeof (pdata_t));
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

void Global_Data::testquad(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;

  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad, *quad2;
  octant_data_t          *qud,*qud2;
  p4est_locidx_t   offset = 0,lpend;
  pdata_t * pad;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      if(qud->poctant && lq>4000){
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

  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
    lpend = qud->lpend;
    for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        pad->flagboundary = (double)qud->flagboundary;
            
      }
       offset = lpend;  
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


void Global_Data::createViewForOctant(){
    p4est_topidx_t      tt;
  
    p4est_locidx_t      lq;
    p4est_locidx_t      offset = 0;
  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  octant_data_t          *qud;
  pdata_t *pads, *padd;
  for (tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
    tree = p8est_tree_array_index (p8est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p8est_quadrant_array_index (&tree->quadrants, lq);
      qud = (octant_data_t *) quad->p.user_data;
      qud->flagboundary = 0;  
      qud->poctant = qud->lpend - offset;
      qud->particle_data_view = sc_array_new_count(sizeof(pdata_t),(size_t)qud->poctant);
      padd = (pdata_t *) sc_array_index_begin(qud->particle_data_view);
      pads = (pdata_t *) sc_array_index(particle_data,(size_t)offset);
      for(size_t i=0;i<(size_t) qud->poctant;i++){
        copyParticle( padd, pads);
        padd++;
        pads++;
      
      } 
          offset = qud->lpend;
    
    }
  }

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
      sc_array_destroy(qud->particle_data_view);
      sc_array_destroy(qud->localneighbourid);
      sc_array_destroy(qud->ghostneighbourid);
    
      lpend = qud->lpend;
      for(int i=offset;i<lpend;i++){
        pad = (pdata_t *)sc_array_index(particle_data,i);
        sc_array_destroy(pad->localneighbour);
        sc_array_destroy(pad->ghostneighbour);
      }
       offset = lpend;  
    
    }
  }
    sc_array_destroy(irecumu);
    sc_array_destroy(irvcumu);


}

void Global_Data::copyParticle(pdata_t *d, pdata_t *s){
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

}

void Global_Data::initParticleNeighbour(){

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
        pad->localneighbour = sc_array_new(sizeof(p4est_locidx_t));
        pad->ghostneighbour = sc_array_new(sizeof(remoteneighbour_t));
      }
       offset = lpend;  
    }
  }


}
