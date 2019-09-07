#include <iostream>
#include <time.h>
#include "initializer.h"
#include <fstream>
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
    int nump;
    Global_Data* g = (Global_Data*) user_data;
     
    Geometry * geom = g->geometry;
    State* state = g->state;
    EOS* eos = g->eos;
    double x,y,z;
    int i,j,k;
    double ls = g->initlocalspacing;
    double initr = g->initperturbation;
    int *lpnum = &g->lpnum; 
    p4est_topidx_t      tt = info->treeid;  /**< the tree containing \a quad */
    p4est_quadrant_t   *quad = info->quad;
    double domain_len = g->domain_len;
    double* bb = g->bb;
    pdata_t *pd;
    octant_data_t *oud = (octant_data_t *)quad->p.user_data;
    p4est_qcoord_t qh;
   


    oud->premain = oud->preceive = 0;
    qh = P4EST_QUADRANT_LEN (quad->level);
    l = qh/(double)P4EST_ROOT_LEN*domain_len;
    p4est_qcoord_to_vertex (g->conn, tt, quad->x, quad->y,
#ifdef P4_TO_P8
                          quad->z,
#endif
                          g->lxyz);
   
    adjustCoordByDomain(domain_len,g->lxyz);
    
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
                pd->xyz[0] = x;
                pd->xyz[1] = y;
#ifdef P4_TO_P8
                pd->xyz[2] = z;
#endif
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
};

Global_Data:: Global_Data(Initializer* init){
    
    lpnum = 0;
    gpnum = 0;
    
    initlevel = init->initlevel;
    maxlevel = init->maxlevel;
    
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
   gplost = 0; 
   srand(time(NULL));   
   P4EST_ASSERT (particle_data == NULL);
   particle_data = sc_array_new(sizeof( pdata_t ));
   p4est_iterate(p4est,NULL,(void *)this,createParticlesInOctant,NULL,NULL,NULL); 
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


