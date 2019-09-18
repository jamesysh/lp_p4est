#ifndef __PARTICLE_DATA_H__
#define __PARTICLE_DATA_H__

#include "initializer.h"
#include "boundary.h"



typedef struct pdata{

    double xyz[3]; //coordinates
    double v[3]; //velocity
    double oldv[3];
    double pressure;
    double soundspeed;
    double temperature;
    double volume;
    double mass;
    double localspacing;
    double flagboundary;
    bool ifhasghostneighbour;
    sc_array_t * ghostneighbour;
    sc_array_t * neighbourparticle;
      
    sc_array_t * neighbourupparticle;
    sc_array_t * neighbourdownparticle;
    sc_array_t * neighbourrightparticle;
    sc_array_t * neighbourleftparticle;
    sc_array_t * neighbourfrontparticle;
    sc_array_t * neighbourbackparticle;
    p4est_gloidx_t      id;

} pdata_t;

typedef struct pdata_copy{

    double xyz[3]; //coordinates
    double v[3]; //velocity
    double pressure;
    double soundspeed;
    double temperature;
    double volume;
    double mass;
    double localspacing;
    double flagboundary;
    

} pdata_copy_t;

/** Data type for payload data inside each quadrant */
typedef struct octant_data
{   
    int flagboundary;  //if true, octant is at cloud boundary;
    sc_array_t *localneighbourid;
    sc_array_t *ghostneighbourid;
   // sc_array_t* particle_data_view; 
    pdata_copy_t localparticle[250];
    int octantid;
    int mpirank;
    p4est_locidx_t    poctant;
  /** Offset into local array of all particles after this quadrant */
    p4est_locidx_t      lpend;
//    p4est_locidx_t ghostneighbourid[300];
  /** counts of local particles remaining on this quadrant and recieved ones */
  p4est_locidx_t      premain, preceive;
}
octant_data_t;

typedef enum pa_mode
{
  PA_MODE_REMAIN,
  PA_MODE_RECEIVE,
  PA_MODE_LOCATE
}
pa_mode_t;


typedef struct neighbour_info{
    
    size_t quadid;
    size_t parid;
    bool ifremote;        // if the particle is in remote processer
    bool ifghost;        // if the particle is a ghost particle       
    double distance;
    double phi;         //right and left
    double theta;          //up and down
    double sigma; //front and back
} neighbour_info_t;       

typedef struct comm_psend
{
  int                 rank;
  sc_array_t          message;     /** Message data to send */
}
comm_psend_t;

/** Array entry for a process that we send messages to */
typedef struct comm_prank
{
  int                 rank;
  comm_psend_t       *psend;        /**< Points to hash table entry */
}
comm_prank_t;

typedef enum comm_tag
{
  COMM_TAG_PART = P4EST_COMM_TAG_LAST,
  COMM_TAG_FIXED,
  COMM_TAG_CUSTOM,
  COMM_TAG_LAST
}
comm_tag_t;

class Global_Data{

    public:
    
        Global_Data(Initializer* init);
        ~Global_Data(); 

        void * sc_array_index_begin (sc_array_t * arr);
        
        void sc_array_paste (sc_array_t * dest, sc_array_t * src);
        void adjustCoordByDomain( double xyz[3]);
        void initFluidParticles();
        void prerun();
        void presearch();
        void packParticles();
        void postsearch();
        void resetOctantData();
        void communicateParticles();
        void regroupParticles();
        void partitionParticles();
        void searchNeighbourOctant(); 
        void searchNeighbourParticle();
        void searchUpwindNeighbourParticle();
        void initParticleNeighbour();
        void copyParticle(pdata_copy_t* d, pdata_t *s);
        void generateGhostParticle();
        void fillArrayWithGhostParticle(sc_array_t * array, pdata_t * pad, int count,int dir);
        void createViewForOctant();
        void cleanForTimeStep();
        void testquad();
        void loopquad (p4est_topidx_t tt, p8est_quadrant_t * quad,double lxyz[3], double hxyz[3], double dxyz[3]);
        

        void split_by_coord ( sc_array_t * in,
                sc_array_t * out[2], pa_mode_t mode, int component,
                const double lxyz[3], const double dxyz[3]);
        void cleanUpArrays();
        void writeVTKFiles();
        void setEOS();

        sc_MPI_Comm mpicomm;
        int mpisize,mpirank;
        int initlevel;
        double timesearchingradius;
        int maxlevel;
        int minlevel; 
        int elem_particles; //max number of particles per octant
        int eoschoice;
        int pelletmaterial;
        double gamma;
        
        double initlocalspacing;
        double initperturbation;
        size_t numrow1st;

        double dt;
        double endt;
        double domain_len = 16; 
        double bb[6]; // bounding box of initial fluid particles  

        double lxyz[3],hxyz[3],dxyz[3]; //boundingbox of octant
        Geometry* geometry;
        State* state;        
        EOS* eos;        
        Boundary *boundary; 
        p4est_locidx_t lpnum; //number of particles on local processor
        p4est_gloidx_t gpnum, gplost; //number of particles on all processor, number of particles on all processers which left domain
        p4est_locidx_t qremain, qreceive;
        int flagrefine, gflagrefine, flagstartrefine;
        octant_data_t *ghost_data; 
        p8est_ghost_t *ghost;
        sc_array_t *particle_data; //local particle data on process
        
        sc_array_t *pfound; //target process of particle
        
        sc_array_t *iremain; /**< locidx_t Index into padata of stay-local particles */

        sc_array_t *ireceive;/**< Index into particle receive buffer */

        sc_array_t *recevs;   /**< comm_prank_t with one entry per receiver, sorted */

        sc_array_t  *recv_req; /**< sc_MPI_Request for receiving */
         
        sc_array_t *send_req; /**< sc_MPI_Request for sending */

        sc_array_t *prebuf;  /**< pdata_t All received particles */

        sc_array_t *cfound;   /**< char Flag for received particles */

        sc_hash_t  *psend;    /**< comm_psend_t with one entry per receiver */
      
        sc_mempool_t *psmem;    /**< comm_psend_t to use as hash table entries */

      p4est_locidx_t      prevlp;   /**< Count local particles in partition callback */
      p4est_locidx_t      qcount;   /**< Count local quadrants in partition callback */
      sc_array_t         *src_fixed;        /**< int Particle counts per quadrant */
      sc_array_t         *dest_fixed;       /**< int Particle counts per quadrant */

        int octantid;  
        p4est_locidx_t      ireindex, ire2, ireindex2;   /**< Running index into iremain */

        p4est_locidx_t      irvindex, irv2, irvindex2;   /**< Running index into ireceive */
        sc_array_t *irecumu;   //cumulative remain particles count
        sc_array_t *irvcumu;
        sc_array_t *ilh[2],*jlh[2],*klh[2];
      //mesh data
        p8est_connectivity_t *conn;
      
        p8est_t            *p8est;

      


};






#endif
