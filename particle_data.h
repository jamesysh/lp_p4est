#ifndef __PARTICLE_DATA_H__
#define __PARTICLE_DATA_H__

#include "initializer.h"



typedef enum pa_mode
{
  PA_MODE_REMAIN,
  PA_MODE_RECEIVE,
  PA_MODE_LOCATE
}
pa_mode_t;


class Global_Data{

    public:
    
        Global_Data(Initializer* init);
        ~Global_Data(); 

        void * sc_array_index_begin (sc_array_t * arr);
        
        void sc_array_paste (sc_array_t * dest, sc_array_t * src);
        void adjustCoordByDomain( double xyz[3]);
        void initFluidParticles();
        void prerun();
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
        int maxlevel;
        int minlevel; 
        int elem_particles; //max number of particles per octant
        int eoschoice;
        int pelletmaterial;
        double gamma;
        
        double initlocalspacing;
        double initperturbation;
        
        double dt;
        double endt;
        double domain_len = 16; 
        double bb[6]; // bounding box of initial fluid particles  

        double lxyz[3],hxyz[3],dxyz[3]; //boundingbox of octant
        Geometry* geometry;
        State* state;        
        EOS* eos;        
        
        p4est_locidx_t lpnum; //number of particles on local processor
        p4est_gloidx_t gpnum, gplost; //number of particles on all processor, number of particles on all processers which left domain
        p4est_locidx_t qremain, qreceive;
        sc_array_t *particle_data; //local particle data on process
        
        sc_array_t *target_proc; //target process of particle
        
        sc_array_t *iremain; /**< locidx_t Index into padata of stay-local particles */

        sc_array_t *ireceive;/**< Index into particle receive buffer */

        sc_array_t *recevs;   /**< comm_prank_t with one entry per receiver, sorted */

        sc_array_t  *recv_req; /**< sc_MPI_Request for receiving */
         
        sc_array_t *send_req; /**< sc_MPI_Request for sending */

        sc_array_t *prebuf;  /**< pdata_t All received particles */

        sc_array_t *iffound;   /**< char Flag for received particles */

        sc_hash_t  *psend;    /**< comm_psend_t with one entry per receiver */
      
        sc_mempool_t *psmem;    /**< comm_psend_t to use as hash table entries */


  
        p4est_locidx_t      ireindex, ire2;   /**< Running index into iremain */

        p4est_locidx_t      irvindex, irv2;   /**< Running index into ireceive */

        sc_array_t *ilh[2],*jlh[2],*klh[2];
      //mesh data
        p8est_connectivity_t *conn;
      
        p8est_t            *p8est;

      


};



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

    p4est_gloidx_t      id;

} pdata_t;



/** Data type for payload data inside each quadrant */
typedef struct octant_data
{
 

  /** Offset into local array of all particles after this quadrant */
    p4est_locidx_t      lpend;

  /** counts of local particles remaining on this quadrant and recieved ones */
  p4est_locidx_t      premain, preceive;
}
octant_data_t;


#endif
