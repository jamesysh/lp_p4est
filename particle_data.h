#ifndef __PARTICLE_DATA_H__
#define __PARTICLE_DATA_H__

#include "initializer.h"
class Global_Data{

    public:
    
        Global_Data(Initializer* init);
        ~Global_Data(); 
        
        void initFluidParticles();
        void cleanUpArrays();
        void writeVTKFiles();

        sc_MPI_Comm mpicomm;
        int mpisize,mpirank;
        int initlevel;
        int maxlevel;
        int elem_particles; //max number of particles per octant
        double initlocalspacing;
        double initperturbation;
        double dt;
        double endt;
        double domain_len = 16; 
        double bb[6]; // bounding box of initial fluid particles  

        double lxyz[3],hxyz[3],dxyz[3]; //boundingbox of octant
        Geometry* geometry;
        State* state;        
        
        
        p4est_locidx_t lpnum; //number of particles on local processor
        p4est_gloidx_t gpnum, gplost; //number of particles on all processor, number of particles on all processers which left domain
        sc_array_t *particle_data; //local particle data on process
        
        sc_array_t *target_proc; //target process of particle
        
        sc_array_t *idremain; /**< locidx_t Index into padata of stay-local particles */

        sc_array_t *idreceive;/**< Index into particle receive buffer */

        sc_array_t *recevs;   /**< comm_prank_t with one entry per receiver, sorted */

        sc_array_t  *recv_req; /**< sc_MPI_Request for receiving */
         
        sc_array_t *send_req; /**< sc_MPI_Request for sending */

        sc_array_t *receivebuffer;  /**< pa_data_t All received particles */

        sc_array_t *iffound;   /**< char Flag for received particles */

        sc_hash_t  *psend;    /**< comm_psend_t with one entry per receiver */
      
        sc_mempool_t *psmem;    /**< comm_psend_t to use as hash table entries */


      //mesh data
        p4est_connectivity_t *conn;
      
        p4est_t            *p4est;

      


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
