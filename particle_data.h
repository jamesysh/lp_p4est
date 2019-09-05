#ifndef __PARTICLE_DATA_H__
#define __PARTICLE_DATA_H__
#include "p4est_to_p8est.h"

#ifndef P4_TO_P8
#include <p4est.h>
#else
#include <p8est.h>
#endif 

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_build.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#else
#include <p4est_to_p8est.h>
#include <p8est_bits.h>
#include <p8est_build.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#endif


#include <sc.h>
#include "initializer.h"
class Global_Data{

    public:
    
        Global_Data(Initializer* init);
        ~Global_Data(); 
        
        sc_MPI_Comm mpicomm;
        int mpisize,mpirank;
        int initlevel;
        int maxlevel;
        int elem_particles; //max number of particles per octant
        double cfl_coefficient;
        double dt;
        double endt;
        double domain_len = 16; 

        double lxyz[3],hxyz[3],dxyz[3]; //boundingbox of octant
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
    double pressure;
    double soundspeed;
    double temperature;

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
