#ifndef __PELLET_SOLVER_H__
#define __PELLET_SOLVER_H__
#include "particle_data.h"
#include "initializer.h"
#include <map>
typedef struct quadrant_data
{   
    
    p4est_locidx_t      lpend;
//    p4est_locidx_t ghostneighbourid[300];
  /** counts of local particles remaining on this quadrant and recieved ones */
  p4est_locidx_t      premain, preceive;
  int quadrantid;
}
quadrant_data_t;

typedef struct integral
{
    p4est_locidx_t id;
    double leftintegral;
    double rightintegral;
    
    }
    integral_t;


class PelletSolver{
    
    public:
    
        PelletSolver(Initializer* init,Global_Data *gdata);
        ~PelletSolver();
        
        void heatingModel(double dt);
        void build_quadtree();
        void prerun();
        void resetQuadrantData();
        void presearch2d(); 
        void packParticles();
        void postsearch2d();
        void regroupParticles2d();
        
        void communicateParticles();
  
        void partitionParticles2d();
        void adaptQuadtree();
        void destoryQuadtree();
  
        void swapXAndZCoordinate();
        void computeDensityIntegral();
        void split_by_coord ( sc_array_t * in,
                sc_array_t * out[2], pa_mode_t mode, int component,
                const double lxyz[3], const double dxyz[3]);
        
        int countNumberinRange(sc_array_t *view, int n, double x, double y); 
        void packParticles_phase2();
        void communicateParticles_phase2();

        void writeIntegralValue();
        
        
        void computeHeatDeposition( double dt);
        
        void computeDensityIntegral1D();
        void getOne_Plus_Zstar(double teinf);
        
       

        int heatingmodel; 
        int elem_particle_box;
        int elem_particle_cell;
        sc_array_t *prebuf;
        sc_array_t *prebuf_integral;
        p4est_t *p4est_heating;
        p4est_connectivity_t *conn;
        sc_array_t *particle_data_copy; //used for pellet problem;
        Global_Data *gdata;
    
        sc_array_t *ilh[2],*jlh[2],*klh[2];
        
       


        
        double mu = 20.1797; 
        double mass = 3.351e-23;
        double Z = 10.;
        double I = 135.5;
        double sublimationenergy = 1363;
        double teinf = 2000.0;
        double neinf = 1.204910e13;
        
        double heatK = 1.602e-18;
        double masse = 9.109e-28;
        double one_plus_Zstar;
        double magneticfield;
    
    };






#endif
