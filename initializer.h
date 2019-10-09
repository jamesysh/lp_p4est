#ifndef __INITIALIZER_H__
#define __INITIALIZER_H__

#include "eos.h"
#include <p4est.h>
#include <p8est.h>

#include <p4est_bits.h>
#include <p4est_build.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p8est_bits.h>
#include <p8est_build.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>

#include "geometry.h"
#include "state.h"
#include <sc.h>
#include <cassert>
class Initializer {

    public:
    
        Initializer();
        ~Initializer(){};

        
        int initlevel ; //init level of octree
        int maxlevel ;
        int minlevel;
        int elem_particles; //max number of particles per octant
        double cfl_coefficient;
        double endt;
        double initlocalspacing;
        double initperturbation;
        int eoschoice;
        int pelletmaterial;
        int gamma;
        double timesearchingradius;



   int getPrintVelocity(){return m_iPrintVelocity;} 

   int getPrintVelocityU(){return m_iPrintVelocityU;}

   int getPrintVelocityV(){return m_iPrintVelocityV;}

   int getPrintVelocityW(){return m_iPrintVelocityW;}

   int getPrintDensity(){return m_iPrintDensity;}

   int getPrintMass(){return m_iPrintMass;}

   int getPrintPressure(){return m_iPrintPressure;}

   int getPrintSoundSpeed(){return m_iPrintSoundSpeed;}

   int getPrintLocalSpacing(){return m_iPrintLocalSpacing;}

   int getPrintTemperature(){return m_iPrintTemperature;}

   int getPrintAllParticle(){return m_iPrintAllParticle;}
    private:

    int m_iPrintAllParticle = 0;
    int m_iPrintVelocity = 0;
    int m_iPrintVelocityU = 0;
    int m_iPrintVelocityV = 0;
    int m_iPrintVelocityW = 0;
    int m_iPrintDensity = 0;
    int m_iPrintMass = 0;
    int m_iPrintPressure = 0;
    int m_iPrintSoundSpeed = 0;
    int m_iPrintLocalSpacing = 0;
    int m_iPrintTemperature = 0;

	double m_fStartTime;///< simulation start time
	double m_fEndTime; ///< simulation end time	
	double m_fWriteTimeInterval;///< write time interval
	std::size_t m_iWriteStep; ///< write step
	double m_fCFLCoeff;///< CFL coeff
	int m_iDimension;///< dimension
	int m_iFluidObjNum;///< number of fluid objects		
	int m_iBoundaryObjNum;///< number of boundary objects
	bool m_iRandomDirSplitOrder;///< if true then the order of directional splitting is randomly set 1:yes 0:no	

	std::size_t m_iNumRow2ndOrder;///< the smallest number of rows of A to solve 2nd order LPF
	std::size_t m_iNumRow1stOrder;///< the smallest number of rows of A to solve 1st order LPF
	std::size_t m_iNumCol2ndOrder;//TODO///< the number of columns of A when solving 2nd order LPF	
	std::size_t  m_iNumCol1stOrder;//TODO///< the number of columns of A when solving 1st order LPF
	double m_iTimeSearchRadius;///< the radius for neighbour search
	double m_fInvalidPressure;///< if p < invalid pressure => invalid state
	double m_fInvalidDensity;///< volume cannot be negative: if volume < invalid volume => invalid state	

    int m_iLPFOrder;///< the order of Local Polynomial Fitting (LPF)  	
	int m_iEOSChoice;///< choice of eos
	double m_fGamma;///< eos parameter gamma
	double m_fPinf;///< eos parameter pinf (stiffened poly gas) 
	double m_fEinf;///<eos parameter einf (stiffened poly gas) 	
	double m_fInitParticleSpacing;///< the initial particle spacing	
	double m_fGravity;///< gravity
	bool m_iUseLimiter;///< if use limiter or not 1:yes 0:no
	double m_fInitialPerturbation;//<amount of maximal initial perturbation in dx
	int m_iPelletDistribution;
    int m_iQuadtreeResolution;
    int m_iBinarytreeResolution;
    int m_iHeatModel;
    double m_iMagneticField;
    



};

#endif // __INITIALIZER_H__
