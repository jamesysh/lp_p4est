#ifndef __BOUNDARY_GRESHO_H__
#define __BOUNDARY_GRESHO_H__

#include "boundary.h"
#include "eos.h"
class Gresho2DSolidBoundary: public Boundary {
public:
	/// constructor
	Gresho2DSolidBoundary();

	/// destructor
	virtual ~Gresho2DSolidBoundary() {}
	
	/**
	 * \brief Get a boundary particle based on a fluid particle      
	 * \param [in] x  The x-coordinate of fluid particle
	 * \param [in] y  The y-coordinate of fluid particle
	 * \param [in] z  The z-coordinate of fluid particle
	 * \param [out] xb  The x-coordinate of boundary particle
	 * \param [out] yb  The y-coordinate of boundary particle
	 * \param [out] zb  The z-coordinate of boundary particle		
	 */
	virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz, 
	double& xb, double& yb, double& zb, 
	double& pressureb, double& vxb, double& vyb, double& vzb);
		
    virtual void generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx);
		
private:
	double radius;	
	double thickness;
	double bo; 
	
};

class Yee2DSolidBoundary: public Boundary {
public:
        /// constructor
        Yee2DSolidBoundary();

        /// destructor
        virtual ~Yee2DSolidBoundary() {}

        /**
         * \brief Get a boundary particle based on a fluid particle      
         * \param [in] x  The x-coordinate of fluid particle
         * \param [in] y  The y-coordinate of fluid particle
         * \param [in] z  The z-coordinate of fluid particle
         * \param [out] xb  The x-coordinate of boundary particle
         * \param [out] yb  The y-coordinate of boundary particle
         * \param [out] zb  The z-coordinate of boundary particle               
         */
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        double& xb, double& yb, double& zb,
        double& pressureb, double& vxb, double& vyb, double& vzb);

        virtual void generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx);

private:
        double radius;
        double thickness;
        double bo;

};

#endif
