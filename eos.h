
#ifndef __EOS_H__
#define __EOS_H__

#include <cassert>
#include <vector>
class EOS {
protected:
    
    double m_fGamma;	
    int m_iEOSChoice; ///< The eos choice: 1=Polytropic gas; 2=Stiffened Polytropic gas; 3=Saha Eos
public:
    
    /// Destructor
	virtual ~EOS() {};

	/// Getter function of the protected data member m_iEOSChoice
	int getEOSChoice() {return m_iEOSChoice;}
    double getGamma() {return m_fGamma;}	
	virtual void getParameters(std::vector<double>& params) = 0;
	virtual double getEnergy(double pressure, double density) = 0;
 	virtual double getTemperature(double pressure, double density) = 0;
      	virtual double getSoundSpeed(double pressure, double density) = 0;
	virtual double getElectricConductivity(double pressure, double density) = 0;
    virtual void diagnosis(double rho0, double rho1, double t0, double t1) = 0;
};





class PolytropicGasEOS : public EOS {
protected:
    
public:
	PolytropicGasEOS(double gamma) {m_iEOSChoice=1; m_fGamma = gamma;}
	
	// Destructor
	virtual ~PolytropicGasEOS() {}
	virtual void getParameters(std::vector<double>& params){params.push_back(m_fGamma);};
	virtual double getEnergy(double pressure, double density);
	virtual double getTemperature(double pressure, double density);
	virtual double getSoundSpeed(double pressure, double density);
	virtual double getElectricConductivity(double pressure, double density);
    virtual void diagnosis(double rho0, double rho1, double t0, double t1){};

};


#endif
