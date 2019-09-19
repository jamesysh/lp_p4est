#include <iostream>
#include "eos.h"
#include <cmath>
#include <cstring>
#include <unistd.h>
////////////////////////////////////////////////////////////////////////////////
// Start of PolytropicGasEOS
////////////////////////////////////////////////////////////////////////////////

double PolytropicGasEOS::getEnergy(double pressure, double density) {
	if(((m_fGamma - 1.) * density) != 0) 
		return ( pressure / ((m_fGamma - 1.) * density) );
	else { // divide by zero
		std::cout<<"Error (Divide by zero)! Computing energy by EOS: "<<std::endl; 
		std::cout<<"gamma = "<<m_fGamma<<", density = "<<density<<std::endl;
		assert(false);
	}
}




double PolytropicGasEOS::getSoundSpeed(double pressure, double density) {
	double cs;
	if(density != 0)
		cs = m_fGamma * pressure / density;

	else {
		std::cout<<"Error (Divide by zero density)! Computing sound speed by EOS: "<<std::endl;
		//std::cout<<"density = "<<density<<std::endl;
		assert(false);
	}
	if(cs > 0) 
		return sqrt(cs);
	else if(cs==0)
		return cs;
	else { // taking square root of a negative number
		std::cout<<"Error (Taking suqre root of a negative number)! Computing sound speed by EOS: "<<std::endl; 
		std::cout<<"gamma = "<<m_fGamma<<", pressure = "<<pressure<<", density = "<<density<<std::endl;
		assert(false);
	}
}

double PolytropicGasEOS::getElectricConductivity(double pressure, double density) {

    double cond = 0.;
    return cond;

}

double PolytropicGasEOS::getTemperature(double pressure, double density) {
  double R,mu;
  
  if(m_iPelletMaterial == 0 ){
    mu = 20.18;
    R = 83.14;
  }
  else if(m_iPelletMaterial == 1){
       R = 83.14;
       mu = 2.014;
      }
  return mu*pressure/(R*density)/11604.525;
}

////////////////////////////////////////////////////////////////////////////////
// End of PolytropicGasEOS
////////////////////////////////////////////////////////////////////////////////

