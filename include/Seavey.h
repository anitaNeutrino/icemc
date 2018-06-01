
#ifndef ICEMC_SEAVEY_H
#define ICEMC_SEAVEY_H

#include "vector.hh"
#include "Detector.h" // for PropagatingSignal

class TCanvas;

namespace icemc {

  class Settings;

  
  /**
   * @class Seavey
   * @brief A lightweight class to store the ANITA antenna information
   */

  class Seavey {
  public:
    Seavey() {;}

    Seavey(const Vector& antennaPositionInPayloadCoordinates)
      : fPosition(antennaPositionInPayloadCoordinates) {;}

    void applyAntennaGain(PropagatingSignal& incomingSignal) const;

    static TCanvas* plotGains();
 
    const icemc::Vector& getPosition() const {return fPosition;}

    ///@todo  make this private when the refactor is complete?
    static void GetEcompHcompEvector(const Vector& n_eplane, const Vector& n_hplane, const Vector& n_pol,
				     double& e_component, double& h_component, double& n_component);

    ///@todo  make this private when the refactor is complete?
    static void GetEcompHcompkvector(const Vector& n_eplane, const Vector& n_hplane, const Vector& n_normal, 
				     const Vector n_exit2bn,
				     double& e_component_kvector, double& h_component_kvector, double& n_component_kvector);

    
    ///@todo make this private when the refactor is complete?
    static void GetHitAngles(double e_component_kvector,double h_component_kvector,double n_component_kvector, 
			     double& hitangle_e, double& hitangle_h); 
    
    
  private:
    icemc::Vector fPosition; ///< Position in payload centered coordinates
    icemc::Vector fEPlaneOrientation;






    
  };
  

}



#endif // ICEMC_SEAVEY_H
