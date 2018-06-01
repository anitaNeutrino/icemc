
#ifndef ICEMC_SEAVEY_H
#define ICEMC_SEAVEY_H

#include "vector.hh"
#include "Detector.h" // for PropagatingSignal

class TCanvas;

namespace icemc {

  class Settings;
  class Balloon;

  /**
   * @class Seavey
   * @brief A lightweight class to store the ANITA antenna information
   */

  class Seavey {
  public:

    /**
     * This is the opposite convention from EventReaderRoot and all downstream things! Be warned!
     * 
     */
    enum Pol {
      sVertical = 0,
      sHorizontal = 1
    };

    Seavey(const Vector& pos, const Vector& ePlane, const Vector& hPlane, const Vector& normal) :
      fPosition(pos), fEPlane(ePlane), fHPlane(hPlane), fNormal(normal)
    {;}

    void applyAntennaGain(PropagatingSignal& incomingSignal) const;
 
    const icemc::Vector& getPosition() const {return fPosition;}
    const icemc::Vector& getEPlane() const {return fEPlane;}
    const icemc::Vector& getHPlane() const {return fHPlane;}
    const icemc::Vector& getNormal() const {return fNormal;}
    
    
    ///@todo  make this private when the refactor is complete?
    static void GetEcompHcompEvector(const Vector& n_eplane, const Vector& n_hplane, const Vector& n_pol,
				     double& e_component, double& h_component, double& n_component);

    ///@todo  make this private when the refactor is complete?
    static void GetEcompHcompkvector(const Vector& n_eplane, const Vector& n_hplane, const Vector& n_normal, const Vector n_exit2bn,
				     double& e_component_kvector, double& h_component_kvector, double& n_component_kvector);
    
    ///@todo make this private when the refactor is complete?
    static void GetHitAngles(double e_component_kvector,double h_component_kvector,double n_component_kvector, 
			     double& hitangle_e, double& hitangle_h); 


    static TCanvas* plotGains();
    
  private:
    const icemc::Vector fPosition; ///< Position in payload centered coordinates
    const icemc::Vector fEPlane; ///< Seavey E-plane
    const icemc::Vector fHPlane; ///< Seavey H-plane
    const icemc::Vector fNormal; ///< Normal to the antenna






    
  };
  

}



#endif // ICEMC_SEAVEY_H
