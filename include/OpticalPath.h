#ifndef ICEMC_OPTICAL_PATH_H
#define ICEMC_OPTICAL_PATH_H

#include "TVector3.h"
#include <vector>

namespace icemc {

  /**
   * @class OpticalPath
   * @brief A class to hold the output of the RayTracing.
   * 
   * Stores the output of the RayTracing in a way that can be used to attenuate the Askaryan signal
   */

  class OpticalPath {
  public:

    double distance() const;
    double attenuation() const;
    // propagate the polarization through all the boundaries and return the result
    TVector3 polarization(const TVector3& initialPolarization) const;

    ///@todo make these private to hide implementation details
    class Step {
    public:
      double distance() const {return direction.Mag();}
      double attenuation() const {return exp(-distance()/attenuationLength);}
      
      TVector3 direction; /// direction
      double n; /// refractive index
      double attenuationLength; /// attenuation length, meters
      TVector3 boundaryNormal; /// normal of the boundary
    };

    void clear(){
      steps.clear();
    }

    double residual = 0;
    std::vector<Step> steps;
    
  };


  
  
};

#endif
