#ifndef ICEMC_OPTICAL_PATH_H
#define ICEMC_OPTICAL_PATH_H

#include "Geoid.h"
#include "TVector3.h"
#include <vector>

namespace icemc {

  /**
   * @class OpticalPath
   * @brief A class to hold the output of the RayTracing.
   * 
   * Stores the output of the RayTracing.
   * Determines the attenuation applied to the Askaryan signal.
   */

  class OpticalPath {
  public:

    double distance() const;
    double attenuation() const;
    double magnification() const;
    // propagate the polarization through all the boundaries and return the result
    TVector3 polarization(const TVector3& initialPolarization) const;

    // Equivalent to icemc's mag1 and mag2, for debugging only
    double magnification1() const;
    double magnification2() const;
    
    ///@todo make these private to hide implementation details
    class Step {
    public:
      double distance() const {return direction().Mag();}
      double attenuation() const {return exp(-distance()/attenuationLength);}
      TVector3 direction() const {return end - start;}

      Geoid::Position start;
      Geoid::Position end;
      double n; /// refractive index
      double attenuationLength; /// attenuation length, meters
      TVector3 boundaryNormal; /// normal of the boundary
    };

    void clear(){
      steps.clear();
    }

    void reset(){
      clear();
      residual = DBL_MAX;
    }

    double residual = DBL_MAX;
    std::vector<Step> steps;
  };


  
  
};

#endif
