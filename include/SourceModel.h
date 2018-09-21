#ifndef ICEMC_SOURCE_MODEL_H
#define ICEMC_SOURCE_MODEL_H


#include "Geoid.h"
#include "TVector3.h"
#include "RNG.h"

namespace icemc {

  /**
   * @class SourceModel
   * @brief Choses what directions neutrinos can come from
   * 
   * During nicemc development, this is just a diffuse source, one can imagine adding fixed sources later
   */

  class SourceModel : public RNG {
  public:
    
    // given an already chosen interaction, 
    virtual TVector3 pickDirection(const Geoid::Position& interaction, double refractiveIndex, const TVector3* rfToDetector = NULL) = 0;

  };


  class DiffuseFlux : public SourceModel {
  public:
    virtual TVector3 pickDirection(const Geoid::Position& interaction, double refractiveIndex, const TVector3* rfToDetector = NULL) override;
  };

}





#endif // ICEMC_SOURCE_MODEL_H
