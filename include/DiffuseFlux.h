#ifndef ICEMC_DIFFUSE_H
#define ICEMC_DIFFUSE_H

#include "AbstractSources.h"

namespace icemc {
  namespace Source {
 
    /**
     * @class DiffuseFlux
     * @brief Diffuse flux, all directions equally likely
     */
    class DiffuseFlux : public DirectionModel, public RNG {
    public:
      virtual TVector3 pickNeutrinoDirection(const Geoid::Position& interaction, const TVector3& rfToDetector) override;
    };

    
  }  
}


#endif
