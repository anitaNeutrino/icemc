#ifndef ICEMC_DIFFUSE_H
#define ICEMC_DIFFUSE_H

#include "AbstractSources.h"
#include "RNG.h"

namespace icemc {
  namespace Source {
 
    /**
     * @class DiffuseFlux
     * @brief Diffuse flux, all directions equally likely
     */
    class DiffuseFlux : public DirectionModel, public RNG {
    public:
      virtual TVector3 pickNeutrinoDirection(const OpticalPath& opticalPath, double dtheta) override;
      virtual double getDirectionWeight() override;
    };

    
  }  
}


#endif
