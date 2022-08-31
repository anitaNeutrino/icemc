#ifndef ICEMC_SOURCE_MODEL_H
#define ICEMC_SOURCE_MODEL_H

#include "OpticalPath.h"
#include "Neutrino.h"
#include "LoopInfo.h"

namespace icemc {
  namespace Source {

    /**
     * @class DirectionModel
     * @brief Abstract base class for chosing directions neutrinos can come from.
     */
    class DirectionModel {
    protected:
      double directionWeight;
    public:
      virtual TVector3 pickNeutrinoDirection(const OpticalPath& opticalPath, const LoopInfo& loop) = 0;
      virtual double getDirectionWeight() = 0;
    };


    /**
     * @class EnergyModel
     * @brief Abstract base class for chosing neutrino energies
     */

    class EnergyModel {
    public:
      virtual Energy pickNeutrinoEnergy() = 0;
    };





   

  }  
}





#endif // ICEMC_SOURCE_MODEL_H
