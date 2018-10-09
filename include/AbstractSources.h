#ifndef ICEMC_SOURCE_MODEL_H
#define ICEMC_SOURCE_MODEL_H

#include "Geoid.h"
#include "TVector3.h"
#include "RNG.h"
#include "Neutrino.h"

namespace icemc {
  namespace Source {

    /**
     * @class DirectionModel
     * @brief Abstract base class for chosing directions neutrinos can come from.
     */
    class DirectionModel {
    public:
      virtual TVector3 pickNeutrinoDirection(const Geoid::Position& interaction, const TVector3& rfToDetector) = 0;
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
