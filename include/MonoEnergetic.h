#ifndef ICEMC_MONOENERGETIC_H
#define ICEMC_MONOENERGETIC_H

#include "AbstractSources.h"

namespace icemc {
  namespace Source {

    /**
     * @class MonoEnergetic
     * @brief A simple class to source monoenergetic neutrinos
     */

    class MonoEnergetic : public EnergyModel {
    public:
      MonoEnergetic(Energy energy) : fEnergy(energy) {;}
      MonoEnergetic(double energy, Energy::Unit unit) : fEnergy(energy, unit) {;}      
      virtual Energy pickNeutrinoEnergy() override {return fEnergy;}
    private:
      const Energy fEnergy;
    };
  }
}


#endif
