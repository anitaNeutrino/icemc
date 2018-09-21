#ifndef MONOENERGETIC_H
#define MONOENERGETIC_H

#include "AbstractSources.h"

namespace icemc {
  namespace Source {

    /**
     * @class MonoEnergetic
     * @brief A simple class to source monoenergetic neutrinos
     */

    class MonoEnergetic : public EnergyModel {
    public:
      MonoEnergetic(double energy) : fEnergy(energy) {;}				     
      virtual double pickNeutrinoEnergy() override {return fEnergy;}
    private:
      const double fEnergy;
    };
  }
}


#endif
