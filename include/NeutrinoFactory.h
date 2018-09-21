#ifndef ICEMC_NEUTRINO_FACTORY_H
#define ICEMC_NEUTRINO_FACTORY_H

#include "Spectra.h"
#include "Primaries.h"
#include "Neutrino.h"

namespace icemc {

  class Settings;

  /**
   * @class NeutrinoFactory
   * @brief Simplify neutrino generation interface to just one place
   */

  class NeutrinoFactory {
  public:

    NeutrinoFactory(const Settings* s);

    Neutrino makeNeutrino();

  private:
    const Settings* fSettings;
    Spectra fSpectra;
    Primaries fPrimaries; // contains an inelasticity distribution thingy
    Interaction fInteraction;    

  };


};

#endif //ICEMC_NEUTRINO_FACTORY_H
