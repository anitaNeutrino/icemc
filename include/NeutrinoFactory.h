#ifndef ICEMC_NEUTRINO_FACTORY_H
#define ICEMC_NEUTRINO_FACTORY_H

#include "Spectra.h"
#include "ConnollyEtAl2011.h"
#include "Neutrino.h"
#include "RNG.h"
#include "Interaction.h"

namespace icemc {

  class Settings;

  /**
   * @class NeutrinoFactory
   * @brief Simplify neutrino generation interface to just one place
   */

  class NeutrinoFactory : public RNG {
  public:

    Neutrino::Flavor pickFlavor();
    
    NeutrinoFactory(const Settings* s);

    Neutrino makeNeutrino();

    

  private:
    const Settings* fSettings;
    Source::Spectra fSpectra;
    ConnollyEtAl2011 fConnollyEtAl2011; // contains an inelasticity distribution thingy
    Interaction fInteraction;    

  };


};

#endif //ICEMC_NEUTRINO_FACTORY_H
