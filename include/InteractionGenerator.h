#ifndef ICEMC_INTERACTION_GENERATOR_H
#define ICEMC_INTERACTION_GENERATOR_H

#include "Neutrino.h"
#include "Detector.h"
#include "Interaction.h"

namespace icemc {

  class InteractionGenerator {
    
  public:
    ///@todo remove the neutrino if possible... does this abstraction even make sense?
    virtual Interaction generateInteraction(const Neutrino& n, const Geoid::Position& detector) = 0;
  };

}



#endif //ICEMC_NEUTRINO_INTERACTION_GENERATOR_H
