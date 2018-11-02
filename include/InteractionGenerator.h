#ifndef ICEMC_INTERACTION_GENERATOR_H
#define ICEMC_INTERACTION_GENERATOR_H

#include "Geoid.h"
#include "RNG.h"
#include "Neutrino.h"

namespace icemc {

  class Settings;
  class WorldModel;
  
  /**
   * @class InteractionGenerator
   * @brief Picks random numbers to model neutrino interactions
   */
  class InteractionGenerator : public RNG {

    const Settings* fSettings = nullptr;
    std::shared_ptr<const WorldModel> fWorldModel = nullptr;
    
  public:    
    InteractionGenerator(const Settings *settings, std::shared_ptr<WorldModel> worldModel); //, int whichray); //, Counting *count1);
    Geoid::Position pickInteractionPosition(const Geoid::Position& detector);
    Neutrino::Interaction::Current pickCurrent();
  };


}



#endif //ICEMC_INTERACTION_GENERATOR_H
