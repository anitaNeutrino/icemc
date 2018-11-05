#ifndef ICEMC_INTERACTION_GENERATOR_H
#define ICEMC_INTERACTION_GENERATOR_H

#include "Geoid.h"
#include "RNG.h"
#include "Interaction.h"
#include "Neutrino.h"

namespace icemc {

  class Settings;
  class WorldModel;
  class CrossSectionModel;
  class YGenerator;
  
  /**
   * @class InteractionGenerator
   * @brief Picks random numbers to model neutrino interactions
   */
  class InteractionGenerator : public RNG {

    const Settings* fSettings = nullptr;
    std::shared_ptr<const WorldModel> fWorldModel = nullptr;
    std::shared_ptr<CrossSectionModel> fCrossSectionModel;
    std::shared_ptr<YGenerator> fYGenerator;

  public:    
    InteractionGenerator(const Settings *settings,
			 std::shared_ptr<WorldModel> worldModel,
			 std::shared_ptr<CrossSectionModel> crossSectionModel,
			 std::shared_ptr<YGenerator> yGenerator);

    Interaction generate(const Neutrino& n, const Geoid::Position& detectorPos);
    Geoid::Position pickInteractionPosition(const Geoid::Position& detector);
    Interaction::Current pickCurrent();
  };


}



#endif //ICEMC_INTERACTION_GENERATOR_H
