#ifndef ICEMC_NEUTRINO_INTERACTION_GENERATOR_H
#define ICEMC_NEUTRINO_INTERACTION_GENERATOR_H

#include "Geoid.h"
#include "RNG.h"
#include "Interaction.h"
#include "InteractionGenerator.h"
#include "Neutrino.h"

namespace icemc {

  class Settings;
  class WorldModel;
  class CrossSectionModel;
  class YGenerator;
  
  /**
   * @class NeutrinoInteractionGenerator
   * @brief Picks random numbers to model neutrino interactions
   */
  class NeutrinoInteractionGenerator : public InteractionGenerator, public RNG {

    const Settings* fSettings = nullptr;
    std::shared_ptr<const WorldModel> fWorldModel = nullptr;
    std::shared_ptr<CrossSectionModel> fCrossSectionModel;
    std::shared_ptr<YGenerator> fYGenerator;

    Geoid::Position fSpecificInteractionCenter; ///< Used in the Settings::SPECIFIC_NU_POSITION case

    /** 
     * Pick a position in the ice near center, reject if outside sphere of rangeMeters 
     */
    Geoid::Position pickInteractionPositionInIce(const Geoid::Position& center, double rangeMeters);
    
    /** 
     * Pick a position inside sphere of rangeMeters around center, reject if not inside ice.
     */
    Geoid::Position pickInteractionPositionInSphere(const Geoid::Position& center, double rangeMeters);
    Geoid::Position pickInteractionPosition(const Geoid::Position& detector);
    Interaction::Current pickCurrent();
    bool fAllowInteractionOutsideIce = true; ///@todo testing, remove me
  public:    
    NeutrinoInteractionGenerator(const Settings *settings,
				 std::shared_ptr<WorldModel>        worldModel,
				 std::shared_ptr<CrossSectionModel> crossSectionModel,
				 std::shared_ptr<YGenerator>        yGenerator);

    virtual Interaction generateInteraction(const Neutrino& n, const Geoid::Position& detector) override;


    void setAllowInteractionOutsideIce(bool allow = true){
      fAllowInteractionOutsideIce = allow;
    }
    bool getAllowInteractionOutsideIce() const {return 	fAllowInteractionOutsideIce;}

  };


}



#endif //ICEMC_NEUTRINO_INTERACTION_GENERATOR_H
