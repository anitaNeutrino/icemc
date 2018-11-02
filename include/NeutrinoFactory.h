#ifndef ICEMC_NEUTRINO_FACTORY_H
#define ICEMC_NEUTRINO_FACTORY_H

#include "AbstractSources.h"
#include "ConnollyEtAl2011.h"
#include "Neutrino.h"
#include "RNG.h"
#include "OpticalPath.h"

namespace icemc {
  namespace Source{
    class EnergyModel;
  }

  class Settings;
  class WorldModel;
  class InteractionGenerator;

  /**
   * @class NeutrinoFactory
   * @brief Simplify neutrino generation interface to just one place
   */

  class NeutrinoFactory : public RNG {
  public:

    Neutrino::Flavor pickFlavor();
    
    NeutrinoFactory(const Settings* settings,
		    std::shared_ptr<Source::EnergyModel> sourceEnergyModel,
		    std::shared_ptr<Source::DirectionModel> sourceDirectionModel,
		    std::shared_ptr<CrossSectionModel> crossSectionModel,
		    std::shared_ptr<YGenerator> yGenerator,
		    std::shared_ptr<WorldModel> worldModel,
		    std::shared_ptr<InteractionGenerator> interactionGenerator);
    

    Neutrino makeNeutrino(const OpticalPath& opticalPath);

    

  private:

    std::pair<double, double> integrateNeutrinoPath(const Geoid::Position& interaction, const TVector3& neutrinoDirection) const;
    
    const Settings* fSettings;
    std::shared_ptr<Source::EnergyModel> fSourceEnergyModel;
    std::shared_ptr<Source::DirectionModel> fSourceDirectionModel; 
    std::shared_ptr<CrossSectionModel> fCrossSectionModel;
    std::shared_ptr<YGenerator> fYGenerator;
    std::shared_ptr<WorldModel> fWorldModel;
    // ConnollyEtAl2011 fConnollyEtAl2011; // contains an inelasticity distribution thingy
    std::shared_ptr<InteractionGenerator> fInteractionGenerator;

  };


};

#endif //ICEMC_NEUTRINO_FACTORY_H
