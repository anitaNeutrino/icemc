#ifndef ICEMC_NEUTRINO_GENERATOR_H
#define ICEMC_NEUTRINO_GENERATOR_H

#include "AbstractSources.h"
#include "ConnollyEtAl2011.h"
#include "Neutrino.h"
#include "RNG.h"
// #include "OpticalPath.h"
// #include "Interaction.h"


namespace icemc {
  namespace Source{
    class EnergyModel;
  }

  class Settings;
  class WorldModel;

  /**
   * @class NeutrinoGenerator
   * @brief Simplify neutrino generation interface to just one place
   */

  class NeutrinoGenerator : public RNG {
  public:

    Neutrino::Flavor pickFlavor();
    Neutrino::L pickMatterType();
    
    NeutrinoGenerator(const Settings* settings,
		    std::shared_ptr<Source::EnergyModel> sourceEnergyModel,
		    std::shared_ptr<Source::DirectionModel> sourceDirectionModel,
		    std::shared_ptr<WorldModel> worldModel);

    // Neutrino makeNeutrino(const Interaction& interaction, const OpticalPath& opticalPath);
    Neutrino generate(); //const Interaction& interaction, const OpticalPath& opticalPath);    

    

  private:

    std::pair<double, double> integrateNeutrinoPath(const Geoid::Position& interaction, const TVector3& neutrinoDirection) const;
    
    const Settings* fSettings;
    std::shared_ptr<Source::EnergyModel> fSourceEnergyModel;
    std::shared_ptr<Source::DirectionModel> fSourceDirectionModel; 
    std::shared_ptr<WorldModel> fWorldModel;
  };


};

#endif //ICEMC_NEUTRINO_FACTORY_H
