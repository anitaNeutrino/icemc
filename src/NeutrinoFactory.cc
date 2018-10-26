#include "NeutrinoFactory.h"
#include "Settings.h"
#include "Report.h"

icemc::NeutrinoFactory::NeutrinoFactory(const Settings* settings,
					std::shared_ptr<Source::EnergyModel> sourceEnergyModel,
					std::shared_ptr<Source::DirectionModel> sourceDirectionModel,
					std::shared_ptr<CrossSectionModel> crossSectionModel,
					std::shared_ptr<YGenerator> yGenerator)
  : fSettings(settings),
    fSourceEnergyModel(sourceEnergyModel),
    fSourceDirectionModel(sourceDirectionModel),
    fCrossSectionModel(crossSectionModel),
    fYGenerator(yGenerator),
    fInteraction(settings)
{

}


icemc::Neutrino::Flavor icemc::NeutrinoFactory::pickFlavor() {
//! pick a neutrino type, flavor ratio 1:1:1
  double r = pickUniform(0, 3);
  if (r < 1){
    return Neutrino::Flavor::e;
  }
  else if(r < 2){
    return Neutrino::Flavor::mu;
  }
  else {
    return Neutrino::Flavor::tau;
  }
}

icemc::Neutrino icemc::NeutrinoFactory::makeNeutrino(const OpticalPath& opticalPath) {

  // neutrino properties
  Neutrino n;
  n.flavor = pickFlavor();
  n.energy = fSourceEnergyModel->pickNeutrinoEnergy();
  n.leptonNumber = Neutrino::L::Matter; ///@todo check

  n.path.direction = fSourceDirectionModel->pickNeutrinoDirection(opticalPath);

  // interaction properties
  n.interaction.position = opticalPath.steps.at(0).start; ///@todo get the *exact* picked position?
  n.interaction.current = fInteraction.pickCurrent();
  // n.interaction.crossSection = fConnollyEtAl2011.getSigma(n.energy, n.leptonNumber,  n.interaction.current); 
  n.interaction.crossSection = fCrossSectionModel->getSigma(n.energy, n.leptonNumber,  n.interaction.current); 
  n.interaction.length = CrossSectionModel::getInteractionLength(n.interaction.crossSection);
  // Energy pnu, Neutrino::L leptonNumber, Neutrino::Interaction::Current current
  // n.interaction.y = fConnollyEtAl2011.pickY(n.energy, n.leptonNumber, n.interaction.current);
  n.interaction.y = fYGenerator->pickY(n.energy, n.leptonNumber, n.interaction.current);  

  return n;
  
}
