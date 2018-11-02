#include "NeutrinoFactory.h"
#include "Settings.h"
#include "Report.h"
#include "WorldModel.h"
#include "InteractionGenerator.h"

icemc::NeutrinoFactory::NeutrinoFactory(const Settings* settings,
					std::shared_ptr<Source::EnergyModel> sourceEnergyModel,
					std::shared_ptr<Source::DirectionModel> sourceDirectionModel,
					std::shared_ptr<CrossSectionModel> crossSectionModel,
					std::shared_ptr<YGenerator> yGenerator,
					std::shared_ptr<WorldModel> worldModel,
					std::shared_ptr<InteractionGenerator> interactionGenerator)
  : fSettings(settings),
    fSourceEnergyModel(sourceEnergyModel),
    fSourceDirectionModel(sourceDirectionModel),
    fCrossSectionModel(crossSectionModel),
    fYGenerator(yGenerator),
    fWorldModel(worldModel),
    fInteractionGenerator(interactionGenerator)    
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

icemc::Neutrino icemc::NeutrinoFactory::makeNeutrino(const Geoid::Position& interaction, const OpticalPath& opticalPath) {

  // neutrino properties
  Neutrino n;
  n.flavor = pickFlavor();
  n.energy = fSourceEnergyModel->pickNeutrinoEnergy();
  n.leptonNumber = Neutrino::L::Matter; ///@todo check

  // interaction properties
  n.setInteractionPosition(interaction);


  // static auto print_geoid = [](const char* words, const Geoid::Position& p){
  // 			      std::cout << words << "\t" << p.Longitude() << "\t" << p.Latitude() << "\t" << p.Altitude() << std::endl;
  // 			    };
  // {
  //   print_geoid("start 0", opticalPath.steps.at(0).start);
  //   print_geoid("end 0", opticalPath.steps.at(0).end);
  //   print_geoid("start 1", opticalPath.steps.at(1).start);
  //   print_geoid("end 1", opticalPath.steps.at(1).end);
  // }
  
  n.interaction.current = fInteractionGenerator->pickCurrent();
  // n.interaction.crossSection = fConnollyEtAl2011.getSigma(n.energy, n.leptonNumber,  n.interaction.current); 
  n.interaction.crossSection = fCrossSectionModel->getSigma(n.energy, n.leptonNumber,  n.interaction.current); 
  n.interaction.length = CrossSectionModel::getInteractionLength(n.interaction.crossSection);  
  // Energy pnu, Neutrino::L leptonNumber, Neutrino::Interaction::Current current
  // n.interaction.y = fConnollyEtAl2011.pickY(n.energy, n.leptonNumber, n.interaction.current);
  n.interaction.y = fYGenerator->pickY(n.energy, n.leptonNumber, n.interaction.current);



  n.path.direction = fSourceDirectionModel->pickNeutrinoDirection(opticalPath);
  std::pair<Geoid::Position, double> entry_columnDepth = fWorldModel->integratePath(n.interaction.position, -n.path.direction);

  std::pair<Geoid::Position, double> exit_columnDepth = fWorldModel->integratePath(n.interaction.position, n.path.direction);
  
  n.path.entry = entry_columnDepth.first;
  n.path.columnDepth = entry_columnDepth.second;

  n.path.exit = exit_columnDepth.first;
  n.path.columnDepthInteractionToExit = exit_columnDepth.second;
  
  // std::cout << "column depth = " << n.path.columnDepth << std::endl;

  return n;
  
}
