#include "NeutrinoGenerator.h"
#include "Settings.h"
#include "Report.h"
#include "WorldModel.h"
#include "InteractionGenerator.h"

icemc::NeutrinoGenerator::NeutrinoGenerator(const Settings* settings,
					std::shared_ptr<Source::EnergyModel> sourceEnergyModel,
					std::shared_ptr<Source::DirectionModel> sourceDirectionModel,
					std::shared_ptr<WorldModel> worldModel)
  : fSettings(settings),
    fSourceEnergyModel(sourceEnergyModel),
    fSourceDirectionModel(sourceDirectionModel),
    fWorldModel(worldModel)
{

}


icemc::Neutrino::Flavor icemc::NeutrinoGenerator::pickFlavor() {
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

icemc::Neutrino::L icemc::NeutrinoGenerator::pickMatterType() {
  // pick matter or antimatter, 1:1 ratio
  double r = pickUniform(0, 2);
  if (r < 1){
    return Neutrino::L::Matter;
  }
  else {
    return Neutrino::L::AntiMatter;
  }
}


icemc::Neutrino icemc::NeutrinoGenerator::generate(){

  // neutrino properties
  Neutrino n;
  n.flavor = pickFlavor();
  n.energy = fSourceEnergyModel->pickNeutrinoEnergy();
  n.leptonNumber = pickMatterType(); ///@todo check
  
  // interaction properties
  // n.path.interaction = interaction.position;


  // static auto print_geoid = [](const char* words, const Geoid::Position& p){
  // 			      std::cout << words << "\t" << p.Longitude() << "\t" << p.Latitude() << "\t" << p.Altitude() << std::endl;
  // 			    };
  // {
  //   print_geoid("start 0", opticalPath.steps.at(0).start);
  //   print_geoid("end 0", opticalPath.steps.at(0).end);
  //   print_geoid("start 1", opticalPath.steps.at(1).start);
  //   print_geoid("end 1", opticalPath.steps.at(1).end);
  // } 

  // n.path.direction = fSourceDirectionModel->pickNeutrinoDirection(opticalPath);

  // std::pair<Geoid::Position, double> entry_columnDepth = fWorldModel->integratePath(n.interaction.position, -n.path.direction);

  // std::pair<Geoid::Position, double> exit_columnDepth = fWorldModel->integratePath(n.interaction.position, n.path.direction);
  
  // n.path.entry = entry_columnDepth.first;
  // n.path.columnDepth = entry_columnDepth.second;

  // n.path.exit = exit_columnDepth.first;
  // n.path.columnDepthInteractionToExit = exit_columnDepth.second;
  
  // std::cout << "column depth = " << n.path.columnDepth << std::endl;

  return n;
  
}
