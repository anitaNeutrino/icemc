#include "Neutrino.h"
#include "WorldModel.h"

//ClassImp(icemc::Neutrino);
//ClassImp(icemc::Neutrino::Path);

int icemc::Neutrino::pdgCode(Neutrino::Flavor flavor) {
  switch(flavor){
  case Neutrino::Flavor::e:   return 12;
  case Neutrino::Flavor::mu:  return 14;
  case Neutrino::Flavor::tau: return 16;
  }
}


int icemc::Neutrino::pdgCode() const {
  return pdgCode(flavor);
}


std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::Flavor& f){
  switch(f){
  case icemc::Neutrino::Flavor::e:
    return os << "Neutrino:Flavor:e";
  case icemc::Neutrino::Flavor::mu:
    return os << "Neutrino:Flavor:mu";
  case icemc::Neutrino::Flavor::tau:
    return os << "Neutrino:Flavor:tau";
  default:
    return os << "Neutrino:Flavor:Unknown";
  }
}


std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::L& l){
  switch(l){
  case icemc::Neutrino::L::Matter:
    return os << "Neutrino:L:Matter";
  case icemc::Neutrino::L::AntiMatter:
    return os << "Neutrino:L:AntiMatter";
  default:
    return os << "Neutrino:L:Unknown";
  }
}


void icemc::Neutrino::Path::integrate(const Geoid::Position& interactionPosition, const std::shared_ptr<WorldModel> world){

  if(weight >= 0){
    std::cerr << "Already performed integral." << std::endl;
    return;
  }
  if(world==nullptr){
    std::cerr << "Can't integrate without world model" << std::endl;
    return;
  }

  /**
   * Here we need to do two integrations:
   * From the interaction position, to the entrance i.e. opposite to the neutrino direction
   * From the interaction position, to the exit i.e. along the neutrino direction
   */

  // structural bindings would be nice here...
  std::pair<Geoid::Position, double> entry_columnDepth = world->integratePath(interactionPosition, -direction);
  std::pair<Geoid::Position, double> exit_columnDepth  = world->integratePath(interactionPosition,  direction);

  entry = entry_columnDepth.first;
  columnDepth = entry_columnDepth.second;

  exit = exit_columnDepth.first;
  columnDepthInteractionToExit = exit_columnDepth.second;

  ///@todo figure out the weight...
  weight = 0;  
}
