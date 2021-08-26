#include "Summary.h"
#include "Report.h"

ClassImp(icemc::Summary)

void icemc::FlavorSummary::addEvent(const EventSummary& event, bool passed){

  nTotal++;

  if(passed){
    nPass++;
    //std::cout << "survival weight=" << event.neutrino.weight() << "  phase weight=" << event.loop.phaseWeight() << std::endl;
    nWeighted += event.neutrino.weight()/event.loop.phaseWeight(); // survival weight / (direction weight * position weight)
    length += event.interaction.length/(1e3); // divide by density of ice 1000 kg/m^3 to get length in meters
  }
  
}

void icemc::FlavorSummary::summarize(double iceVolume){

  if(nPass>0){
    calculateEffectiveVolume(iceVolume);
    length = length/nPass; // average now that we're done simulating
    calculateEffectiveArea();
  }
  reportSummary();
}

void icemc::FlavorSummary::reportSummary(){
  
  icemc::report() << "\n~~~~~ Summary for ";  
  switch(flavor){
  case Flavor::e:
    icemc::report() << "electron"; break;
  case Flavor::mu:
    icemc::report() << "mu"; break;
  case Flavor::tau:
    icemc::report() << "tau"; break;
  case Flavor::all:
    icemc::report() << "all";
  }
    
  icemc::report() << " neutrinos  ~~~~~~~~~~~~~~~ " << std::endl;
  icemc::report() << "\tNumber simulated: " << nTotal << std::endl;
  icemc::report() << "\tNumber passed (unweighted): " << nPass << std::endl;
  icemc::report() << "\tNumber passed (weighted): " << nWeighted << std::endl;
  icemc::report() << "\tEffective volume: " << effectiveVolume << " km^3 sr" << std::endl;
  icemc::report() << "\tEffective area: " << effectiveArea << " km^2 sr" << std::endl;
  icemc::report() << "\n";
}

void icemc::FlavorSummary::calculateEffectiveVolume(double iceVolume){
  effectiveVolume = (nWeighted * (iceVolume/1e9) * 4 * TMath::Pi()) / nTotal; // Eq. 8.1 in Cremonesi 2019, convert iceVolume in m^3 to km^3
}

void icemc::FlavorSummary::calculateEffectiveArea(){
  effectiveArea = effectiveVolume / length;
}

void icemc::Summary::addEvent(const EventSummary& event, bool passed){
  // Add relevant event properties into summary
  switch(event.neutrino.flavor){
  case icemc::Neutrino::Flavor::e:
    eSummary.addEvent(event, passed); break;
  case icemc::Neutrino::Flavor::mu:
    muSummary.addEvent(event, passed); break;
  case icemc::Neutrino::Flavor::tau:
    tauSummary.addEvent(event, passed); break;
  default:
    icemc::report() << severity::error << "Attempting to add unknown neutrino flavor to summary!" << std::endl;
  }
}

void icemc::Summary::summarize(){

  icemc::report() << "Volume of antarctic ice is " << (iceVolume/1e9) << " km^3" << std::endl;

  eSummary.summarize(iceVolume);
  muSummary.summarize(iceVolume);
  tauSummary.summarize(iceVolume);

  nTotal = eSummary.nTotal + muSummary.nTotal + tauSummary.nTotal;
  nPass = eSummary.nPass + muSummary.nPass + tauSummary.nPass;
  nWeighted = eSummary.nWeighted + muSummary.nWeighted + tauSummary.nWeighted;
  length = eSummary.length + muSummary.length + tauSummary.length;
    
  calculateEffectiveVolume(iceVolume);
  length = length/nPass; // average now that we're done simulating
  calculateEffectiveArea();

  reportSummary();
}

