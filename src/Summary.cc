#include "Summary.h"
#include "Report.h"

ClassImp(icemc::Summary)

void icemc::FlavorSummary::addEvent(const Event& event){
 
  nPass++;
  std::cout << "Adding event, phaseWeight = " << event.loop.phaseWeight() << "  survival weight = " << event.neutrino.weight() << std::endl;
  nWeighted += event.loop.phaseWeight()*event.neutrino.weight(); // direction*position*survival weights
  length += event.interaction.length;
  
}

void icemc::FlavorSummary::summarize(int n, double iceVolume){
  nTotal = n;

  calculateEffectiveVolume(iceVolume);
  length = length/nPass; // average now that we're done simulating
  calculateEffectiveArea();

  reportSummary();
}

void icemc::FlavorSummary::reportSummary(){
  
  std::cout << "\n~~~~~ Summary for ";  
  switch(flavor){
  case Flavor::e:
    std::cout << "electron"; break;
  case Flavor::mu:
    std::cout << "mu"; break;
  case Flavor::tau:
    std::cout << "tau"; break;
  case Flavor::all:
    std::cout << "all";
  }
    
  std::cout << " neutrinos  ~~~~~~~~~~~~~~~ " << std::endl;
  std::cout << "\tNumber simulated: " << nTotal << std::endl;
  std::cout << "\tNumber passed (unweighted): " << nPass << std::endl;
  std::cout << "\tNumber passed (weighted): " << nWeighted << std::endl;
  std::cout << "\tEffective volume: " << effectiveVolume << " km^3 sr" << std::endl;
  std::cout << "\tEffectiv area: " << effectiveArea << " km^2 sr" << std::endl;
  std::cout << "\n";
}


void icemc::FlavorSummary::calculateEffectiveVolume(double iceVolume){
  effectiveVolume = (nPass * iceVolume * 4 * TMath::Pi()) / nTotal; // Eq. 8.1 in Linda's paper
}

void icemc::FlavorSummary::calculateEffectiveArea(){
  effectiveArea = effectiveVolume / length;
}

void icemc::Summary::addEvent(const Event& event){
  // Add relevant event properties into summary
  switch(event.neutrino.flavor){
  case icemc::Neutrino::Flavor::e:
    eSummary.addEvent(event); break;
  case icemc::Neutrino::Flavor::mu:
    muSummary.addEvent(event); break;
  case icemc::Neutrino::Flavor::tau:
    tauSummary.addEvent(event); break;
  default:
    icemc::report() << severity::error << "Attempting to add unknown neutrino flavor to summary!" << std::endl;
  }
}


void icemc::Summary::summarize(int nE, int nMu, int nTau){

  eSummary.summarize(nE, iceVolume);
  muSummary.summarize(nMu, iceVolume);
  tauSummary.summarize(nTau, iceVolume);

  nTotal = nE + nMu + nTau;
    
  calculateEffectiveVolume(iceVolume);
  length = length/nPass; // average now that we're done simulating
  calculateEffectiveArea();

  reportSummary();
}

