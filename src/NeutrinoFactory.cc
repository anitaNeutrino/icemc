#include "NeutrinoFactory.h"
#include "Settings.h"
#include "Report.h"

icemc::NeutrinoFactory::NeutrinoFactory(const Settings* s)
  : fSettings(s),
    fSpectra(s),
    fConnollyEtAl2011(s),
    fInteraction(s)
{

}


icemc::Neutrino::Flavor icemc::NeutrinoFactory::pickFlavor() {
//! pick a neutrino type, flavor ratio 1:1:1
  double r = pickUniform(0, 3);
  if (r <= 1){
    return Neutrino::Flavor::e;
  }
  else if(r <= 2){
    return Neutrino::Flavor::mu;
  }
  else {
    return Neutrino::Flavor::tau;
  }
}

icemc::Neutrino icemc::NeutrinoFactory::makeNeutrino() {

    // ///@todo finalize
    // int someKindOfError = primary1->GetSigma(pnu, sigma, interactionLength_kgm2, &settings1, xsecParam_nutype, xsecParam_nuint);
    // /// now generate some askaryan rf
  
  Neutrino n;
  n.flavor = pickFlavor();

  n.energy = fSpectra.pickNeutrinoEnergy();
  n.leptonNumber = Neutrino::L::Matter; ///@todo check
  n.interactionCurrent = fInteraction.pickCurrent();
  n.crossSection = fConnollyEtAl2011.getSigma(n.energy, n.leptonNumber,  n.interactionCurrent);
  n.interactionLength = CrossSectionModel::getInteractionLength(n.crossSection);

  return n;
  
}
