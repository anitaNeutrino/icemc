#include "NeutrinoFactory.h"
#include "Settings.h"

icemc::NeutrinoFactory::NeutrinoFactory(const Settings* s)
  : fSettings(s),
    fSpectra(s), ///@todo check cast to int is ok here!
    fConnollyEtAl2011(s),
    fInteraction(&fConnollyEtAl2011, fSettings)
{

}

icemc::Neutrino icemc::NeutrinoFactory::makeNeutrino() {

    // ///@todo finalize
    // int someKindOfError = primary1->GetSigma(pnu, sigma, interactionLength_kgm2, &settings1, xsecParam_nutype, xsecParam_nuint);
    // /// now generate some askaryan rf
  
  Neutrino n;

  n.energy = fSpectra.pickNeutrinoEnergy();
  n.flavor = fConnollyEtAl2011.pickFlavor();
  n.leptonNumber = Neutrino::L::Matter; ///@todo check
  n.interactionCurrent = fInteraction.pickCurrent();
  n.crossSection = fConnollyEtAl2011.getSigma(n.energy, n.leptonNumber,  n.interactionCurrent);
  n.interactionLength = CrossSectionModel::getInteractionLength(n.crossSection);

  // pick flavour!

  return n;
  
}
