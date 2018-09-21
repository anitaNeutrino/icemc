#include "NeutrinoFactory.h"
#include "Settings.h"

icemc::NeutrinoFactory::NeutrinoFactory(const Settings* s)
  : fSettings(s),
    fSpectra(s), ///@todo check cast to int is ok here!
    fPrimaries(),
    fInteraction(&fPrimaries, fSettings)
{

}

icemc::Neutrino icemc::NeutrinoFactory::makeNeutrino() {

    // ///@todo finalize
    // int someKindOfError = primary1->GetSigma(pnu, sigma, interactionLength_kgm2, &settings1, xsecParam_nutype, xsecParam_nuint);
    // /// now generate some askaryan rf
  
  Neutrino n;

  n.energy = fSpectra.pickNeutrinoEnergy();
  n.flavor = fPrimaries.GetNuFlavor();
  n.leptonNumber = Neutrino::L::Matter; ///@todo check
  n.interactionCurrent = fInteraction.GetCurrent();
  fPrimaries.GetSigma(n.energy, n.crossSection, n.interactionLength,
		      fSettings,
		      n.leptonNumber,  n.interactionCurrent);
		      // xsecParam_nutype, xsecParam_nuint);
  // pick flavour!

  return n;
  
}
