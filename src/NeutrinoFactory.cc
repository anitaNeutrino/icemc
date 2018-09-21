#include "NeutrinoFactory.h"
#include "Settings.h"

icemc::NeutrinoFactory::NeutrinoFactory(const Settings* s)
  : fSettings(s),
    fSpectra((int)s->EXPONENT), ///@todo check cast to int is ok here!
    fPrimaries(),
    fInteraction(&fPrimaries, fSettings)
{

}

icemc::Neutrino icemc::NeutrinoFactory::makeNeutrino() {


    
    // ///@todo finalize
    // int someKindOfError = primary1->GetSigma(pnu, sigma, interactionLength_kgm2, &settings1, xsecParam_nutype, xsecParam_nuint);
    // /// now generate some askaryan rf
  
  Neutrino n;

  ///@todo put this inside Spectra and have a single function call
  if(fSettings->USEDARTBOARD){
    n.energy = fSpectra.GetNuEnergy();
  }
  else {
    n.energy = fSpectra.GetCDFEnergy();
  }

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
