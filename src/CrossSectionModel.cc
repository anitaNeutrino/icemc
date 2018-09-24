#include "CrossSectionModel.h"
#include "RNG.h"
#include "Report.h"
#include "Settings.h"

#include "TCanvas.h"
#include "TH2D.h"


TCanvas* icemc::CrossSectionModel::plotSigma(Neutrino::L l, Neutrino::Current c, int nSamples) const {

  RNG r;
  TH2D* h = new TH2D("hSigma","Cross-section model;log_{10}(Energy)(GeV);log_{10}(#sigma) (units?)", 600, 7., 12., 600, -40., -30.);

  for(int i=0; i < nSamples; i++){

    double energy_eV = r.pickUniform(fMinEnergy_eV, fMaxEnergy_eV);
    double sigma = getSigma(energy_eV, l, c);
    
    double energy_GeV=energy_eV/1.E9;//Convert eV to GeV.
    double epsilon=log10(energy_GeV);

    h->Fill(epsilon, log10(sigma));
  }
  h->SetBit(kCanDelete);

  auto can = new TCanvas();
  h->Draw("colz");
  
  return can;
}





double icemc::Reno::getSigma(double energy_eV, Neutrino::L leptonNumber, Neutrino::Current current) const {

  if(!validEnergy(energy_eV)){
    icemc::report() << severity::error << "Need a parameterization for this energy region, energy (eV) = "
		    << energy_eV << " but min = " << fMinEnergy_eV << ", max = " << fMaxEnergy_eV << std::endl;
    return -1;
  }

  // fit to cross sections calculated by M.H. Reno using the same method as Gandhi et al, but with the CTEQ6-DIS parton distribution functions instead of the CTEQ4-DIS distribution functions
  double sigma=(2.501E-39)*pow(energy_eV/1.E9,0.3076)*fSettings->SIGMA_FACTOR; // 10^18 eV - 10^21 eV(use this one for ANITA)
  //sigma=(1.2873E-39)*pow(pnu/1.E9,0.33646)*SIGMA_FACTOR; // 10^17 eV - 10^20 eV (use this one for SalSA)

  return sigma;  
}

















