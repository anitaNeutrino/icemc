#include "AskaryanFreqs.h"
#include "IcemcLog.h"
#include "anita.hh"
#include "TGraph.h"
#include "secondaries.hh"

ClassImp(icemc::AskaryanFreqs)


icemc::AskaryanFreqs::AskaryanFreqs()
{

}





icemc::AskaryanFreqs::AskaryanFreqs(int nf, double minFreqHz, double deltaFreqHz, double cherenkovAngleRad, const ShowerProperties* sp, const double* vmmhz_input)
  : vmmhz(vmmhz_input, vmmhz_input + nf),  fMinFreqHz(minFreqHz < 0 ? 0 : minFreqHz), fDeltaFreqHz(deltaFreqHz), fCherenkovAngleRad(cherenkovAngleRad) 
{

  fSpreadTestFreqHz = fMinFreqHz <= 0 ? deltaFreqHz : fMinFreqHz;
  fEmFrac = sp->emFrac;
  fHadFrac = sp->hadFrac;  
}




void icemc::AskaryanFreqs::taperAmplitudesForOffConeViewing(double viewAngleRadians, double testFreqHz, double deltaThetaEmTest, double deltaThetaHadTest){  

  constexpr double epsilon = 1e-10;
  const double sin_viewAngle = sin(viewAngleRadians);

  if(testFreqHz <= 0){
    testFreqHz = fSpreadTestFreqHz;
    deltaThetaEmTest = fDeltaThetaEmTest;
    deltaThetaHadTest = fDeltaThetaHadTest;
  }

  double freq = fMinFreqHz;
  for(double& amp : vmmhz){
    // the angular spread (deltaTheta) type variables follow a 1/freq dependence    
    // so it can just be calculated once and scaled.

    if(freq > 0){
      // there should be 0 power in a DC offset bin, so skip
      // (it would cause a division by 0...)

      // here we encode the 1/f  dependence... for the freq in the j-th element
      double deltaThetaEm = deltaThetaEmTest*testFreqHz/freq;
      double deltaThetaHad = deltaThetaHadTest*testFreqHz/freq;

      constexpr double maxSigma = 20; // ignore anything more than 20 sigma away from cone
    
      // V/m/MHz at 1m due to EM component of shower        
      double nSigmaEm = (viewAngleRadians-fCherenkovAngleRad)*(viewAngleRadians-fCherenkovAngleRad)/(deltaThetaEm*deltaThetaEm);
      double vmmhz1m_em = (fEmFrac > epsilon && nSigmaEm < maxSigma) ? amp*exp(-nSigmaEm) : 0;

      // V/m/MHz at 1m due to HAD component of shower
      double nSigmaHad = (viewAngleRadians-fCherenkovAngleRad)*(viewAngleRadians-fCherenkovAngleRad)/(deltaThetaHad*deltaThetaHad);
      double vmmhz1m_had = (fHadFrac != 0 && nSigmaHad) < maxSigma ? amp*exp(-nSigmaHad) : 0;

      // Sum the EM and hadronic components
      amp = sin_viewAngle*(fEmFrac*vmmhz1m_em + fHadFrac*vmmhz1m_had);
    }

    freq += fDeltaFreqHz;
  }

  fTotalPowerDirty = true;
}


double icemc::AskaryanFreqs::operator[](int i) const {
  if(i >= 0 && i < vmmhz.size()){
    return vmmhz[i];
  }
  else {
    Log() << icemc::error << "Attempt to access icemc::AskaryanFreqs at index "
	  << i << ", which is out of bounds. (min=0, max="
	  << vmmhz.size() <<"). Returning 0.\n";
    return 0;
  }
}




TGraph icemc::AskaryanFreqs::makeGraph() const {
  std::vector<double> freqs;
  freqs.reserve(vmmhz.size());
  for(int i=0; i < vmmhz.size(); i++){
    freqs.push_back(fMinFreqHz + i*fDeltaFreqHz);
  }
  TGraph gr(vmmhz.size(), &freqs[0], &vmmhz[0]);
  gr.SetTitle(";Frequency (Hz); Askaryan E-field (V/m/MHz)");
  return gr;
}

