#include "AskaryanFreqs.h"
#include "IcemcLog.h"
#include "anita.hh"
#include "TGraph.h"

ClassImp(icemc::AskaryanFreqs)

icemc::AskaryanFreqs::AskaryanFreqs()
: vmmhz(Anita::NFREQ, 0), vmmhz_em(Anita::NFREQ, 0), minFreqHz(0), maxFreqHz(0)
{

}


icemc::AskaryanFreqs::AskaryanFreqs(int nf, double fMinHz, double fMaxHz, const double* vmmhz_input, const double* vmmhz_em_input)
  : vmmhz(vmmhz_input, vmmhz_input + nf),  minFreqHz(fMinHz), maxFreqHz(fMaxHz)
{
  if(vmmhz_em_input){
    vmmhz_em.assign(vmmhz_em_input, vmmhz_em_input + nf);
  }
  else{
    vmmhz_em.assign(nf, 0);
  }
  // std::cout << minFreqHz  << "\t" << maxFreqHz << std::endl;
}



void icemc::AskaryanFreqs::applyTapering(double viewAngleRadians){

  for(UInt_t k=0; k < vmmhz.size(); k++){
    // deltheta_em[k] = deltheta_em_max*anita1->FREQ_LOW/anita1->freq[k];
    // deltheta_had[k] = deltheta_had_max*anita1->FREQ_LOW/anita1->freq[k];
  }
}


double icemc::AskaryanFreqs::operator[](int i) const {
  if(i >= 0 && i < vmmhz.size()){
    return vmmhz[i];
  }
  else {
    Log() << icemc::error << "Attempt to access icemc::AskaryanFreqs at index " << i << ", which is out of bounds. (min=0, max=" << vmmhz.size() <<"). Returning 0.\n";
    return 0;
  }
}


TGraph icemc::AskaryanFreqs::makeGraph() const {
  const double df = (maxFreqHz - minFreqHz)/vmmhz.size();
  std::vector<double> freqs;
  freqs.reserve(vmmhz.size());
  for(int i=0; i < vmmhz.size(); i++){
    freqs.push_back(i*df);
  }
  TGraph gr(vmmhz.size(), &freqs[0], &vmmhz[0]);
  gr.SetTitle(";Frequency (Hz); Askaryan E-field (V/m/MHz)");
  return gr;
}

