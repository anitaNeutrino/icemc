#include "VoltsRX.h"


icemc::VoltsRX::VoltsRX(int nRX) : fNRX(nRX) {

  max = 0;
  ave = 0;
  sum = 0;

  max_highband = 0; // max voltage seen on an antenna - just for debugging purposes
  max_lowband = 0; // max voltage seen on an antenna - just for debugging purposes

  rfcm_lab_e_all.reserve(nRX);
  rfcm_lab_h_all.reserve(nRX);
  for(int i=0; i < fNRX; i++){
    rfcm_lab_e_all.emplace_back(std::array<double, Anita::HALFNFOUR>());
    rfcm_lab_h_all.emplace_back(std::array<double, Anita::HALFNFOUR>());
  }
  reset();
}


void icemc::VoltsRX::reset() {
  
  for(auto& rfcm_lab_e : rfcm_lab_e_all){
    rfcm_lab_e.fill(0);
  }
  for(auto& rfcm_lab_h : rfcm_lab_h_all){
    rfcm_lab_h.fill(0);
  }
}
