#include "VoltsRX.h"

icemc::VoltsRX::VoltsRX(){

  max = 0;
  ave = 0;
  sum = 0;

  max_highband = 0; // max voltage seen on an antenna - just for debugging purposes
  max_lowband = 0; // max voltage seen on an antenna - just for debugging purposes

  memset(rfcm_lab_e_all, 0, sizeof(double)*nRX*nSamp);
  memset(rfcm_lab_h_all, 0, sizeof(double)*nRX*nSamp);  
  
}
