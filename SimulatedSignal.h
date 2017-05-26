#ifndef SIMULATEDSIGNAL_H
#define SIMULATEDSIGNAL_H
#include "RFSignal.h"

//!  This is a wrapper class for an Simulated Signal
/*!
  It inheriths from RFSignal
*/
class SimulatedSignal : public RFSignal {

 public:

  SimulatedSignal(); ///<Default constructor
  
  SimulatedSignal(int nfreq0, double *freq0, double *freqAmp0); ///< Constructor from icemc askaryan field (vmmhz) defined from 200MHz to 1200MHz
  
  SimulatedSignal(Int_t numPoints,Double_t *tVals,Double_t *vVals,Int_t mvNs); ///< Constructor from time domain waveform

    ~SimulatedSignal(); ///<Destructor
  
 private:
  
  // The final number of frequencies is 256 between 0 and 1300 MHz
  Int_t    nfreqs = 256; //Anita::HALFNFOUR/2;
  Double_t newdf  = 1300e6/nfreqs;
      
};

#endif // SIMULATEDSIGNAL_H
