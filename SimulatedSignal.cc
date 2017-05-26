#include <iostream>
//#include "anita.hh"
#include "SimulatedSignal.h"

///< SimulatedSignal inherits from RFSignal

SimulatedSignal::SimulatedSignal()
 :RFSignal()
{
  // Default constructor
}


///< Define an SimulatedSignal from an icemc askaryan field (vmmhz)
///< Icemc askaryan field is defined between 200 MHz and 1200 MHz
///< So we interpolate it from 0 to 200 MHz
///< 0 pad it until 1300 MHz
///< And renormalise it taking into account the number of non-zero
///< amplitudes before and after
///< For the moment the phase is always -90 degrees

SimulatedSignal::SimulatedSignal(int nfreqs0, double *freqs0, double *freqAmp0)
  :RFSignal()
{

  // First implementation of Askaryan field approximate all phases to -90 degrees
  double phase = (TMath::Pi()/2)*(-1);
  
  Double_t tempFreqVals[500]   ;
  Double_t tempMags[500]       ;
  Double_t tempPhases[500]     ;

  tempFreqVals[0] = 0;
  tempMags[0]     = 0;
  tempPhases[0]   = phase;

  int ndt0 = ((nfreqs0)*2+1);
  int ndt  = floor(1000e6/newdf)*2+1;
  // Renormalise because we are adding power when interpolating between 0 and 200 MHz
  double norm = TMath::Sqrt(ndt0*1./ndt);

  // Graph to interpolate between points properly
  TGraph *gfreq0 = new TGraph(nfreqs0, freqs0, freqAmp0);

  for (int ifreq=1;ifreq<=nfreqs;ifreq++){
    tempFreqVals[ifreq]=newdf*ifreq;
    if (tempFreqVals[ifreq]<freqs0[0]) {
      // Points between 0 and 200 MHz are interpolated with a straight line
      tempMags[ifreq] = freqAmp0[0]*tempFreqVals[ifreq]*norm/0.2e9;
    } else if (tempFreqVals[ifreq]<=(freqs0[nfreqs0-1]+newdf)) {
      // Points between 200 and 1200 MHz are interpolated using the TGraph
      tempMags[ifreq] = gfreq0->Eval(tempFreqVals[ifreq])*norm;
    } else {
      // Points between 1200 and 1300 MHz are set to 0
      tempMags[ifreq] = 0;
    }
    tempPhases[ifreq]=phase;
   
  }
  
  // clean up
  delete gfreq0;

  // Set frequencies, magnitudes and phases of the RFSignal
  setFreqs(nfreqs, tempFreqVals);
  setMagsPhases(tempMags, tempPhases);

  // Update the time domain
  updateTimeDomain();
}


///< Constructor from time domain values

SimulatedSignal::SimulatedSignal(Int_t numPoints,Double_t *tVals,Double_t *vVals,Int_t mvNs)
  :RFSignal(numPoints,tVals,vVals,mvNs)
{

}

///< Default destructor

SimulatedSignal::~SimulatedSignal() 
{
//Default destructor
}



///< Add CW to simulated signal
///< frequency: is the CW frequency (from power sprectum?)
///< phase: is the phase, coming from the dt of propagation
///<        from the source to the antenna face
///< amplitude: is the CW amplitude (from power spectrum?) // good is 0.01

void SimulatedSignal::AddCW(double frequency, double phase, double amplitude){

  double deltaT = (1/2.6)*1e-9;
  double omega;
  double volts_cw[512];
  
  for (int itime=0; itime<fNpoints; itime++){
    omega=TMath::Pi()*2*frequency;
    volts_cw[itime]=amplitude*TMath::Sin(omega*itime*deltaT + phase);
  }

  RFSignal *tmp = new RFSignal(fNpoints, fX, volts_cw);

  addToSignal(tmp);

  delete tmp;
}
