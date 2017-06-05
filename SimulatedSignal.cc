#include <iostream>
#include <array>
#include "vector.hh"
#include "Tools.h"
#include "anita.hh"
#include "SimulatedSignal.h"

///< SimulatedSignal inherits from RFSignal

SimulatedSignal::SimulatedSignal()
 :RFSignal()
{
  // Default constructor
  nfreqs = 257; //Anita::HALFNFOUR/2;
  newdf  = 1300e6/nfreqs;
  
}


///< Constructor from time domain values

SimulatedSignal::SimulatedSignal(Int_t numPoints,Double_t *tVals,Double_t *vVals,Int_t mvNs)
  :RFSignal(numPoints,tVals,vVals,mvNs)
{
  nfreqs = 257; //Anita::HALFNFOUR/2;
  newdf  = 1300e6/nfreqs;
}

///< Default destructor

SimulatedSignal::~SimulatedSignal() 
{
//Default destructor
}




///< Define a SimulatedSignal from an icemc askaryan field (vmmhz)
///< Icemc askaryan field is defined between 200 MHz and 1200 MHz
///< So we interpolate it from 0 to 200 MHz
///< 0 pad it until 1300 MHz
///< And renormalise it taking into account the number of non-zero
///< amplitudes before and after
///< For the moment the phase is always -90 degrees

void SimulatedSignal::updateSimSignalFromVmmhz(int nfreqs0, double *freqs0, double *freqAmp0)
{

  // First implementation of Askaryan field approximate all phases to -90 degrees
  double phase = (TMath::Pi()/2)*(-1);
  
  Double_t tempFreqVals[500]   ;
  Double_t tempMags[500]       ;
  Double_t tempPhases[500]     ;

  tempFreqVals[0] = 0;
  tempMags[0]     = 0;
  tempPhases[0]   = phase;

  // Renormalise because we are adding power when interpolating between 0 and freqs0[0] (200 MHz)
  int firstNonZero = Tools::Getifreq(freqs0[0],0,newdf*nfreqs,nfreqs);
  int lastNonZero  = Tools::Getifreq(freqs0[nfreqs0-1],0,newdf*nfreqs,nfreqs);
  double norm=TMath::Sqrt(double(nfreqs0)/double(lastNonZero-firstNonZero));
  // cout << firstNonZero << " " << lastNonZero << " " << lastNonZero-firstNonZero << " " << norm << endl;

  
  //  cout << norm << endl;
  // Graph to interpolate between points properly
  TGraph *gfreq0 = new TGraph(nfreqs0, freqs0, freqAmp0);

  for (int ifreq=1;ifreq<nfreqs;ifreq++){
    tempFreqVals[ifreq]=newdf*ifreq;
    if (tempFreqVals[ifreq]<freqs0[0]) {
      // Points between 0 and 200 MHz are interpolated with a straight line
      tempMags[ifreq] = freqAmp0[0]*tempFreqVals[ifreq]*norm/freqs0[0];
    } else if (tempFreqVals[ifreq]<=(freqs0[nfreqs0-1]+newdf)) {
      // Points between 200 and 1200 MHz are interpolated using the TGraph
      tempMags[ifreq] = gfreq0->Eval(tempFreqVals[ifreq])*norm;
    } else {
      // Points between 1200 and 1300 MHz are set to 0
      tempMags[ifreq] = 0;
    }
    tempPhases[ifreq]=phase;
    // cout << "Set frequency " << ifreq << " " << tempFreqVals[ifreq] << endl;
  }

  // clean up
  delete gfreq0;

  // Set frequencies, magnitudes and phases of the RFSignal
  setFreqs(nfreqs-1, tempFreqVals);
  setMagsPhases(tempMags, tempPhases);

  // Update the time domain
  updateTimeDomain();
}




///< Add CW to simulated signal
///< frequency: is the CW frequency (from power sprectum?)
///< phase: is the phase, coming from the dt of propagation
///<        from the source to the antenna face
///< amplitude: is the CW amplitude (from power spectrum?) // good is 0.01

void SimulatedSignal::addCW(double frequency, double phase, double amplitude){

  // double deltaT = (1/2.6)*1e-9;
  double omega;
  double volts_cw[512];
  double *times = GetX();
  
  for (int itime=0; itime<fNpoints; itime++){
    omega=TMath::Pi()*2*frequency;
    volts_cw[itime]=amplitude*TMath::Sin(omega*times[itime] + phase);
  }

  RFSignal *tmp = new RFSignal(fNpoints, times, volts_cw);

  addToSignal(tmp);

  delete tmp;
}


void SimulatedSignal::getVmmhz(Anita *anita1, double *vmmhz){


  double *freqs = getFreqs();
  double *mags  = getMags();
  int numFreqs  = nfreqs-1;
  double freqMin = freqs[0];
  double freqMax = freqs[numFreqs-1];//+freqs[1]-freqs[0];
  double freqMin1 = anita1->freq[0];
  double freqMax1 = anita1->freq[anita1->NFREQ-1];

  int firstNonZero = Tools::Getifreq(freqMin1,freqMin,freqMax,numFreqs);
  int lastNonZero  = Tools::Getifreq(freqMax1,freqMin,freqMax,numFreqs);
  double norm=TMath::Sqrt(double(lastNonZero-firstNonZero)/double(anita1->NFREQ));
  // cout << firstNonZero << " " << lastNonZero << " " << lastNonZero-firstNonZero << " " << norm << endl;

  
  // Graph to interpolate between points properly
  TGraph *gfreq0 = new TGraph(numFreqs, freqs, mags);
  for (int ifreq=0; ifreq<anita1->NFREQ; ifreq++){
    
    // int ifour=Tools::Getifreq(anita1->freq[ifreq],freqMin,freqMax,numFreqs);
    //  vmmhz[ifreq]=mags[ifour]*norm;
    // cout << freqs[ifour] << " " << anita1->freq[ifreq] << " " << vmmhz[ifreq] << " " << endl;
    vmmhz[ifreq]=gfreq0->Eval(anita1->freq[ifreq])*norm;
    //    cout << vmmhz[ifreq] << " " << endl;
  }

  delete gfreq0;
  
}
