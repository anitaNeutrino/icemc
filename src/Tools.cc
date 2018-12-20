#include "Tools.h"
#include <iostream>
#include <cmath>
#include "TSpline.h"
#include "TH2F.h"
#include <fstream>
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Constants.h"
#include "TStyle.h"

void icemc::Tools::ShiftLeft(double *x,const int n,int ishift) {
    
  double x_temp[n];
  // shift the x array to the left by ishift bins and fill the gap with zeroes
  for (int i=0;i<n;i++) {
    x_temp[i]=x[i];
  }

  for (int i=0;i<n-ishift;i++) {
    if(i+ishift >= n || i+ishift < 0 || i < 0 || i >= n){
      std::cerr << "It's retard 1 with " << i+ishift << ", " << i << std::endl;
    }
    x[i]=x_temp[i+ishift];
  }
  for (int i=n-ishift;i<n;i++) {
    if(i >= n || i < 0){
      std::cerr << "It's retard 2 with " << i << std::endl;
    }
    x[i]=0.;
  }
}


void icemc::Tools::ShiftRight(double *x,const int n,int ishift) {
    
  double x_temp[n];
  // shift the x array to the right by ishift bins and fill the gap with zeroes
  for (int i=0;i<n;i++) {
    x_temp[i]=x[i];
  }
  for (int i=ishift;i<n;i++) {
    x[i]=x_temp[i-ishift];
  }
  for (int i=0;i<ishift;i++) {
    x[i]=0.;
  }
}



double icemc::Tools::GetFWHM(TH1 *h1) {
  int imax=h1->GetMaximumBin();
  double max=h1->GetMaximum();
    
  //  std::cout << "imax, max are " << imax << " " << max << "\n";
  int ibin_plus=0;
  int ibin_minus=0;
  // now step to the right until it's half
  for (int ibin=imax;ibin<=h1->GetNbinsX();ibin++) {
    if (h1->GetBinContent(ibin)<max/2.) {
      ibin_plus=ibin;
      ibin=h1->GetNbinsX()+1;
      //  std::cout << "ibin_plus is " << ibin_plus << "\n";
    }
  }
  // now step to the left
  for (int ibin=imax;ibin>=1;ibin--) {
    if (h1->GetBinContent(ibin)<max/2.) {
      ibin_minus=ibin;
      ibin=0;
      //  std::cout << "ibin_minus is " << ibin_minus << "\n";
    }
  }

  if (ibin_plus>0 && ibin_minus==0) {
    ibin_minus=1;
    //std::cout << "bin_minus is " << ibin_minus << "\n";
  }
    
  if (ibin_plus==0 && ibin_minus==0) {
    std::cout << "Found 0 FWHM.\n";
    return 0.;
  }
    
  return (h1->GetBinCenter(ibin_plus)-h1->GetBinCenter(ibin_minus))/2.;
}




int icemc::Tools::Getifreq(double freq, double freq_low, double freq_high, int n) {

  if (freq>=freq_high || freq < freq_low){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << " requested out of bounds frequency " << freq
	      << " when bounds are " << freq_low << "->" << freq_high << "\n";
    return -1;
  }
  return (int)((freq-freq_low)/(freq_high-freq_low)*(double)n);
} //Getifreq



void icemc::Tools::NormalTimeOrdering(const int n, double* volts) {
  double volts_temp[n];
  for (int i=0;i<n/2;i++) {
    volts_temp[i]=volts[i+n/2];
    volts_temp[i+n/2]=volts[i];
  }
  for (int i=0;i<n;i++) {
    volts[i]=volts_temp[i];
  }
}

void icemc::Tools::reverseTimeOrdering(const int n, int* bitsin,int *bitsout) {
  int bits_temp[n];
  for (int i=0;i<n;i++) {
    bits_temp[i]=bitsin[n-i-1];
  }
  for (int i=0;i<n;i++) {
    bitsout[i]=bitsin[i];
  }
}


int icemc::Tools::findIndex(double *freqlab,double freq,int npoints,double min,double max) {
    
  //  std::cout << "inside findIndex, freq, min are " << freq << " " << min << "\n";
    
  if (freq<min)
    return -1;
  if (freq>=min && freq<freqlab[1])
    return 0;
  if (freq>=freqlab[npoints-1] && freq<max)
    return npoints-1;
  return ((int)((freq-freqlab[0])/(freqlab[1]-freqlab[0])))+1;

  return -1;
}




TGraph *icemc::Tools::getInterpolatedGraph(TGraph *grIn, Double_t deltaT)
{
  //Will use the ROOT::Math::Interpolator function to do this.
  std::vector<double> tVec;
  std::vector<double> vVec;
   
  Int_t numIn=grIn->GetN();
  Double_t tIn,vIn;

  Double_t startTime=0;
  Double_t lastTime=0;
  for (int samp=0;samp<numIn;samp++) {
    grIn->GetPoint(samp,tIn,vIn);
    tVec.push_back(tIn);
    vVec.push_back(vIn);
    //std::cout << "samp " << samp << " t " << tIn << " v " << vIn << " this-last " << tIn-tVec[tVec.size()-2] << std::endl;
    if(samp==0)
      startTime=tIn;
    lastTime=tIn;
  }
  if(tVec.size()<1) {
    std::cout << "Insufficent points for interpolation\n";
    return NULL;
  }

  //Bastards
  ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);
   
  Int_t roughPoints=Int_t((lastTime-startTime)/deltaT);
   

  Double_t *newTimes = new Double_t[roughPoints+100]; //Will change this at some point, but for now
  Double_t *newVolts = new Double_t[roughPoints+100]; //Will change this at some point, but for now
  Int_t numPoints=0;
  for(Double_t time=startTime;time<=lastTime;time+=deltaT) {
    newTimes[numPoints]=time;
    newVolts[numPoints]=chanInterp.Eval(time);
    //      std::cout << numPoints << "\t" << newTimes[numPoints]
    //      		<< "\t" << newVolts[numPoints] << std::endl;
	       
    numPoints++;
  }

  TGraph *grInt = new TGraph(numPoints,newTimes,newVolts);
  delete [] newTimes;
  delete [] newVolts;
  return grInt;

}
