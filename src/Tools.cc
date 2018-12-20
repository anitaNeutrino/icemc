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




void icemc::Tools::Zero(int *anarray,int n) {
  for (int i=0;i<n;i++) {
    anarray[i]=0;
  } //for
} //Zero (int*,int)


void icemc::Tools::Zero(double *anarray,int n) {
  for (int i=0;i<n;i++) {
    anarray[i]=0.;
  } //for
} //Zero (int*,int)



double icemc::Tools::dMax(const double *x,int n) {
    
  double max=x[0];
  for (int k=1;k<n;k++) {
    if (x[k]>max){
      max=x[k];
    }
  }
  return max;
} //dMax(double*, int)



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


void icemc::Tools::get_random_rician(double signal_amplitude, double signal_phase, double sigma, double &amplitude, double &phase){
  double rand_gauss_a, rand_gauss_b;
  get_circular_bivariate_normal_random_variable(rand_gauss_a, rand_gauss_b);
    
  // Gives the gaussian-distributed random variables a standard deviation of sigma
  rand_gauss_a *= sigma;
  rand_gauss_b *= sigma;
    
  // Gives the gaussian-distributed random variables a mean of (v*cos(theta), v*sin(theta)) when v is the mean of the desired rician distribution
  rand_gauss_a += signal_amplitude * cos(signal_phase);
  rand_gauss_b += signal_amplitude * sin(signal_phase);
    
  // The Rician Distribution produces the probability of the the absolute value (radius) of a circular bivariate normal random variable:
  amplitude = sqrt(rand_gauss_a * rand_gauss_a + rand_gauss_b * rand_gauss_b);
  // Thus, the descriptor other than amplitude for the circular bivariate is given by a phase:
  phase = atan2(rand_gauss_b, rand_gauss_a);
  return;
}


void icemc::Tools::get_circular_bivariate_normal_random_variable(double& rand_gauss_a, double& rand_gauss_b){
  double rand_uni_a = gRandom->Rndm(); //gRandom->Rndm() produces uniformly-distributed floating points in ]0,1]
  double rand_uni_b = gRandom->Rndm();
  double first_term = sqrt(-2. * log(rand_uni_a));
  // Box-Muller transform from a bivariate uniform distribution from 0 to 1 to a gaussian with mean = 0 and sigma = 1
  rand_gauss_a = first_term * cos(2. * M_PI * rand_uni_b);
  rand_gauss_b = first_term * sin(2. * M_PI * rand_uni_b);
  return;
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



TStyle* icemc::Tools::RootStyle() {

  TStyle *RootStyle = new TStyle("Root-Style", "The Perfect Style for Plots ;-)");

#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                  // save the global style reference
  gStyle = RootStyle;
#endif
  // otherwise you need to call TROOT::SetStyle("Root-Style")

  // Paper size
  RootStyle->SetPaperSize(TStyle::kUSLetter);

  // Canvas
  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderSize(10);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH      (600);
  RootStyle->SetCanvasDefW      (600);
  RootStyle->SetCanvasDefX      (10);
  RootStyle->SetCanvasDefY      (10);

  // Pads
  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderSize  (10);
  RootStyle->SetPadBorderMode  (0);
  //  RootStyle->SetPadBottomMargin(0.13);
  RootStyle->SetPadBottomMargin(0.16);
  RootStyle->SetPadTopMargin   (0.08);
  RootStyle->SetPadLeftMargin  (0.18);
  RootStyle->SetPadRightMargin (0.05);
  RootStyle->SetPadGridX       (0);
  RootStyle->SetPadGridY       (0);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);

  // Frames
  RootStyle->SetFrameFillStyle ( 0);
  RootStyle->SetFrameFillColor ( 0);
  RootStyle->SetFrameLineColor ( 1);
  RootStyle->SetFrameLineStyle ( 0);
  RootStyle->SetFrameLineWidth ( 2);
  RootStyle->SetFrameBorderSize(10);
  RootStyle->SetFrameBorderMode( 0);


  // Histograms
  RootStyle->SetHistFillColor(0);
  RootStyle->SetHistFillStyle(1);
  RootStyle->SetHistLineColor(1);
  RootStyle->SetHistLineStyle(0);
  RootStyle->SetHistLineWidth(2);

  // Functions
  RootStyle->SetFuncColor(1);
  RootStyle->SetFuncStyle(0);
  RootStyle->SetFuncWidth(1);

  //Legends
  RootStyle->SetStatBorderSize(2);
  RootStyle->SetStatFont      (42);
  //  RootStyle->SetOptStat       (111111);
  RootStyle->SetOptStat       (0);
  RootStyle->SetStatColor     (0);
  RootStyle->SetStatX         (0.93);
  RootStyle->SetStatY         (0.90);
  RootStyle->SetStatFontSize  (0.07);
  //  RootStyle->SetStatW         (0.2);
  //  RootStyle->SetStatH         (0.15);

  // Labels,  Ticks,  and Titles
  RootStyle->SetTickLength ( 0.015, "X");
  RootStyle->SetTitleSize  ( 0.10, "X");
  RootStyle->SetTitleOffset( 1.20, "X");
  RootStyle->SetTitleBorderSize(0);
  //  RootStyle->SetTitleFontSize((double)3.);
  RootStyle->SetLabelOffset( 0.015, "X");
  RootStyle->SetLabelSize  ( 0.050, "X");
  RootStyle->SetLabelFont  ( 42   , "X");
  RootStyle->SetTickLength ( 0.015, "Y");
  RootStyle->SetTitleSize  ( 0.10, "Y");
  RootStyle->SetTitleOffset( 0.600, "Y");
  RootStyle->SetLabelOffset( 0.015, "Y");
  RootStyle->SetLabelSize  ( 0.050, "Y");
  RootStyle->SetLabelFont  ( 42   , "Y");
  RootStyle->SetTitleFont  (42, "XY");
  RootStyle->SetTitleColor  (1);

  // Options
  RootStyle->SetOptFit     (0);
  RootStyle->SetMarkerStyle(20);
  RootStyle->SetMarkerSize(0.4);

  //  std::cout << ">> Style initialized with the Root Style!" << std::endl;
  //  std::cout << ">> " << modified << std::endl << std::endl;
  return RootStyle;
}
//end RootStyle()
