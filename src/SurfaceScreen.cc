#include <math.h>
#include <iostream>

#include "TVector3.h"
#include "Geoid.h"
#include "SurfaceScreen.h"
#include "Detector.h"
#include "FTPair.h"
#include "Tools.h"
#include "RayTracer.h" ///@todo remove when test complete

icemc::SurfaceScreen::SurfaceScreen(int a){
  if( !(a%2) )
    a++;                    // force 'a' odd to set screen steps correctly

  //std::cerr << "Generating default screen" << std::endl;
  fedgeLength=1.;
  fcentralPoint = Geoid::Position(0, 0, 0);
  fnormal = TVector3(1.,1.,1.);
  funit_x = TVector3(1.,1.,1.);
  funit_y = TVector3(1.,1.,1.);

  fNsamples = a;
};


void icemc::SurfaceScreen::SetNsamples(int i){
  fNsamples = i;
};


int icemc::SurfaceScreen::GetNsamples() const{
  return fNsamples;
};


void icemc::SurfaceScreen::SetEdgeLength(double a){
  fedgeLength = a;
};


void icemc::SurfaceScreen::SetCentralPoint(Geoid::Position a){
  fcentralPoint = a;
};

void icemc::SurfaceScreen::SetCosineProjectionFactor(double a){
  fcosineProjectionFactor = a;
};


double icemc::SurfaceScreen::GetCosineProjectionFactor() const {
  return fcosineProjectionFactor;
};

void icemc::SurfaceScreen::SetNormal(TVector3 a){
  fnormal = a.Unit();
};

void icemc::SurfaceScreen::SetUnitX(TVector3 a){
  funit_x = a;
};

void icemc::SurfaceScreen::SetUnitY(TVector3 a){
  funit_y = a;
};

double icemc::SurfaceScreen::GetEdgeLength() const {
  return fedgeLength;
};


Geoid::Position icemc::SurfaceScreen::GetCentralPoint() const {
  return fcentralPoint;
};


TVector3 icemc::SurfaceScreen::GetNormal() const {
  return fnormal;
};


TVector3 icemc::SurfaceScreen::GetUnitX() const {
  return funit_x;
};


TVector3 icemc::SurfaceScreen::GetUnitY() const {
  return funit_y;
};


double icemc::SurfaceScreen::CalcXindex(int i) const {
  return (double) (i % (fNsamples));
};


double icemc::SurfaceScreen::CalcYindex(int i) const {
  return (double) floor(i / (fNsamples));
};


Geoid::Position icemc::SurfaceScreen::GetPosition(int i, int j) const {
  Geoid::Position pos;


  // this picks points that are NOT on the edge
  pos = fcentralPoint                                       // base
        //- 0.5*fedgeLength*funit_x - 0.5*fedgeLength*fcosineProjectionFactor*funit_y   // shift to a corner
        //+ (1./(2.*(double)(fNsamples)))*fedgeLength*(funit_x + fcosineProjectionFactor*funit_y)  // move off the edge
        //+ (xindex/((double)(fNsamples)))*fedgeLength*funit_x   // move by x-increment
        //+ (yindex/((double)(fNsamples)))*fedgeLength*fcosineProjectionFactor*funit_y;  // move by y-increment with the cosine projection correction
        + i * fedgeLength*funit_x
        + j * fedgeLength*funit_y;

  return pos;
};


void icemc::SurfaceScreen::AddVmmhz_freq(double A){
  //here i refers to the 'screen' loop
  fVmmhz_freq.push_back(A);
};


double icemc::SurfaceScreen::GetVmmhz_freq(int i) const {
  return fVmmhz_freq[i];
};


void icemc::SurfaceScreen::AddVmmhz0(double A){
  fVmmhz0.push_back(A);
};

double icemc::SurfaceScreen::GetVmmhz0(int i) const {
  return fVmmhz0[i];
};

void icemc::SurfaceScreen::AddViewangle(double A){
  fViewangle.push_back(A);
};

double icemc::SurfaceScreen::GetViewangle(int i) const {
  return fViewangle[i];
};


void icemc::SurfaceScreen::AddDelay(double A){
  fDelays.push_back(A);
};


double icemc::SurfaceScreen::GetDelay(int i) const {
  return fDelays[i];
};


void icemc::SurfaceScreen::SetNvalidPoints(int i){
  fNvalidpoints = i;
};


double icemc::SurfaceScreen::GetNvalidPoints() const {
  return fNvalidpoints;
};


void icemc::SurfaceScreen::AddVec2bln(TVector3 v){
  fVec2blns.push_back(v);
};


TVector3 icemc::SurfaceScreen::GetVec2bln(int i) const {
  return fVec2blns[i];
};


void icemc::SurfaceScreen::AddPol(TVector3 v){
  fPols.push_back(v);
};


TVector3 icemc::SurfaceScreen::GetPol(int i) const {
  return fPols[i];
};


void icemc::SurfaceScreen::AddImpactPt(Geoid::Position p){
  fImpactPt.push_back(p);
};


Geoid::Position icemc::SurfaceScreen::GetImpactPt(int i) const {
  return fImpactPt[i];
};


void icemc::SurfaceScreen::AddWeight(double a){
  fWeight.push_back(a);
};

double icemc::SurfaceScreen::GetWeight(int i) const {
  return fWeight[i];
};

void icemc::SurfaceScreen::SetWeightNorm(double a){
  fWeightNorm = a;
};

double icemc::SurfaceScreen::GetWeightNorm() const {
  return fWeightNorm;
};

void icemc::SurfaceScreen::AddIncidenceAngle(double A){
  fIncAngles.push_back(A);
};

double icemc::SurfaceScreen::GetIncidenceAngle(int i) const {
  return fIncAngles[i];
};

void icemc::SurfaceScreen::AddTransmissionAngle(double A){
  fTransAngles.push_back(A);
};

double icemc::SurfaceScreen::GetTransmissionAngle(int i) const {
  return fTransAngles[i];
};

void icemc::SurfaceScreen::AddFacetLength(double A){
  fFacetLength.push_back(A);
};

double icemc::SurfaceScreen::GetFacetLength(int i) const {
  return fFacetLength[i];
};

void icemc::SurfaceScreen::AddTparallel_polParallel(double A){
  fTcoeff_parl_polparl.push_back(A);
};

double icemc::SurfaceScreen::GetTparallel_polParallel(int i) const {
  return fTcoeff_parl_polparl[i];
};

void icemc::SurfaceScreen::AddTperpendicular_polParallel(double A){
  fTcoeff_perp_polparl.push_back(A);
};

double icemc::SurfaceScreen::GetTperpendicular_polParallel(int i) const {
  return fTcoeff_perp_polparl[i];
};

void icemc::SurfaceScreen::AddTparallel_polPerpendicular(double A){
  fTcoeff_parl_polperp.push_back(A);
};

double icemc::SurfaceScreen::GetTparallel_polPerpendicular(int i) const {
  return fTcoeff_parl_polperp[i];
};

void icemc::SurfaceScreen::AddTperpendicular_polPerpendicular(double A){
  fTcoeff_perp_polperp.push_back(A);
};

double icemc::SurfaceScreen::GetTperpendicular_polPerpendicular(int i) const {
  return fTcoeff_perp_polperp[i];
};

void icemc::SurfaceScreen::ResetParameters(){
  // reset these in icemc:
  // Nsamples
  // edge length
  // central point
  // normal
  // cosine projection factor
  // unit x /y vectors

  //need to reset everything else, like the vectors...
  fNvalidpoints = 0;

  fVmmhz_freq.clear();
  fVmmhz0.clear();
  fViewangle.clear();
  fDelays.clear();
  fVec2blns.clear();
  fPols.clear();
  fImpactPt.clear();
  fWeight.clear();
  fWeightNorm = 1.;
  fIncAngles.clear();
  fTransAngles.clear();
  fFacetLength.clear();
  fTcoeff_parl_polparl.clear();
  fTcoeff_perp_polparl.clear();
  fTcoeff_parl_polperp.clear();
  fTcoeff_perp_polperp.clear();
};


// void icemc::SurfaceScreen::PropagateSignalsToDetector(const Settings* settings1, ANITA* d, int inu) const {

//   //@todo remove hardcoded number of frequencies!
//   std::vector<std::complex<double> > tmp_vhz(Anita::NFREQ, 0);
//   const double df = d->freq[1] - d->freq[0];


//   // std::cout << " in screen df = " << df << std::endl;
//   // int n;
//   // double dt;
//   // d->getDesiredNDt(n, dt);

//   ///@todo remove hardcoding here!
//   double TIMESTEP=(1./2.6)*1.E-9; // time step between samples

//   const Geoid::Position detPos = d->getPosition();
//   double firstDelay = 0;

//   for (int jpt=0; jpt<GetNvalidPoints(); jpt++){
//     const Geoid::Position&rfExit = GetImpactPt(jpt);
//     const double nominalTimeOfFlightSeconds = (detPos - rfExit).Mag()/constants::CLIGHT;

//     for(int rx = 0; rx < d->getNumRX(); rx++){
//       TVector3 rxPos = d->getPositionRX(rx);
//       const double rxTimeOfFlightSeconds = (rxPos - rfExit).Mag()/constants::CLIGHT;

//       ///@todo remove hardcoded number of frequencies!
//       for (int k=0;k<Anita::NFREQ;k++) {

// 	/// @todo put this scaling inside ChanTrigger because it's a change in units the ChanTrigger
// 	// cares about, but it has nothing to do with the Screen

// 	double amp = GetVmmhz_freq(jpt*Anita::NFREQ + k)/sqrt(2)/(TIMESTEP*1.E6);
// 	// here we mimic the "putting everything in imag" as done in MakeArrayForFFT
// 	// the normal ordering nonsense (which hugely affects the phase) is done in FTPair
// 	tmp_vhz[k].real(0);
// 	tmp_vhz[k].imag(amp);

// 	/**
// 	 * @todo this is the most disgusting fudge factor hack and needs to be dealt with comprehensively in the future.
// 	 * When migrating to FTPair from stupid raw frequency arrays icemc pads 128 frequency bins to 256.
// 	 * i.e. It has 128 frequency bins (what would be 129 if the nyquist was properly handled), as if there were 256 time samples.
// 	 * MakeArrayForFFT pads things for there to be 256 frequency bins (257 w/ nyquist) i.e. 512 time samples.
// 	 * The inverse FFT (correctly) scales things down by a factor of N.
// 	 * But is the correct N 256, or 512, or something else?
// 	 * The FTPair class handles the general case of this normalization, but I'm not sure icemc does.
// 	 * (It must  be one of the random factors of 2 or sqrt(2) littered around?)
// 	 * Anyway I need to apply this factor to make sure the PSDs agree before and after refactoring.
// 	 * (Where the PSD accounts for the difference in length correctly.)
// 	 * But I need to fudge the amplitudes here so the PSDs agree during the refactor.
// 	 * Before doing this the new PSD (which squares the voltages) had twice the amplitude of the old chanTrigger, as you would expect.
// 	 */
// 	tmp_vhz[k] /= sqrt(2);
//       }

//       FTPair signal(tmp_vhz, df, true);

//       double relativeDelaySeconds = rxTimeOfFlightSeconds - nominalTimeOfFlightSeconds;
//       if(rx==0){
// 	firstDelay = relativeDelaySeconds;
//       }
//       relativeDelaySeconds -= firstDelay;
      
//       //@todo set up delay to antennas
//       // signal.applyConstantGroupDelay(relativeDelaySeconds + 20e-9, false);
//       signal.applyConstantGroupDelay(relativeDelaySeconds, false);      

//       // TGraph& gr = signal.changeTimeDomain();
//       // // add a point to force  up to next power of 2...
//       // gr.SetPoint(gr.GetN(), gr.GetX()[gr.GetN()-1] + gr.GetX()[1] - gr.GetX()[0], 0);

//       // add a point to force  up to next power of 2...
//       // auto& cs = signal.changeFreqDomain();
//       // for(int i=0; i < 10; i++){
//       // 	cs.push_back(0);
//       // }

//       // is that the correct geometry?
//       // @todo find out where we get vec2bln from... and see if Screen really needs to know
//       PropagatingSignal s(signal, GetVec2bln(jpt), GetPol(jpt));

//       d->addSignalToRX(s, rx, inu);
//     }
//   }
// }
