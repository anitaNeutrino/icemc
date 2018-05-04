#include <math.h>
#include <iostream>

#include "vector.hh"
#include "position.hh"
#include "screen.hh"
#include "Detector.h"
#include "ANITA.h"
#include "anita.hh"
#include "FTPair.h"
#include "Tools.h"

icemc::Screen::Screen(int a){
  if( !(a%2) )
    a++;                    // force 'a' odd to set screen steps correctly

  //std::cerr << "Generating default screen" << std::endl;
  fedgeLength=1.;
  fcentralPoint = Position(1.,1.,1.);
  fnormal = Vector(1.,1.,1.);
  funit_x = Vector(1.,1.,1.);
  funit_y = Vector(1.,1.,1.);

  fNsamples = a;
};


void icemc::Screen::SetNsamples(int i){
  fNsamples = i;
};


int icemc::Screen::GetNsamples() const{
  return fNsamples;
};


void icemc::Screen::SetEdgeLength(double a){
  fedgeLength = a;
};


void icemc::Screen::SetCentralPoint(Position a){
  fcentralPoint = a;
};

void icemc::Screen::SetCosineProjectionFactor(double a){
  fcosineProjectionFactor = a;
};


double icemc::Screen::GetCosineProjectionFactor() const {
  return fcosineProjectionFactor;
};

void icemc::Screen::SetNormal(Vector a){
  fnormal = a.Unit();
};

void icemc::Screen::SetUnitX(Vector a){
  funit_x = a;
};

void icemc::Screen::SetUnitY(Vector a){
  funit_y = a;
};

double icemc::Screen::GetEdgeLength() const {
  return fedgeLength;
};


icemc::Position icemc::Screen::GetCentralPoint() const {
  return fcentralPoint;
};


icemc::Vector icemc::Screen::GetNormal() const {
  return fnormal;
};


icemc::Vector icemc::Screen::GetUnitX() const {
  return funit_x;
};


icemc::Vector icemc::Screen::GetUnitY() const {
  return funit_y;
};


double icemc::Screen::CalcXindex(int i) const {
  return (double) (i % (fNsamples));
};


double icemc::Screen::CalcYindex(int i) const {
  return (double) floor(i / (fNsamples));
};


icemc::Position icemc::Screen::GetPosition(int i, int j) const {
  Position pos;


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


void icemc::Screen::AddVmmhz_freq(double A){
  //here i refers to the 'screen' loop
  fVmmhz_freq.push_back(A);
};


double icemc::Screen::GetVmmhz_freq(int i) const {
  return fVmmhz_freq[i];
};


void icemc::Screen::AddVmmhz0(double A){
  fVmmhz0.push_back(A);
};

double icemc::Screen::GetVmmhz0(int i) const {
  return fVmmhz0[i];
};

void icemc::Screen::AddViewangle(double A){
  fViewangle.push_back(A);
};

double icemc::Screen::GetViewangle(int i) const {
  return fViewangle[i];
};


void icemc::Screen::AddDelay(double A){
  fDelays.push_back(A);
};


double icemc::Screen::GetDelay(int i) const {
  return fDelays[i];
};


void icemc::Screen::SetNvalidPoints(int i){
  fNvalidpoints = i;
};


double icemc::Screen::GetNvalidPoints() const {
  return fNvalidpoints;
};


void icemc::Screen::AddVec2bln(Vector v){
  fVec2blns.push_back(v);
};


icemc::Vector icemc::Screen::GetVec2bln(int i) const {
  return fVec2blns[i];
};


void icemc::Screen::AddPol(Vector v){
  fPols.push_back(v);
};


icemc::Vector icemc::Screen::GetPol(int i) const {
  return fPols[i];
};


void icemc::Screen::AddImpactPt(Position p){
  fImpactPt.push_back(p);
};


icemc::Position icemc::Screen::GetImpactPt(int i) const {
  return fImpactPt[i];
};


void icemc::Screen::AddWeight(double a){
  fWeight.push_back(a);
};

double icemc::Screen::GetWeight(int i) const {
  return fWeight[i];
};

void icemc::Screen::SetWeightNorm(double a){
  fWeightNorm = a;
};

double icemc::Screen::GetWeightNorm() const {
  return fWeightNorm;
};

void icemc::Screen::AddIncidenceAngle(double A){
  fIncAngles.push_back(A);
};

double icemc::Screen::GetIncidenceAngle(int i) const {
  return fIncAngles[i];
};

void icemc::Screen::AddTransmissionAngle(double A){
  fTransAngles.push_back(A);
};

double icemc::Screen::GetTransmissionAngle(int i) const {
  return fTransAngles[i];
};

void icemc::Screen::AddFacetLength(double A){
  fFacetLength.push_back(A);
};

double icemc::Screen::GetFacetLength(int i) const {
  return fFacetLength[i];
};

void icemc::Screen::AddTparallel_polParallel(double A){
  fTcoeff_parl_polparl.push_back(A);
};

double icemc::Screen::GetTparallel_polParallel(int i) const {
  return fTcoeff_parl_polparl[i];
};

void icemc::Screen::AddTperpendicular_polParallel(double A){
  fTcoeff_perp_polparl.push_back(A);
};

double icemc::Screen::GetTperpendicular_polParallel(int i) const {
  return fTcoeff_perp_polparl[i];
};

void icemc::Screen::AddTparallel_polPerpendicular(double A){
  fTcoeff_parl_polperp.push_back(A);
};

double icemc::Screen::GetTparallel_polPerpendicular(int i) const {
  return fTcoeff_parl_polperp[i];
};

void icemc::Screen::AddTperpendicular_polPerpendicular(double A){
  fTcoeff_perp_polperp.push_back(A);
};

double icemc::Screen::GetTperpendicular_polPerpendicular(int i) const {
  return fTcoeff_perp_polperp[i];
};

void icemc::Screen::ResetParameters(){
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


void icemc::Screen::PropagateSignalsToDetector(const Settings* settings1, ANITA* d) const {

  double tmp_vhz[Anita::NFREQ];
  double tmp_volts[Anita::NFOUR/2];

  // int n;
  // double dt;
  // d->getDesiredNDt(n, dt);

  ///@todo remove hardcoding here!
  double TIMESTEP=(1./2.6)*1.E-9; // time step between samples

  for(int rx = 0; rx < d->getNumRX(); rx++){  
    for (int jpt=0; jpt<GetNvalidPoints(); jpt++){
      for (int k=0;k<Anita::NFREQ;k++) {

	//Copy frequency amplitude to screen point
	// tmp_vhz[0][k] = tmp_vhz[1][k] = GetVmmhz_freq(jpt*Anita::NFREQ + k)/sqrt(2)/(anita1->TIMESTEP*1.E6);
	tmp_vhz[k] = GetVmmhz_freq(jpt*Anita::NFREQ + k)/sqrt(2)/(TIMESTEP*1.E6);

      } // end looping over frequencies.

      d->MakeArrayforFFT(tmp_vhz,tmp_volts, 90., true);

      // and now it is really in units of V
      FTPair::realft(tmp_volts,-1,Anita::NFOUR/2);

      // put it in normal time ording -T to T
      // instead of 0 to T, -T to 0
      Tools::NormalTimeOrdering(Anita::NFOUR/2,tmp_volts);

      int numBinShift = int(this->GetDelay(jpt) / TIMESTEP);
      if(fabs(numBinShift) >= Anita::HALFNFOUR){
	//cout<<"skipping"<<"\n";
	//don't bother adding it to the total since it's shifted out of range
      }
      else{
	if( GetDelay(jpt)>0 ){
	  Tools::ShiftLeft(tmp_volts, Anita::NFOUR/2, numBinShift );
	}
	else if( GetDelay(jpt)<0 ){
	  Tools::ShiftRight(tmp_volts, Anita::NFOUR/2, -1*numBinShift );
	}
      }

      AskaryanSignal s;
      //@todo verify what is a complete guess at the geometry
      for(int i=0; i < Anita::NFOUR/2; i++){
	s.waveform.SetPoint(i, i*TIMESTEP, tmp_volts[i]);
      }
      s.polarization = GetPol(jpt);
      s.poynting = d->getPositionRX(rx) - GetImpactPt(jpt);

      d->addSignalToRX(s, rx);
    }
  }
}
