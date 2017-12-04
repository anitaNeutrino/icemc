#include <math.h>
#include <iostream>

#include "vector.hh"
#include "position.hh"
#include "screen.hh"

Screen::Screen(int a){
  if( !(a%2) )
    a++;                    // force 'a' odd to set screen steps correctly

  std::cerr << "Generating default screen" << std::endl;
  fedgeLength=1.;
  fcentralPoint = Position(1.,1.,1.);
  fnormal = Vector(1.,1.,1.);
  funit_x = Vector(1.,1.,1.);
  funit_y = Vector(1.,1.,1.);

  fNsamples = a;
  fpositionindex = 0;
};


void Screen::SetNsamples(int i){
  fNsamples = i;
};


int Screen::GetNsamples(){
  return fNsamples;
};


void Screen::SetEdgeLength(double a){
  fedgeLength = a;
};


void Screen::SetCentralPoint(Position a){
  fcentralPoint = a;
};

void Screen::SetCosineProjectionFactor(double a){
  fcosineProjectionFactor = a;
};


double Screen::GetCosineProjectionFactor(){
  return fcosineProjectionFactor;
};

void Screen::SetNormal(Vector a){
  fnormal = a.Unit();
};

void Screen::SetUnitX(Vector a){
  funit_x = a;
};

void Screen::SetUnitY(Vector a){
  funit_y = a;
};

double Screen::GetEdgeLength(){
  return fedgeLength;
};


Position Screen::GetCentralPoint(){
  return fcentralPoint;
};


Vector Screen::GetNormal(){
  return fnormal;
};


Vector Screen::GetUnitX(){
  return funit_x;
};


Vector Screen::GetUnitY(){
  return funit_y;
};


void Screen::ResetPositionIndex(){
  fpositionindex = 0;
};


double Screen::CalcXindex(int i){
  return (double) (i % (fNsamples));
};


double Screen::CalcYindex(int i){
  return (double) floor(i / (fNsamples));
};


Position Screen::GetNextPosition(int i){
  Position pos;

  double yindex = CalcYindex(i);
  double xindex = CalcXindex(i);

  // this picks points that are NOT on the edge
  pos = fcentralPoint                                       // base
        - 0.5*fedgeLength*funit_x - 0.5*fedgeLength*fcosineProjectionFactor*funit_y   // shift to a corner
        + (1./(2.*(double)(fNsamples)))*fedgeLength*(funit_x + fcosineProjectionFactor*funit_y)  // move off the edge
        + (xindex/((double)(fNsamples)))*fedgeLength*funit_x   // move by x-increment
        + (yindex/((double)(fNsamples)))*fedgeLength*fcosineProjectionFactor*funit_y;  // move by y-increment with the cosine projection correction

  fpositionindex++;
  //std::cerr<<fpositionindex<<"  "<<yindex<<"  "<<xindex<<std::endl;
  return pos;
};


void Screen::AddVmmhz_freq(double A){
  //here i refers to the 'screen' loop
  fVmmhz_freq.push_back(A);
};


double Screen::GetVmmhz_freq(int i){
  return fVmmhz_freq[i];
};


void Screen::AddVmmhz0(double A){
  fVmmhz0.push_back(A);
};

double Screen::GetVmmhz0(int i){
  return fVmmhz0[i];
};

void Screen::AddViewangle(double A){
  fViewangle.push_back(A);
};

double Screen::GetViewangle(int i){
  return fViewangle[i];
};


void Screen::AddDelay(double A){
  fDelays.push_back(A);
};


double Screen::GetDelay(int i){
  return fDelays[i];
};


void Screen::SetNvalidPoints(int i){
  fNvalidpoints = i;
};


double Screen::GetNvalidPoints(){
  return fNvalidpoints;
};


void Screen::AddVec2bln(Vector v){
  fVec2blns.push_back(v);
};


Vector Screen::GetVec2bln(int i){
  return fVec2blns[i];
};


void Screen::AddPol(Vector v){
  fPols.push_back(v);
};


Vector Screen::GetPol(int i){
  return fPols[i];
};


void Screen::AddImpactPt(Position p){
  fImpactPt.push_back(p);
};


Position Screen::GetImpactPt(int i){
  return fImpactPt[i];
};


void Screen::AddWeight(double a){
  fWeight.push_back(a);
};

double Screen::GetWeight(int i){
  return fWeight[i];
};

void Screen::SetWeightNorm(double a){
  fWeightNorm = a;
};

double Screen::GetWeightNorm(){
  return fWeightNorm;
};

void Screen::AddIncidenceAngle(double A){
  fIncAngles.push_back(A);
};

double Screen::GetIncidenceAngle(int i){
  return fIncAngles[i];
};

void Screen::AddTransmissionAngle(double A){
  fTransAngles.push_back(A);
};

double Screen::GetTransmissionAngle(int i){
  return fTransAngles[i];
};

void Screen::ResetParameters(){
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
};

