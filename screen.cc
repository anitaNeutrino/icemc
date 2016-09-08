#include <math.h>
#include <iostream>

#include "vector.hh"
#include "position.hh"
#include "screen.hh"

Screen::Screen(int a){
  if( !(a%2) )
    a++;                    // force 'a' odd to set screen steps correctly

  std::cerr << "Generating default screen with " << a << " sampling points" << std::endl;
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


Vector Screen::CalculateTransmittedPolarization(const Vector &nnu, Vector vec_specularnormal, Vector vec_localnormal, Vector vec_pos_current_to_balloon, Vector vec_nnu_to_impactPoint, Vector npol_local_inc){
  Vector temp_a = nnu.Cross(vec_specularnormal).Unit();
  Vector temp_b = nnu.Cross(vec_localnormal).Unit();
  double mtrx_cos_inc = temp_a.Dot(temp_b) / temp_a.Mag() / temp_b.Mag();
  double mtrx_sin_inc = sqrt(1. - mtrx_cos_inc*mtrx_cos_inc);

  temp_a = vec_pos_current_to_balloon.Cross(vec_specularnormal).Unit();
  temp_b = vec_pos_current_to_balloon.Cross(vec_localnormal).Unit();
  double mtrx_cos_trans = temp_a.Dot(temp_b) / temp_a.Mag() / temp_b.Mag();
  double mtrx_sin_trans = sqrt(1. - mtrx_cos_trans*mtrx_cos_trans);

  // now define the incidence plane using int.point-imp.point vector and local surface normal
  Vector vec_incplane_normal = vec_nnu_to_impactPoint.Cross( (const Vector)vec_localnormal ).Unit();

  // now define the scattering plane using imp.point-balloon vector and local surface normal
  Vector vec_scatplane_normal = vec_pos_current_to_balloon.Cross( (const Vector)vec_localnormal ).Unit();

  // determine incident vertical and horizontal field components (vertical = in-plane)
  Vector E_local_h_inc = npol_local_inc.Dot(vec_incplane_normal) * vec_incplane_normal;
  Vector E_local_v_inc = npol_local_inc - E_local_h_inc;

  // recast the values with the transformation matrix
  double E_local_v_trans_mag = mtrx_cos_trans * (mtrx_cos_inc*E_local_v_inc.Mag() - mtrx_sin_inc*E_local_h_inc.Mag())
          - mtrx_sin_trans * (mtrx_sin_inc*E_local_v_inc.Mag() + mtrx_cos_inc*E_local_h_inc.Mag());
  double E_local_h_trans_mag = mtrx_sin_trans * (mtrx_cos_inc*E_local_v_inc.Mag() - mtrx_sin_inc*E_local_h_inc.Mag())
          + mtrx_cos_trans * (mtrx_sin_inc*E_local_v_inc.Mag() + mtrx_cos_inc*E_local_h_inc.Mag());

  // transmitted polarization needs to be perpendicular to to-balloon vector, and the horizontal component is 'set', so need to find appropriate vector for the vertical component to ensure perpendicularity
  Vector npol_local_trans = (E_local_h_trans_mag*vec_scatplane_normal
          + E_local_v_trans_mag* vec_pos_current_to_balloon.Cross( (const Vector)vec_scatplane_normal ).Unit() ).Unit();

  return npol_local_trans;
};
