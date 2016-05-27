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


void Screen::SetEdgeLength(double a){
  fedgeLength = a;
};


void Screen::SetCentralPoint(Position a){
  fcentralPoint = a;
};

void Screen::SetNormal(Vector a){
  fnormal = a.Unit();
  if( fabs(fnormal.GetX())>1.){
    funit_x = Vector( -1.*(fnormal.GetY()/fnormal.GetX()), 1., 0. );
    funit_x = funit_x.Unit();
    funit_y = fnormal.Cross(funit_x);
    funit_y = funit_y.Unit();
  }
  else{
    funit_x = Vector( 1., -1.*(fnormal.GetX()/fnormal.GetY()), 0. );
    funit_x = funit_x.Unit();
    funit_y = fnormal.Cross(funit_x);
    funit_y = funit_y.Unit();
  }
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


Position Screen::GetNextPosition(){
  Position pos;

  float yindex = (float) floor(fpositionindex / fNsamples);
  float xindex = (float) (fpositionindex % fNsamples);

  pos = fcentralPoint                                       // base
        - (fedgeLength/2.)*(funit_x + funit_y)              // shift to a corner
        + (xindex/((float)fNsamples))*fedgeLength*funit_x   // move by x-increment
        + (yindex/((float)fNsamples))*fedgeLength*funit_y;  // move by y-increment

  fpositionindex++;
  return pos;
};






