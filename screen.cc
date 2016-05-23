#include <math.h>
#include <iostream>

#include "vector.hh"
#include "screen.hh"

Screen::Screen(int a){
  if( !(a%2) )
    a++;                    // force 'a' odd to set screen steps correctly

  std::cerr << "Generating default screen with " << a << " sampling points" << std::endl;
  fedgeLength=1.;
  fcentralPoint = Vector(1.,1.,1.);
  fnormal = Vector(1.,1.,1.);
  funit_x = Vector(1.,1.,1.);
  funit_y = Vector(1.,1.,1.);

  fNsamples = a;
  fpositionindex = 0;
};


void Screen::SetEdgeLength(double a){
  fedgeLength = a;
};


void Screen::SetCentralPoint(Vector a){
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

void Screen::SetCornerPosition(){
  fcornerPosition = fcentralPoint - fedgeLength/2. * (funit_x + funit_y);
};


double Screen::GetEdgeLength(){
  return fedgeLength;
};


Vector Screen::GetCentralPoint(){
  return fcentralPoint;
};


Vector Screen::GetNormal(){
  return fnormal;
};


Vector Screen::GetCornerPosition(){
  return fcornerPosition;
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


Vector Screen::GetNextPosition(){
  Vector pos;

  int yindex = floor(fpositionindex / fNsamples);
  int xindex = fpositionindex % fNsamples;

  pos = fcornerPosition + xindex*funit_x + yindex*funit_y + fedgeLength/4. * (funit_x + funit_y);
        // base         //delta due to moving along screen     // shift back relative to central point
  
  fpositionindex++;
  return pos;
};






