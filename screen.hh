#ifndef SCREEN_H_
#define SCREEN_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdexcept>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <vector>

#include "vector.hh"
#include "position.hh"

class Vector;
class Position;

class Screen {
private:
  double fedgeLength;        // the full length of one side

  Position fcentralPoint;      // coordinates of screen center

  double fcosineProjectionFactor;    // cosine projection factor of the screen onto the ground,
                                    // corrects for the long extension so sampling is faster
                                    // = cos(angle between local normal at RF exit and vector to balloon)
  
  Vector fnormal;            // screen orientation, '+' = pointing back to balloon
  Vector funit_x;            // X unit vector in screen (parallel to ground surface, perp. to screen normal)
  Vector funit_y;            // Y unit vector in screen (~ perp. to ground surface, perp. to screen normal)

  int fNsamples;          // number of samples in X-direction (and Y-, assuming symmetry)
  int fpositionindex;

  int fNvalidpoints;

  std::vector<double> fVmmhz_freq; //container for the valid screen points giving the frequency dependence magnitude for each point; every anita::NFREQ will be each screen point; final size will be (anita::NFREQ * fNsamples)
  std::vector<double> fDelays;     //container for the relative propagation phase delays; final size will be anita::NFREQ after the push_backs

  

public:
  Screen(int a);

  void SetEdgeLength(double a);
  
  void SetCentralPoint(Position a);

  void SetCosineProjectionFactor(double a);

  double GetCosineProjectionFactor();

  void SetNormal(Vector a);

  void SetUnitX(Vector a);

  void SetUnitY(Vector a);

  double GetEdgeLength();
  
  Position GetCentralPoint();
  
  Vector GetNormal();

  Vector GetUnitX();

  Vector GetUnitY();

  void ResetPositionIndex();

  double CalcXindex(int i);

  double CalcYindex(int i);

  Position GetNextPosition(int i);

  void SetVmmhz_freq(double A);

  double GetVmmhz_freq(int i);

  void SetDelay(double A);

  double GetDelay(int i);

  void SetNvalidPoints(int i);

  double GetNvalidPoints();
};
#endif