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
  std::vector<double> fVmmhz0;     // container for vmmhz[0]
  std::vector<double> fViewangle;
  std::vector<double> fDelays;     //container for the relative propagation phase delays; final size will be anita::NFREQ after the push_backs
  std::vector<Vector> fVec2blns;  //container of 'vector to balloon'
  std::vector<Vector> fPols;      //container of transmitted polarizations
  std::vector<Position> fImpactPt;//container of ground impact points

  

public:
  Screen(int a);

  void SetNsamples(int i);

  int GetNsamples();

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

  void AddVmmhz_freq(double A);

  double GetVmmhz_freq(int i);

  void AddVmmhz0(double A);

  double GetVmmhz0(int i);

  void AddViewangle(double A);

  double GetViewangle(int i);

  void AddDelay(double A);

  double GetDelay(int i);

  void SetNvalidPoints(int i);

  double GetNvalidPoints();

  void AddVec2bln(Vector v);

  Vector GetVec2bln(int i);

  void AddPol(Vector v);

  Vector GetPol(int i);

  void AddImpactPt(Position p);

  Position GetImpactPt(int i);

  Vector CalculateTransmittedPolarization(const Vector &nnu, Vector vec_specularnormal, Vector vec_localnormal, Vector vec_pos_current_to_balloon, Vector vec_nnu_to_impactPoint, Vector npol_local_inc);

};
#endif