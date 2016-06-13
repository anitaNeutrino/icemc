#ifndef SCREEN_H_
#define SCREEN_H_

#include "vector.hh"
#include "position.hh"

class Vector;
class Position;

class Screen {
private:
  double fedgeLength;        // the full length of one side

  Position fcentralPoint;      // coordinates of screen center
  
  Vector fnormal;            // screen orientation, '+' = pointing back to balloon
  Vector funit_x;            // X unit vector in screen (parallel to ground surface, perp. to screen normal)
  Vector funit_y;            // Y unit vector in screen (~ perp. to ground surface, perp. to screen normal)

  int fNsamples;          // number of samples in X-direction (and Y-, assuming symmetry)
  int fpositionindex;

public:
  Screen(int a);

  void SetEdgeLength(double a);
  
  void SetCentralPoint(Position a);

  void SetNormal(Vector a);

  void SetUnitX(Vector a);

  void SetUnitY(Vector a);

  double GetEdgeLength();
  
  Position GetCentralPoint();
  
  Vector GetNormal();

  Vector GetUnitX();

  Vector GetUnitY();

  void ResetPositionIndex();

  Position GetNextPosition();
};
#endif