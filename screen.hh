#ifndef SCREEN_H_
#define SCREEN_H_

#include "vector.hh"

class Vector;

class Screen {
private:
  double fedgeLength;        // the full length of one side

  Vector fcentralPoint;      // coordinates of screen center
  
  Vector fnormal;            // screen orientation, '+' = pointing back to balloon
  Vector funit_x;            // X unit vector in screen (parallel to ground surface, perp. to screen normal)
  Vector funit_y;            // Y unit vector in screen (~ perp. to ground surface, perp. to screen normal)

  int fNsamples;          // number of samples in X-direction (and Y-, assuming symmetry)
  int fpositionindex;

  Vector fcornerPosition;

public:
  Screen(int a);

  void SetEdgeLength(double a);
  
  void SetCentralPoint(Vector a);

  void SetNormal(Vector a);

  void SetCornerPosition();

  double GetEdgeLength();
  
  Vector GetCentralPoint();
  
  Vector GetNormal();

  Vector GetCornerPosition();

  Vector GetUnitX();

  Vector GetUnitY();

  void ResetPositionIndex();

  Vector GetNextPosition();
};
#endif