#ifndef ROUGHNESS_H_
#define ROUGHNESS_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdexcept>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <vector>

#include "Constants.h"



class TF2;
class Ray;
class Vector;
class Settings;
class IceModel;
class Balloon;
class Interaction;
class Screen;

class Roughness {

private:

  std::string rough_dir_str;


public:

  Roughness();

  //! Interpolates the power value for the specified angles 
  /**
  * @param tcoeff_perp - perpendicular transmission coefficient
  * @param tcoeff_parl - parallel transmission coefficient
  * @param T0 - incident angle [degrees]
  * @param T - transmitted polar angle [degrees]
  * @param A - transmitted azimuthal angle [degrees]
  * @return double  - fractional transmitted power
  */
  void InterpolatePowerValue(double &tcoeff_perp, double &tcoeff_parl, double T0, double T, double A);

  //! Calculates the angular distance between two points on the sphere
  /**
  * @param long1 - first longitude [rad]
  * @param lat1 - first latitude [rad]
  * @param long2 - second longitude [rad]
  * @param lat2 - second latitude [rad]
  * @return double  - arclength between points [radians]
  */
  double greatCircleDist(double long1, double lat1, double long2, double lat2);

};

#endif
