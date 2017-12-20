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

#ifdef USE_HEALPIX
#include "healpix_base.h"
#endif

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
#ifdef USE_HEALPIX
  int order;
  Healpix_Base H;
#endif

  std::string roughscale_str;
public:

  Roughness();

  void SetRoughScale(double a);

#ifdef USE_HEALPIX
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
#endif

};

#endif
