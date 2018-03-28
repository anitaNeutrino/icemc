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
#include "Settings.h"
#include "vector.hh"

#ifdef USE_HEALPIX
#include "healpix_base.h"
#endif

class TF2;

namespace icemc{
  class Ray;
  // class Vector;
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
    std::string roughnsims_str;
    std::string roughmaterial_str;
    double NINDEX;
  public:

    Roughness(const Settings *settings1);

    void SetRoughScale(double a);

    std::string incAngle_asString(double T0);

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
    void InterpolatePowerValue(double &tcoeff_perp_polperp, double &tcoeff_parl_polperp, double &tcoeff_perp_polparl, double &tcoeff_parl_polparl, double T0, double T, double A);
#endif

  };
}

#endif
