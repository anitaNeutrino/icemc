#include <cmath>
#include <algorithm>
#include <cstdio>
#include <fstream>

#include "vector.hh"
#include "position.hh"
#include "roughness.hh"
#include "Constants.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "signal.hh"
#include "Settings.h"
#include "Primaries.h"
#include "anita.hh"
#include "balloon.hh"
#include "earthmodel.hh"
#include "icemodel.hh"


#include "ray.hh"
//#include "healpix_base.h"

Roughness::Roughness(){

  rough_dir_str=std::getenv("ICEMC_SRC_DIR");

};



void Roughness::InterpolatePowerValue(double &tcoeff_perp, double &tcoeff_parl, double T0, double T, double A){
  // tcoeff_* are the transmission coefficients which will be copied by address for use in icemc.cc's main function
  // T0 [degrees] is the incident angle of the ray-in-ice with respect the the surface normal pointed into the air
  // T [degrees] is the exiting angle from the surface towards the balloon, measured with respect to the surface normal pointing into the air
  // A [degrees] is the azimuthal ground angle towards the balloon, measured with respect to the incidence plane (negative azimuth lies to the right, if facing parallel to the direction of incidence)

  // in the new implementation, we will use HEALPix look-up tables and average the values for the incidence angle tables that bracket the T0 value (floor and ceil)

  int thispixel;

  int pixel;
  double junk, Tparl, Tperp;

  double Tparl_down, Tperp_down, Tparl_up, Tperp_up;

  Tparl_down = Tperp_down = Tparl_up = Tperp_up = 0.;

  char base_rough_file_down[100];
  char base_rough_file_up[100];
  std::string base_rough_file_str="";
  ifstream ifs;
  std::string full_rough_file;

  // here, calculate the healpix pixel corresponding to (T, A)



  // "lower" table filename
  std::sprintf(base_rough_file_down, "/data/roughness_tables/combined_inc%dp0_nsims10000_hp2048_beckmann.hpx", floor(T0));
  base_rough_file_str = base_rough_file_down;
  full_rough_file = rough_dir_str + base_rough_file_str;
  // open and read table
  ifs.open (full_rough_file, std::ifstream::in);
  while (ifs.good()) {
    ifs >> pixel >> junk >> junk >> Tparl >> Tperp;
    if (pixel == thispixel){
      Tparl_down = Tparl;
      Tperp_down = Tperp_down;
      break;
    }
  }
  ifs.close();

  // "upper" table filename
  base_rough_file_str="";
  std::sprintf(base_rough_file_up, "/data/roughness_tables/combined_inc%dp0_nsims10000_hp2048_beckmann.hpx", ceil(T0));
  base_rough_file_str = base_rough_file_up;
  full_rough_file = rough_dir_str + base_rough_file_str;
  // open and read table
  ifs.open (full_rough_file, std::ifstream::in);
  while (ifs.good()) {
    ifs >> pixel >> junk >> junk >> Tparl >> Tperp;
    if (pixel == thispixel){
      Tparl_up = Tparl;
      Tperp_up = Tperp_down;
      break;
    }
  }
  ifs.close();

  tcoeff_perp = (Tperp_down + Tperp_up)/2.;
  tcoeff_parl = (Tparl_down + Tparl_up)/2.;

};




