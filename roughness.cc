#include <cmath>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <math.h>
#include <iostream>
#include <string>

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
#ifdef USE_HEALPIX
#include "healpix_base.h"
#include "pointing.h"
#endif

Roughness::Roughness(){

  rough_dir_str = std::getenv("ICEMC_SRC_DIR");
#ifdef USE_HEALPIX
  order = 6;
  H = Healpix_Base(order, RING);  // RING is an enum, and is the default used to make the maps with the microfacet python simulation
#endif

};


#ifdef USE_HEALPIX
void Roughness::InterpolatePowerValue(double &tcoeff_perp, double &tcoeff_parl, double T0, double T, double A){
  // tcoeff_* are the transmission coefficients which will be copied by address for use in icemc.cc's main function
  // T0 [degrees] is the incident angle of the ray-in-ice with respect to the surface normal
  // T [degrees] is the exiting angle from the surface towards the balloon, measured with respect to the surface normal
  // A [degrees] is the azimuthal ground angle towards the balloon, measured with respect to the incidence plane (negative azimuth lies to the right, if facing parallel to the direction of incidence)

  // in the new implementation, we will use HEALPix look-up tables and average the values for the incidence angle tables that bracket the T0 value (floor and ceil)

  // determine which pixel corresponds to (T, A)
  pointing ptg = pointing(T*PI/180., A*PI/180.);
  int thispixel = H.ang2pix( ptg );

  int pixel;
  double ptheta, ptheta_g, pphi, Tparl, Tperp;  // temporary values while reading file

  double Tparl_down, Tperp_down, Tparl_up, Tperp_up;
  Tparl_down = Tperp_down = Tparl_up = Tperp_up = 0.; //default to zeros in case entry isn't present in file

  std::string header;
  std::string base_rough_file_str="";
  ifstream ifs;
  std::string full_rough_file;

  // "lower" table: read through table looking for specific pixel
  // open and read table, discard header
  base_rough_file_str = "/data/roughness_tables/combined_inc"+(std::string)Form("%i",int(floor(T0)))+"p0_nsims10000_hp"+Form("%i",H.Nside())+"_beckmann.hpx";
  full_rough_file = rough_dir_str + base_rough_file_str;
  ifs.open (full_rough_file, std::ifstream::in);
  std::getline(ifs, header);
  while (ifs.good()) {
    ifs >> pixel >> pphi >> ptheta >> Tparl >> Tperp;
    if(pixel == thispixel){
      Tparl_down = Tparl;
      Tperp_down = Tperp;
    }
  }
  ifs.close();
  // "upper" table filename: same procedure
  // open and read table, discard header
  base_rough_file_str = base_rough_file_str = "/data/roughness_tables/combined_inc"+(std::string)Form("%i",int(ceil(T0)))+"p0_nsims10000_hp"+Form("%i",H.Nside())+"_beckmann.hpx";;
  full_rough_file = rough_dir_str + base_rough_file_str;
  ifs.open (full_rough_file, std::ifstream::in);
  std::getline(ifs, header);
  while (ifs.good()) {
    ifs >> pixel >> pphi >> ptheta >> Tparl >> Tperp;
    if(pixel == thispixel){
      Tparl_up = Tparl;
      Tperp_up = Tperp;
    }
  }
  ifs.close();

  /*std::cerr<<"Inter[: "<<T0<<"  "
  <<T<<"  "
  <<A<<"  "
  <<thispixel<<"  "
  <<Tparl_down<<"  "
  <<Tparl_up<<"  "
  <<Tperp_down<<"  "
  <<Tperp_up<<std::endl;*/

  // now, average
  tcoeff_perp = (Tperp_down + Tperp_up)/2.;
  tcoeff_parl = (Tparl_down + Tparl_up)/2.;
};
#endif



