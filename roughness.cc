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



double Roughness::greatCircleDist(double long1, double lat1, double long2, double lat2){
  // based on 'geographical' longitude and latitude (0 lat is the equator, poles +-90), i.e. NOT Healpix-coords
  return acos( cos(lat1)*cos(lat2)*cos(long1-long2) + sin(lat1)*sin(lat2) );
}



void Roughness::InterpolatePowerValue(double &tcoeff_perp, double &tcoeff_parl, double T0, double T, double A){
  // tcoeff_* are the transmission coefficients which will be copied by address for use in icemc.cc's main function
  // T0 [degrees] is the incident angle of the ray-in-ice with respect to the surface normal
  // T [degrees] is the exiting angle from the surface towards the balloon, measured with respect to the surface normal
  // A [degrees] is the azimuthal ground angle towards the balloon, measured with respect to the incidence plane (negative azimuth lies to the right, if facing parallel to the direction of incidence)

  // in the new implementation, we will use HEALPix look-up tables and average the values for the incidence angle tables that bracket the T0 value (floor and ceil)

  double T_g = 90. - T;  //convert to geographical latitude
  // azimuth A is ok

  int pixel, thispixel;
  double ptheta, ptheta_g, pphi, Tparl, Tperp;  // temporary values while reading file
  double arcdist, arcdist_min;

  double Tparl_down, Tperp_down, Tparl_up, Tperp_up;

  Tparl_down = Tperp_down = Tparl_up = Tperp_up = 0.; //default to zero in case entry isn't present in file

  char base_rough_file_down[100];
  char base_rough_file_up[100];
  std::string base_rough_file_str="";
  ifstream ifs;
  std::string full_rough_file;

  // "lower" table: read through table and find 'nearest pixel center'
  std::sprintf(base_rough_file_down, "/data/roughness_tables/combined_inc%dp0_nsims10000_hp2048_beckmann.hpx", floor(T0));
  base_rough_file_str = base_rough_file_down;
  full_rough_file = rough_dir_str + base_rough_file_str;
  // open and read table
  ifs.open (full_rough_file, std::ifstream::in);
  //
  while (ifs.good()) {
    arcdist_min = 180.;
    //
    ifs >> pixel >> ptheta >> pphi >> Tparl >> Tperp;
    ptheta_g = 90. - ptheta;
    //
    arcdist = greatCircleDist(ptheta_g*PI/180., pphi*PI/180., T_g*PI/180., A*PI/180.);
    if(arcdist < arcdist_min){
      arcdist_min = arcdist;
      Tparl_down = Tparl;
      Tperp_down = Tperp;
    }
  }
  ifs.close();
  //
  // check if 'nearest' pixel was actually within ~ 1 pixel radius of the point,
  // discard coefficients if too far
  if(arcdist_min > 0.0125*PI/180.){
    Tparl_down = 0.;
    Tperp_down = 0.;
  }

  // "upper" table filename: same procedure
  base_rough_file_str="";
  std::sprintf(base_rough_file_up, "/data/roughness_tables/combined_inc%dp0_nsims10000_hp2048_beckmann.hpx", ceil(T0));
  base_rough_file_str = base_rough_file_up;
  full_rough_file = rough_dir_str + base_rough_file_str;
  // open and read table
  ifs.open (full_rough_file, std::ifstream::in);
  while (ifs.good()) {
    arcdist_min = 180.;
    //
    ifs >> pixel >> ptheta >> pphi >> Tparl >> Tperp;
    ptheta_g = 90. - ptheta;
    //
    arcdist = greatCircleDist(ptheta_g*PI/180., pphi*PI/180., T_g*PI/180., A*PI/180.);
    if(arcdist < arcdist_min){
      arcdist_min = arcdist;
      Tparl_up = Tparl;
      Tperp_up = Tperp;
    }
  }
  ifs.close();
  //
  if(arcdist_min > 0.0125*PI/180.){
    Tparl_up = 0.;
    Tperp_up = 0.;
  }

  // now, average
  tcoeff_perp = (Tperp_down + Tperp_up)/2.;
  tcoeff_parl = (Tparl_down + Tparl_up)/2.;
};




