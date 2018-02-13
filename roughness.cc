#include <cmath>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>

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

Roughness::Roughness(Settings *settings1){

  rough_dir_str = std::getenv("ICEMC_SRC_DIR");
#ifdef USE_HEALPIX
  order = 6;  // order=6 corresponds to nside = 64
  H = Healpix_Base(order, RING);  // RING is an enum, and is the default used to make the maps with the microfacet python simulation
#endif

  roughscale_str = "0p10";
  roughnsims_str = "10000000";

  if (settings1->FIRN){
    roughmaterial_str="firn";
    NINDEX=NFIRN;
  }
  else{
    roughmaterial_str="ice";
    NINDEX=NICE;
  }
};


void Roughness::SetRoughScale(double a){
  std::ostringstream strs;
  strs.precision(2);
  strs.setf( std::ios::fixed, std:: ios::floatfield );
  strs << a;
  roughscale_str = strs.str();
  roughscale_str[1] = 'p';

  if(roughscale_str=="0p00")
    roughnsims_str = "1";
  //std::cerr<<"Roughness srms: "<<roughscale_str<<std::endl;
};


#ifdef USE_HEALPIX


//lower cache 

//  map from path to map to pixel to TParl, Tperp 
static std::map<std::string, std::map<int, std::pair<double,double> > >lower_cache; 
static std::map<std::string, std::map<int, std::pair<double,double> > >upper_cache; 


void Roughness::InterpolatePowerValue(double &tcoeff_perp, double &tcoeff_parl, double T0, double T, double A){
  // tcoeff_* are the transmission coefficients which will be copied by address for use in icemc.cc's main function
  // T0 [degrees] is the incident angle of the ray-in-ice with respect to the surface normal
  // T [degrees] is the exiting angle from the surface towards the balloon, measured with respect to the surface normal
  // A [degrees] is the azimuthal ground angle towards the balloon, measured with respect to the incidence plane (negative azimuth lies to the right, if facing parallel to the direction of incidence)

  // in the new implementation, we will use HEALPix look-up tables and average the values for the incidence angle tables that bracket the T0 value (floor and ceil)

  pointing ptg;
  int pixel, thispixel_up, thispixel_low;
  double ptheta, pphi, Tparl, Tperp;  // temporary values while reading file

  double Tparl_down, Tperp_down, Tparl_up, Tperp_up;
  Tparl_down = Tperp_down = Tparl_up = Tperp_up = 0.; //default to zeros in case entry isn't present in file

  std::string header;
  std::string base_rough_file_str="";
  ifstream ifs;
  std::string full_rough_file_lower;
  std::string full_rough_file_upper;

  // "lower" table: read through table looking for specific pixel
  // open and read table, discard header
  base_rough_file_str = "/data/roughness_tables/"+roughmaterial_str+"/"+roughscale_str+"/combined_inc"+(std::string)Form("%i",int(floor(T0)))+"p0_nsims"+roughnsims_str+"_hp"+Form("%i",H.Nside())+"_beckmann.hpx";
  full_rough_file_lower = rough_dir_str + base_rough_file_str;
  //std::cerr<<full_rough_file_lower<<"  :  "<<lower_cache.count(full_rough_file_lower)<<std::endl;
  // determine which pixel corresponds to (T, A)
  T = asin(NINDEX * sin( floor(T0)*PI/180. ));
  //std::cerr<<"low: "<<T0<<"  "<<floor(T0)*PI/180.<<"  "<<NINDEX<<"  "<<T<<std::endl;
  if( !isnan(T)){
    ptg = pointing(T, A*PI/180.);
    thispixel_low = H.ang2pix( ptg );
    //
    if (!lower_cache.count(full_rough_file_lower)){
      //std::cerr<<"Not in cache"<<std::endl;
      ifs.open (full_rough_file_lower, std::ifstream::in);
      //std::cerr<<ifs.good()<<std::endl;
      if(ifs.good())
      {
        //would be more efficient to use std::emplace, but meh
        //also, could use a std::array for inner part but that requires a b tmore logic 
        std::map<int, std::pair<double,double> > this_lower; 
        std::getline(ifs, header);
        while (ifs.good()) {
          ifs >> pixel >> pphi >> ptheta >> Tparl >> Tperp;
          //std::cerr<<pphi<<"  "<<ptheta<<"  "<<Tparl<<"  "<<Tperp<<std::endl;
          if (Tparl < 0) Tparl = 0; 
          if (Tperp < 0) Tperp = 0; 
          this_lower[pixel]=std::pair<double,double>(Tparl,Tperp); 
        }
        lower_cache[full_rough_file_lower] = this_lower; 
      }
      ifs.close();
    }
  }

  // "upper" table filename: same procedure
  // open and read table, discard header
  base_rough_file_str = base_rough_file_str = "/data/roughness_tables/"+roughmaterial_str+"/"+roughscale_str+"/combined_inc"+(std::string)Form("%i",int(ceil(T0)))+"p0_nsims"+roughnsims_str+"_hp"+Form("%i",H.Nside())+"_beckmann.hpx";;
  full_rough_file_upper = rough_dir_str + base_rough_file_str;
  // determine which pixel corresponds to (T, A)
  T = asin(NINDEX * sin( ceil(T0)*PI/180. ));
  //std::cerr<<"low: "<<T0<<"  "<<ceil(T0)*PI/180.<<"  "<<NINDEX<<"  "<<T<<std::endl;
  if( !isnan(T)){
    ptg = pointing(T, A*PI/180.);
    thispixel_up = H.ang2pix( ptg );
    //
    if (!upper_cache.count(full_rough_file_upper)){
      ifs.open (full_rough_file_upper, std::ifstream::in);
      if(ifs.good())
      {
        std::map<int, std::pair<double,double> > this_upper; 
        std::getline(ifs, header);
        while (ifs.good()) {
          ifs >> pixel >> pphi >> ptheta >> Tparl >> Tperp;
          if (Tparl < 0) Tparl = 0; 
          if (Tperp < 0) Tperp = 0; 
          this_upper[pixel]=std::pair<double,double>(Tparl,Tperp); 
        }
        upper_cache[full_rough_file_upper] = this_upper; 
      }
      ifs.close();
    }
  }
  Tparl_down = lower_cache[full_rough_file_lower][thispixel_low].first; 
  Tperp_down = lower_cache[full_rough_file_lower][thispixel_low].second; 
  Tparl_up = upper_cache[full_rough_file_upper][thispixel_up].first; 
  Tperp_up = upper_cache[full_rough_file_upper][thispixel_up].second; 
  


  //std::cerr<<"Inter[: "<<T0<<"  "
  //<<T<<"  "
  //<<A<<"  "
  //<<thispixel<<"  "
  //<<Tparl_down<<"  "
  //<<Tparl_up<<"  "
  //<<Tperp_down<<"  "
  //<<Tperp_up<<std::endl;

  // now, average
  tcoeff_perp = (Tperp_down + Tperp_up)/2.;
  tcoeff_parl = (Tparl_down + Tparl_up)/2.;
};
#endif



