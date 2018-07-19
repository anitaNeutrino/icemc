#include <cmath>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>

#include "TVector3.h"
#include "GeoidModel.h"
#include "Roughness.h"
#include "Constants.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "Settings.h"
#include "Primaries.h"
#include "anita.hh"
#include "balloon.hh"
#include "Earth.h"
#include "Antarctica.h"
#include "spline.h"
#include "RayTracer.h"

#ifdef USE_HEALPIX
#include "healpix_base.h"
#include "pointing.h"
#endif

icemc::Roughness::Roughness(const Settings *settings1){

  rough_dir_str = std::getenv("ICEMC_SRC_DIR");
#ifdef USE_HEALPIX
  order = 6;  // order=6 corresponds to nside = 64
  H = Healpix_Base(order, RING);  // RING is an enum, and is the default used to make the maps with the microfacet python simulation
#endif

  roughscale_str = "0p10";      //just some defaults
  roughnsims_str = "1000000";

  if (settings1->FIRN){
    roughmaterial_str="firn";
    NINDEX=constants::NFIRN;
  }
  else{
    roughmaterial_str="ice";
    NINDEX=constants::NICE;
  }
};


void icemc::Roughness::SetRoughScale(double a){
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


std::string icemc::Roughness::incAngle_asString(double T0){
  char s_f[10];
  sprintf(s_f, "%.1f", T0);
  std::string s = std::string() + s_f;
  s.replace(s.find("."), 1, "p");
  return s;
};

#ifdef USE_HEALPIX

//  map from path to map to pixel to TParl, Tperp for each polarization
static std::map<std::string, std::map<int, std::pair<double,double> > >lower_cache_polperp;
//
static std::map<std::string, std::map<int, std::pair<double,double> > >lower_cache_polparl;


void icemc::Roughness::InterpolatePowerValue(double &tcoeff_perp_polperp, double &tcoeff_parl_polperp, double &tcoeff_perp_polparl, double &tcoeff_parl_polparl, double T0, double T, double A){
  // tcoeff_* are the transmission coefficients which will be copied by address for use in icemc.cc's main function
  // T0 [degrees] is the incident angle of the ray-in-ice with respect to the surface normal
  // T [degrees] is the exiting angle from the surface towards the balloon, measured with respect to the surface normal
  // A [degrees] is the azimuthal ground angle towards the balloon, measured with respect to the incidence plane (negative azimuth lies to the right, if facing parallel to the direction of incidence)

  // in the new implementation, we will use HEALPix look-up tables and average the values for the incidence angle tables that bracket the T0 value (floor and ceil)

  pointing ptg;
  int pixel, thispixel_low;
  // temporary values while reading file
  double ptheta, pphi;
  double Tparl_polperp, Tperp_polperp;
  double Tparl_polparl, Tperp_polparl;
  double T0_i;

  std::string header;
  std::string base_rough_file_str="";
  ifstream ifs;
  std::string full_rough_file_lower;
  std::string full_rough_file_upper;

  std::vector<double> X, Yperp_polparl, Yparl_polparl, Yperp_polperp, Yparl_polperp; // here X-> transmitted angle, Y-> T**2 entry

  tk::spline spl_parl_polparl, spl_perp_polparl, spl_parl_polperp, spl_perp_polperp;

  // Here is the header info from the table generation code:
  // ## pixel  azimuth   hpx_theta   P_polparl_parl   P_polparl_perp   P_polperp_parl   P_polperp_perp

  for (double i=0; i<90; i+=.1){         // loop over incidence angle tables to read entries
    base_rough_file_str = "/data/roughness_tables/"+roughmaterial_str+"/"+roughscale_str+"/combined_inc"+incAngle_asString((double)i)+"_nsims"+roughnsims_str+"_hp"+Form("%i",H.Nside())+"_beckmann.hpx";
    full_rough_file_lower = rough_dir_str + base_rough_file_str;

    T0_i = asin(NINDEX * sin(i*PI/180.))*180./PI;

    if( !isnan(T0_i)){
      ptg = pointing(T0_i*PI/180., A*PI/180.);
      thispixel_low = H.ang2pix( ptg );
      //std::cerr<<T<<"  "<<thispixel_low<< std::endl;
      //
      if (!lower_cache_polperp.count(full_rough_file_lower)){
        //std::cerr<<"Not in cache"<<std::endl;
        ifs.open (full_rough_file_lower, std::ifstream::in);
        //std::cerr<<ifs.good()<<std::endl;
        if(ifs.good())
        {
          //would be more efficient to use std::emplace, but meh
          //also, could use a std::array for inner part but that requires a b tmore logic 
          std::map<int, std::pair<double,double> > this_lower_polperp;
          std::map<int, std::pair<double,double> > this_lower_polparl;
          std::getline(ifs, header);
          while (ifs.good()) {
            ifs >> pixel >> pphi >> ptheta >> Tparl_polparl >> Tperp_polparl >> Tparl_polperp >> Tperp_polperp;
            //std::cerr<<pphi<<"  "<<ptheta<<"  "<<Tparl<<"  "<<Tperp<<std::endl;
            if (Tparl_polparl < 0) Tparl_polparl = 0; 
            if (Tperp_polparl < 0) Tperp_polparl = 0; 
            if (Tparl_polperp < 0) Tparl_polperp = 0; 
            if (Tperp_polperp < 0) Tperp_polperp = 0; 
            //std::cerr<<Tparl<<"  "<<Tperp<<std::endl;
            this_lower_polparl[pixel]=std::pair<double,double>(Tparl_polparl,Tperp_polparl);
            this_lower_polperp[pixel]=std::pair<double,double>(Tparl_polperp,Tperp_polperp);
          }
          lower_cache_polparl[full_rough_file_lower] = this_lower_polparl;
          lower_cache_polperp[full_rough_file_lower] = this_lower_polperp;
        }
        ifs.close();
      }//end !lower_cache.count()

      //std::cerr<<i<<"  "<<lower_cache[full_rough_file_lower][thispixel_low].first<<"  "<<lower_cache[full_rough_file_lower][thispixel_low].second<<std::endl;
      X.push_back(T0_i);
      Yparl_polparl.push_back( lower_cache_polparl[full_rough_file_lower][thispixel_low].first );
      Yperp_polparl.push_back( lower_cache_polparl[full_rough_file_lower][thispixel_low].second );
      Yparl_polperp.push_back( lower_cache_polperp[full_rough_file_lower][thispixel_low].first );
      Yperp_polperp.push_back( lower_cache_polperp[full_rough_file_lower][thispixel_low].second );
    }//end !isnan(T)
  }//end for i loop

  X.push_back(90.);
  Yparl_polparl.push_back( 0. );
  Yperp_polparl.push_back( 0. );
  Yparl_polperp.push_back( 0. );
  Yperp_polperp.push_back( 0. );

  //std::cerr<<X.size()<<"  "<<Yparl.size()<<"  "<<Yperp.size()<<std::endl;

  if(X.size()>0){
    spl_parl_polparl.set_points(X,Yparl_polparl);
    spl_perp_polparl.set_points(X,Yperp_polparl);
    spl_parl_polperp.set_points(X,Yparl_polperp);
    spl_perp_polperp.set_points(X,Yperp_polperp);

    tcoeff_perp_polperp = spl_perp_polperp(T);
    tcoeff_parl_polperp = spl_parl_polperp(T);
    tcoeff_perp_polparl = spl_perp_polparl(T);
    tcoeff_parl_polparl = spl_parl_polparl(T);
  }

};
#endif



