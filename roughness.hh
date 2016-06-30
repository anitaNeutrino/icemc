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

#include "Constants.h"
#include "spline.h"


class TF2;
class Ray;
class Vector;
class Settings;
class IceModel;
class Balloon;
class Interaction;
class Screen;

////////////////////////////////////////////////////////////////////////////////////////////////
//class Roughness:
////////////////////////////////////////////////////////////////////////////////////////////////
// #include "vector.hh"

//!  Some roughness work that was attempted
class Roughness {

private:

  std::string file_roughness;

  int froughsetting;  // default to flat glass
  // other options 400 grit (1), 1000 grit (2), 1500 grit (3), otherwise switch to default
  int gritvalue;  // [-1, 400, 1000, 1500]
  int Ntheta0;

  double laserPower;

  std::vector<double> theta_0;
  std::vector<double> theta_0_unique;
  std::vector<double> theta;
  std::vector<double> power;

  std::vector<tk::spline*> spline_theta0;

  std::vector<double> corrfactor_thetas;      // both fresnel and loss taken directly from Table 4 ELOG #077
  std::vector<double> corrfactor_fresnel;
  std::vector<double> corrfactor_loss;

  tk::spline *spl_cf_fresnel_ptr;
  tk::spline *spl_cf_loss_ptr;

  tk::spline *spl_ag2ga_ptr;   // spline ptr for air-glass incidence angle value to glass-air value
  tk::spline *spl_ga2ag_ptr;   // spline ptr for opposite

  std::vector<double> theta_g2a;

  //parameters for 2d gaussian fits from astropy
  double amplitude;                             // in [muW]
  double x_mean, y_mean, x_stddev, y_stddev;    // in [rad]
  double gaustheta;                                 // in [rad]
  double maxmeaspower, fitfuncmax;

  void ReadDataFile(void);
  
  void ConstructSplines(void);


public:

  Roughness(int a);

  double GetLaserPower();

  double InterpolatePowerValue(double T0, double T);

  double evaluate2dGaussian(double T0, double T);

  double GetFresnelCorrectionFactor(double T0);

  double GetLossCorrectionFactor(double T0);

  double ConvertTheta0AirGlass_to_GlassAir(double T0);

  double ConvertTheta0GlassAir_to_AirGlass(double T1);
  
  void GetFresnel(const Vector &surface_normal, const Vector &air_rf, const Vector &ice_rf, Vector &pol, double efield, double emfrac, double hadfrac, double deltheta_em_max, double deltheta_had_max, double &fresnel, double &mag);

  








  //
  //  Below here are the old roughness members, need to decide what to do with them
  //
  TF2 *fpokey2;
  TF2 *fslappy2;
  void GetExitforRoughness(Settings *settings1,IceModel *antarctica,double emfrac,double hadfrac,double deltheta_em_max,double deltheta_had_max,Ray *ray1,Vector &nnu,Vector &r_bn,Vector &posnu);
  double GetPokey(double incident_angle,double transmitted_angle, double emfrac,double hadfrac,double deltheta_em,double deltheta_had);
  double GetSlappy(double incident_angle,double transmitted_angle, double emfrac,double hadfrac,double deltheta_em,double deltheta_had);
  double GetCombinedDeltheta(double emfrac,double hadrac,double deltheta_em_max,double deltheta_had_max);
    // need to find balloon location in coordinates where 
  // the surface normal is the +z axis
  // the x-z plane is defined by the neutrino trajectory
  // these are set in getballoonlocation
  double balloonphi; 
  double balloontheta;
  double balloondist;
  double rough_sigma;
  void GetBalloonLocation(Interaction *interaction1,Ray *ray1,Balloon *bn1,IceModel *antarctica);

  Vector nrf_iceside_specular; // ice side "specular" ray - for now, this is just a line from the interaction to the specular exit point, since we are still working on the formula to trace the ray back to the interaction through the firn.  For interactions within the firn, this will be equal to nrf_iceside[3]

};
//double rough_sigma_0=0.002;
#endif
