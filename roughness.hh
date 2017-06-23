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
#include "spline.h"


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

  std::string file_roughness;

  int froughsetting;                       ///< default to 400 grit (1), 1000 grit (2), 1500 grit (3)
  int gritvalue;                           ///< [-1, 400, 1000, 1500]
  int Ntheta0;                             ///< number of UCLA incidence angles for a given grit value
  std::vector<double> theta_0;             ///< vector of UCLA theta0 values (incidence)
  std::vector<double> theta_0_unique;      ///< vector of unique UCLA theta0 values
  std::vector<double> theta;               ///< vector of UCLA theta values (transmitted)
  std::vector<double> power;               ///< vector of UCLA power measurements for each (theta0, theta) pair

  std::vector<tk::spline*> spline_theta0;  ///< spline parameterization for theta0

  std::vector<double> corrfactor_thetas;   ///< vector of angles used for Fresnel and loss correction factors, from Table 4 ELOG #077
  std::vector<double> corrfactor_fresnel;  ///< vector of Fresnel correction factors
  std::vector<double> corrfactor_loss;     ///< vector of loss correction factors

  tk::spline *spl_cf_fresnel_ptr;          ///< spline parameterization of the Fresnel correction factors
  tk::spline *spl_cf_loss_ptr;             ///< spline parameterization of the loss correction factors

  tk::spline *spl_ag2ga_ptr;               ///< spline parameterization for air-glass incidence angle value to glass-air value
  tk::spline *spl_ga2ag_ptr;               ///< spline parameterization for glass-air to air-glass

  std::vector<double> theta_g2a;           ///< vector of glass-air incidence angles

  //parameters for 2d gaussian fits from astropy
  double amplitude;                             ///< amplitude of the 2-d gaussian fit [\muW]
  double x_mean, y_mean, x_stddev, y_stddev;    ///< means and sigmas of the 2-d gaussian fit in [rad]
  double gaustheta;                             ///< rotation angle of the 2-d gaussian fit [rad]
  double maxmeaspower, fitfuncmax;

  void ReadDataFile(void);
  
  void ConstructSplines(void);


public:

  Roughness(int a);

  /// Gets the total laser power
  /// @return returns 580. [muW]
  double GetLaserPower();

  /// Interpolates the power value using a two-dimensional gaussian for the specified angles 
  /// @param T0 - incident angle [degrees]
  /// @param T - transmitted angle [degrees]
  /// @return double [muW]
  double InterpolatePowerValue(double T0, double T);

  /// Evaluate the two-dimensional gaussian for the specified angle
  /// @param T0 - incident angle [degrees]
  /// @param T - transmitted angle [degrees]
  /// @return double [muW]
  double evaluate2dGaussian(double T0, double T);

  /// Get the Fresnel correction factor
  /// @param T0 - incident angle [degrees]
  /// @return double
  double GetFresnelCorrectionFactor(double T0);

  /// Get the loss correction factor through the glass
  /// @param T0 - incident angle [degrees]
  /// @return double
  double GetLossCorrectionFactor(double T0);

  /// Calculate the corresponding incidence angle if the interface is reversed
  /// @param T0 - incident angle [degrees]
  /// @return double
  double ConvertTheta0AirGlass_to_GlassAir(double T0);

  /// Calculate the corresponding incidence angle if the interface is reversed
  /// @param T0 - incident angle [degrees]
  /// @return double
  double ConvertTheta0GlassAir_to_AirGlass(double T1);

  /// Calculate the transmitted polarization vector
  /// @param nnu - neutrino direction
  /// @param vec_specularnormal - surface normal vector at the specular RF exit point
  /// @param vec_localnormal - surface normal vector at the screen point's projected impact position
  /// @param vec_pos_current_to_balloon - vector from the projected impact position to the balloon
  /// @param vec_nnu_to_impactPoint - vector from neutrino interaction position to the projected impact point
  /// @param npol_local_inc - incident polarization vector
  /// @return transmitted polarization vector
  Vector CalculateTransmittedPolarization(const Vector &nnu, Vector vec_specularnormal, Vector vec_localnormal, Vector vec_pos_current_to_balloon, Vector vec_nnu_to_impactPoint, Vector npol_local_inc);

  /// For a given interface between n0 and n1_old, re-calculate the transmitted angle based on switching n1_old to n1_new
  /// @param n1_old - old index of refraction
  /// @param n1_new - new index of refraction
  /// @param trans_angle_old - old angle of transmission
  double AdjustTransmissionAngle(double n1_old, double n1_new, double trans_angle_old);
};

#endif
