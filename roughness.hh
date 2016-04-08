#ifndef ROUGHNESS_H_
#define ROUGHNESS_H_

class Signal;
class TH2D;
class TF2;
class Ray;
class Vector;
class Settings;
class IceModel;
class Vector;
class Balloon;
class Interaction;

////////////////////////////////////////////////////////////////////////////////////////////////
//class Roughness:
////////////////////////////////////////////////////////////////////////////////////////////////
// #include "vector.hh"

//!  Some roughness work that was attempted
class Roughness {

private:

  TF2 *fpokey2;
  TF2 *fslappy2;

public:
  // need to find balloon location in coordinates where 
  // the surface normal is the +z axis
  // the x-z plane is defined by the neutrino trajectory
  // these are set in getballoonlocation
  double balloonphi; 
  double balloontheta;
  double balloondist;

  double rough_sigma;
  Roughness();
  void GetExitforRoughness(Settings *settings1,IceModel *antarctica,double emfrac,double hadfrac,double deltheta_em_max,double deltheta_had_max,Ray *ray1,Vector &nnu,Vector &r_bn,Vector &posnu);
  
  double GetPokey(double incident_angle,double transmitted_angle, double emfrac,double hadfrac,double deltheta_em,double deltheta_had);
  
  double GetSlappy(double incident_angle,double transmitted_angle, double emfrac,double hadfrac,double deltheta_em,double deltheta_had);
  
  double GetCombinedDeltheta(double emfrac,double hadrac,double deltheta_em_max,double deltheta_had_max);

  void GetBalloonLocation(Interaction *interaction1,Ray *ray1,Balloon *bn1,IceModel *antarctica);

  Vector nrf_iceside_specular; // ice side "specular" ray - for now, this is just a line from the interaction to the specular exit point, since we are still working on the formula to trace the ray back to the interaction through the firn.  For interactions within the firn, this will be equal to nrf_iceside[3]

};
//double rough_sigma_0=0.002;
#endif
