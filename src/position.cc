#include "vector.hh"
#include "TRandom3.h"
#include "Settings.h"
#include "position.hh"
#include "earthmodel.hh"
#include <cmath>
#include "Constants.h"
 icemc::Position::Position() : Vector() 
{
  //This method intentionally left blank.
} //Position default constructor

 icemc::Position::Position(Vector vec) : Vector(vec[0],vec[1],vec[2]) 
{
  //This method intentionally left blank.
} //Position constructor from Vector

 icemc::Position::Position(double longitude, double latitude, double altitude) {
  Vector location = z_axis;
  theta = latitude * RADDEG;

  phi=EarthModel::LongtoPhi_0isPrimeMeridian(longitude); // convert longitude (-180 to 180) to phi (0 to 2pi wrt 90E, counter-clockwise)

  location = location.RotateY(theta);
  location = location.RotateZ(phi);
  location = altitude * location;

  x = location[0];
  y = location[1];
  z = location[2];
} //Position constructor from longitude and latitude

 icemc::Position::Position(double theta_inp, double phi_inp) : Vector(theta_inp,phi_inp) 
{
  //This method intentionally left blank.
} //Constructor Position(theta,phi)

 double icemc::Position::Distance(const Position &second) const {
  return sqrt((x - second.x)*(x-second.x) 
	      + (y - second.y)*(y-second.y) 
	      + (z - second.z)*(z-second.z)); // it saves time to multiply them like this rather than use pow
} //icemc::Position::Distance

 double icemc::Position::SurfaceDistance(const Position &second, double local_surface) const {
  return  this->Angle(second) * local_surface;
} //icemc::Position::SurfaceDistance

 double icemc::Position::Lat() const {
  return theta*DEGRAD;
} //icemc::Position::Lat

 double icemc::Position::Lon() const {
  double phi_deg = phi*DEGRAD;

  if (phi_deg > 270.)
    phi_deg = phi_deg-360.;
  
  return (270. - phi_deg);
} //icemc::Position::Lon
