#include "vector.hh"
#include "TRandom3.h"
#include "Settings.h"
#include "position.hh"
#include "Earth.h"
#include <cmath>
#include "Constants.h"

icemc::Position::Position() : Vector() 
{
  //This method intentionally left blank.
} //Position default constructor

icemc::Position::Position(const Vector& vec) : Vector(vec[0],vec[1],vec[2]) 
{
  //This method intentionally left blank.
} //Position constructor from Vector

icemc::Position::Position(double longitude, double latitude, double altitude) {
  Vector location = z_axis;
  double theta = latitude * constants::RADDEG;
  double phi = Earth::LongtoPhi_0isPrimeMeridian(longitude); // convert longitude (-180 to 180) to phi (0 to 2pi wrt 90E, counter-clockwise)
  location = location.RotateY(theta);
  location = location.RotateZ(phi);
  location = location*altitude;
  SetXYZ(location[0], location[1], location[2]);
} //Position constructor from longitude and latitude

icemc::Position::Position(double theta, double phi)
  : Vector(theta,phi)
{
}

double icemc::Position::Distance(const Position &second) const {
  double dx = GetX()-second.GetX();
  double dy = GetY()-second.GetY();
  double dz = GetZ()-second.GetZ();  
  return sqrt(dx*dx + dy*dy + dz*dz);
}

double icemc::Position::SurfaceDistance(const Position &second, double local_surface) const {
  return  this->Angle(second) * local_surface;
} //icemc::Position::SurfaceDistance

double icemc::Position::Lat() const {
  return Theta()*constants::DEGRAD;
} //icemc::Position::Lat

double icemc::Position::Lon() const {
  double phi_deg = Phi()*constants::DEGRAD;

  if (phi_deg > 270.){
    phi_deg = phi_deg-360.;
  }  
  return (270. - phi_deg);
}
