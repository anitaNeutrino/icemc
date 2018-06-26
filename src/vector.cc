#include "vector.hh"
#include "Constants.h"
#include <cmath>
#include "TVector3.h"
#include "TRotation.h"
#include "IcemcLog.h"



icemc::Vector::Vector(double theta, double phi) {
	
  if (theta < 0.0 || theta >= constants::PI) {
		
    icemcLog () << icemc::warning << __PRETTY_FUNCTION__
		<< " attempted to construct Vector with invalid theta. "
		<< "Expected 0 < theta < pi, but got theta = " << theta << "!\n";
  }  
  // if we really trust the input theta/phi we could
  // use them to set fTheta/fPhi... but for now

  double sinTheta = sin(theta);
  fX = sinTheta * cos(phi);
  fY = sinTheta * sin(phi);
  fZ = cos(theta);
}







icemc::Vector icemc::Vector::Cross(const Vector &vec) const {
  return Vector(fY * vec.fZ - fZ * vec.fY,
		-fX * vec.fZ + fZ * vec.fX,
		fX * vec.fY - fY * vec.fX);
}

double icemc::Vector::Dot(const Vector &vec) const {
  return fX*vec.fX + fY*vec.fY + fZ*vec.fZ ;
}


//Takes the dot product  this x vec.

double icemc::Vector::Mag() const {
  return sqrt(fX*fX + fY*fY + fZ*fZ);
} //icemc::Vector::Mag


double icemc::Vector::Angle(const Vector &vec) const {
  return acos(((*this).Dot(vec)) / (this->Mag() * vec.Mag()));
} //icemc::Vector::Angle


icemc::Vector icemc::Vector::Unit() const {
  return (*this) / this->Mag();
} //icemc::Vector::Unit


















/**
 * Rotation/transformation functions, 
 * None of these alter this vector, but construct new vectors 
 */


icemc::Vector icemc::Vector::ChangeCoord(const Vector &new_x_axis,const Vector &new_y_axis) const {
	
	
  TVector3 temp;
  temp.SetX(this->GetX());
  temp.SetY(this->GetY());
  temp.SetZ(this->GetZ());

  TVector3 tnew_x_axis;
  tnew_x_axis.SetX(new_x_axis.GetX());
  tnew_x_axis.SetY(new_x_axis.GetY());
  tnew_x_axis.SetZ(new_x_axis.GetZ());
	
  TVector3 tnew_y_axis;
  tnew_y_axis.SetX(new_y_axis.GetX());
  tnew_y_axis.SetY(new_y_axis.GetY());
  tnew_y_axis.SetZ(new_y_axis.GetZ());

  //cout << "before trotation.\n";	
  TRotation r;
  r.SetXAxis(tnew_x_axis);
  r.SetYAxis(tnew_y_axis);
  r.SetZAxis(tnew_x_axis.Cross(tnew_y_axis));
  //cout << "tnew_x_axis is ";
  //tnew_x_axis.Print();
  //cout << "tnew_y_axis is ";
  //tnew_y_axis.Print();
  //cout << "tnew_z_axis is ";
  //(tnew_x_axis.Cross(tnew_y_axis)).Print();

  //	cout << "after trotation.\n";	
	
  //r.SetToIdentity();
  //r.RotateAxes(newX,newY,newZ);
  //r.Invert();

  temp.Transform(r);

  Vector new_vector;
  new_vector.SetX(temp.X());
  new_vector.SetY(temp.Y());
  new_vector.SetZ(temp.Z());

  //Vector new_vector = this->RotateY(new_z_axis.theta);
  //new_vector = new_vector.RotateZ(new_z_axis.phi);

  return new_vector;
} //icemc::Vector::ChangeCoord

icemc::Vector icemc::Vector::ChangeCoord(const Vector &new_z_axis) const {
  // double theta = Theta();
  // double phi = Phi();
  Vector new_vector = this->RotateY(new_z_axis.Theta());
  new_vector = new_vector.RotateZ(new_z_axis.Phi());

  return new_vector;

}


icemc::Vector icemc::Vector::RotateX(double angle) const {

  double x = fX;
  double y = cos(angle)*fY - sin(angle)*fZ;
  double z = sin(angle)*fY + cos(angle)*fZ;
  Vector rotated_vector(x,y,z);
  return rotated_vector;
} //RotateX

icemc::Vector icemc::Vector::RotateY(double angle) const {
  double x = cos(angle)*fX + sin(angle)*fZ;
  double y = fY;
  double z = -sin(angle)*fX + cos(angle)*fZ;
  Vector rotated_vector(x,y,z);
  return rotated_vector;
} //RotateY

icemc::Vector icemc::Vector::RotateZ(double angle) const {
  // double new_x = cos(angle)*x - sin(angle)*y;
  // double new_y = sin(angle)*x + cos(angle)*y;
  // double new_z = z;
  double cosangle = cos(angle);
  double sinangle = sin(angle);
  Vector rotated_vector(cosangle*fX - sinangle*fY,sinangle*fX + cosangle*fY,fZ);
  return rotated_vector;
} //RotateZ

icemc::Vector icemc::Vector::Rotate(const double angle,const Vector& axis) const {
  //Code blatently stolen from Root's TRotation::Rotate method
  //Example: If you rotate the vector (0,0,1) by 90 degrees around the vector (0,1,0), the result is (1,0,0).
  double length = axis.Mag();
	
  double s = sin(angle);
  double c = cos(angle);
  double dx = axis.fX / length;
  double dy = axis.fY / length;
  double dz = axis.fZ / length;

  double newx = (c+(1-c)*dx*dx) * fX + ((1-c)*dx*dy-s*dz) * fY + ((1-c)*dx*dz+s*dy) * fZ;
  double newy = ((1-c)*dy*dx+s*dz) * fX + (c+(1-c)*dy*dy) * fY + ((1-c)*dy*dz-s*dx) * fZ;
  double newz = ((1-c)*dz*dx-s*dy) * fX + ((1-c)*dz*dy+s*dx) * fY + (c+(1-c)*dz*dz) * fZ;
	
  return Vector(newx,newy,newz);
} //icemc::Vector::Rotate









/**
 * Setter functions
 * @note All of these should set fAnglesDirty = true
 */


void icemc::Vector::SetX(double x) {
  if(fX!=x){
    fAnglesDirty = true;
  }
  fX = x;
  
}



void icemc::Vector::SetY(double y) {
  if(fY!=y){
    fAnglesDirty = true;
  }
  fY = y;
}



void icemc::Vector::SetZ(double z) {
  if(fZ!=z){
    fAnglesDirty = true;
  }
  fZ = z;
}


void icemc::Vector::SetXYZ(double x,double y,double z) {
  if(fX!=x || fY!=y || fZ!=z){
    fAnglesDirty = true;
  }
  fX = x;
  fY = y;
  fZ = z;
}



void icemc::Vector::updateThetaPhi() const {
  if(fAnglesDirty){
    double rho = sqrt(fX*fX+fY*fY);
    fTheta = atan2(rho,fZ);
    fPhi = atan2(fY,fX);

    if (fPhi<0){
      fPhi += 2*constants::PI;
    }
    fAnglesDirty = false;
  }
}


std::ostream& operator <<(std::ostream& outs, const icemc::Vector& vec) {
  outs << "(" << vec.GetX() << "," << vec.GetY() << "," << vec.GetZ() << ")";
  return outs;
} //Vector overloaded operator <<
