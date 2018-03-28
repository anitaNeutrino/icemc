#include "vector.hh"
#include "Constants.h"
#include <cmath>
#include "TVector3.h"
#include "TRotation.h"

using std::cout;

icemc::Vector::Vector(double x_inp,double y_inp,double z_inp) {
  x = x_inp;
  y = y_inp;
  z = z_inp;
  UpdateThetaPhi();
} //Constructor Vector(double,double,double)
icemc::Vector::Vector(double* xarray) {
  x=xarray[0];
  y=xarray[1];
  z=xarray[2];
}
icemc::Vector::Vector(double theta_inp, double phi_inp) {
	
  if (theta_inp < 0.0 || theta_inp > constants::PI) {
		
    cout<<"Error!  Attempt to construct Vector from invalid theta!\n";
		
    x = 0.;
    y = 0.;
    z = 1.;
    UpdateThetaPhi();
  } //check to see if theta is valid
	
  x = sin(theta_inp) * cos(phi_inp);
  y = sin(theta_inp) * sin(phi_inp);
  z = cos(theta_inp);
	
  UpdateThetaPhi();
} //Constructor Vector(theta,phi)

icemc::Vector::Vector() {
  x = 0.;
  y = 0.;
  z = 1.;
  UpdateThetaPhi();
} //Vector default constructor

icemc::Vector icemc::Vector::Cross(const Vector &vec) const {
  return Vector(y * vec.z - z * vec.y,
		-x * vec.z + z * vec.x,
		x * vec.y - y * vec.x);
} //icemc::Vector::Cross

double icemc::Vector::Dot(const Vector &vec) const
{return x * vec.x +y * vec.y + z * vec.z ;}
//Takes the dot product  this x vec.

double icemc::Vector::Mag() const {
  return sqrt(x*x + y*y + z*z);
} //icemc::Vector::Mag

double icemc::Vector::Angle(const Vector &vec) const {
  return acos(((*this).Dot(vec)) / (this->Mag() * vec.Mag()));
} //icemc::Vector::Angle

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
	
  Vector new_vector = this->RotateY(new_z_axis.theta);
  new_vector = new_vector.RotateZ(new_z_axis.phi);
	
  return new_vector;
	
}

icemc::Vector icemc::Vector::Unit() const {
  return (*this) / this->Mag();
} //icemc::Vector::Unit

void icemc::Vector::Print() const {
  cout << x << " " << y << " " << z << "\n";
}
double icemc::Vector::GetX() const {
  return x;
} //icemc::Vector::GetX

double icemc::Vector::GetY() const {
  return y;
} //icemc::Vector::GetY

double icemc::Vector::GetZ() const {
  return z;
} //icemc::Vector::GetZ

double icemc::Vector::Theta() const {
  return theta;
} //icemc::Vector::Theta

double icemc::Vector::Phi() const {
  return phi;
} //icemc::Vector::Phi

void icemc::Vector::SetX(double inp) {
  x = inp;
  UpdateThetaPhi();
} //icemc::Vector::SetX

void icemc::Vector::SetY(double inp) {
  y = inp;
  UpdateThetaPhi();
} //icemc::Vector::SetY

void icemc::Vector::SetZ(double inp) {
  z = inp;
  UpdateThetaPhi();
} //icemc::Vector::SetZ
void icemc::Vector::SetXYZ(double inpx,double inpy,double inpz) {
  x = inpx;
  y = inpy;
  z = inpz;
  UpdateThetaPhi();
} //icemc::Vector::SetXYZ

void icemc::Vector::Reset(double x_inp, double y_inp, double z_inp) {
  x = x_inp;
  y = y_inp;
  z = z_inp;
  UpdateThetaPhi();
} //icemc::Vector::Reset

icemc::Vector icemc::Vector::RotateX(double angle) const {
  double new_x = x;
  double new_y = cos(angle)*y - sin(angle)*z;
  double new_z = sin(angle)*y + cos(angle)*z;
  Vector rotated_vector(new_x,new_y,new_z);
  return rotated_vector;
} //RotateX

icemc::Vector icemc::Vector::RotateY(double angle) const {
  double new_x = cos(angle)*x + sin(angle)*z;
  double new_y = y;
  double new_z = -sin(angle)*x + cos(angle)*z;
  Vector rotated_vector(new_x,new_y,new_z);
  return rotated_vector;
} //RotateY

icemc::Vector icemc::Vector::RotateZ(double angle) const {
  // double new_x = cos(angle)*x - sin(angle)*y;
  // double new_y = sin(angle)*x + cos(angle)*y;
  // double new_z = z;
  double cosangle = cos(angle);
  double sinangle = sin(angle);
  Vector rotated_vector(cosangle*x - sinangle*y,sinangle*x + cosangle*y,z);
  return rotated_vector;
} //RotateZ

icemc::Vector icemc::Vector::Rotate(const double angle,const Vector& axis) const {
  //Code blatently stolen from Root's TRotation::Rotate method
  //Example: If you rotate the vector (0,0,1) by 90 degrees around the vector (0,1,0), the result is (1,0,0).
  double length = axis.Mag();
	
  double s = sin(angle);
  double c = cos(angle);
  double dx = axis.x / length;
  double dy = axis.y / length;
  double dz = axis.z / length;
	
  double newx = (c+(1-c)*dx*dx) * x + ((1-c)*dx*dy-s*dz) * y + ((1-c)*dx*dz+s*dy) * z;
  double newy = ((1-c)*dy*dx+s*dz) * x + (c+(1-c)*dy*dy) * y + ((1-c)*dy*dz-s*dx) * z;
  double newz = ((1-c)*dz*dx-s*dy) * x + ((1-c)*dz*dy+s*dx) * y + (c+(1-c)*dz*dz) * z;
	
  return Vector(newx,newy,newz);
} //icemc::Vector::Rotate

void icemc::Vector::UpdateThetaPhi() {
  //This is a private method that will calculate values of theta and phi from x,y,z coordinates,
  //and store the results in the class variables.
  //double transverse = hypot(x,y);
  double transverse = sqrt(x*x+y*y);
	
  // atan2 outputs in the range -pi to pi
  theta = atan2(transverse,z);
	
  phi=atan2(y,x);
	
  if (phi<0){
    phi+=2*constants::PI;
  }
  // phi is now from 0 to 2*pi wrt +x
	
} //UpdateThetaPhi

icemc::Vector icemc::Vector::Zero() {
  //Zero the vector
	
  x=0;
  y=0;
  z=0;
  return Vector(x,y,z);
} // Zero


std::ostream& operator <<(std::ostream& outs, const icemc::Vector& vec) {
  outs << vec.GetX() << "," << vec.GetY() << "," << vec.GetZ();
  return outs;
} //Vector overloaded operator <<
