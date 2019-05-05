#include "vector.hh"
#include "Constants.h"
#include <cmath>
#include "TVector3.h"
#include "TRotation.h"

using std::cout;

Vector::Vector(double x_inp,double y_inp,double z_inp) {
	x = x_inp;
	y = y_inp;
	z = z_inp;
  angles_need_updating = true; 
} //Constructor Vector(double,double,double)
Vector::Vector(double* xarray) {
	x=xarray[0];
	y=xarray[1];
	z=xarray[2];
  angles_need_updating = true; 
}
Vector::Vector(double theta_inp, double phi_inp) {
	
	if (theta_inp < 0.0 || theta_inp > PI) {
		
		cout<<"Error!  Attempt to construct Vector from invalid theta!\n";
		
		x = 0.;
		y = 0.;
		z = 1.;
		UpdateThetaPhi();
	} //check to see if theta is valid
	
	x = sin(theta_inp) * cos(phi_inp);
	y = sin(theta_inp) * sin(phi_inp);
	z = cos(theta_inp);
	
  theta = theta_inp; 
  phi = phi_inp; 
  angles_need_updating = false;
} //Constructor Vector(theta,phi)

Vector::Vector() {
	x = 0.;
	y = 0.;
	z = 1.;
  theta = 0; 
  phi = 0; 
  angles_need_updating = false;
} //Vector default constructor

Vector Vector::Cross(const Vector &vec) const {
	return Vector(y * vec.z - z * vec.y,
				  -x * vec.z + z * vec.x,
				  x * vec.y - y * vec.x);
} //Vector::Cross

double Vector::Dot(const Vector &vec) const
{return x * vec.x +y * vec.y + z * vec.z ;}
//Takes the dot product  this x vec.

double Vector::Mag() const {
	return sqrt(x*x + y*y + z*z);
} //Vector::Mag

double Vector::Angle(const Vector &vec) const {
	return acos(((*this)*vec) / (this->Mag() * vec.Mag()));
} //Vector::Angle

Vector Vector::ChangeCoord(const Vector &new_x_axis,const Vector &new_y_axis) const {
	
	
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
} //Vector::ChangeCoord

Vector Vector::ChangeCoord(const Vector &new_z_axis) const {
	
	Vector new_vector = this->RotateY(new_z_axis.Theta());
	new_vector = new_vector.RotateZ(new_z_axis.Phi());
	
	return new_vector;
	
}

Vector Vector::Unit() const {
	return (*this) / this->Mag();
} //Vector::Unit

void Vector::Print() const {
	cout << x << " " << y << " " << z << "\n";
}
double Vector::Theta() const {
        if (angles_need_updating) UpdateThetaPhi(); 
	return theta;
} //Vector::Theta

double Vector::Phi() const {
  if (angles_need_updating) UpdateThetaPhi(); 
  return phi;
} //Vector::Phi

Vector Vector::RotateX(double angle) const {
	double new_x = x;
	double new_y = cos(angle)*y - sin(angle)*z;
	double new_z = sin(angle)*y + cos(angle)*z;
	Vector rotated_vector(new_x,new_y,new_z);
	return rotated_vector;
} //RotateX

Vector Vector::RotateY(double angle) const {
	double new_x = cos(angle)*x + sin(angle)*z;
	double new_y = y;
	double new_z = -sin(angle)*x + cos(angle)*z;
	Vector rotated_vector(new_x,new_y,new_z);
	return rotated_vector;
} //RotateY

Vector Vector::RotateZ(double angle) const {
	// double new_x = cos(angle)*x - sin(angle)*y;
	// double new_y = sin(angle)*x + cos(angle)*y;
	// double new_z = z;
        double cosangle = cos(angle);
        double sinangle = sin(angle);
   	Vector rotated_vector(cosangle*x - sinangle*y,sinangle*x + cosangle*y,z);
	return rotated_vector;
} //RotateZ

Vector Vector::Rotate(const double angle,const Vector& axis) const {
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
} //Vector::Rotate

void Vector::UpdateThetaPhi() const {
	//This is a private method that will calculate values of theta and phi from x,y,z coordinates,
	//and store the results in the class variables.
	//double transverse = hypot(x,y);
	double transverse = sqrt(x*x+y*y);
	
	// atan2 outputs in the range -pi to pi
	theta = atan2(transverse,z);
	
	phi=atan2(y,x);
	
	if (phi<0)
		phi+=2*PI;
	// phi is now from 0 to 2*pi wrt +x
        angles_need_updating = false; 
	
} //UpdateThetaPhi
Vector Vector::Zero() {
	//Zero the vector
	
	x=0;
	y=0;
	z=0;
	return Vector(x,y,z);
} // Zero







ostream& operator <<(ostream& outs, const Vector& vec) {
	cout<<vec.x<<","<<vec.y<<","<<vec.z;
	return outs;
} //Vector overloaded operator <<


Vector Vector::Orthogonal() const 
{
  //this is based on ROOT's implementation which might be based on Geant4's? 

     Double_t xx = fabs(x); 
     Double_t yy = fabs(y); 
     Double_t zz = fabs(z); 
     if (xx < yy) {
        return xx < zz ? Vector(0,z,-y) : Vector(y,-x,0);
     } else {
        return yy < zz ? Vector(-z,0,x) : Vector(y,-x,0);
     }
}
