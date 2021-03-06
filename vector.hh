#ifndef VECTOR_H_
#define VECTOR_H_
////////////////////////////////////////////////////////////////////////////////////////////////
//class Vector:
//This class represents a three-vector.  Operators are overloaded to provide for the
//familiar operations of vector addition, subtraction, scalar multiplication and division,
//and the dot product.  The x,y,z components can be output with the << operator onto an
//output stream.
//Methods are provided to take the vector product, find the magnitude, find a three-dimensional
//angle, and perform various rotations.
//
//The only methods that will alter an established vector are Reset and the three Set functions.
//All other methods (cross product, rotations, etc.) return a new Vector object.
////////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>

//#include <fstream>
//#include <sstream>
//#include <math.h>
//#include <string>
//#include <stdio.h>
//#include <stdlib.h>

using std::ostream;

//! This class represents a three-vector.  Operators are overloaded to provide for the familiar operations of vector addition, subtraction, scalar multiplication and division, and the dot product.
class Vector {
public:
  friend Vector operator +(const Vector& vector1, const Vector& vector2)
  {
	return Vector(vector1.x + vector2.x , vector1.y + vector2.y , vector1.z + vector2.z);
  }
  //Add two vectors component by component

  friend void operator +=(Vector& vector1, const Vector& vector2)
  {
    vector1 = vector1 + vector2;
  }

  friend Vector operator -(const Vector& vector1, const Vector& vector2)
  {
    return Vector(vector1.x - vector2.x , vector1.y - vector2.y , vector1.z - vector2.z);
  }
  //Subtract two vectors component by component

  friend void operator -=(Vector& vector1, const Vector& vector2)
  {
    vector1 = vector1 - vector2;
  }

  friend Vector operator -(const Vector& vec)
  {
    return Vector(-vec.x,-vec.y,-vec.z);
  }
  //Gives the negative of a vector

  friend double operator *(const Vector& vector1, const Vector& vector2)
  {
    return ((vector1.x * vector2.x) + (vector1.y * vector2.y) + (vector1.z * vector2.z));
  }

  //Take the dot product of two vectors

  friend Vector operator /(const Vector &v, const double& a)
  {
    return Vector(v.x / a, v.y / a, v.z / a);
  }
  //Divide a vector by a scalar

  friend Vector operator *(const double& a, const Vector& v)
  {
    return Vector(a*v.x , a*v.y , a*v.z);
  }

  friend Vector operator *(const Vector& v, const double& a)
  {
    return Vector(a*v.x , a*v.y , a*v.z);
  }
 
  //Multiply a vector by a scalar

  friend ostream& operator <<(ostream& outs, const Vector& vec);
  //Print the three rectangular components of the vector to the screen

  double operator [](int i) const
  {
    switch(i) 
    { 
      case 0: return x; 
      case 1: return y; 
      case 2: return z; 
      case 3: return phi; 
      case 4: return theta; 
      default: 
        return 0; 
    }
  }
  //Returns an element of the vector.  vector[0] is the x component,
  //vector[1] is the y component, and vector[2] is the z component.

  Vector(double x_inp,double y_inp,double z_inp);
  //Constructor: Initialize a new vector with given values of x, y, and z.

  Vector(double *xarray);
  //Constructor: Initialize a new vector with elements of xarray giving values 
  // of x,y,z.

  Vector(double theta, double phi);
  //Constructor: Initialize a new vector with unit length and in the
  //theta, phi direction.
  //theta and phi must be in RADIANS!
  //Accepts theta from 0 to PI, and any phi.

  Vector();
  //Default constructor: Initialize a unit vector in the z direction.

  Vector RotateX(double angle) const;
  //Returns the vector rotated counterclockwise (right handed coordinates) 
  //by "angle" radians about the X axis.
  //N.B. : Returns a new Vector object.  Does not change the vector it is called from.

  Vector RotateY(double angle) const;
  //Returns the vector rotated counterclockwise (right handed coordinates) 
  //by "angle" radians about the Y axis.
  //N.B. : Returns a new Vector object.  Does not change the vector it is called from.

  Vector RotateZ(double angle) const;
  //Returns the vector rotated counterclockwise (right handed coordinates) 
  //by "angle" radians about the Z axis.
  //N.B. : Returns a new Vector object.  Does not change the vector it is called from.

  Vector Cross(const Vector &vec) const;
  //Takes the cross product  this x vec.

  double Dot(const Vector &vec) const;
//Takes the dot product  this x vec.

  Vector Rotate(double angle, const Vector &axis) const;
  //Returns the vector that is this vector rotated around the vector "axis" by angle (in radians) "angle".

  Vector Zero();
  //zero the vector

  double Mag2() const { return X()*X()* + Y()*Y() + Z()*Z(); }

  double Mag() const;
  //Returns the magnitude of this vector.

  double Angle(const Vector &vec) const;
  //Returns the 3-dimensional angle between this vector and the argument.

  Vector ChangeCoord(const Vector &new_x_axis,const Vector &new_y_axis) const;
  Vector ChangeCoord(const Vector &new_z_axis) const;
  //Returns this vector, rotated to a new coordinate system with the argument as the z-axis.
  //The vector is rotated in such a way that a vector pointing in the z direction will be 
  //rotated onto the new axis.

  Vector Orthogonal() const; 
  // Returns a vector orthogonal to this one. 


  Vector Unit() const;
  //Returns a unit vector in the same direction as this vector.

  //Accessor functions
  double GetX() const { return x; }
  double GetY() const { return y; }
  double GetZ() const { return z; } 
  double X() const { return x; }
  double Y() const { return y; }
  double Z() const { return z; } 



  double Theta() const;
  double Phi() const;
  void Print() const;

  //Mutator functions
  void SetX(double inp) { x = inp; angles_need_updating = true ;}
  void SetY(double inp) { y = inp; angles_need_updating = true ;}
  void SetZ(double inp) { z = inp; angles_need_updating = true ;}
  void SetXYZ(double inpx,double inpy,double inpz)
  {
    x = inpx; 
    y = inpy; 
    z = inpz; 
    angles_need_updating = true; 
  }
  void Reset(double x_inp, double y_inp, double z_inp) { SetXYZ(x_inp, y_inp, z_inp); } 

protected:
  //Class variables
  double x; //x component of vector
  double y; //y component of vector
  double z; //z component of vector
  mutable double theta;  //theta component of vector in radians
  mutable double phi; //phi component of vector in radians
  mutable bool angles_need_updating; 

  //Class private functions
  void UpdateThetaPhi() const; 
  //This method finds theta and phi from the x,y,z Cartesian coordinates.  
  //It should be called at any time that the x,y,z components are modified,
  //so that the theta and phi components are current at all times.

}; //class Vector

////////////////////////////////////////////////////////////////////////////////////////////////
// Vector Constants

static Vector x_axis = Vector(1,0,0);
static Vector y_axis = Vector(0,1,0);
static Vector z_axis = Vector(0,0,1);
#endif
