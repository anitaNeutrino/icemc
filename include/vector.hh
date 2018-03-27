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

namespace icemc {

  //! This class represents a three-vector.  Operators are overloaded to provide for the familiar operations of vector addition, subtraction, scalar multiplication and division, and the dot product.
  class Vector {
  public:

    inline double operator[](int i) const {
      //code taken from ROOT's TVector3 class
      switch(i) {
      case 0:
	return x;
      case 1:
	return y;
      case 2:
	return z;
      default:
	printf("icemc::Vector::operator[](i) has been given a bad index: %d.  Returning zero.",i);
      } //end switch
      return 0.;
    }
    inline void operator +=(const Vector& v2) {
      x += v2.x;
      y += v2.y;
      z += v2.z;
    }
    inline void operator -=(const Vector& v2) {
      x -= v2.x;
      y -= v2.y;
      z -= v2.z;      
    }
    inline void operator /= (const double a){
      x /= a;
      y /= a;
      z /= a;
    }
    inline void operator*= (const double a){
      x *= a;
      y *= a;
      z *= a;
    }
    friend Vector operator +(Vector lhs, const Vector& rhs){
      lhs += rhs;
      return lhs;
    }
    friend Vector operator -(const Vector& rhs){
      return -1*rhs;
    }    
    friend Vector operator -(Vector lhs, const Vector& rhs){
      lhs -= rhs;
      return lhs;
    }
    friend Vector operator *(Vector lhs, double a){
      lhs *= a;
      return lhs;
    }
    friend Vector operator *(double a,  Vector rhs){
      return rhs*a;
    }    
    friend Vector operator /(Vector lhs, double a){
      lhs/=a;
      return lhs;
    }
    friend std::ostream& operator <<(std::ostream& outs, const Vector& vec);

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

    double Mag() const;
    //Returns the magnitude of this vector.

    double Angle(const Vector &vec) const;
    //Returns the 3-dimensional angle between this vector and the argument.

    Vector ChangeCoord(const Vector &new_x_axis,const Vector &new_y_axis) const;
    Vector ChangeCoord(const Vector &new_z_axis) const;
    //Returns this vector, rotated to a new coordinate system with the argument as the z-axis.
    //The vector is rotated in such a way that a vector pointing in the z direction will be 
    //rotated onto the new axis.

    Vector Unit() const;
    //Returns a unit vector in the same direction as this vector.

    //Accessor functions
    double GetX() const;
    double GetY() const;
    double GetZ() const;
    double Theta() const;
    double Phi() const;
    void Print() const;

    //Mutator functions
    void SetX(double inp);
    void SetY(double inp);
    void SetZ(double inp);
    void SetXYZ(double inpx,double inpy,double inpz);
    void Reset(double x_inp, double y_inp, double z_inp);

  protected:
    //Class variables
    double x; //x component of vector
    double y; //y component of vector
    double z; //z component of vector
    double theta;  //theta component of vector in radians
    double phi; //phi component of vector in radians

    //Class private functions
    void UpdateThetaPhi(); 
    //This method finds theta and phi from the x,y,z Cartesian coordinates.  
    //It should be called at any time that the x,y,z components are modified,
    //so that the theta and phi components are current at all times.

  }; //class Vector

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Vector Constants

  static Vector x_axis = Vector(1,0,0);
  static Vector y_axis = Vector(0,1,0);
  static Vector z_axis = Vector(0,0,1);

}
#endif
