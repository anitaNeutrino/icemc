#ifndef VECTOR_H_
#define VECTOR_H_

#include <iostream>

namespace icemc {

  /**
   * @class Vector
   * @brief A three dimensional vector
   *
   * Operators are overloaded to provide for the familiar operations: vector addition, subtraction,
   * scalar multiplication and division, the dot product and cross products.
   * The x,y,z components can be inserted into an ostream with the << operator.
   *
   * Methods are provided to perform various rotations.
   *
   * The only methods that will alter an established vector are Set functions.
   * All other methods (cross product, rotations, etc.) return a new Vector object.
   */
  class Vector {

  private:
    // use getter functions GetX(), GetY(), GetZ() to access these
    // use setter functions SetX(x), SetY(y), SetZ(z) to change them
    double fX = 0; ///< x-component of vector
    double fY = 0; ///< y-component of vector
    double fZ = 1; ///< z-component of vector

  public:
    /**
     * Default constructor: initializes a unit vector along the z-axis
     */
    Vector() {}

    /**
     * Constructor: Initialize a new vector with given values of x, y, and z.
     *
     * @param x x-component
     * @param y y-component
     * @param z z-component
     */
    Vector(double x, double y, double z)
      : fX(x), fY(y), fZ(z)
    {
    }



    /**
     * Constructor: Initialize a new unit vector in the theta, phi direction.
     *
     * theta and phi must be in RADIANS!
     * Accepts theta from 0 to PI, and any phi.
     *
     * @param theta polar angle
     * @param phi radial angle
     */
    Vector(double theta, double phi);








    /**
     * Functions that make new vectors, but don't change this vector
     *
     */

    /**
     * Vector addition
     *
     * @param v1 first of the vectors to be summed
     * @param v2 second of the vectors to be summed
     *
     * @return v1 + v2
     */
    friend Vector operator +(Vector v1, const Vector& v2){
      v1 += v2;
      return v1;
    }

    /**
     * Vector reversal
     *
     * @param v the vector to reverse
     *
     * @return a vector of equal magnitude but opposite direction
     */
    friend Vector operator -(const Vector& v){
      return -1*v;
    }

    /**
     * Vector subtraction
     *
     * @param v1 is the initial vector
     * @param v2 is subtracted from the initial vector
     *
     * @return v1 - v2
     */
    friend Vector operator -(Vector v1, const Vector& v2){
      v1 -= v2;
      return v1;
    }

    /**
     * Scalar multiplication
     *
     * @param v1 is the vector to scale
     * @param a is the scale factor to scale v1 up by
     *
     * @return a copy of v1 scaled up by the factor a
     */
    friend Vector operator *(Vector v1, double a){
      v1 *= a;
      return v1;
    }

    /**
     * Scalar multiplication
     *
     * @param a is the scale factor
     * @param v1 is the vector to scale
     *
     * @return a copy of v1 scaled by the factor a
     */
    friend Vector operator *(double a,  Vector v1){
      return v1*a;
    }

    /**
     * Scalar division
     *
     * @param v1 is the vector to scale
     * @param a is the factor to scale v1 down by
     *
     * @return a copy of v1 scaled down by the factor a
     */
    friend Vector operator /(Vector v1, double a){
      v1/=a;
      return v1;
    }




    /**
     * Returns the vector rotated counterclockwise (right handed coordinates)
     * by angle radians around the x-axis.
     * N.B. : Returns a new Vector object.  Does not change the vector it is called from.
     *
     * @param angle of rotation in radians
     *
     * @return a copy of this vector, rotated angle radians around the x-axis
     */
    Vector RotateX(double angle) const;

    /**
     * Returns the vector rotated counterclockwise (right handed coordinates)
     * by "angle" radians about the Y axis.
     * N.B. : Returns a new Vector object.  Does not change the vector it is called from.
     *
     * @param angle of rotation in radians
     *
     * @return a copy of this vector, rotated angle radians around the y-axis
     */
    Vector RotateY(double angle) const;


    /**
     * Returns the vector rotated counterclockwise (right handed coordinates)
     * by angle radians about the Z axis.
     * N.B. : Returns a new Vector object.  Does not change the vector it is called from.
     *
     * @param angle of rotation in radians
     *
     * @return a copy of this vector, rotated angle radians around the z-axis
     */
    Vector RotateZ(double angle) const;


    /**
     * Takes the cross product this vector with another
     *
     * @param vec The vector to cross with this vector
     *
     * @return this vector crossed with vec
     */
    Vector Cross(const Vector &vec) const;


    /**
     * Take the dot product of this vector with vec
     *
     * @param vec is the vector with which to find the dot product
     *
     * @return this vector dotted with vec
     */
    double Dot(const Vector &vec) const;


    Vector Rotate(double angle, const Vector &axis) const;
    //Returns the vector that is this vector rotated around the vector "axis" by angle (in radians) "angle".

    double Mag() const;
    //Returns the magnitude of this vector.

    double Square() const;
    //Returns the squared magnitude of this vector.

    double Angle(const Vector &vec) const;
    //Returns the 3-dimensional angle between this and vec.

    Vector ChangeCoord(const Vector &new_x_axis,const Vector &new_y_axis) const;
    Vector ChangeCoord(const Vector &new_z_axis) const;
    //Returns this vector, rotated to a new coordinate system with the argument as the z-axis.
    //The vector is rotated in such a way that a vector pointing in the z direction will be
    //rotated onto the new axis.

    /**
     * Create a normalized copy of this vector.
     * @warning, if a this is called on a vector of length 0, then a vector of length 0 is also returned.
     *
     * @return A unit vector in the same direction as this vector.
     */
    Vector Unit() const;






    /****************************************************************
     * Getter functions (don't change anything or make anything new)
     ***************************************************************/


    /**
     * Get the x-component
     *
     * @return #fX
     */
    double GetX() const {return fX;}

    /**
     * Get the y-component
     *
     * @return #fY
     */
    double GetY() const {return fY;}

    /**
     * Get the z-component
     *
     * @return #fZ
     */
    double GetZ() const {return fZ;}


    double Theta() const {updateThetaPhi(); return fTheta;}
    double Phi() const {updateThetaPhi(); return fPhi;}

    inline double operator[](int i) const {
      //code taken from ROOT's TVector3 class
      switch(i) {
      case 0:
	return fX;
      case 1:
	return fY;
      case 2:
	return fZ;
      default:
	std::cerr << "icemc::Vector::operator[" << i << "] is only valid for 0, 1, 2, returning zero." << std::endl;
      }
      return 0.;
    }






    /****************************************************************
     * Setter functions, which DO modify the vector contents
     ***************************************************************/

    void SetX(double x);
    void SetY(double y);
    void SetZ(double z);
    void SetXYZ(double x, double y, double z);

    /**
     * Vector addition
     *
     * @param v2 is added to this vector
     */
    inline void operator +=(const Vector& v2) {
      fX += v2.fX;
      fY += v2.fY;
      fZ += v2.fZ;
      if(v2.fX != 0 || v2.fY != 0 || v2.fZ != 0){
	fAnglesDirty = true;
      }
    }

    /**
     * Vector subtraction
     *
     * @param v2 is subtracted form this vector
     */
    inline void operator -=(const Vector& v2) {
      fX -= v2.fX;
      fY -= v2.fY;
      fZ -= v2.fZ;
      if(v2.fX != 0 || v2.fY != 0 || v2.fZ != 0){
	fAnglesDirty = true;
      }
    }

    /**
     * Scalar division
     * (@note this does not dirty fTheta/fPhi)
     *
     * @param a is what all components are divided by
     */
    inline void operator /= (const double a){
      fX /= a;
      fY /= a;
      fZ /= a;
    }

    /**
     * Scalar multiplication
     * (@note this does not dirty fTheta/fPhi)
     *
     * @param a is what all components are multiplied by
     */
    inline void operator*= (const double a){
      fX *= a;
      fY *= a;
      fZ *= a;
    }





  private:
    /**
     * Instead of calculating fTheta and fPhi at every change,
     * we now just mark them dirty (fAnglesDirty=true) if fX, fY, or fZ change.
     * They are then lazily calculated with updateThetaPhi() if accessed with Theta() or Phi().
     * fAnglesDirty is then set to false.
     * Since this can potentially happen if the vector is const,  we mark them mutable
     */
    mutable double fTheta = 0;  ///< theta component of vector in radians
    mutable double fPhi = 0; ///< phi component of vector in radians
    mutable bool fAnglesDirty = true; ///< Book keep whether or not we need to recalculate fTheta/fPhi
    /**
     * Calculates fTheta and fPhi from the fX, fY, fZ cartesian coordinates.
     * This is done on demand when Theta() or Phi() is called if fAnglesDirty is true.
     * fAnglesDirty is set to true if any of fX, fY, fZ is modified.
     */
    void updateThetaPhi() const;

  };

  /**
   * Useful vector Constants
   */
  static Vector x_axis = Vector(1,0,0);
  static Vector y_axis = Vector(0,1,0);
  static Vector z_axis = Vector(0,0,1);

}

std::ostream& operator <<(std::ostream& outs, const icemc::Vector& vec);

#endif
