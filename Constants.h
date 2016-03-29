#ifndef CONSTANTS_H
#define CONSTANTS_H
#include "vector.hh"


// constants in math
const double TWOPI=6.2831852;
const double PI=3.141592654;
const double ALOG2=0.693147;        // natural log of 2
const double INV_E=0.36787944;      // 1/e
const double sr=4*PI;


// conversion constants
const double CMINCH=2.54;          // inches to cm
const double RADDEG=0.017453292;   // radians/degree  
const double DEGRAD=57.2957795;    // degree/rad

// physical constants
const double CLIGHT=3.0E8;            // speed of light m/s
const double KBOLTZ=1.38E-23;         // Boltzmann constant J/K
const double Z0=377.;                 // resistivity of free space
const double Zr=50.;  // radiation resistance (50 Ohms?)
const double M_NUCL=1.66E-27;         // amu mass in kg
const double MTAU=1.777E9;            // mass of the tau
const double TAUDECAY_TIME=290.6E-15; // lifetime of tau

// properties of water
const double X0H20=0.361;          // radiation length of water (meters)


// properties of air


const double Z_AIR=377;            // resistance of air = sqrt(epsilon/mu)
  const double RHOAIR=1.25;          // density of air (kg/m**3)
// // properties of ice


const double NFIRN=1.3250;                   // index of refraction at the very surface - Peter
const double NICE=1.79;                      // index of refraction of ice



// constant vectors used in balloon class - oindree 

const Vector const_z(0,0,1);
const Vector const_y(0,1,0);
const Vector const_x(1,0,0);



#endif //CONSTANTS_H
