#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TVector3.h"

namespace icemc {

  /**
   * @namespace constants
   * @brief Well, it's for things that are constants :P
   */
  namespace constants {
  
    // Mathematical constants     
    constexpr double TWOPI = 6.2831852;
    constexpr double PI    = 3.141592654;
    constexpr double ALOG2 = 0.693147;        // natural log of 2
    constexpr double INV_E = 0.36787944;      // 1/e
    constexpr double sr    = 4*PI;
    
    /** 
     * Get the Poisson error associated from n events, currently valid up to 20
     * 
     * @param n the number of observed events,  must be less than #maxPoissonStats
     * @param poissonErrorPlus the error in the positive direction
     * @param poissonErrorMinus the error in the negative direction
     */
    void getPoissonError(int n, double& poissonErrorPlus, double& poissonErrorMinus);
    const int maxPoissonStats = 20;

    
    // conversion constants
    const double CMINCH=2.54;          // inches to cm
    const double RADDEG=0.017453292;   // radians/degree  
    const double DEGRAD=57.2957795;    // degree/rad

    // physical constants
    constexpr double CLIGHT=3.0E8;            // speed of light m/s
    constexpr double KBOLTZ=1.38E-23;         // Boltzmann constant J/K
    constexpr double Z0=377.;                 // resistivity of free space
    constexpr double Zr=50.;  // radiation resistance (50 Ohms?)
    constexpr double M_NUCL=1.66E-27;         // amu mass in kg
    constexpr double MTAU=1.777E9;            // mass of the tau
    constexpr double TAUDECAY_TIME=290.6E-15; // lifetime of tau
    constexpr double VACUUM_PERMEABILITY = 1e-7*sr; //mu_0 
    constexpr double VACUUM_PERMITTIVITY = 1./(CLIGHT*CLIGHT*VACUUM_PERMEABILITY); //epsilon_0

    // properties of water
    const double X0H20=0.361;          // radiation length of water (meters)
    const double RHOH2O = 1000;        // density of water (kg/m^3)

    // properties of air
    const double Z_AIR=377;            // resistance of air = sqrt(epsilon/mu)
    const double RHOAIR=1.25;          // density of air (kg/m**3)

    // // properties of ice
    const double NFIRN=1.3250;                   // index of refraction at the very surface - Peter
    const double NICE=1.79;                      // index of refraction of ice 

    //Constants relating to all ice models
    const double FIRNDEPTH=150.0;                // depth of the firn, in meters: currently a constant over all ice
   
  
    // constant vectors used in balloon class - oindree 
    const TVector3 const_z(0,0,1);
    const TVector3 const_y(0,1,0);
    const TVector3 const_x(1,0,0);

    // TUFF configuration switching times in ChanTrigger.cc
    const    int TUFFconfig_B_end_1 = 1480713195;
    const    int TUFFconfig_P_end_1 = 1480730719;
    const    int TUFFconfig_C_end_1 = 1480731802;
    const    int TUFFconfig_P_end_2 = 1480807284;
    const    int TUFFconfig_G_end_1 = 1481013795;
    const    int TUFFconfig_O_end_1 = 1481100915;
    const    int TUFFconfig_G_end_2 = 1481173515;
    const    int TUFFconfig_O_end_2 = 1481490208;
    const    int TUFFconfig_P_end_3 = 1481642754;
    const    int TUFFconfig_B_end_2 = 1482121239;
    const    int TUFFconfig_P_end_4 = 1482168627;
    const    int TUFFconfig_B_end_3 = 1482205359;
    const    int TUFFconfig_A_end_1 = 1482206201;
    const    int TUFFconfig_B_end_4 = 1482286948;
    const    int TUFFconfig_P_end_5 = 1482347440;
    const    int TUFFconfig_B_end_5 = 1482445716;
    const    int TUFFconfig_P_end_6 = 1482465408;
    const    int TUFFconfig_B_end_6 = 1482964570;
    const    int TUFFconfig_P_end_7 = 1482987942;
  }
}
#endif //CONSTANTS_H
