#ifndef ICEMC_ENERGY_H
#define ICEMC_ENERGY_H

#include <iostream>

namespace icemc{

  /**
   * @class Energy
   * @brief A class to represent particle energy and do unit conversion
   */
  class Energy {
  public:

    ///@todo joules and the like?
    enum class Unit {eV,
		     keV,
		     MeV,
		     GeV,
		     TeV,
		     PeV,
		     EeV
    };

    /** 
     * Default constructor, 0 energy
     * 
     */
    Energy() : feV(0) {}

    /** 
     * Useful constructor, represent an energy with given units.
     * 
     * @param value 
     * @param unit 
     */
    Energy(double value, Unit unit){
      set(value, unit);
    }

    /** 
     * Make this object represent an energy with given units.
     * 
     * @param value 
     * @param unit 
     */
    void set(double value, Unit unit){
      feV = value/to(unit);
    }

    /** 
     * Access the energy value as a double ONLY when the unit is specified.
     * 
     * @param unit what units do we want this number in?
     * 
     * @return value of the energy in the desired units
     */
    double in(Unit unit) const {return feV*to(unit);}    
    

    /**
     * Add or subtract energy only with another Energy object
     */
    Energy& operator+=(const Energy& e){feV += e.feV; return *this;}
    Energy& operator-=(const Energy& e){feV += e.feV; return *this;}
    Energy  operator+ (const Energy& e) const {Energy e2 = *this; e2 += e; return e2;}
    Energy  operator- (const Energy& e) const {Energy e2 = *this; e2 -= e; return e2;}

    
    /**
     * Scale energy only by a unitless number
     */
    Energy& operator*=(double s){feV *= s; return *this;}
    Energy& operator/=(double s){feV /= s; return *this;}
    Energy  operator* (double s) const {Energy e; e*=s; return e;}    
    Energy  operator/ (double s) const {Energy e; e/=s; return e;}
    /* Or take a ratio of two energies */
    double  operator/(const Energy& e){return feV/e.feV;}

    /**
     * Relational operators
     * 
     */
    bool operator==(const Energy& e) const {return feV == e.feV;}
    bool operator!=(const Energy& e) const {return feV != e.feV;}    
    bool operator< (const Energy& e) const {return feV <  e.feV;}
    bool operator<=(const Energy& e) const {return feV <= e.feV;}    
    bool operator> (const Energy& e) const {return feV >  e.feV;}
    bool operator>=(const Energy& e) const {return feV >= e.feV;}    

    
  private:
    double feV; // internally we store the energy in electron_volts
    /** 
     * Internal function, gets conversion factor from eV to desired unit.
     * 
     * @param unit is the desired unit
     * 
     * @return converstion factor from eV to desired unit
     */
    double to(Unit unit) const {
      switch(unit){
      case Unit::eV:  return 1;
      case Unit::keV: return 1e-3;
      case Unit::MeV: return 1e-6;
      case Unit::GeV: return 1e-9;
      case Unit::TeV: return 1e-12;
      case Unit::PeV: return 1e-15;
      case Unit::EeV: return 1e-18;
      default: return 0; ///@todo warn?
      }
    }
  };
}


/** 
 * For scaling with left multiplication
 * 
 * @param lhs is a double
 * @param rhs is an Energy
 * 
 * @return a scaled energy value
 */
icemc::Energy operator* (double lhs, icemc::Energy rhs);


/** 
 * For a nice cout/cerr/logging experience
 * 
 * @param os is a output string stream
 * @param e is the Energy class
 * 
 * @return the updated output string stream
 */
std::ostream& operator<<(std::ostream& os, const icemc::Energy& e);



#endif
