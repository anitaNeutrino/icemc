#ifndef ICEMC_NEUTRINO_H
#define ICEMC_NEUTRINO_H

#include "Geoid.h"

namespace icemc {

  /**
   * @class Neutrino
   * @brief All the randomly selected (and derived) numbers that describe a simulated neutrino
   */
  class Neutrino {
  public:

    /**
     * @class Flavor
     * @brief enum for neutrino flavour
     */
    enum class Flavor : int {
			     e = 1,
			     mu = 2,
			     tau = 3
    };

    /**
     * @class CurrentType
     * @brief enum for type of interaction
     */
    enum class CurrentType : int {
				  Charged = 0,
				  Neutral = 1
    };

  public:

    double Energy() const {return fEnergy;}

  private:
    double fEnergy; ///< (eV )electron volts
    Geoid::Position fInteraction; ///< Where the interaction happens    
    Flavor fFlavor; ///< Neutrino flavor
    CurrentType fCurrent; ///< Interaction current
    
  };
  




};

/** 
 * For a nice cout/cerr/logging experience
 * 
 * @param os is a output string stream
 * @param f is the Flavor class enum
 * 
 * @return the updated output string stream
 */
std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::Flavor& f);




/** 
 * For a nice cout/cerr/logging experience
 * 
 * @param os is a output string stream
 * @param c is the Current class enum
 * 
 * @return the updated output string stream
 */
std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::CurrentType& c);

#endif
