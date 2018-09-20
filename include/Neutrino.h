#ifndef ICEMC_NEUTRINO_H
#define ICEMC_NEUTRINO_H

#include <iostream>

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
     * @class Current
     * @brief enum for type of interaction
     */
    enum class Current : int {
			      Charged = 0,
			      Neutral = 1
    };

    /**
     * @class L
     * @brief For lepton number, are we matter (neutrinos) or anti-matter (anti-neutrinos)
     */

    enum class L : int { ///@todo it would make more sense to newcomers if lepton number of 1 was for matter, -1 for anti-matter
			Matter = 0,
			AntiMatter = 1
    };

    double energy; ///< (eV )electron volts
    double crossSection; /// @todo units
    double interactionLength; ///@todo units
    
    Flavor flavor; ///< Neutrino flavor
    Current interactionCurrent; ///< Interaction current
    L leptonNumber = L::Matter;
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
std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::Current& c);


/** 
 * For a nice cout/cerr/logging experience
 * 
 * @param os is a output string stream
 * @param l is the L class enum
 * 
 * @return the updated output string stream
 */
std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::L& l);

#endif


