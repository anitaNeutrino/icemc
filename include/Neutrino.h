#ifndef ICEMC_NEUTRINO_H
#define ICEMC_NEUTRINO_H

#include <iostream>

#include "Geoid.h"
#include "Energy.h"

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
     * @class L
     * @brief For lepton number, are we matter (neutrinos) or anti-matter (anti-neutrinos)
     */

    enum class L : int { ///@todo would it make more sense to newcomers if lepton number of 1 was for matter, -1 for anti-matter?
			Matter = 0,
			AntiMatter = 1
    };

    class Interaction {
    public:

      /**
       * @class Current
       * @brief enum for type of interaction
       */
      enum class Current : int {
				Charged = 0,
				Neutral = 1
      };
    
      Geoid::Position position;
      double crossSection;
      double length;
      Current current;
      double y;
    };

    class Path {
    public:
      TVector3 direction;
      Geoid::Position entry;
      Geoid::Position exit;
      double weight;      
    };
    
    Energy energy;
    Flavor flavor = Flavor::e; ///< Neutrino flavor
    L leptonNumber = L::Matter;
    Interaction interaction;
    Path path;
  };
}



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
std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::Interaction::Current& c);


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


