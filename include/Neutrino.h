#ifndef ICEMC_NEUTRINO_H
#define ICEMC_NEUTRINO_H

#include <iostream>

#include "Geoid.h"
#include "Energy.h"
#include "TObject.h"

namespace icemc {

  /**
   * @class Neutrino
   * @brief All the randomly selected (and derived) numbers that describe a simulated neutrino
   */
  class Neutrino {
    friend class NeutrinoFactory;
    
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

    static int pdgCode(Flavor f);
    int pdgCode() const;


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
      double crossSection = -1;
      double length = -1;
      Current current = Current::Charged;
      double y = -1;
      ClassDef(Interaction, 1);
    };

    class Path {
    public:
      TVector3 direction; ///< unit vector direction of neutrino travel
      Geoid::Position entry; ///< Where the neutrino entered the Earth.
      Geoid::Position interaction; ///< Same as interaction.position ///@todo deduplicate?
      Geoid::Position exit; ///< Where it would have left the Earth if it didn't interact
      double columnDepth = -1; ///< Line integral of Earth density from entry to interaction
      double columnDepthInteractionToExit = -1; ///< Line integral of Earth density from interaction to exit
      double weight = -1; ///@todo

      double length(bool includeDistanceToExit = false) const {
	const Geoid::Position& end = includeDistanceToExit ? exit : interaction;
	const TVector3 d = end - entry;
	return d.Mag();
      }

      ClassDef(Path, 3);
    };
    
    Energy energy;
    Flavor flavor = Flavor::e; ///< Neutrino flavor
    L leptonNumber = L::Matter;
    Interaction interaction;
    Path path;
    ClassDef(Neutrino, 1);




    double weight(){
      return 0; ///@todo product of all weights;
    };

  private:
    void setInteractionPosition(const Geoid::Position& p){
      // keep these the same
      interaction.position = p;
      path.interaction = p;
    }
    
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


