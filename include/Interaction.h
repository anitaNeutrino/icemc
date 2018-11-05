#ifndef ICEMC_INTERACTION_H
#define ICEMC_INTERACTION_H

#include "Geoid.h"
#include "TObject.h"

namespace icemc {

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
}



/** 
 * For a nice cout/cerr/logging experience
 * 
 * @param os is a output string stream
 * @param c is the Current class enum
 * 
 * @return the updated output string stream
 */
std::ostream& operator<<(std::ostream& os, const icemc::Interaction::Current& c);

#endif // ICEMC_INTERACTION_H
