#include "TObject.h"


#include "NeutrinoPath.h"

namespace icemc {

  /**
   * @class Event
   * @brief A class which contains all the info about an icemc event
   */

  class Weight {
    double position;
    double direction;
    operator double(){return position*direction;}
  };

  class Event : public TObject {

  public:
    Event(UInt_t eventNumber_, int run_) :
      eventNumber(eventNumber_), run(run_)
    {
      
    }
    
    UInt_t realTime;
    UInt_t eventNumber; 
    Int_t run;
    Geoid::Position detector; ///< detector pos
    Geoid::Position interaction; ///
    NeutrinoPath path;
    Weight weight;
  };

}
