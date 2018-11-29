#ifndef ICEMC_EVENT_H
#define ICEMC_EVENT_H

#include "TObject.h"
#include "TMath.h"

#include "Geoid.h"

#include "LoopInfo.h"
#include "ShowerModel.h"
#include "Neutrino.h"
#include "PropagatingSignal.h"

namespace icemc {

  /**
   * @class Event
   * @brief A class which contains all the info about an icemc event.
   */
  class Event {
    // class Weights {
    //   double position;
    //   double direction;
    //   operator double(){return position*direction;}
    // };

  public:
    Event(Int_t run=-1, UInt_t eventNumber=0, double eventTime=TMath::QuietNaN());
    
    LoopInfo loop;
    Neutrino neutrino;
    Interaction interaction;
    Shower shower;
    Geoid::Position detector;
    SignalSummary signalAtInteraction;
    SignalSummary signalAtDetector;

    ClassDef(Event, 1)
  };
}


#endif // ICEMC_EVENT_H
