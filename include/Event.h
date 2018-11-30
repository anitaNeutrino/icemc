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
   * @class EventSummary
   * @brief A class which contains sumary info about an icemc event.
   */
  class EventSummary {

  public:
    EventSummary(Int_t run=-1, UInt_t eventNumber=0, double eventTime=TMath::QuietNaN());
    
    LoopInfo loop;
    Neutrino neutrino;
    Interaction interaction;
    Shower shower;
    Geoid::Position detector;
    SignalSummary signalSummaryAt1m;
    SignalSummary signalSummaryAtDetector;

    ClassDef(EventSummary, 3)
  };



  /**
   * @class Event
   * @brief A class which contains all the info about an icemc event.
   */
  class Event : public EventSummary {
  public:
    Event() : EventSummary() {}
    void copy(const EventSummary& other);
    
    TGraph signalAt1m;
    TGraph signalAtDetector;    
    ClassDef(Event, 3)
  };
  
}


#endif // ICEMC_EVENT_H
