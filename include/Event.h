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
   * @brief Saved for each generated event, a small summary.
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

    // Phase weights
    double positionWeight;
    double directionWeight;
    double weight() const { return positionWeight*directionWeight; };

    ClassDef(EventSummary, 3)
  };



  /**
   * @class Event
   * @brief Saved for each passing event, more data.
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
