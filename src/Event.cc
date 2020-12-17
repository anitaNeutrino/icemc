#include "Event.h"

ClassImp(icemc::Event)


icemc::EventSummary::EventSummary(Int_t run, UInt_t eventNumber, double eventTime){
  loop.run = run;
  loop.next(eventNumber, eventTime);
}

void icemc::Event::copy(const EventSummary& summary){
  loop                       = summary.loop;
  neutrino                   = summary.neutrino;
  interaction                = summary.interaction;
  shower                     = summary.shower;
  detector                   = summary.detector;
  signalSummaryAt1m          = summary.signalSummaryAt1m;
  signalSummaryAtDetector    = summary.signalSummaryAtDetector;
}
