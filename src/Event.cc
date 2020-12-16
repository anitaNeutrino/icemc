#include "Event.h"

ClassImp(icemc::Event)


icemc::EventSummary::EventSummary(Int_t run, UInt_t eventNumber, double eventTime){
  loop.run = run;
  loop.next(eventNumber, eventTime);
}

void icemc::EventSummary::setWeights(double volumeFraction, double dtheta){

  positionWeight = volumeFraction > 0 ? 1/volumeFraction : DBL_MAX; //@reasonable if no ice visible?
  
  if (dtheta > TMath::Pi()/2)
    dtheta = TMath::Pi()/2;

  if (dtheta == -1 || dtheta == 0) //@todo is this reasonable?
    directionWeight = 0;
  else
    directionWeight = 2/(1-cos(dtheta));
}


void icemc::Event::copy(const EventSummary& summary){
  loop                       = summary.loop;
  neutrino                   = summary.neutrino;
  interaction                = summary.interaction;
  shower                     = summary.shower;
  detector                   = summary.detector;
  signalSummaryAt1m          = summary.signalSummaryAt1m;
  signalSummaryAtDetector    = summary.signalSummaryAtDetector;
  positionWeight             = summary.positionWeight;
  directionWeight             = summary.directionWeight;
}
