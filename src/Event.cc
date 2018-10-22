#include "Event.h"

ClassImp(icemc::Event)


icemc::Event::Event(Int_t run, UInt_t eventNumber, double eventTime){
  loop.run = run;
  loop.next(eventNumber, eventTime);
}
