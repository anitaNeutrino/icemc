#include "Detector.h"
#include "Report.h"
#include "RNG.h"

/**
 * icemc::Detector us an abstract class with pure virtual functions...
 * so there's nothing here! This class serves purely to enforce a 
 * standard interface for icemc to interact with any detector.
 * 
 * It's up to you to implement the pure virtual functions in your derived class.
 */




double icemc::Detector::pickEventTime(){
  static RNG random;
  return random.pickUniform(getStartTime(), getEndTime());
}
