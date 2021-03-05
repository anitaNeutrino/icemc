#ifndef ICEMC_LOOP_INFO_H
#define ICEMC_LOOP_INFO_H

#include "TObject.h"
#include "TMath.h"

namespace icemc {

  
  /**
   * @class LoopInfo
   * @brief Data to track how the neutrino faired in the main neutrino generation loop, including weights
   */
  class LoopInfo {
  public:    
    
    /**
     * @class Step
     * @brief A glorifed pair of boolians. Did we evaluate this bit of the loop? And what was the result?
     */
    class Step {
    public:
      Step (){;}
      Step(bool b);
      Step& operator=(bool b);
      bool operator==(bool b) const;
      bool operator!=(bool b) const;
      bool wasEvaluated() const;
      bool getResult() const;
    private:
      bool evaluated = false;
      bool result = false;
      ClassDef(Step, 1)
    };

    Int_t run;
    UInt_t eventNumber;
    double eventTime;
    Step rayTracingSolution;
    Step chanceInHell;
    Step passesTrigger;

    void setPositionWeight(double volumeFraction);
    double positionWeight = 0;
    double directionWeight = 0;
    double dTheta;
    double weight() const { return directionWeight*positionWeight; };
    
    void next(UInt_t nextEventNumber, double nextEventTime){
      eventNumber = nextEventNumber;
      eventTime = nextEventTime;
      rayTracingSolution = false;
      chanceInHell = false;
      passesTrigger = false;
    }
    
    ClassDef(LoopInfo, 1);
  };


  
}


bool operator==(bool b, icemc::LoopInfo::Step s);
bool operator!=(bool b, icemc::LoopInfo::Step s);
std::ostream& operator<<(std::ostream& os, icemc::LoopInfo::Step s);

#endif
