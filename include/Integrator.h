#ifndef ICEMC_INTEGRATOR_H
#define ICEMC_INTEGRATOR_H


#include "TRandom3.h"


namespace icemc {

  /**
   * @class Integrator
   * @brief Any class which samples random numbers to evaluate part of the icemc integral should inherit from this class
   */
  
  class Integrator : public TObject {
  public:
    Integrator(int run);
    void reset();

  private:
    UInt_t hash;
  protected:
    TRandom3 fRandom;
  };
}



#endif // ICEMC_INTEGRATOR_H
