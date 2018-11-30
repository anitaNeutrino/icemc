#ifndef WAIS_PULSER_H
#define WAIS_PULSER_H

#include "ImpulsiveRadioGenerator.h"

namespace icemc {

  class WaisPulser : public ImpulsiveRadioGenerator : public InteractionGenerator {
    void readModel();
    std::vector<double> fVoltsPerMeter;
    std::vector<double> fTimes;
    
  public :
    WaisPulser();
    virtual ~WaisPulser(){}

    PropagatingSignal generate(const OpticalPath& opticalPath, const Neutrino& nu, const Shower& shower) const override;
    
  };

}


#endif // WAIS_PULSER_H
