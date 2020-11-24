#ifndef ICEMC_IMPULSIVE_RADIO_GENERATOR_H
#define ICEMC_IMPULSIVE_RADIO_GENERATOR_H

#include "PropagatingSignal.h"
#include "OpticalPath.h"
#include "Neutrino.h"
#include "ShowerModel.h"
#include "Detector.h"

namespace icemc {

  class ImpulsiveRadioGenerator {
  public:

    ///@todo some of these aren't necessary for pulsers e.g. WAIS, can this be made a little more elegant?
    virtual PropagatingSignal generateImpulse(const OpticalPath& opticalPath, const Neutrino& nu, const Shower& shower) const = 0;

    virtual double getThetaRange(const Detector& detector, const Neutrino& nu, const Shower& dummyShower, const OpticalPath& opticalPath) const = 0;
  };
  
};


#endif
