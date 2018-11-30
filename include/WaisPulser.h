#ifndef WAIS_PULSER_H
#define WAIS_PULSER_H

#include "ImpulsiveRadioGenerator.h"
#include "InteractionGenerator.h"

namespace icemc {

  class WorldModel;
  
  class WaisPulser : public ImpulsiveRadioGenerator, public InteractionGenerator {

    const Settings* fSettings;    
    std::vector<double> fVoltsPerMeter;
    std::vector<double> fTimes;
    Geoid::Position fWaisPulserPosition;
    TVector3 fSurfaceNormal;

    void readModel();
    
  public :
    WaisPulser(const Settings* settings,  std::shared_ptr<const WorldModel> worldModel);
    virtual ~WaisPulser();

    
    virtual PropagatingSignal generateImpulse(const OpticalPath& opticalPath, const Neutrino& nu, const Shower& shower) const override;
    virtual Interaction generateInteraction(const Neutrino& n, const Geoid::Position& detector) override;    
    
  };

}


#endif // WAIS_PULSER_H
