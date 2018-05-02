#ifndef ICEMC_ANITA_FULL_H
#define ICEMC_ANITA_FULL_H

#include "Detector.h"
#include "anita.hh"
#include "balloon.hh"



namespace icemc {

  class Ray;
  class Counting;
  class Screen;

  /**
   * @class ANITA
   * @brief Implements the Detector virtual functions and combines the different aspects of the simulated ANITA payload into a single class.
   */

  class ANITA : public Detector, public Anita, public Balloon {
  public:
    ANITA(const Settings* settings, Counting* retardedClass, Ray* sillyRay, Screen* sillyPanel);
    virtual ~ANITA();
    
    virtual GeographicCoordinate getCenterOfDetector();
    virtual const std::vector<TVector>& getFieldCalcLocationsRelativeToAveDetPos();
    virtual bool applyTrigger(const std::vector<AskaryanSignal>& signals);    
    virtual void getNDtForTimeDomainAskaryanSignals(int& n, double& dt) const;


    double GetAverageVoltageFromAntennasHit(const Settings *settings1,  int *nchannels_perrx_triggered,  double *voltagearray,  double& volts_rx_sum) const;
  private:
    std::vector<TVector> testVecNotRealYet;
    const Settings* fSettingsPtrIDontOwn;
    Counting* fCountingPtrIDontOwn;
    Ray* fRayPtrIDontOwn;
    Screen* fScreenPtrIDontOwn;    










    // Complete junk from EventGenerator.h


    
  };




};


#endif //ICEMC_ANITA_FULL_H
