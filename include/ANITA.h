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

    virtual int getNumRX() const {return 96;} ///@todo make this proper
    virtual const icemc::Vector& getPositionRX(int i) const{
      return testVecNotRealYet[i];
    }
    virtual GeographicCoordinate getCenterOfDetector(){
      GeographicCoordinate gc;
      return gc;
    }
    virtual bool applyTrigger();
    virtual void getDesiredNDt(int& n, double& dt) const {
      n = 128;
      dt = 0.1;
    }

    virtual void addSignalToRX(const AskaryanSignal& signal, int rx);


    double GetAverageVoltageFromAntennasHit(const Settings *settings1,  int *nchannels_perrx_triggered,  double *voltagearray,  double& volts_rx_sum) const;
  private:
    std::vector<icemc::Vector> testVecNotRealYet;
    const Settings* fSettingsPtrIDontOwn;
    Counting* fCountingPtrIDontOwn;
    Ray* fRayPtrIDontOwn;
    Screen* fScreenPtrIDontOwn;


    // Complete junk from EventGenerator.h


    
  };




};


#endif //ICEMC_ANITA_FULL_H
