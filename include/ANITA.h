#ifndef ICEMC_ANITA_FULL_H
#define ICEMC_ANITA_FULL_H

#include "Detector.h"
#include "anita.hh"
#include "balloon.hh"



namespace icemc {

  class RayTracer;
  class Screen;

  /**
   * @class ANITA
   * @brief Implements the Detector virtual functions and combines the different aspects of the simulated ANITA payload into a single class.
   */

  class ANITA : public Detector, public Anita, public Balloon {
  public:
    ANITA(const Settings* settings, RayTracer* sillyRay, Screen* sillyPanel);
    virtual ~ANITA();

    virtual int getNumRX() const {return 96;} ///@todo make this proper
    virtual icemc::Vector getPositionRX(int i) const;

    virtual icemc::Position getCenterOfDetector(UInt_t unixTime = 0);
    virtual bool applyTrigger(){return  applyTrigger(-1);}
    virtual bool applyTrigger(int inu);
    virtual void getDesiredNDt(int& n, double& dt) const {
      n = 1024;
      dt = 1e-9*1./2.6;
    }

    virtual void addSignalToRX(const PropagatingSignal& signal, int rx){
      addSignalToRX(signal, rx, -1);
    }
    virtual void addSignalToRX(const PropagatingSignal& signal, int rx, int inu); // just for debugging


    double GetAverageVoltageFromAntennasHit(const Settings *settings1,  int *nchannels_perrx_triggered,  double *voltagearray,  double& volts_rx_sum) const;


    // Indices... uuurrrggghh
    // @todo make these private again after debugging complete.
    void getAntPolFromRX(int rx, int&ant, int& pol) const;
    void getLayerFoldFromRX(int rx, int& ilayer, int& ifold) const;

  private:

    const Settings* fSettingsPtrIDontOwn;
    RayTracer* fRayPtrIDontOwn;
    Screen* fScreenPtrIDontOwn;


  };
}

#endif //ICEMC_ANITA_FULL_H
