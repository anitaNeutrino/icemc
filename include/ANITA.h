#ifndef ICEMC_ANITA_FULL_H
#define ICEMC_ANITA_FULL_H

#include "Detector.h"
#include "anita.hh"
#include "balloon.hh"



namespace icemc {

  class Ray;
  class Screen;

  /**
   * @class ANITA
   * @brief Implements the Detector virtual functions and combines the different aspects of the simulated ANITA payload into a single class.
   */

  class ANITA : public Detector, public Anita, public Balloon {
  public:
    ANITA(const Settings* settings, Ray* sillyRay, Screen* sillyPanel);
    virtual ~ANITA();

    virtual int getNumRX() const {return 96;} ///@todo make this proper
    virtual icemc::Vector getPositionRX(int i) const;

    virtual icemc::Position getCenterOfDetector(UInt_t* unixTime = NULL);
    virtual bool applyTrigger();
    virtual void getDesiredNDt(int& n, double& dt) const {
      n = 1024;
      dt = 1e-9*1./2.6;
    }

    virtual void addSignalToRX(const PropagatingSignal& signal, int rx);


    double GetAverageVoltageFromAntennasHit(const Settings *settings1,  int *nchannels_perrx_triggered,  double *voltagearray,  double& volts_rx_sum) const;
  private:

    const Settings* fSettingsPtrIDontOwn;
    Ray* fRayPtrIDontOwn;
    Screen* fScreenPtrIDontOwn;

    // Indices... uuurrrggghh
    void getAntPolFromRX(int rx, int&ant, int& pol) const;
    void getLayerFoldFromRX(int rx, int& ilayer, int& ifold) const;

  };
}

#endif //ICEMC_ANITA_FULL_H
