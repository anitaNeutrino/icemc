#ifndef ICEMC_ANITA_FULL_H
#define ICEMC_ANITA_FULL_H

#include "Detector.h"
#include "anita.hh"
#include "balloon.hh"
#include "Seavey.h"
#include "AnitaSimOutput.h"

namespace icemc {

  class RayTracer;
  class Screen;
  class RootOutput;

  /**
   * @class ANITA
   * @brief Implements the Detector virtual functions and combines the different aspects of the simulated ANITA payload into a single class.
   */

  class ANITA : public Detector, public Balloon, public Anita {
  public:
    ANITA(const Settings* settings, const RayTracer* sillyRay, const RootOutput* ro);
    virtual ~ANITA();

    virtual int getNumRX() const override {return fNumRX;}
    virtual icemc::Vector getPositionRX(int i) const override;

    virtual icemc::Position getCenterOfDetector(UInt_t unixTime = 0) override;
    virtual bool applyTrigger() override {return  applyTrigger(-1);}
    virtual bool applyTrigger(int inu);
    virtual void getDesiredNDt(int& n, double& dt) const override {
      n = 1024;
      dt = 1e-9*1./2.6;
    }

    virtual void addSignalToRX(const PropagatingSignal& signal, int rx) override {
      addSignalToRX(signal, rx, -1);
    }
    virtual void addSignalToRX(const PropagatingSignal& signal, int rx, int inu); // just for debugging

    double GetAverageVoltageFromAntennasHit(const Settings *settings1,  int *nchannels_perrx_triggered,  const double *voltagearray,  double& volts_rx_sum) const;

    UInt_t getLastEventNumber() const {return fEventNumber;}

    const RayTracer* fRayPtrIDontOwn; /// @todo temp public for test
    /** 
     * What's the ilayer/ifold of given RX?
     * 
     * @param rx index of the fSeaveys
     * @param ilayer layer of ANITA 
     * @param ifold index in phi, maybe...
     */
    void getLayerFoldFromRX(int rx, int& ilayer, int& ifold) const; ///@todo temp public for test

    /** 
     * What's the ilayer/ifold of given trigger RX?
     * 
     * @param rx index of the fSeaveys
     * @param ilayer layer of ANITA 
     * @param ifold index in phi, maybe...
     */
    void getLayerFoldFromTriggerRX(int rx, int& ilayer, int& ifold) const; ///@todo temp public for test
    
  private:

    const Settings* fSettings;
    int fNumRX;

    std::vector<Seavey> fSeaveys; ///< The set of Seavey antennas on the payload
    UInt_t fEventNumber = 0;

    ///@todo maybe move or remove these things...? Although they are used in the ROOTified ANITA style data trees...
    VoltsRX fVoltsRX;
    int fL3trig[Anita::NPOL] = {0};  // 16 bit number which says which phi sectors pass L3 V-POL
    // For each trigger layer,  which "clumps" pass L2.  16 bit,  16 bit and 8 bit for layers 1 & 2 and nadirs
    int fL2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX] = {{0}};
    //For each trigger layer,  which antennas pass L1.  16 bit,  16 bit and 8 bit and layers 1,  2 and nadirs
    int fL1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX] = {{0}};

    friend class AnitaSimOutput; ///@todo Can I do this and respect privacy with getters?
    AnitaSimOutput fAnitaOutput; ///< Handles converting the MC output into the same format as real ANITA data



    void initSeaveys(const Settings *settings1, const Anita *anita1);
  };
}

#endif //ICEMC_ANITA_FULL_H
