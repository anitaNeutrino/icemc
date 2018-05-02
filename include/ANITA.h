#ifndef ICEMC_ANITA_FULL_H
#define ICEMC_ANITA_FULL_H

#include "Detector.h"
#include "anita.hh"
#include "balloon.hh"



namespace icemc {

  /**
   * @class ANITA
   * @brief Implements the Detector virtual functions and combines the different aspects of the simulated ANITA payload into a single class.
   */

  class ANITA : public Detector, public Anita, public Balloon {
  public:
    ANITA();
    virtual ~ANITA();
    
    virtual GeographicCoordinate getCenterOfDetector();
    virtual const std::vector<TVector>& getFieldCalcLocationsRelativeToAveDetPos();
    virtual bool applyTrigger(const std::vector<TGraph>& pureSignalVoltageTimeGraphs, const TVector& poyntingVector, const TVector& polarizationVector);
    virtual void getNDtForTimeDomainAskaryanSignals(int& n, double& dt) const;

  private:
    std::vector<TVector> testVecNotRealYet;
    
  };




};


#endif //ICEMC_ANITA_FULL_H
