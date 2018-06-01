
#ifndef ICEMC_SEAVEY_H
#define ICEMC_SEAVEY_H

#include "vector.hh"
#include "FTPair.h"


class TCanvas;

namespace icemc {

  /**
   * @class Seavey
   * @brief A lightweight class to store the ANITA antenna information
   */

  class Seavey {
  public:
    Seavey() {;}

    Seavey(const Vector& antennaPositionInPayloadCoordinates)
      : fPosition(antennaPositionInPayloadCoordinates) {;}

    // void applyAntennaGain(FTPair& incomingSignal) const;

    TCanvas* plotGains() const;

    const icemc::Vector& getPosition() const {return fPosition;}
    
  private:
    icemc::Vector fPosition; ///< Position in payload centered coordinates
    icemc::Vector fEPlaneOrientation;
  };
  

}



#endif // ICEMC_SEAVEY_H
