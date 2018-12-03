#ifndef ICEMC_PROPAGATING_SIGNAL_H
#define ICEMC_PROPAGATING_SIGNAL_H

#include "OpticalPath.h"
#include "TVector3.h"
#include "FTPair.h"
#include "Energy.h"

namespace icemc {

  class PropagatingSignal;
  /**
   * @class SignalSummary
   * @brief Summary class for plotting
   */

  class SignalSummary {
  public:
    SignalSummary(const PropagatingSignal* s = nullptr) {set(s);}
    void set(const PropagatingSignal* s = nullptr);

    double maxEField;
    double energy; ///@todo make an instance of the Energy class
  };





  /**
   * @class PropagatingSignal
   * @brief A waveform traveling along a particular direction with a particular polarization
   */
  class PropagatingSignal {
  public:
    PropagatingSignal (const FTPair& sig, const TVector3& direction,  const TVector3& pol)
      : waveform(sig), poynting(direction), polarization(pol) {;}
    
    icemc::FTPair waveform; ///< E-field (Volts/m) vs. time (seconds).
    TVector3 poynting; /// < Direction of signal travel (in the icemc coordinate system).
    TVector3 polarization; ///< Polarization vector (in the icemc coordinate system).

    /** 
     * Apply the effects of attenuation, fresnel, free-space path loss, from propagation along the optical path
     * 
     * @param opticalPath
     */
    void propagate(const icemc::OpticalPath& opticalPath);


    /** 
     * Create a summary object, to describe the signal so we don't have to store the whole thing.
     */
    SignalSummary summarize() const {return SignalSummary(this);}

    double energy() const;
    double maxEField() const;
  };



  
}
#endif
