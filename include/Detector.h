#ifndef ICEMC_DETECTOR_H
#define ICEMC_DETECTOR_H

#include "vector.hh" ///@todo use TVector or icemc::Vector?
#include "position.hh" ///@todo use TVector or icemc::Vector?
#include "FTPair.h"

namespace icemc {


  /**
   * @class PropagatingSignal
   * @brief A waveform traveling along a particular direction with a particular polarization
   */
  class PropagatingSignal {
  public:
    PropagatingSignal (const FTPair& sig, const Vector& direction,  const Vector& pol)
      : waveform(sig), poynting(direction), polarization(pol) {;}
    
    icemc::FTPair waveform; ///< E-field (Volts/m) vs. time (seconds).
    icemc::Vector poynting; /// < Direction of signal travel (in the icemc coordinate system).
    icemc::Vector polarization; ///< Polarization vector (in the icemc coordinate system).
  };  

  

  /**
   * @class Detector
   * @brief Abstract detector class 

   * Enforces separation of detector simulation from UHEN and Askaryan RF simulation.
   * All detectors that interact with icemc must inherit from this class 
   * and implement the pure virtual functions.
   * 
   */
  class Detector {
  public:    

    
    /** 
     * @brief Where is the detector?
     */
    virtual icemc::Position getCenterOfDetector(UInt_t unixTime = 0) = 0;


    /** 
     * How many antennas does the detector have? For looping.
     * 
     * @return The total number of receivers (one per polarization per antenna).
     */
    virtual int getNumRX() const = 0;


    /** 
     * Get the antennna position in the icemc coordinate system
     * 
     * @param rx is the index of the receiver
     * 
     * @return a vector containing the position in the icemc coordinate system
     */
    virtual icemc::Vector getPositionRX(int rx) const = 0;
    



    /** 
     * @brief Add a time domain signal to a receiver
     * 
     * This function should be called at least once per receiver per trigger.
     * Any application of off-axis gain (and delay) should be handled internally by the detector.
     * This function can be called an arbitrary number of times before the trigger is run.
     * This allows for the implementation of multipath signals from the Screen class.
     * @see icemc::Screen
     * 
     * @param signal is the signal to add.
     * @param i is the index of the receiver 
     */    
    virtual void addSignalToRX(const PropagatingSignal& signal, int i) = 0;
    

    /**
     * @brief Run the propagate the accumulated signals at each receiver through the detector trigger
     * 
     * @return true if event passes trigger, false otherwise.
     */
    virtual bool applyTrigger() = 0;



    /** 
     * @brief Tell icemc how you like your Askaryan signals.
     * 
     * @param n number of samples for the generated time-domain Askaryan signal
     * @param dt time step (@todo units?) for the generated time-domain Askaryan signal
     */
    virtual void getDesiredNDt(int& n, double& dt) const = 0;

    
  protected:
    std::vector<TGraph> fWaveformsRX; ///< Time domain signals at each receiver
    
  };
  
  
}


#endif // ICEMC_DETECTOR_H
