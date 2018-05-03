#ifndef ICEMC_DETECTOR_H
#define ICEMC_DETECTOR_H

#include "TGraph.h"
#include "TVector.h" ///@todo use TVector or icemc::Vector?
#include "vector.hh" ///@todo use TVector or icemc::Vector?

namespace icemc {

  /**
   * @todo Move this to the EarthModel class?
   * 
   */
  struct GeographicCoordinate {
    double longitude;
    double latitude;
    double altitude;
    UInt_t unixTime;
  };
  


  /**
   * @class AskaryanSignal
   * @brief A waveform traveling along a particular direction with a particular polarization
   */

  class AskaryanSignal {
  public:
    TGraph waveform; ///< E-field (Volts/m) vs. time (seconds)
    icemc::Vector poynting; /// < Direction of signal travel (in the icemc coordinate system).
    icemc::Vector polarization; ///< Polarization vector (in the icemc coordinate system).
  };  


  /**
   * @class Detector
   * @brief Abstract detector class 

   * Enforces separation of detector simulation from UHEN and Askaryan RF simulation.
   * All detectors that interact with icemc *MUST* inherit from this class 
   * and implement the pure virtual functions.
   * 
   */
  class Detector {
  public:    

    
    /** 
     * @brief Where is the detector?
     * 
     * icemc will try to generate neutrinos inside the horizon of your detector.
     * That depends on where your detector is.
     * For in-ice detectors this should be trivial to implement.
     * For moving detectors, like ANITA, you can return the GPS position (and time) to generate the neutrino.
     * 
     * @return lon/lat/alt/unixTime 
     */
    virtual GeographicCoordinate getCenterOfDetector() = 0;


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
     * @return a reference to the relative position
     */
    virtual const icemc::Vector& getPositionRX(int rx) const = 0;
    



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
    virtual void addSignalToRX(const AskaryanSignal& signal, int i) = 0;
    

    /**
     * @brief Run the propagate the accumulated signals at each receiver through the detector trigger
     * 
     * @return true if event passes trigger, false otherwise.
     */
    virtual bool applyTrigger() = 0;



    /** 
     * @brief Tell icemc how you like your time-domain Askaryan signals.
     * 
     * icemc will generate time domain Askaryan signals for your detector of any length and time-step.
     * Tell it what you want here.
     * @todo specify units of dt.
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
