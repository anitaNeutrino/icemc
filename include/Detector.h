#ifndef ICEMC_DETECTOR_H
#define ICEMC_DETECTOR_H

#include "Geoid.h"
#include "PropagatingSignal.h"

namespace icemc {


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
    virtual ~Detector(){};
    
    virtual double getStartTime() const = 0;
    virtual double getEndTime() const = 0;
    
    /** 
     * @brief Where is the detector?
     */
    virtual const Geoid::Position& getPosition(double time = TMath::QuietNaN()) = 0;


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
    virtual TVector3 getPositionRX(int rx) const = 0;
    


    /** 
     * Function which can be optionally implemented in a derived detector class.
     * The idea is to do a very simple test on the signal to see if it is worth
     * doing the full trigger simulation. For example if there is basically no
     * power in the signal or the polarization is nowhere near on axis, etc.
     * then you can return false.
     * 
     * @param signal is a signal at the payload.
     * 
     * @return true if there is some 
     */
    virtual bool chanceInHell(const PropagatingSignal& signal){
      (void) signal;
      return true;
    }
    


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

  };
  
  
}


#endif // ICEMC_DETECTOR_H
