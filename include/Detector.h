#ifndef ICEMC_DETECTOR_H
#define ICEMC_DETECTOR_H

#include "TGraph.h"
#include "TVector.h" ///@todo use TVector or icemc::Vector?

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
     * @brief Where and how many places should the Askaryan field be calculated relative to the "average" detector position?
     * 
     * You might want to calculate the Askaryan E-field separately at each antenna, or just the center of the detector.
     * 
     * @return A const reference to 
     */
    virtual const std::vector<TVector>& getFieldCalcLocationsRelativeToAveDetPos() = 0;

    

    /**
     * @brief Detect (or not) the Askaryan radiation arriving at the payload.
     * 
     * @todo figure out what else is required, polarization and direction?
     * 
     * @param pureSignalVoltageTimeGraphs A set of voltage/time Askaryan signals arriving at the locations specified by getFieldCalcLocationsRelativeToAveDetPos.
     * The size of the parameter should be equal to the size of the what was returned by getFieldCalcLocationsRelativeToAveDetPos.
     *
     * @param poyntingVector The direction of travel of the waveform (@todo in what coordinate system?)
     * @param polarizationVector The polarization vector of the waveform (@todo in what coordinate system?)
     * 
     * @return true if the waveform triggered the instrument, false if it did not.
     */
    virtual bool applyTrigger(const std::vector<TGraph>& pureSignalVoltageTimeGraphs, const TVector& poyntingVector, const TVector& polarizationVector) = 0;

    /** 
     * @brief Tell icemc how you like your Askaryan signals.
     * 
     * icemc will generate time domain Askaryan signals for your detector of any length and time-step.
     * Tell it what you want here.
     * @todo specify units of dt.
     * 
     * @param n number of samples for the generated time-domain Askaryan signal
     * @param dt time step (@todo units?) for the generated time-domain Askaryan signal
     */
    virtual void getNDtForTimeDomainAskaryanSignals(int& n, double& dt) const = 0;
    
  };  

}


#endif // ICEMC_DETECTOR_H
