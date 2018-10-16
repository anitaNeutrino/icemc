#ifndef BALLOON_H
#define BALLOON_H


////////////////////////////////////////////////////////////////////////////////////////////////
//class Balloon:
////////////////////////////////////////////////////////////////////////////////////////////////
#include "anita.hh"
#include "Geoid.h"

#include "TVector3.h"

#include <iostream>

#ifdef ANITA_UTIL_EXISTS
#include "Adu5Pat.h"
#endif

class TChain;
class TTreeIndex;

namespace icemc {
  class RayTracer;
  class Interaction;
  class Antarctica;
  class FancyTTreeInterpolator;

  /**
   * @class FlightPath
   * @brief Make a little enum for the flight path for easier human reading of code
   */
  enum class FlightPath {FixedPosition        = 0,
			 Circle80DegreesSouth = 1,
			 AnitaLite            = 2,
			 Custom               = 5,
			 Anita1               = 6,
			 Anita2               = 7,
			 Anita3               = 8,
			 Anita4               = 9
  };

  /**
   * @class Balloon
   * @brief Handles everything related to balloon positions, payload orientation over the course of a flight.
   */
  class Balloon {
 
  public:
    Balloon(const Settings* settings);
    Balloon(FlightPath path, const Settings* settings = nullptr);
    virtual ~Balloon();

    
  protected:
    void getBalloonPosition(double eventTime,Anita *anita1 = nullptr);
  public:
 
    ///< This function sets the default balloon position
    void SetDefaultBalloonPosition();

    TVector3 RotatePayload(TVector3 ant_pos) const;
    
    ///< This function UN-rotates the payload
    /**
     * Rotate from earth to payload coord. (undoes RotatePayload)
     *
     *
     * @param  ant_pos -
     * @return returns vector
     */
    TVector3 unRotatePayload(TVector3 ant_pos) const;

    inline FlightPath whichPath() const {return WHICHPATH;}
    inline const Geoid::Position& position() const {return fPosition;}
    inline double getLatitude() const {return latitude;}
    inline double getLongitude() const {return longitude;}
    inline double getAltitude() const {return altitude;}
    inline double getHeading() const {return heading;}
    inline UInt_t getRealTime() const {return realTime;}
    double getPitch() const; ///< Converts to constant for A2 and A3 (A4?)
    double getRoll() const; ///< Converts to constant for A2 and A3 (A4?)
    // double getSurfaceUnderBalloon() const {return surface_under_balloon;}

    virtual double getStartTime() const {return fFirstRealTime;} 
    virtual double getEndTime() const {return fLastRealTime;}

    
#ifdef ANITA_UTIL_EXISTS
    /** 
     * Construct an ANITA data style Adu5Pat from the GPS info in #fChain
     * @return the constructed Adu5Pat
     */
    Adu5Pat pat() const;
#endif
    
    TChain *fChain = nullptr;
    double BN_ALTITUDE;                                         ///< pick balloon altitude
    double MAXHORIZON;                                          ///< pick the interaction within this distance from the balloon so that it is within the horizon
    double BN_LONGITUDE;                                        ///< balloon longitude for fixed balloon location
    double BN_LATITUDE;                                         ///< balloon latitude for fixed balloon location
    unsigned short surfTrigBandMask[9][2];                      ///< Ryan's 16 bit masks for 9 surfs.  2x16 bit masks gives 32 channels per surf
    // int NPOINTS;                                                ///< number of GPS positions we're picking from.
    int REDUCEBALLOONPOSITIONS;                                 ///< only take every 100th entry in the flight data file
    double altitude_bn;
    // double dtryingposition;                                     ///< weighting factor: how many equivalent tries each neutrino counts for after having reduced possible interaction positions to within horizon
    std::vector<double> latitude_bn_anitalite;                  ///< latitude at times along flightpath, equally distributed among gps data. This is filled with anita or anita-lite data, depending on which the user specifies
    std::vector<double> longitude_bn_anitalite;                 ///< same for longitude
    std::vector<double> altitude_bn_anitalite;                  ///< same for altitude
    std::vector<double> heading_bn_anitalite;                   ///< same for heading of the balloon
    std::vector<double> realtime_bn_anitalite;                  ///< same for real life time

    float flatitude,flongitude,faltitude,fheading,froll, fpitch;

  private:
    // double surface_under_balloon;                               ///< distance between center of the earth and the surface of earth under balloon
    
    void InitializeBalloon(const Settings* settings);    
    void ReadAnitaliteFlight();
    void setr_bn(double latitude,double longitude);
    
    const Settings* fSettings;
    const FlightPath WHICHPATH;
    int igps;                                                   ///< which balloon position do we use out of the 25000 anitalite GPS positions.
    int ibnposition;
    unsigned int realTime;                           ///< realtime from the flight data file

    double latitude,longitude,altitude,heading,roll,pitch;
    
    double MINALTITUDE;                                         ///< minimum altitude balloon needs to be before we consider it a good event to read from the flight data file
    int igps_previous;                                          ///< which entry from the flight data file the previous event was so we can just take the next one.
    int RANDOMIZE_BN_ORIENTATION;                               ///< 0=fixed balloon orientation,1=randomized
    float powerthresh[9][32];                                   ///< power threshold in Watts
    float meanp[9][32];                                         ///< mean power in Watts
    Geoid::Position fPosition; ///< Balloon position
    FancyTTreeInterpolator* fInterp = nullptr;
    // Geoid::Position r_bn_shadow;                                       ///< position of the balloon projected on earth surface - point just below balloon at surface of the earth

  protected:
    double fFirstRealTime = -1;
    double fLastRealTime = -1;

  }; //class Balloon
}

/** 
 * For a nice cout/cerr/logging experience
 * 
 * @param os is a output string stream
 * @param fp is the FlightPath class enum
 * 
 * @return the updated output string stream
 */
std::ostream& operator<<(std::ostream& os, const icemc::FlightPath& fp);

#endif


