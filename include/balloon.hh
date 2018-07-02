#ifndef BALLOON_H
#define BALLOON_H


////////////////////////////////////////////////////////////////////////////////////////////////
//class Balloon:
////////////////////////////////////////////////////////////////////////////////////////////////
#include "anita.hh"
#include "position.hh"
#include "vector.hh"

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

  /**
   * @class FlightPath
   * @brief Make a little enum for the flight path for easier human reading of code
   */
  enum class FlightPath {FixedPosition        = 0,
			 Circle80DegreesSouth = 1,
			 AnitaLite            = 2,
			 BananaPlot           = 3,
			 PeterEvent           = 4,
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
    virtual ~Balloon() {;}

  


    void setObservationLocation(Interaction *interaction1,int inu, const Antarctica *antarctic, const Settings *settings1);
    // void GetBoresights(const Settings *settings1, Anita *anita1, Position r_bn, double phi_spin, Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX]);
    void PickDownwardInteractionPoint(Interaction *interaction1,Anita *anita1,const Settings *settings1, const Antarctica *antarctica1,
				      RayTracer *ray1, int &beyondhorizon); 
  
    

    /**
     * This function centers the payload
     * 
     * @param  hitangle_e 
     */
    void CenterPayload(double& hitangle_e);

    
    /**
     * This function picks the balloon position
     *
     * @param  straightup -
     * @param  antarctica -
     * @param  settings1 -
     * @param  anita1 -
     */
    void PickBalloonPosition(Vector straightup, const Antarctica *antarctica, const Settings *settings1, Anita *anita1);

    
    /**
     * Position of spot under balloon
     *
     * @param  antarctica1 -
     * @param  settings1 -
     * @param  inu -
     * @param  anita1 -
     * @param  randomNumber -
     */
    void PickBalloonPosition(const Antarctica *antarctica1, const Settings *settings1,int inu,Anita *anita1, double randomNumber);

    /**
     * Gets ith balloon position
     *
     *
     * @return returns int
     */
    int Getibnposition();

    ///< This function gets the spin of the balloon whatever that means
    /**
     * Get the azimuth of the balloon
     *
     *
     * @param  heading -
     * @return returns double
     */
    double GetBalloonSpin(double heading) const;

    
    ///< This function gets the antenna orientation
    /**
     * Long description
     *
     *
     * @param  settings1 -
     * @param  anita1 -
     * @param  ilayer -
     * @param  ifold
     * @param  n_eplane -
     * @param  n_hplane -
     * @param  n_normal -
     * @return returns void
     */
    void GetAntennaOrientation(const Settings *settings1,
			       Anita *anita1, 
			       int ilayer, int ifold, 
			       Vector& n_eplane,
			       Vector& n_hplane, 
			       Vector& n_normal) const;

    

 
    ///< This function sets the default balloon position
    /**
     * Long description
     *
     *
     * @param  antarctica1 -
     * @return returns void
     */
    void SetDefaultBalloonPosition(const Antarctica *antarctica1);


    // ///< This function adjusts the slac balloon position
    // /**
    //  * move payload around like we did at slac
    //  *
    //  *
    //  * @param  inu -
    //  * @return returns void
    //  */
    // void AdjustSlacBalloonPosition(int inu);
    
    // ///< This function gets the slac balloon positions
    // /**
    //  * Long description
    //  *
    //  *
    //  * @param  anita1 -
    //  * @return returns void
    //  */
    // void GetSlacPositions(Anita *anita1);

    ///< This function gets boresights
    /**
     * Long description
     *
     *
     * @param  settings1 -
     * @param  anita1 -
     * @return returns void
     */
    void GetBoresights(const Settings *settings1, const Anita *anita1);
    
    ///< This function calculates antenna positions
    /**
     * Long description
     *
     *
     * @param  settings1 -
     * @param  anita1 -
     * @return returns void
     */
    void calculate_antenna_positions(const Settings *settings1,Anita *anita1) const ;
    
    ///< This function rotates the payload
    /**
     * Rotate from payload coord to earth coord
     *
     *
     * @param  ant_pos -
     * @return returns vector
     */
    Vector RotatePayload(Vector ant_pos) const;
    
    ///< This function UN-rotates the payload
    /**
     * Rotate from earth to payload coord. (undoes RotatePayload)
     *
     *
     * @param  ant_pos -
     * @return returns vector
     */
    Vector unRotatePayload(Vector ant_pos) const;

    inline FlightPath whichPath() const {return WHICHPATH;}

    inline const Position& position() const {return r_bn;}

    unsigned int realTime() const {return realTime_flightdata;}

    inline double getLatitude() const {return latitude;}
    inline double getLongitude() const {return longitude;}
    inline double getHeading() const {return heading;}
    double getPitch() const;
    double getRoll() const;

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
    int NPOINTS;                                                ///< number of GPS positions we're picking from.
    int REDUCEBALLOONPOSITIONS;                                 ///< only take every 100th entry in the flight data file
    double theta_bn;
    double phi_bn;                                              ///< theta,phi of balloon wrt south pole
    double altitude_bn;
    double dtryingposition;                                     ///< weighting factor: how many equivalent tries each neutrino counts for after having reduced possible interaction positions to within horizon
    Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX]; ///< position of antenna boresights
    std::vector<double> latitude_bn_anitalite;                  ///< latitude at times along flightpath, equally distributed among gps data. This is filled with anita or anita-lite data, depending on which the user specifies
    std::vector<double> longitude_bn_anitalite;                 ///< same for longitude
    std::vector<double> altitude_bn_anitalite;                  ///< same for altitude
    std::vector<double> heading_bn_anitalite;                   ///< same for heading of the balloon
    std::vector<double> realtime_bn_anitalite;                  ///< same for real life time

    float flatitude,flongitude,faltitude,fheading,froll, fpitch;

    double surface_under_balloon;                               ///< distance between center of the earth and the surface of earth under balloon        
    
  private:

    void InitializeBalloon();
    void ReadAnitaliteFlight();
    void setr_bn(double latitude,double longitude);
    
    const FlightPath WHICHPATH;                                 ///< 0=fixed balloon position,1=randomized,2=ANITA-lite GPS data,3=banana plot
    // GPS positions of Anita-lite balloon flight
    int igps;                                                   ///< which balloon position do we use out of the 25000 anitalite GPS positions.
    int ibnposition;
    // static const int MAX_POSITIONS=50;                          ///< for the slac beam test
    // Vector slacpositions[MAX_POSITIONS];
    // std::string sslacpositions[MAX_POSITIONS];
    // int islacposition;
    unsigned int realTime_flightdata_temp;                      ///< realtime from the flight data file
    unsigned int realTime_flightdata;                           ///< realtime from the flight data file
    double latitude,longitude,altitude,heading,roll,pitch;
    
    double MINALTITUDE;                                         ///< minimum altitude balloon needs to be before we consider it a good event to read from the flight data file
    int igps_previous;                                          ///< which entry from the flight data file the previous event was so we can just take the next one.
    int RANDOMIZE_BN_ORIENTATION;                               ///< 0=fixed balloon orientation,1=randomized
    float powerthresh[9][32];                                   ///< power threshold in Watts
    float meanp[9][32];                                         ///< mean power in Watts
    Position r_bn;                                              ///< position of balloon
    double horizcoord_bn;                                       ///< x component of balloon position
    double vertcoord_bn;                                        ///< y component of balloon position

    // Vector x_axis_rot;
    // Vector y_axis_rot;
    // Vector z_axis_rot;
    // Vector n_bn;                                                ///< normalized r_bn
    Vector n_east;                                              ///< east, as seen from the balloon position
    Vector n_north;                                             ///< north, as seen from the balloon position

    Position r_bn_shadow;                                       ///< position of the balloon projected on earth surface - point just below balloon at surface of the earth
    double phi_spin;                                            ///< orientation of the balloon

    int NPOINTS_MIN;                                            ///< min and max index for gps positions we want to include in the simulation (to exclude launch and fall).  These are set in ReadFlight
    int NPOINTS_MAX;

    

   
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


