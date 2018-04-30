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
  class Ray;
  class Interaction;


  /**
   * @class BalloonInfo
   * @brief Where is the payload? A very minimalistic  class for RootOutput::allTree
   * 
   * Just like an Adu5Pat, but if this is stand alone icemc, that won't exist
   */
  class BalloonInfo {
  public:
    BalloonInfo();
    UInt_t realTime;    
    float heading;
    float pitch;
    float roll;
    float latitude;
    float longitude;
    float altitude;
    float horizcoord_bn;
    float vertcoord_bn;

    #ifdef ANITA_UTIL_EXISTS
    operator Adu5Pat() const; /// Type conversion to a Adu5Pat
    #endif

    ClassDef(BalloonInfo, 1);    
  };
  

  ///< Handles everything related to balloon positions, payload orientation over the course of a flight.
  class Balloon {    

  private:
    //  std::string anitaliteflight; // the gps path of the anita-lite flight
    //std::string anitaflight;// gps path of anita flight
 
  public:
    Balloon();

  
    // GPS positions of Anita-lite balloon flight
    int igps;                                                   ///< which balloon position do we use out of the 25000 anitalite GPS positions.
    int ibnposition;
    double dtryingposition;                                     ///< weighting factor: how many equivalent tries each neutrino counts for after having reduced possible interaction positions to within horizon
    static const int MAX_POSITIONS=50;                          ///< for the slac beam test
    Vector slacpositions[MAX_POSITIONS];
    std::string sslacpositions[MAX_POSITIONS];
    int islacposition;
    TChain *flightdatachain;
    TTreeIndex *tindex;
    unsigned int realTime_flightdata_temp;                      ///< realtime from the flight data file
    unsigned int realTime_flightdata;                           ///< realtime from the flight data file
    float flatitude,flongitude,faltitude,fheading,froll, fpitch;
    double latitude,longitude,altitude,heading,roll,pitch;
    
    double MINALTITUDE;                                         ///< minimum altitude balloon needs to be before we consider it a good event to read from the flight data file
    int igps_previous;                                          ///< which entry from the flight data file the previous event was so we can just take the next one.
    int REDUCEBALLOONPOSITIONS;                                 ///< only take every 100th entry in the flight data file
    int WHICHPATH;                                              ///< 0=fixed balloon position,1=randomized,2=ANITA-lite GPS data,3=banana plot
    int RANDOMIZE_BN_ORIENTATION;                               ///< 0=fixed balloon orientation,1=randomized
    double BN_ALTITUDE;                                         ///< pick balloon altitude
    unsigned short surfTrigBandMask[9][2];                      ///< Ryan's 16 bit masks for 9 surfs.  2x16 bit masks gives 32 channels per surf
    float powerthresh[9][32];                                   ///< power threshold in Watts
    float meanp[9][32];                                         ///< mean power in Watts
    double altitude_bn;
    double theta_bn;
    double phi_bn;                                              ///< theta,phi of balloon wrt south pole
    Position r_bn;                                              ///< position of balloon
    double horizcoord_bn;                                       ///< x component of balloon position
    double vertcoord_bn;                                        ///< y component of balloon position
    Position r_boresights[icemc::Anita::NLAYERS_MAX][icemc::Anita::NPHI_MAX]; ///< position of antenna boresights
    Vector x_axis_rot;
    Vector y_axis_rot;
    Vector z_axis_rot;
    Vector n_bn;                                                ///< normalized r_bn
    Vector n_east;                                              ///< east, as seen from the balloon position
    Vector n_north;                                             ///< north, as seen from the balloon position
    double surface_under_balloon;                               ///< distance between center of the earth and the surface of earth under balloon
    Position r_bn_shadow;                                       ///< position of the balloon projected on earth surface - point just below balloon at surface of the earth
    double MAXHORIZON;                                          ///< pick the interaction within this distance from the balloon so that it is within the horizon
    double phi_spin;                                            ///< orientation of the balloon
    int NPOINTS;                                                ///< number of GPS positions we're picking from.
    int NPOINTS_MIN;                                            ///< min and max index for gps positions we want to include in the simulation (to exclude launch and fall).  These are set in ReadFlight
    int NPOINTS_MAX;
    double latitude_bn_anitalite[100000];                       ///< latitude at times along flightpath, equally distributed among gps data. This is filled with anita or anita-lite data, depending on which the user specifies
    double longitude_bn_anitalite[100000];                      ///< same for longitude
    double altitude_bn_anitalite[100000];                       ///< same for altitude
    double heading_bn_anitalite[100000];                        ///< same for heading of the balloon
    double realtime_bn_anitalite[100000];                       ///< same for real life time
    double BN_LONGITUDE;                                        ///< balloon longitude for fixed balloon location
    double BN_LATITUDE;                                         ///< balloon latitude for fixed balloon location


    ///<This function sets the observation location
    /**
     * This is a long description that I dont know yet
     *
     *
     * @param  interaction1 -
     * @param  inu -
     * @param  antarctic -
     * @param  settings1 -
     * @return returns void
     */
    void setObservationLocation(Interaction *interaction1,int inu,IceModel *antarctic,const Settings *settings1);
    
    ///< This function gets the boresights
    /**
     * This is a long description that I dont know yet
     *
     *
     * @param  settings1 -
     * @param  anita1 -
     * @param  r_boresights - [NLAYERS_MAX][NPHI_MAX] 
     * @return returns void
     */
    // void GetBoresights(const Settings *settings1,Anita *anita1,Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX]);
    void GetBoresights(const Settings *settings1, Anita *anita1, Position r_bn, double phi_spin, Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX]);
    
    ///< This function picks downward interaction point
    /**
     * This is a long description that I dont know yet
     *
     *
     * @param  interaction1 -
     * @param  anita1 -
     * @param  settings1 -
     * @param  antarctica1 -
     * @param  ray1 -
     * @param  beyondhorizon -
     * @return returns void
     */
    void PickDownwardInteractionPoint(Interaction *interaction1,Anita *anita1,const Settings *settings1,IceModel *antarctica1,
				      Ray *ray1, int &beyondhorizon); 
  
    ///< This function initializes the balloon or the specific flight
    /**
     * This is a long description that I dont know yet
     *
     *
     * @return returns void
     */
    void InitializeBalloon();
    
    ///< This function reads in the ANITA LITE flight
    /**
     * ANITA Lite is the first prototype ANITA flight
     *
     *
     * @return returns void
     */
    void ReadAnitaliteFlight();
    
    ///< This function centers the payload
    /**
     * Long description
     *
     *
     * @param  hitangle_e -
     * @return returns void
     */
    void CenterPayload(double& hitangle_e);

    
    ///< This function picks the balloon position
    /**
     * Long description
     *
     *
     * @param  straightup -
     * @param  antarctica -
     * @param  settings1 -
     * @param  anita1 -
     * @return returns void
     */
    void PickBalloonPosition(Vector straightup,IceModel *antarctica,const Settings *settings1, Anita *anita1);
    
    ///< This function picks the balloon position
    /**
     * Position of spot under balloon
     *
     *
     * @param  antarctica1 -
     * @param  settings1 -
     * @param  inu -
     * @param  anita1 -
     * @param  randomNumber -
     * @return returns void
     */
    /// @todo Default NULL ptr for BalloonInfo won't be required when we're not supporting all the silly "copy programs" and just have simple programs calling class functions
    void PickBalloonPosition(IceModel *antarctica1,const Settings *settings1,int inu,Anita *anita1, double randomNumber, BalloonInfo* bi = NULL);

    ///< This function gets ith balloon position
    /**
     * Long description
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
    double GetBalloonSpin(double heading);

    
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
			       Vector& n_normal); 

    ///< This function gets the e-component, h-component and e-vector
    /**
     * Long description
     *
     *
     * @param  settings1 -
     * @param  n_eplane -
     * @param  n_hplane -
     * @param  n_pol
     * @param  e_component -
     * @param  h_component -
     * @param  n_component -
     * @return returns void
     */
    void GetEcompHcompEvector(const Settings *settings1,
			      Vector n_eplane, 
			      Vector n_hplane, 
			      const Vector n_pol, 
			      double& e_component, 
			      double& h_component, 
			      double& n_component);
		
    ///< This function gets the e-component, h-component and k-vector
    /**
     * Long description
     *
     *
     * @param  n_eplane -
     * @param  n_hplane -
     * @param  n_normal -
     * @param  n_exit2bn -
     * @param  e_component_kvector -
     * @param  h_component_kvector - 
     * @param  n_component_kvector -
     * @return returns void
     */
    void GetEcompHcompkvector(Vector n_eplane,
			      Vector n_hplane, 
			      Vector n_normal, 
			      const Vector n_exit2bn, 
			      double& e_component_kvector, 
			      double& h_component_kvector, 
			      double& n_component_kvector);
    
    ///< This function gets the hit angles
    void GetHitAngles(
		      double e_component_kvector,
		      double h_component_kvector,
		      double n_component_kvector, 
		      double& hitangle_e,
		      double& hitangle_h); 
 
    ///< This function sets the default balloon position
    /**
     * Long description
     *
     *
     * @param  antarctica1 -
     * @return returns void
     */
    void SetDefaultBalloonPosition(IceModel *antarctica1);

    ///< This function sets r of the balloon
    /**
     * Long description
     *
     *
     * @param  latitude -
     * @param  longitude -
     * @return returns void
     */
    void setr_bn(double latitude,double longitude);

    ///< This function adjusts the slac balloon position
    /**
     * move payload around like we did at slac
     *
     *
     * @param  inu -
     * @return returns void
     */
    void AdjustSlacBalloonPosition(int inu);
    
    ///< This function gets the slac balloon positions
    /**
     * Long description
     *
     *
     * @param  anita1 -
     * @return returns void
     */
    void GetSlacPositions(Anita *anita1);

    ///< This function gets boresights
    /**
     * Long description
     *
     *
     * @param  settings1 -
     * @param  anita1 -
     * @return returns void
     */
    void GetBoresights(const Settings *settings1,Anita *anita1);
    
    ///< This function calculates antenna positions
    /**
     * Long description
     *
     *
     * @param  settings1 -
     * @param  anita1 -
     * @return returns void
     */
    void calculate_antenna_positions(const Settings *settings1,Anita *anita1);
    
    ///< This function rotates the payload
    /**
     * Rotate from payload coord to earth coord
     *
     *
     * @param  ant_pos -
     * @return returns vector
     */
    Vector RotatePayload(Vector ant_pos);
    
    ///< This function UN-rotates the payload
    /**
     * Rotate from earth to payload coord. (undoes RotatePayload)
     *
     *
     * @param  ant_pos -
     * @return returns vector
     */
    Vector unRotatePayload(Vector ant_pos);



    /** 
     * Pass on the balloon info for storage in allTree, 
     * to be used after the Balloon position has been picked.
     * @see PickBalloonPosition
     * 
     * @param bi 
     */
    void fillBalloonInfo(BalloonInfo& bi) const ;
   
  }; //class Balloon
}

#endif


