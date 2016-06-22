#ifndef BALLOON_H
#define BALLOON_H


////////////////////////////////////////////////////////////////////////////////////////////////
//class Balloon:
////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>



class Position;
class Vector;
class TChain;
class TTree;
class TFile;
class TTreeIndex;
class IceModel;
class Anita;
class Settings;
class TRandom3;
class Interaction;
class Ray;
class TH1F;
class TGraph;

using std::string;
using std::cout;

//! Handles everything related to balloon positions, payload orientation over the course of a flight.
class Balloon {

private:

  //  string anitaliteflight; // the gps path of the anita-lite flight
  //string anitaflight;// gps path of anita flight
 
  TRandom3 Rand3;
public:
  Balloon();

  void setObservationLocation(Interaction *interaction1,int inu,IceModel *antarctic,Settings *settings1);
  void GetBoresights(Settings *settings1,Anita *anita1,Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX]);
  // GPS positions of Anita-lite balloon flight
  int igps; // which balloon position do we use out of the 25000 anitalite GPS positions.
  int ibnposition;
  double dtryingposition; // weighting factor: how many equivalent tries each neutrino counts for after having reduced possible interaction positions to within horizon
  
  void PickDownwardInteractionPoint(Interaction *interaction1,Anita *anita1,Settings *settings1,IceModel *antarctica1,
				    Ray *ray1, int &beyondhorizon); 
  
  void InitializeBalloon();
  void ReadAnitaliteFlight();
  void CenterPayload(double& hitangle_e);

  void PickBalloonPosition(Vector straightup,IceModel *antarctica,Settings *settings1, Anita *anita1);
  void PickBalloonPosition(IceModel *antarctica1,Settings *settings1,int inu,Anita *anita1, double randomNumber); // position of spot under balloon

  int Getibnposition();

  double GetBalloonSpin(double heading); // get the azimuth of the balloon

  void GetAntennaOrientation(Settings *settings1, 
			     Anita *anita1, 
			     int ilayer, int ifold, 
			     Vector& n_eplane, 
			     Vector& n_hplane, 
			     Vector& n_normal); 


  void GetEcompHcompEvector(Settings *settings1, 
			    Vector n_eplane, 
			    Vector n_hplane, 
			    const Vector n_pol, 
			    double& e_component, 
			    double& h_component, 
			    double& n_component);

		
  void GetEcompHcompkvector(Vector n_eplane, 
			    Vector n_hplane, 
			    Vector n_normal, 
			    const Vector n_exit2bn, 
			    double& e_component_kvector, 
			    double& h_component_kvector, 
			    double& n_component_kvector);

  void GetHitAngles(
		    double e_component_kvector,
		    double h_component_kvector,
		    double n_component_kvector, 
		    double& hitangle_e,
		    double& hitangle_h); 
 
  void SetDefaultBalloonPosition(IceModel *antarctica1);

  void setr_bn(double latitude,double longitude);

  void AdjustSlacBalloonPosition(int inu); // move payload around like we did at slac
  void GetSlacPositions(Anita *anita1); 




  // for the slac beam test
  static const int MAX_POSITIONS=50;
  Vector slacpositions[MAX_POSITIONS];
  string sslacpositions[MAX_POSITIONS];
  int islacposition;

  TChain *flightdatachain;

  TTreeIndex *tindex;

  unsigned int realTime_flightdata_temp; // realtime from the flight data file
  unsigned int realTime_flightdata; // realtime from the flight data file
   
  float flatitude,flongitude,faltitude,fheading,froll, fpitch;
  double latitude,longitude,altitude,heading,roll,pitch;
  double MINALTITUDE; // minimum altitude balloon needs to be before we consider it a good event to read from the flight data file
  int igps_previous;  // which entry from the flight data file the previous event was so we can just take the next one.
  int REDUCEBALLOONPOSITIONS; // only take every 100th entry in the flight data file
  int WHICHPATH; // 0=fixed balloon position,1=randomized,2=ANITA-lite GPS data,3=banana plot
  int RANDOMIZE_BN_ORIENTATION;// 0=fixed balloon orientation,1=randomized
  double BN_ALTITUDE; // pick balloon altitude
  unsigned short surfTrigBandMask[9][2]; // Ryan's 16 bit masks for 9 surfs.  2x16 bit masks gives 32 channels per surf
  float powerthresh[9][32]; // power threshold in Watts
  float meanp[9][32]; // mean power in Watts

  int CENTER; // whether or not to center one phi sector of the payload on the incoming signal (for making signal efficiency curves)

  double altitude_bn;
  double theta_bn;
  double phi_bn;// theta,phi of balloon wrt south pole
  Position r_bn; // position of balloon
  double horizcoord_bn; // x component of balloon position
  double vertcoord_bn; // y component of balloon position
  Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //position of antenna boresights

  void GetBoresights(Settings *settings1,Anita *anita1);
  Vector RotatePayload(Vector ant_pos);
  void calculate_antenna_positions(Settings *settings1,Anita *anita1);// this calculates the above 


  Vector n_bn; // normalized r_bn
  Vector n_east; // east, as seen from the balloon position
  Vector n_north; // north, as seen from the balloon position
  double surface_under_balloon; // distance between center of the earth and the surface of earth under balloon
  Position r_bn_shadow;//position of the balloon projected on earth surface - point just below balloon at surface of the earth
  double MAXHORIZON;                // pick the interaction within this distance from the balloon so that it is within the horizon
  double phi_spin; // orientation of the balloon

  int NPOINTS; // number of GPS positions we're picking from.
  int NPOINTS_MIN; // min and max index for gps positions we want to include in the simulation (to exclude launch and fall).  These are set in ReadFlight
  int NPOINTS_MAX;

  double latitude_bn_anitalite[100000]; // latitude at times along flightpath, equally distributed among gps data. This is filled with anita or anita-lite data, depending on which the user specifies
  double longitude_bn_anitalite[100000]; // same for longitude
  double altitude_bn_anitalite[100000]; // same for altitude
  double heading_bn_anitalite[100000]; // same for heading of the balloon
  double realtime_bn_anitalite[100000]; // same for real life time

    
  double BN_LONGITUDE; //balloon longitude for fixed balloon location
  double BN_LATITUDE; //balloon latitude for fixed balloon location
}; //class Balloon




const string anitaliteflight="data/BalloonGPS.txt"; // the gps path of the anita-lite flight
const string anitaflight="data/anitagps.txt";// gps path of anita flight


#endif


