////////////////////////////////////////////////////////////////////////////////////////////////
//class Tools:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <fstream>
#include <vector>

class Anita;
class Secondaries;
class Signal;
class Balloon;
class Ray;

using std::string;
using std::ifstream;
using std::ofstream;
using std::vector;

#include "TString.h"
#include <TObject.h>
#include <map>

// from RVersion.h
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
#if ROOT_VERSION_CODE < ROOT_VERSION(6,22,0)
#include "TClingRuntime.h"
#endif
#else
#include "TCint.h"
#endif

//! Reads in and stores input settings for the run

class Settings : public TObject {

  /* protected:  */

 public:

  Settings();
  ~Settings(); 
  void Initialize();
  void printAllKeyValuePairStrings();

  void getSetting(const char* key, int& value, bool nonag=false);
  void getSetting(const char* key, float& value, bool nonag=false);
  void getSetting(const char* key, double& value, bool nonag=false);

  void getSetting(const char * key, std::string & value, bool nonag=false); 
  void getSetting(const char* key, vector<int>& valueArray, bool nonag=false);
  void getSetting(const char* key, vector<float>& valueArray, bool nonag=false);
  void getSetting(const char* key, vector<double>& valueArray, bool nonag=false);

  void ReadInputs(const char* fileName , ofstream &foutput,
		  // Anita* anita1, Secondaries* sec1, Signal* sig1, Balloon* bn1, Ray* ray1,
		  int& NNU, double& RANDOMISEPOL);

  void ApplyInputs(Anita* anita1, Secondaries* sec1, Signal* sig1, Balloon* bn1, Ray* ray1);


  int UNBIASED_SELECTION;
  double UNBIASED_PS_MAX_DISTANCE_KM; 
  double UNBIASED_CHORD_STEP_M; 
  int WHICH; // which payload to use 0=Anita-lite,1=Ross,2=Smex,3=make your own
  int ANITAVERSION;
  int CYLINDRICALSYMMETRY; // is it cylindrically symmetric =1 if which=1,2, =0 if which=0
  // if which=3 then 0 or 1
  double SIGMA_FACTOR; // factor to multiply cross section by for error analysis
  int SIGMAPARAM; // 0=Reno, 1=Connolly et al. 2011 for cross section parametrization
  int YPARAM; // 0=Reno, 1=Connolly et al. 2011 for cross section parametrization
  int SIGNAL_FLUCT;  // 1=add noise fluctuation to signal or 0=do not
  int TRIGGERSCHEME;  // frequency domain voltage, frequency domain energy, time domain diode integration
  int ZEROSIGNAL;  // zero the signal to see how many of our hits are noise hits
  int REMOVEPOLARIZATION; //Disable polarizations

  double INCLINE_TOPTHREE;
  double INCLINE_NADIR;
  int USEDARTBOARD;
  int GAINS;
  int BANDING;
  int NBANDS;
  int PERCENTBW; 
  int trigRequirements[4];//  0th element - L1 - how many channels per antenna should pass
  // 1st element- L2 - how many antennas on a layer
  // 2nd element - L3 - how many L2 triggers should be coincident
  int REQUIRE_CENTRE; // require centre antenna in clump to be one of those hit
  double INCLUDE_NADIRONLY; // cant angle of nadir (bottom) layer of antennas
  int PULSER;
  double SIGMA_THETA; // resolution on the polar angle of the signal
  double FREQ_LOW;       ///< lowest frequency
  double FREQ_HIGH;
  
  int SECONDARIES;
  int TAUDECAY; // is tau decay counted as a secondary interaction

  int trigEffScanPhi;                      // central phi sector of trigger efficiency scan

  
  int WHICHPATH;
  int BN_LATITUDE;
  int BN_LONGITUDE;
  int BN_ALTITUDE;
  int RANDOMIZE_BN_ORIENTATION;
  int CENTER;                                                                ///< whether or not to center one phi sector of the payload on the incoming signal (for making signal efficiency curves)
  double MAXHORIZON;
  
  int EVENTSMAP;//whether draw the events distribution map

  int WHICHRAYS;  // how many rays to look at (1) direct only (2) direct and down-going.
  int MAKEVERTICAL; // option in the input file to force the signal to hit the payload with completely vertical polarisation.  For making signal efficiency curves.

  // trigger
  int LCPRCP; // 1 for circular polarization trigger, 0 for V and H
  int JUSTVPOL; // 0 for both polarizations, 1 for just V polarization
  // doesn't allow for both LCPRCP=1 and JUSTVPOL=1
  //int FIFTHBAND; // 1 to include 0.2-1.2 GHz as a frequency band if JUSTVPOL==1
  //int NFOLD=3;  // how many channels must pass the trigger - in old mechanism - only used for anita-lite
  int NFOLD;  // how many channels must pass the trigger - in old mechanism - only used for anita-lite


  //int CHMASKING=1; // whether or not to include channel masking
  //int PHIMASKING=1; // whether or not to include phi masking
  int CHMASKING; // whether or not to include channel masking
  int PHIMASKING; // whether or not to include phi masking

  //int NLAYERS=0;
  //int NANTENNAS=0;

  int NLAYERS;
  int NANTENNAS;

  /* int ONLYFINAL=1; // only write to final histogram */
  /* int HIST_MAX_ENTRIES=10000; //maximum number of events to put in histograms */
  /* int HIST=1;          //write to histograms  */

  int ONLYFINAL; // only write to final histogram
  int HIST_MAX_ENTRIES; //maximum number of events to put in histograms
  int HIST;          //write to histograms
  double BW; // BANDWIDTH
  //int DISCONES=1; // whether or not to use discones
  int DISCONES; // whether or not to use discones

  //double NDISCONES_PASS=3; // number of discones needed to pass
  double NDISCONES_PASS; // number of discones needed to pass

  int BORESIGHTS; // whether to loop over boresights
  int SLAC; // whether or not we are simulating the slac run
  double SLACSLOPE; // slope of the ice
  double SLACICELENGTH;  // length of the block of ice
  double SLAC_HORIZDIST; // horizontal distance from interaction to center of payload at slac beam test
  double SLAC_DEPTH; // vertical depth of interaction at slac beam test
  double SLAC_HORIZ_DEPTH; // horizontal depth of interaction at slac

  std::string SOURCE;  // the source option: FAVA (blazars), GRB (gamma ray bursts), SN (supernovae),  see source.hh for more info
  std::string WHICH_SOURCES; // which sources? All, or just a specific one
  std::string WHICH_SUBTYPE; // which subtype? All, or just a specific one
  std::string WHICH_START_TIME; // which start time? 0, or an anita flight option
  std::string WHICH_END_TIME; // which end time? 0 or an anita flight option
  
  int SOURCE_USE_EXPONENT; //Use the exponent for a source (if exponent <= 21) 
  double SOURCE_MIN_E;  // log10 of minimum energy for sources 
  double SOURCE_MAX_E;   // log10 of maximinum energy for sources 

  int SOURCE_SKIP_WHEN_NONE; //Whether or not to reroll position if no source is available! 
  
  int ROUGHNESS; // include effects of surface roughness
  int FIRN; // whether or not to include the firn

  //int SLOPEY=1; // 1=slopeyness on, 0=slopeyness off
  //double SLOPEYSIZE=0.012; // This determines size of the slopeyness (0.10=5.4, 0.20=7.4 deg mean)

  int SLOPEY; // 1=slopeyness on, 0=slopeyness off
  double SLOPEYSIZE; // This determines size of the slopeyness (0.10=5.4, 0.20=7.4 deg mean)

  bool DEBUG;

  //double THERMALNOISE_FACTOR=1.0; // factor to multiply thermal noise for error analysis
  double THERMALNOISE_FACTOR; // factor to multiply thermal noise for error analysis

  //double FREQ_LOW_SEAVEYS=200.E6; // min frequency for seaveys
  //const double FREQ_HIGH_SEAVEYS=1200.E6; // max frequency for seaveys

  double FREQ_LOW_SEAVEYS; // min frequency for seaveys
  double FREQ_HIGH_SEAVEYS; // max frequency for seaveys
  double BW_SEAVEYS;
  //int FORSECKEL=1; // Make array of strength of signal across frequencies for different viewing angles.
  int FORSECKEL; // Make array of strength of signal across frequencies for different viewing angles.

  double ROUGHSIZE; // roughness size
  double SCREENEDGELENGTH;        // edge length of screen used if there is roughness
  double SCREENSTEPSIZE;        // step size of screen grid if there is roughness

  int ICE_MODEL; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP.
  int NOFZ; // 1=depth dependent index of refraction,0=off
  int CONSTANTCRUST; // set crust density and thickness to constant values.
  int CONSTANTICETHICKNESS; // set ice thickness to constant value
  int FIXEDELEVATION; // fix the elevation to the thickness of ice.
  int MOOREBAY; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data
  int USEPOSITIONWEIGHTS;// whether or not to restrict the neutrino position so it is within the horizon of the balloon
  int WRITE_FILE; //Select whether or not to write a new input file for CreateHorizons

  int MINRAY;
  int MAXRAY;

  int horizontal_banana_points;
  int vertical_banana_points;
  double EXPONENT; //Select neutrino flux exponent value or flux model. Detail : READ_EXPONENT


  // Bunch of variables which were global in icemc.cc but are settings:
  int FILLRAYTREES; // fill tree for each ray in roughness simulation
  int SEED;      // random number seed.
  double THETA_TH_FACTOR; // factor to multiply theta_th to check code is working properly
  double CHANCEINHELL_FACTOR; // loosen chance in hell cuts to check code is working properly
  int WEIGHTABSORPTION; // whether or not to weight for earth absorption
  int CONSTANTY; // whether or not to set y to a constant=0.2
  int taumodes; //whether to choose a flat distribution for y  (tau made in rock) or not
  int VARIABLE_ATTEN; //  0=depth dependent attenuation length, 1=fixed
  int TRIGTYPE; //1=Trigger scheme as in the SMEX proposal  or 0= Just at least 8 channels pass with 2.3 sigma signal
  int ATMOSPHERE;// include atmosphere
  int SCALEDOWNLCPRX1; // scale down lcp voltage of antenna 1 by sqrt(2)
  int SCALEDOWNEPOLRX1; // scale down power of e pol. of antenna 1 by factor of 2
  int SCALEDOWNHPOLRX1; // scale down power of h pol. of antenna 1 by factor of 2
  int SCALEDOWNEPOLRX2; // scale down power of e pol. of antenna 2 by factor of user's choice
  double SCALEFACTOREPOLRX2; // scale power of e pol. of antenna 2 by this factor
  int SCALEDOWNHPOLRX2; // scale down power of h pol. of antenna 2 by factor of 2
  int EPOLRX2ZERO; // lcp channel on anita-lite is not considered for triggering.
  int HPOLRX2ZERO; // h pol. of 2nd antenna set to zero.
  int RCPRX2ZERO; // rcp of 2nd antenna set to zero.
  int LCPRX2ZERO; // lcp of 2nd antenna set to zero.
  int FLATSURFACE; // Normals of all positions on the surface are straight up.
  int WRITEPOSFILE; //Write neutrino position information to file
  int SKIPCUTS; //See every neutrino through to the end - don't make any of the various cuts designed to speed up the program.  (For checking distributions.)
  int USEDIRECTIONWEIGHTS;// whether or not to restrict the neutrino angle so that the viewing angle is near the cerenkov cone
  int SHOWERTYPE; // Type of shower for previous option
  int antennaclump; //number of antenna in clump (L2)
  // End of the once-global varibles.
  double COHERENT_THRESHOLD;
  int APPLYIMPULSERESPONSEDIGITIZER;       // apply impulse response in the digitizer path
  int APPLYIMPULSERESPONSETRIGGER;         // apply impulse response in the trigger path
  int USETIMEDEPENDENTTHRESHOLDS;          // use time-dependent thresholds
  int USEDEADTIME;                         // use dead time from flight
  int NOISEFROMFLIGHTTRIGGER;              // use thermal noise from flight in trigger path
  int NOISEFROMFLIGHTDIGITIZER;            // use thermal noise from flight in digitizer path
  int MINBIAS;                             // generate minimum bias sample
  int TRIGGEREFFSCAN;                      // do a trigger efficiency scan
  int TRIGGEREFFSCAPULSE;                  // Apply pulse at AMPA (0) or at SURF (1)

  int TUFFSTATUS;                             // Are the TUFFs on for the whole flight?

  int ADDCW;                               // Add CW
  
  int PAYLOAD_USE_SPECIFIC_TIME;           //Instead of using the entire flight path, only generate neutrinos for a specific time for the paylaod (0 to disable). 
  int PAYLOAD_USE_SPECIFIC_TIME_DELTA;     //How much before and after the specific time can we use payload locations? 
  int SPECIFIC_NU_POSITION;                //Use a specific interaction position 
  double SPECIFIC_NU_POSITION_LATITUDE, SPECIFIC_NU_POSITION_LONGITUDE, SPECIFIC_NU_POSITION_ALTITUDE; //the specific interaction position 
  double SPECIFIC_NU_POSITION_DISTANCE; //Max distance from place
  int IGNORE_CROSSPOL; //Ignore the crosspol polarization component
  int POL_SIGN_HACK; // patch up the sign of e/h 
  double CUTONWEIGHTS;
  double CUTONWEIGHTPROBS;
  double DEC_CUT;
  int ALL_SKY_MAP;
  int WRITE_WAVEFORMS; 
  
  // custom sources
  std::string CUSTOM_NAME;
  double CUSTOM_RA; // in decimal degrees
  double CUSTOM_DEC; // in decimal degrees
  double CUSTOM_GAMMA;
  
  double HORIZON_OFFSET; 
  int useLPM;

  // In-header intialization is to old gcc as Domino's pizza is to real Italians
  double jamieFactor;// = 0;
  int medium;// = 0;
  int askaryanParameterization;// = 0;
  int SAVE_TRUTH_NU_TREE; 
    
  //  TString outputdir; // directory where outputs go

  ClassDef(Settings,2);
  
 private:
  typedef std::map<TString, TString> kvpMap;

  kvpMap keyValuePairStrings; //< The raw key value pairs as string, from parsing the config file
  Bool_t newKvpPassesSanityChecks(const TString& key, const TString& value, const char* fileName, int lineNum);
  void complainAboutNotFindingKey(const TString& key);
  void parseValueArray(const char* valueString, vector<int>& values);
  void parseValueArray(const char* valueString, vector<float>& values);
  void parseValueArray(const char* valueString, vector<double>& values);
  void parseSettingsFile(const char* fileName, ofstream& outputFile);

  vector<double> efficiencyScanOffAxisAttenuations;
  vector<double> efficiencyScanPhiSectorDelay;
  vector<double> efficiencyScanRingDelay;
  vector<int> efficiencyScanRingsUsed;
  vector<int> efficiencyScanApplyRingDelay;
  vector<int> whichTUFFsON;
  vector<double> tempThresholds;  
  vector<double> bandLowEdgesMHz;
  vector<double> bandHighEdgesMHz;
  vector<int> requiredBands;
  vector<int> allowedBands;
  vector<double> notchFilterLimitsMHz;
  vector<int> channelRequirePol;
  vector<int> channelAllowedPol;

    
};
#endif
