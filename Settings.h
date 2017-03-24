////////////////////////////////////////////////////////////////////////////////////////////////
//class Tools:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <fstream>

class Anita;
class Secondaries;
class Signal;
class Balloon;
class Ray;

using std::string;
using std::ifstream;
using std::ofstream;

#include "TString.h"
#include <map>


//! Reads in and stores input settings for the run


class Settings {

protected:


public:

  Settings();
  void Initialize();
  void printAllKeyValuePairStrings();

  void getSetting(const char* key, int& value);
  void getSetting(const char* key, float& value);
  void getSetting(const char* key, double& value);

  void getSetting(const char* key, std::vector<int>& valueArray);
  void getSetting(const char* key, std::vector<float>& valueArray);
  void getSetting(const char* key, std::vector<double>& valueArray);

  void ReadInputs(const char* fileName , ofstream &foutput, Anita* anita1, Secondaries* sec1, Signal* sig1, Balloon* bn1, Ray* ray1, int& NNU, double& RANDOMISEPOL);



  int UNBIASED_SELECTION;
  int WHICH; // which payload to use 0=Anita-lite,1=Ross,2=Smex,3=make your own
  int CYLINDRICALSYMMETRY; // is it cylindrically symmetric =1 if which=1,2, =0 if which=0
  // if which=3 then 0 or 1
  double SIGMA_FACTOR; // factor to multiply cross section by for error analysis
  int SIGMAPARAM; // 0=Reno, 1=Connolly et al. 2011 for cross section parametrization
  int YPARAM; // 0=Reno, 1=Connolly et al. 2011 for cross section parametrization
  int SIGNAL_FLUCT;  // 1=add noise fluctuation to signal or 0=do not
  int TRIGGERSCHEME;  // frequency domain voltage, frequency domain energy, time domain diode integration
  int ZEROSIGNAL;  // zero the signal to see how many of our hits are noise hits
  int REMOVEPOLARIZATION; //Disable polarizations

  int EVENTSMAP;//whether draw the events distribution map

  int WHICHRAYS;  // how many rays to look at (1) direct only (2) direct and down-going.

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

  int ROUGHNESS; // include effects of surface roughness
  int FIRN; // whether or not to include the firn

  //int SLOPEY=1; // 1=slopeyness on, 0=slopeyness off
  //double SLOPEYSIZE=0.012; // This determines size of the slopeyness (0.10=5.4, 0.20=7.4 deg mean)

  int SLOPEY; // 1=slopeyness on, 0=slopeyness off
  double SLOPEYSIZE; // This determines size of the slopeyness (0.10=5.4, 0.20=7.4 deg mean)

  bool DEBUG;
  string outputdir; // directory where outputs go

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
  double ROUGH_INTPOS_SHIFT;      // furthest distance to shift the neutrino interaction position from the balloon if roughness
  int ROUGHSCREENDIV_BASE;        // (N x N) grid for the base screen (to preselect)
  int ROUGHSCREENDIV_SUB;         // (n x n) subgrids for the preselected regions
  int ROUGHSCREENFRAC_BASE;       // fraction threshold for ratio of minimum Efield to maximum Efield magnitude for the base screen
  int ROUGHSCREENFRAC_SUB;        // fraction threshold for ratio of minimum Efield to maximum Efield magnitude for the subgrids
  int ROUGHMAXGEN;                // number of maximum generations (inclusive)

  /* int ICE_MODEL=0; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP. */
  /* int NOFZ=1; // 1=depth dependent index of refraction,0=off */
  /* int CONSTANTCRUST=0; // set crust density and thickness to constant values. */
  /* int CONSTANTICETHICKNESS=0; // set ice thickness to constant value */
  /* int FIXEDELEVATION=0; // fix the elevation to the thickness of ice. */
  /* int MOOREBAY=0; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data */
  /* int USEPOSITIONWEIGHTS=1;// whether or not to restrict the neutrino position so it is within the horizon of the balloon */
  /* int WRITE_FILE=0; //Select whether or not to write a new input file for CreateHorizons */

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
  int USETIMEDEPENDENTTHRESHOLDS; // use time-dependent thresholds
  int NOISEFROMFLIGHTTRIGGER;            // use thermal noise from flight in trigger path
  int NOISEFROMFLIGHTDIGITIZER;          // use thermal noise from flight in digitizer path
  int MINBIAS;                    // generate minimum bias sample
  int TRIGGEREFFSCAN;                      // do a trigger efficiency scan
  int TRIGGEREFFSCAPULSE;                  // Apply pulse at AMPA (0) or at SURF (1)

private:
  typedef std::map<TString, TString> kvpMap;

  kvpMap keyValuePairStrings; //< The raw key value pairs as string, from parsing the config file
  Bool_t newKvpPassesSanityChecks(const TString& key, const TString& value, const char* fileName, int lineNum);
  void complainAboutNotFindingKey(const TString& key);
  void parseValueArray(const char* valueString, std::vector<int>& values);
  void parseValueArray(const char* valueString, std::vector<float>& values);
  void parseValueArray(const char* valueString, std::vector<double>& values);
  void parseSettingsFile(const char* fileName, std::ofstream& outputFile);

};
#endif
