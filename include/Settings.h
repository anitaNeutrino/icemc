#ifndef ICEMC_SETTINGS_H
#define ICEMC_SETTINGS_H

#include <fstream>
#include <vector>

#include "TString.h"
#include <TObject.h>
#include <map>

// from RVersion.h
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
#include "TClingRuntime.h"
#else
#include "TCint.h"
#endif

class TNamed;

namespace icemc {
  
  class ShowerGenerator;
  class AskaryanRadiationModel;
  class RayTracer;

  /**
   * @class Settings
   * @brief Reads in and stores input settings for the run
   */

  class Settings : public TObject {

    friend class CommandLineOptions;

  public:    

    Settings();
    virtual ~Settings();
    void Initialize();
    void printAllKeyValuePairStrings() const;

    void getSetting(const char* key, int& value) const;
    void getSetting(const char* key, float& value) const;
    void getSetting(const char* key, double& value) const;
    void getSetting(const char* key, std::string& value) const;

    void getSetting(const char* key, std::vector<int>& valueArray) const;
    void getSetting(const char* key, std::vector<float>& valueArray) const;
    void getSetting(const char* key, std::vector<double>& valueArray) const;
    void getSetting(const char* key, std::vector<std::string>& valueArray) const;

    virtual void ReadInputs(const char* fileName , std::ofstream &foutput);
		    // Anita* anita1, ShowerGenerator* sec1, AskaryanRadiationModel* askFreqGen, Balloon* bn1, Ray* ray1,
		    // int& NNU, double& RANDOMISEPOL);

    // void ApplyInputs(Anita* anita1) const;

    const char* getOutputDir() const {return fOutputDir.c_str();}
    int getRun() const {return fRun;}
    int getStartNu() const {return fStartNu;}

    int NNU; ///< The number of neutrinos
    double RANDOMISEPOL; ///< Randomize the polarity?
    
    int UNBIASED_SELECTION;
    double SIGMA_FACTOR; // factor to multiply cross section by for error analysis
    int SIGMAPARAM; // 0=Reno, 1=Connolly et al. 2011 for cross section parametrization
    int YPARAM; // 0=Reno, 1=Connolly et al. 2011 for cross section parametrization
    int SIGNAL_FLUCT;  // 1=add noise fluctuation to signal or 0=do not
    int ZEROSIGNAL;  // zero the signal to see how many of our hits are noise hits
    int REMOVEPOLARIZATION; //Disable polarizations

    int USEDARTBOARD;
    int PERCENTBW; 
    int PULSER;
  
    int SECONDARIES;
    int TAUDECAY; // is tau decay counted as a secondary interaction


  
  
    int EVENTSMAP;//whether draw the events distribution map

    int WHICHRAYS;  // how many rays to look at (1) direct only (2) direct and down-going.
    int MAKEVERTICAL; // option in the input file to force the signal to hit the payload with completely vertical polarisation.  For making signal efficiency curves.

    // trigger

    //int NLAYERS=0;
    //int NANTENNAS=0;


    /* int ONLYFINAL=1; // only write to final histogram */
    /* int HIST_MAX_ENTRIES=10000; //maximum number of events to put in histograms */
    /* int HIST=1;          //write to histograms  */

    int ONLYFINAL; // only write to final histogram
    int HIST_MAX_ENTRIES; //maximum number of events to put in histograms
    int HIST;          //write to histograms
    double BW; // BANDWIDTH

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
    int USEPOSITIONWEIGHTS;// whether or not to restrict the neutrino position so it is within the horizon of the detector
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
    int MINBIAS;                             // generate minimum bias sample
    int PAYLOAD_USE_SPECIFIC_TIME;           //Instead of using the entire flight path, only generate neutrinos for a specific time for the paylaod (0 to disable). 
    int PAYLOAD_USE_SPECIFIC_TIME_DELTA;     //How much before and after the specific time can we use payload locations? 
    int SPECIFIC_NU_POSITION;                //Use a specific interaction position 
    double SPECIFIC_NU_POSITION_LATITUDE, SPECIFIC_NU_POSITION_LONGITUDE, SPECIFIC_NU_POSITION_ALTITUDE; //the specific interaction position 
    double SPECIFIC_NU_POSITION_DISTANCE; //Max distance from place

    int useLPM;

    // In-header intialization is to old gcc as Domino's pizza is to real Italians
    double jamieFactor;// = 0;
    int medium;// = 0;
    int askaryanParameterization;// = 0;
    
    //  TString outputdir; // directory where outputs go

    TNamed* makeRootSaveableSettings() const;

    ClassDef(Settings,2);
  
  private:

    
    typedef std::map<TString, TString> kvpMap;

    kvpMap keyValuePairStrings; //< The raw key value pairs as string, from parsing the config file
    Bool_t newKvpPassesSanityChecks(const TString& key, const TString& value, const char* fileName, int lineNum) const;
    void complainAboutNotFindingKey(const TString& key) const;
    void parseValueArray(const char* valueString, std::vector<int>& values) const;
    void parseValueArray(const char* valueString, std::vector<float>& values) const;
    void parseValueArray(const char* valueString, std::vector<double>& values) const;
    void parseSettingsFile(const char* fileName, std::ofstream& outputFile);
    void processStrings(const std::string& raw, std::vector<std::string>& processed) const;
    void processLine(const std::string& thisLine, std::ofstream& outputFile, const char* fileName, int lineNum);


    TString wholeSettingsFile;
    std::string fOutputDir;
    int fRun;
    int fStartNu;
  };
}



#endif
