#include "TF1.h"
#include <array>
#include "position.hh"
#include "Constants.h"
#include "Settings.h"
#include "Tools.h"
#include "vector.hh"
#include "roughness.hh"
#include "anita.hh"
#include "balloon.hh"
#include "icemodel.hh"
#include "Spectra.h"
#include "signal.hh"
#include "secondaries.hh"
#include "ray.hh"
#include "counting.hh"
#include "Primaries.h"


#include <string.h>

#include "TString.h"
#include "TRegexp.h"
#include "TObjString.h"



// Prettify warnings because, why not?
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

ClassImp(Settings);

/**
 * Default constructor
 *
 */
Settings::Settings() : jamieFactor(0), medium(0), askaryanParameterization(0)
{
  Initialize();
}


/**
 * Default destructor
 *
 */
Settings::~Settings() {

}






/**
 * Copies the contents of a yaml settings file into internal memory as strings
 * Later on these strings get turned into ints, floats, doubles...
 *
 * @param fileName the name of the settings file to read
 */
void Settings::parseSettingsFile(const char* fileName, std::ofstream& outputFile){

  std::ifstream settingsFile(fileName);

  // Print error message if I can't read the file
  if(!settingsFile.is_open()){
    std::cerr << "Error in " << ANSI_COLOR_BLUE << __FILE__ << ANSI_COLOR_RESET
	      << ", could not open file " << ANSI_COLOR_RED << fileName << ANSI_COLOR_RESET << std::endl;
    exit(1);
  }
  else{

    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    outputFile << "Current date and time are: " << asctime(timeinfo) << std::endl;


    int lineNum = 1;

    // Read every line in the file...
    while(!settingsFile.eof()){

      std::string thisLine;
      std::getline(settingsFile, thisLine);

      // Copy to output file
      outputFile << thisLine << std::endl;

      // First cut out the comment, which is all characters after the first #, including the #
      std::size_t found = thisLine.find("#");


      // Here we switch to TString because of its lovely tokenization methods
      TString thisLineCommentsRemoved(thisLine.substr(0, found));

      // Now have a TObjArray of TObjStrings split by out delimeter :
      TObjArray* tokens = thisLineCommentsRemoved.Tokenize(":");

      int nTokens = tokens->GetEntries();

      // If there are two tokens, then there was one delimeter
      if(nTokens == 2){

	TString key = ((TObjString*) tokens->At(0))->GetString();
	TString value = ((TObjString*) tokens->At(1))->GetString();

	Bool_t addVariable = newKvpPassesSanityChecks(key, value, fileName, lineNum);

	if(addVariable){
	  keyValuePairStrings[key.Data()] = value.Data();
	}

      }
      else{

	TRegexp reggie("[a-zA-Z0-9]");
	Ssiz_t len = thisLineCommentsRemoved.Length();

	Bool_t isAlphaNumeric = reggie.Index(thisLineCommentsRemoved, &len) != -1;

	if(nTokens > 2 || isAlphaNumeric){
	  // complain if there are more than two tokens, i.e. more that one colon.
	  // complain if there are non whitespace characters but no colon.
	  std::cerr << "Warning in " ANSI_COLOR_RED << __FILE__ << ANSI_COLOR_RESET
		    << ". I couldn't parse line " << ANSI_COLOR_RED << lineNum << ANSI_COLOR_RESET
		    << " in " << fileName << ". It said: " << std::endl;
	  std::cerr << ANSI_COLOR_BLUE << thisLine << ANSI_COLOR_RESET << std::endl;
	}
      }
      delete tokens;

      lineNum++;
    }
  }


  outputFile << std::endl << std::endl;
  outputFile << __FILE__ << " has finished parsing " << fileName << std::endl;
  outputFile << std::endl << std::endl;

}



/**
 * Perform some basic checks on the key value pair parsed on this line
 * Prints an appropriate warning message if there was a problem
 * Gets it own function as the warnings are a little verbose
 * @param key is the key (Setting name)
 * @param value is the value (Setting value)
 * @param fileName is the name of the input.conf file
 * @param lineNum is the line being parsed.
 *
 * @return true if we should insert the key into the kvp map, false if there was a problem.
 */
Bool_t Settings::newKvpPassesSanityChecks(const TString& key, const TString& value, const char* fileName, int lineNum){

  Bool_t isGood = true;

  if(key.Length()==0){
    std::cerr << "Warning in " << ANSI_COLOR_BLUE << __FILE__ << ANSI_COLOR_RESET << ", "
	      << ANSI_COLOR_BLUE << fileName << ANSI_COLOR_RESET
	      << " has a variable with no name at " << ANSI_COLOR_RED << "line " << lineNum
	      << ANSI_COLOR_RESET << "." << std::endl;
    isGood = false;
  }

  else if(value.Length()==0){
    std::cerr << "Warning in " << ANSI_COLOR_BLUE << __FILE__ << ANSI_COLOR_RESET << ", "
	      << ANSI_COLOR_BLUE << fileName << ANSI_COLOR_RESET
	      << " has a variable with no value at " << ANSI_COLOR_RED << "line " << lineNum
	      << ANSI_COLOR_RESET << "." << std::endl;
    isGood = false;
  }
  else{
    kvpMap::iterator it = keyValuePairStrings.find(key);

    if(it!=keyValuePairStrings.end()){
    std::cerr << "Warning in " << ANSI_COLOR_BLUE << __FILE__ << ANSI_COLOR_RESET << ", "
	      << ANSI_COLOR_BLUE << fileName << ANSI_COLOR_RESET
	      << " already has a variable named " << ANSI_COLOR_RED << key.Data()
	      << ANSI_COLOR_RESET << "." << std::endl;

      isGood = false;
    }
  }

  if(!isGood){
    std::cerr << "Will ignore " << ANSI_COLOR_RED << "line " << lineNum << ANSI_COLOR_RESET << std::endl;
  }

  return isGood;
}



/**
 * Print all entries in the config key/value pair string map.
 *
 * For debugging and testing
 *
 */
void Settings::printAllKeyValuePairStrings(){

  kvpMap::iterator it;
  for(it = keyValuePairStrings.begin(); it!=keyValuePairStrings.end(); ++it){
    std::cout << it->first << "\t" << it->second << std::endl;
  }
}



/**
 * Set member variables to default values
 *
 */
void Settings::Initialize() {
  NDISCONES_PASS=3;
  DEBUG=false;                   // debugging option
  // outputdir="outputs"; // directory where outputs go
  FREQ_LOW_SEAVEYS=200.E6;
  FREQ_HIGH_SEAVEYS=1200.E6;
  BW_SEAVEYS=FREQ_HIGH_SEAVEYS-FREQ_LOW_SEAVEYS;
  SIGMAPARAM=1;  // Connolly et al. 2011 default cross section parametrization
  YPARAM=1;  // Connolly et al. 2011 default y parametrization
  UNBIASED_SELECTION=1.; // (0) pick neutrino interaction in the ice and neutrino from any direction or (1) choose neutrino interaction point in the horizon on the balloon in the ice and neutrino direction on the cerenkov cone
  SIGMA_FACTOR=1;
  USEDARTBOARD=0;

  // Bunch of variables which were global in icemc.cc but are settings:
  SEED=65540;      // random number seed.
  THETA_TH_FACTOR=1.0; // factor to multiply theta_th to check code is working properly
  CHANCEINHELL_FACTOR=1.0; // loosen chance in hell cuts to check code is working properly
  CONSTANTY=0; // whether or not to set y to a constant=0.2
  VARIABLE_ATTEN=0; //  0=depth dependent attenuation length, 1=fixed
  TRIGTYPE=1; //1=Trigger scheme as in the SMEX proposal  or 0= Just at least 8 channels pass with 2.3 sigma signal
  ATMOSPHERE=1;// include atmosphere
  SCALEDOWNLCPRX1=1; // scale down lcp voltage of antenna 1 by sqrt(2)
  SCALEDOWNEPOLRX1=1; // scale down power of e pol. of antenna 1 by factor of 2
  SCALEDOWNHPOLRX1=1; // scale down power of h pol. of antenna 1 by factor of 2
  SCALEDOWNEPOLRX2=1; // scale down power of e pol. of antenna 2 by factor of user's choice
  SCALEFACTOREPOLRX2=0.05; // scale power of e pol. of antenna 2 by this factor
  SCALEDOWNHPOLRX2=1; // scale down power of h pol. of antenna 2 by factor of 2
  EPOLRX2ZERO=1; // lcp channel on anita-lite is not considered for triggering.
  HPOLRX2ZERO=1; // h pol. of 2nd antenna set to zero.
  RCPRX2ZERO=1; // rcp of 2nd antenna set to zero.
  LCPRX2ZERO=1; // lcp of 2nd antenna set to zero.
  FLATSURFACE=0; // Normals of all positions on the surface are straight up.
  WRITEPOSFILE=0; //Write neutrino position information to file
  SKIPCUTS=0; //See every neutrino through to the end - don't make any of the various cuts designed to speed up the program.  (For checking distributions.)
  USEDIRECTIONWEIGHTS=1;// whether or not to restrict the neutrino angle so that the viewing angle is near the cerenkov cone
  SHOWERTYPE=0; // Type of shower for previous option
  // End of the once-global varibles.
  taumodes = 1; //Taumodes =1, taucreated in the rock.
  SCREENEDGELENGTH=25.;
  ROUGH_INTPOS_SHIFT=100000.;
  TUFFSON=0;
}









void Settings::ReadInputs(const char* inputFileName, std::ofstream &foutput,
			  // Anita* anita1, Secondaries* sec1, Signal* sig1,
			  // Balloon* bn1, Ray* ray1,
			  int& NNU, double& RANDOMISEPOL) {

  parseSettingsFile(inputFileName, foutput);

  getSetting("Number of neutrinos", NNU);
  getSetting("Energy exponent", EXPONENT);
  getSetting("Energy CDF or dartboard", USEDARTBOARD);

  getSetting("Neutrino position", UNBIASED_SELECTION);
  getSetting("Write hists and trees", HIST);
  getSetting("Write ray", FILLRAYTREES);

  getSetting("Only final tree", ONLYFINAL);
  getSetting("Max histogram entries", HIST_MAX_ENTRIES);
  getSetting("Random seed", SEED);
  std::cout << "INITIAL SEED is " << SEED << std::endl;
  gRandom->SetSeed(SEED);

  getSetting("Write neutrino position", WRITEPOSFILE);

  if (WRITEPOSFILE==1){
    std::cout << "Non-default setting: WRITEPOSFILE= " << WRITEPOSFILE << std::endl;
  }
  getSetting("Events map", EVENTSMAP);



// ################################################################################################
// # Balloon and payload
// ################################################################################################



  getSetting("Which payload", WHICH);
  getSetting("Antenna layers", NLAYERS);

  if (WHICH==9)       ANITAVERSION=3;
  else if (WHICH==10) ANITAVERSION=4;
  else                ANITAVERSION=0;
  
  if(((WHICH==1 || WHICH==6) && NLAYERS!=4) || (WHICH==0 && NLAYERS!=1) || (WHICH==7 && NLAYERS!=1)){
    std::cout << "Non-default setting: WHICH = " << WHICH << " and NLAYERS= " << NLAYERS << std::endl;
  }


  getSetting("Inclination top three layers", INCLINE_TOPTHREE);

  if(INCLINE_TOPTHREE!=10){
    std::cout << "Non-default setting: INCLINE_TOPTHREE= " << INCLINE_TOPTHREE << std::endl;
  }

  getSetting("Inclination fourth layer", INCLINE_NADIR);

  if(INCLINE_NADIR!=10){
    std::cout << "Non-default setting: INCLINE_NADIR= " << INCLINE_NADIR << std::endl;
  }

  getSetting("Flight path", WHICHPATH);

  if((WHICH==0 && WHICHPATH!=2) || (WHICH==2 && WHICHPATH!=6)){
    std::cout << "Non-default setting: WHICHPATH = " << WHICHPATH << " and WHICH = "
	      << WHICH << std::endl;
  }

  getSetting("Balloon latitude", BN_LATITUDE);

  getSetting("Balloon longitude", BN_LATITUDE);

  if((BN_LONGITUDE!=999 || BN_LATITUDE!=999) && WHICHPATH==0){
    std::cout << "BN_LATITUDE: "<< BN_LATITUDE << ", BN_LONGITUDE: " << BN_LONGITUDE << std::endl;
  }

  if (BN_LONGITUDE>180. && BN_LONGITUDE!=999){
    std::cout << "Entered balloon longitude wrong!  Should be between -180 and 180 degrees." << std::endl;
  }
  getSetting("Balloon orientation", RANDOMIZE_BN_ORIENTATION);


  if (RANDOMIZE_BN_ORIENTATION==1 && (WHICHPATH==2 || WHICHPATH==6 ||
					   WHICHPATH==7 || WHICHPATH==8 || WHICHPATH==9)){
    std::cout << "Warning:: Strangely you asked for a real flight path but a randomized balloon orientation.  WILL BE OVERRIDDEN." << std::endl;
  }

  getSetting("Balloon altitude", BN_ALTITUDE);


  // whether to use constant gains as entered in GetBeamWidths (0) or to use Ped's measurements as entered in ReadGains (1)
  // GAINS is actually an int, not a double...
  getSetting("Gain setting", GAINS);



  getSetting("Trigger scheme", TRIGGERSCHEME);

  // Ben S, leaving this here... should it get its own config file entry?
  TRIGTYPE=1; // ANITA, not anita-lite.  But the anita-lite code back in later

  getSetting("Band thresholds", tempThresholds);

  getSetting("Banding", BANDING);


  if(BANDING !=0 && BANDING!= 1 && BANDING!=2 && BANDING!=4) {
    std::cout << "Banding should be set to 0 (Anita 1), 1 (custum), 2 (Anita 2), 3 (Satellite) or 4 (Anita 3)."
	      << std::endl;
    exit(1);
  }

  if ((TRIGGERSCHEME==0 || TRIGGERSCHEME==1) && BANDING!=1) {
    std::cout << "Frequency domain trigger schemes can only be used with user-set sub-bands." << std::endl;
    exit(1);
  }

  if (TRIGGERSCHEME==2 && BANDING==1) {
    std::cout << "Time domain trigger scheme only works with Anita 1, Anita 2 or Anita 3 banding data, you can't set your own bands." << std::endl;
    exit(1);
  }


  getSetting("Lower band edges", bandLowEdgesMHz);
  getSetting("Upper band edges", bandHighEdgesMHz);
  getSetting("Required bands", requiredBands);
  getSetting("Allowed bands", allowedBands);
  getSetting("Number of bands", NBANDS);
  getSetting("Percent bandwidth", PERCENTBW);

  getSetting("Notch filter limits", notchFilterLimitsMHz);

  getSetting("Num antenna channels for L1 trigger", trigRequirements[0]);
  getSetting("Num L1 hits to pass L2", trigRequirements[1]);
  getSetting("Num antenna for L2 trigger", antennaclump);
  getSetting("Require centre antenna", REQUIRE_CENTRE);
  getSetting("L3 trigger requirement", trigRequirements[2]);
  getSetting("LCP/RCP or V/H", LCPRCP);
  getSetting("Channels required polarization", channelRequirePol);
  getSetting("Channels allowed polarization", channelAllowedPol);

  getSetting("Exta antennas in trigger", DISCONES);
  getSetting("Nadir only trigger", INCLUDE_NADIRONLY);

  if (INCLUDE_NADIRONLY!=0){
    std::cout << "Non-default setting:  INCLUDE_NADIRONLY = " << INCLUDE_NADIRONLY << std::endl;
  }

  getSetting("ANITA-1 channel masking", CHMASKING);

  if (WHICHPATH!=6 && CHMASKING==1) {
    std::cout << "Cannot include masking for flights other than the ANITA-1 flight." << std::endl;
    std::cout << "For the ANITA-3 channel masking, it is implemented together with the phi masking and it's turned on whenever the PHIMASKING is ON." << std::endl;
    std::cout << "CHMASKING set to 0." << std::endl;
    CHMASKING=0;
  }

  getSetting("ANITA-2 channel masking", PHIMASKING);



  getSetting("Scale down LCP voltage 1st ant", SCALEDOWNLCPRX1);
  getSetting("Scale down E pol 1st ant", SCALEDOWNEPOLRX1);
  getSetting("Scale down H pol 1st ant", SCALEDOWNHPOLRX1);
  getSetting("Scale down E pol 2nd ant", SCALEDOWNEPOLRX2);
  getSetting("E pol scale down factor 2nd ant", SCALEFACTOREPOLRX2);
  getSetting("H pol scale down factor 2nd ant", SCALEDOWNHPOLRX2);
  getSetting("E pol 2nd ant dead", EPOLRX2ZERO);
  getSetting("H pol 2nd ant dead", HPOLRX2ZERO);
  getSetting("RCP 2nd ant dead", RCPRX2ZERO);
  getSetting("LCP 2nd ant dead", LCPRX2ZERO);

  if (WHICH==0 && !(SCALEDOWNEPOLRX1==1 && RCPRX2ZERO==1)){
    std::cout << "Non-default setting:  WHICH= " << WHICH << " and EPOLRX2ZERO= " << EPOLRX2ZERO << std::endl;
  }













  getSetting("Add noise to signal", SIGNAL_FLUCT);

  if (SIGNAL_FLUCT!=1){
    std::cout << "Non-default setting:  SIGNAL_FLUCT= " << SIGNAL_FLUCT << std::endl;
  }

  getSetting("Zero signal", ZEROSIGNAL); 
  if (ZEROSIGNAL!=0){
    std::cout << "Non-default setting:  ZEROSIGNAL= " << ZEROSIGNAL << std::endl;
  }
  
  getSetting("Random rotation polarization", RANDOMISEPOL);

  getSetting("LPM effect", useLPM);
  getSetting("E-field factor", jamieFactor);
  getSetting("Thermal noise factor", THERMALNOISE_FACTOR);

  if (THERMALNOISE_FACTOR!=1){
    std::cout << "Non-default setting:  THERMALNOISE_FACTOR= " << THERMALNOISE_FACTOR << std::endl;
  }

  getSetting("Disable polarization vectors", REMOVEPOLARIZATION);

  if (REMOVEPOLARIZATION==1){
    std::cout << "Non-default setting:  Polarizations turned off!" << std::endl;
  }
  getSetting("Use pulser spectrum", PULSER);
  if (PULSER!=0){
    std::cout << "Warning! Injecting a pulser spectrum- not simulating neutrinos!  PULSER = "
	      << PULSER << std::endl;
  }

  getSetting("Centre one phi-sector", CENTER);

  if (CENTER!=0){
    std::cout << "WARNING!!  Rotating payload to center one phi sector on the incoming signal for each event."
	      << std::endl;
  }

  getSetting("Force vertical polarization", MAKEVERTICAL);

  if (MAKEVERTICAL!=0){
    std::cout << "WARNING!!  Rotating polarization so it is always vertical approaching the payload" << std::endl;
  }









  getSetting("Slopeyness", SLOPEYSIZE);
  if (SLOPEYSIZE!=0.012){
    std::cout << "Non-default setting:  SLOPEYSIZE= " << SLOPEYSIZE << std::endl;
  }
  getSetting("Enable slopeyness", SLOPEY);
  if (SLOPEY!=1){
    std::cout << "Non-default setting:  SLOPEY= " << SLOPEY << std::endl;
  }
  getSetting("Depth dependent refractive index", NOFZ);
  if (NOFZ!=1){
    std::cout << "Non-default setting:  NOFZ= " << NOFZ << std::endl;
  }
  getSetting("Variable attenuation length", VARIABLE_ATTEN);
  if (VARIABLE_ATTEN!=0){
    std::cout << "Non-default setting:  VARIABLE_ATTEN= " << VARIABLE_ATTEN << std::endl;
  }
  getSetting("Constant ice thickness", CONSTANTICETHICKNESS);

  if (CONSTANTICETHICKNESS==1){
    std::cout << "Non-default setting:  CONSTANTICETHICKNESS= " << CONSTANTICETHICKNESS << std::endl;
  }

  getSetting("Antarctic ice model", ICE_MODEL);
  if ((CONSTANTICETHICKNESS || FIXEDELEVATION) && ICE_MODEL != 0) {
    ICE_MODEL=0;
    std::cout << "Constant ice thickness and/or fixed elevation requested.  Using Crust 2.0 ice model." << std::endl;
  } //use the Crust 2.0 data if set to constant icethickness or ground elevation

  if (ICE_MODEL==0){
    std::cout << "Using Crust 2.0 ice model." << std::endl;
  }
  else if (ICE_MODEL==1){
    std::cout << "Using BEDMAP ice model." << std::endl;
  }

  getSetting("Flat surface", FLATSURFACE);
  if (FLATSURFACE==1){
    std::cout << "Non-default setting: all surface segments are flat." << std::endl;
  }

  getSetting("Fixed ice elevation", FIXEDELEVATION);
  if (FIXEDELEVATION==1){
    std::cout << "Non-default setting:  FIXEDELEVATION= " << FIXEDELEVATION << std::endl;
  }

  getSetting("Medium", medium);

  getSetting("Enable surface roughness", ROUGHNESS);
  getSetting("Surface roughness", ROUGHSIZE);
  getSetting("Screen edge length [meters]", SCREENEDGELENGTH);
  getSetting("Interaction position shift [meters]", ROUGH_INTPOS_SHIFT);

  getSetting("Base screen divisions", ROUGHSCREENDIV_BASE);
  getSetting("Subgrid divisions", ROUGHSCREENDIV_SUB);
  getSetting("Base screen fraction threshold", ROUGHSCREENFRAC_BASE);
  getSetting("Subgrid fraction threshold", ROUGHSCREENFRAC_SUB);
  getSetting("Number of subgrid generations", ROUGHMAXGEN);

  getSetting("FIRN", FIRN);
  if (FIRN==0){
    std::cout << "Warning!  Non-standard parameter setting.  FIRN = " << FIRN << std::endl;
  }
  getSetting("Which attenuation length", MOOREBAY);








  getSetting("Cross-section factor", SIGMA_FACTOR);
  if (SIGMA_FACTOR!=1){
    std::cout << "Non-default setting:  settings->SIGMA_FACTOR= " << SIGMA_FACTOR << std::endl;
  }
  getSetting("Theta_th factor", THETA_TH_FACTOR);
  if (THETA_TH_FACTOR!=1){
    std::cout << "Non-default setting:  THETA_TH_FACTOR= " << THETA_TH_FACTOR << std::endl;
  }
  getSetting("Chance in hell factor", CHANCEINHELL_FACTOR);
  if (CHANCEINHELL_FACTOR!=1){
    std::cout << "Non-default setting:  CHANCEINHELL_FACTOR= " << CHANCEINHELL_FACTOR << std::endl;
  }
  getSetting("Skip neutrinos", SKIPCUTS);
  if (SKIPCUTS==1){
    std::cout << "Non-default setting:  Skipping all cuts!" << std::endl;
  }
  getSetting("Restrict neutrino directions", USEDIRECTIONWEIGHTS);
  getSetting("Restrict neutrino positions", USEPOSITIONWEIGHTS);
  if (USEPOSITIONWEIGHTS==0){
    std::cout << "Non-default setting:  Not selecting events within the horizon." << std::endl;
  }
  getSetting("Weight on absorption", WEIGHTABSORPTION);
  getSetting("Phi points banana", horizontal_banana_points);
  getSetting("Theta points banana", vertical_banana_points);
  getSetting("Signal across frequencies", FORSECKEL);
  getSetting("Shower type", SHOWERTYPE);
  getSetting("Loop over boresights", BORESIGHTS);
  if (BORESIGHTS==0){
    std::cout << "Warning!  Non-standard parameter setting.  BORESIGHTS = " << BORESIGHTS << std::endl;
  }







  getSetting("Askaryan parameterization", askaryanParameterization);

  getSetting("Cross-section parameterization", SIGMAPARAM);
  getSetting("Inelasticity parameterization", YPARAM);
  getSetting("Secondary interactions", SECONDARIES);
  if (SECONDARIES!=1){
    std::cout << "Non-default setting:  SECONDARIES= " << SECONDARIES << std::endl;
  }
  getSetting("Tau decay as secondary interaction", TAUDECAY);
  if (TAUDECAY!=1){
    std::cout << "Non-default setting:  TAUDECAY= " << TAUDECAY << std::endl;
  }
  getSetting("Include atmosphere", ATMOSPHERE);
  if (ATMOSPHERE!=1){
    std::cout << "Non-default setting:  ATMOSPHERE= " << ATMOSPHERE << std::endl;
  }
  getSetting("Constant crust density", CONSTANTCRUST);
  if (CONSTANTCRUST==1){
    std::cout << "Non-default setting:  CONSTANTCRUST= " << CONSTANTCRUST << std::endl;
  }
  getSetting("Constant y", CONSTANTY);
  getSetting("Max interaction distance", MAXHORIZON);
  getSetting("Set tau modes", taumodes);














  getSetting("Which rays", WHICHRAYS);
  if (WHICHRAYS!=1){
    std::cout << "Non-default setting:  WHICHRAYS= " << WHICHRAYS << std::endl;
  }
  getSetting("CreateHorizons file", WRITE_FILE);
  if (WRITE_FILE){
    std::cout<<"Writing CreateHorizons input file." << std::endl;
  }
  getSetting("Theta resolution", SIGMA_THETA);
  if (SIGMA_THETA==1){
    std::cout << "Non-default setting:  SIGMA_THETA = 1" << std::endl;
  }
  SIGMA_THETA*=RADDEG; // immediately convert degrees to radians
  getSetting("Low frequency", FREQ_LOW);

  // if (FREQ_LOW_SEAVEYS>anita1->FREQ_LOW){
  //   FREQ_LOW_SEAVEYS=anita1->FREQ_LOW;
  // }
  getSetting("High frequency", FREQ_HIGH);
  BW = FREQ_HIGH - FREQ_LOW; // total bandwidth of simulation






  getSetting("SLAC run", SLAC);
  // this rotates surface slope 10 degrees away from south pole
  if (SLAC==1){
    std::cout << "Warning!  Non-standard parameter setting.  SLAC = " << SLAC << std::endl;
  }
  if (SLAC) {
    foutput << "!!!!SLAC setting causes some settings to be superseded:" << std::endl;
    FIRN=0; // no firn
    foutput << "FIRN=0" << std::endl;
    SLOPEY=0; // slopeyness off
    foutput << "SLOPEY=0" << std::endl;
    BORESIGHTS=1; // loop over boresights
    foutput << "BORESIGHTS=1" << std::endl;
    BN_ALTITUDE=4.22/0.3; // balloon altitude in ft.!!
    foutput << "BN_ALTITUDE=4.22/0.3" << std::endl;
    RANDOMIZE_BN_ORIENTATION=0; // don't randomize the balloon orientation
    foutput << "RANDOMIZE_BN_ORIENTATION=0" << std::endl;
    SKIPCUTS=1; // don't make chance in hell cuts
    foutput << "SKIPCUTS=1" << std::endl;
    SLACSLOPE=5.8; // slope of the ice in degrees
    foutput << "SLACSLOPE=5.8" << std::endl;
    SLACICELENGTH=5.02; // length of the block of ice
    foutput << "SLACICELENGTH=5.02" << std::endl;
  }
  getSetting("SLAC horizontal distance", SLAC_HORIZDIST);
  getSetting("SLAC ice slope", SLACSLOPE);
  getSetting("SLAC block length", SLACICELENGTH);
  getSetting("SLAC interaction depth", SLAC_HORIZ_DEPTH);
  SLAC_DEPTH=tan(SLACSLOPE*RADDEG)*(SLACICELENGTH-SLAC_HORIZ_DEPTH) // height from lowest point of ice
    +21.375*CMINCH/100.; // height from beam to lowest point of ice








  getSetting("Coherent power threshold", COHERENT_THRESHOLD );



  // default values are 0
  APPLYIMPULSERESPONSEDIGITIZER=0;
  APPLYIMPULSERESPONSETRIGGER=0;
  USETIMEDEPENDENTTHRESHOLDS=0;
  USEDEADTIME=0;
  getSetting("Digitizer path impulse response", APPLYIMPULSERESPONSEDIGITIZER);
  std::cout << "Apply impulse response to digitizer path: " << APPLYIMPULSERESPONSEDIGITIZER << std::endl;
  getSetting("Trigger path impulse response", APPLYIMPULSERESPONSETRIGGER);
  std::cout << "Apply impulse response to trigger path: " << APPLYIMPULSERESPONSETRIGGER << std::endl;

#ifdef ANITA_UTIL_EXISTS
  if ( (APPLYIMPULSERESPONSEDIGITIZER || APPLYIMPULSERESPONSETRIGGER) && WHICH!=8 && WHICH!=9 && WHICH!=10) {
    std::cout << "Signal chain impulse response is only available for anita-2 and anita-3." << std::endl;
    exit(1);
  }
#endif
#ifndef ANITA_UTIL_EXISTS
  if (APPLYIMPULSERESPONSEDIGITIZER || APPLYIMPULSERESPONSETRIGGER){
    std::cout << "Signal chain impulse response can only be applied when the Anita tools are sourced." << std::endl;
    exit(1);
  }
#endif
  getSetting("Time dependent thresholds", USETIMEDEPENDENTTHRESHOLDS);
  std::cout << "Use time-dependent thresholds: " << USETIMEDEPENDENTTHRESHOLDS << std::endl;
  getSetting("Dead time", USEDEADTIME);
  std::cout << "Use dead time from flight: " << USEDEADTIME << std::endl;
  
  if ( (USETIMEDEPENDENTTHRESHOLDS || USEDEADTIME) && WHICH!=9 && WHICH!=10) {
    std::cout << "Time-dependent thresholds are only available for anita-3." << std::endl;
    exit(1);
  }


  getSetting("Digitizer noise from flight", NOISEFROMFLIGHTDIGITIZER);
  std::cout << "Use noise from flight for digitizer path: " << NOISEFROMFLIGHTDIGITIZER << std::endl;

  getSetting("Trigger noise from flight", NOISEFROMFLIGHTTRIGGER);
  std::cout << "Use noise from flight for trigger path: " << NOISEFROMFLIGHTTRIGGER << std::endl;

#ifdef ANITA3_EVENTREADER
  if ( (NOISEFROMFLIGHTDIGITIZER || NOISEFROMFLIGHTTRIGGER) && (WHICH!=9 && WHICH!=10)) {
    std::cout << "Noise from flight only available for anita-3." << std::endl;
    exit(1);
  }
  if (!APPLYIMPULSERESPONSETRIGGER && NOISEFROMFLIGHTTRIGGER ){
    std::cout << "Noise from flight can only be applied to trigger path if impulse reponse is also used " << std::endl;
    exit(1);
  }
#endif

#ifndef ANITA_UTIL_EXISTS
  if (NOISEFROMFLIGHTDIGITIZER || NOISEFROMFLIGHTTRIGGER){
    std::cout << "Noise from flight can only be applied when the Anita tools are sourced." << std::endl;
    exit(1);
  }
#endif

  getSetting("Min bias", MINBIAS);
  if (MINBIAS){
    std::cout << "Generate Minimum Bias sample: " << MINBIAS << std::endl;
  }


  getSetting("Efficiency scan",      TRIGGEREFFSCAN                       );
  if (TRIGGEREFFSCAN==1){
    getSetting("Central phi-sector",   trigEffScanPhi                     );
    getSetting("Apply pulse at surf",  TRIGGEREFFSCAPULSE                 );
    getSetting("Off-axis attenuation", efficiencyScanOffAxisAttenuations  );
    getSetting("Rings used",           efficiencyScanRingsUsed            );
    getSetting("Phi-sectors delays",   efficiencyScanPhiSectorDelay       );
    getSetting("Ring delays",          efficiencyScanRingDelay            );
    getSetting("Ring delays to phi",   efficiencyScanApplyRingDelay       );
  }
  
  getSetting("Simulate TUFFs", TUFFSON);
  getSetting("Which TUFFs are on", whichTUFFsON);

  
} //method ReadInputs



  

void Settings::ApplyInputs(Anita* anita1, Secondaries* sec1, Signal* sig1,
			   Balloon* bn1, Ray* ray1){
  
   //When you look at the Anita payload there are 4 layers, with 8,8,16 and 8 antennas each.  But in the trigger, the top two become one layer of 16 antennas. 
  if (WHICH==2 || WHICH==6 || WHICH==8 || WHICH==9 || WHICH==10){
    anita1->NTRIGGERLAYERS = NLAYERS - 1;
  }
  else{
    anita1->NTRIGGERLAYERS=NLAYERS;
  }


  anita1->INCLINE_TOPTHREE=INCLINE_TOPTHREE;
  anita1->INCLINE_NADIR=INCLINE_NADIR;

  bn1->WHICHPATH=WHICHPATH;

  
  if(bn1->WHICHPATH==2){
    anita1->LIVETIME=45.*24.*3600.*0.75; // 45 days for anita-lite
  }
  else if (bn1->WHICHPATH==0){
    anita1->LIVETIME=6.02*24.*3600.; // anita-lite
  } else if (bn1->WHICHPATH==6){
    // kim's livetime for anita
    anita1->LIVETIME=17.*24.*3600.; // for anita, take 34.78 days * 0.75 efficiency
  }
  else if (bn1->WHICHPATH==7){
    anita1->LIVETIME=28.5*24*3600;  // Anita-2 livetime taken from paper
  }
  else if (bn1->WHICHPATH==8){
    anita1->LIVETIME=17.4*24*3600;  // Anita-3 livetime taken from Ben Strutt's thesis (elog note 698)
  } else {
    anita1->LIVETIME=14.*24.*3600.; // otherwise use 2 weeks by default
  }
  
  if (WHICH==7){
    // EeVEX
    anita1->LIVETIME=100.*24.*3600.; // ultra-long duration balloon flight of 100 days
  }

  bn1->BN_LATITUDE              = BN_LATITUDE;
  bn1->BN_LONGITUDE             = BN_LONGITUDE;
  bn1->BN_ALTITUDE              = BN_ALTITUDE;
  bn1->RANDOMIZE_BN_ORIENTATION = RANDOMIZE_BN_ORIENTATION;
  bn1->MAXHORIZON               = MAXHORIZON;

  anita1->GAINS   = GAINS;
  anita1->BANDING = BANDING;
  anita1->NBANDS  = NBANDS;
  anita1->PERCENTBW = PERCENTBW;
  anita1->PULSER    = PULSER;
  anita1->SIGMA_THETA = SIGMA_THETA;
  anita1->FREQ_LOW    = FREQ_LOW;
  anita1->FREQ_HIGH   = FREQ_HIGH;

  for (unsigned int i=0;i<tempThresholds.size();i++) {
    anita1->bwslice_thresholds[i] = tempThresholds.at(i);
  }


  for (unsigned int i=0; i < bandLowEdgesMHz.size(); i++) {
    anita1->bwslice_min[i] = 1e6*bandLowEdgesMHz.at(i);
    anita1->bwslice_max[i] = 1e6*bandHighEdgesMHz.at(i);
    anita1->bwslice_center[i] = 0.5*(anita1->bwslice_min[i] + anita1->bwslice_max[i]);
    anita1->bwslice_width[i] = anita1->bwslice_max[i] - anita1->bwslice_min[i];
  }

for(unsigned int i=0; i < requiredBands.size(); i++){
    anita1->bwslice_required[i] = requiredBands.at(i);
  }

  for(unsigned int i=0; i < allowedBands.size(); i++){
    anita1->bwslice_allowed[i] = allowedBands.at(i);
  }

  anita1->maxthreshold=0.;
  anita1->bwmin=1.E10;
  if (anita1->BANDING!=1){
    anita1->bwmin=200.E6;
  }

  const int numBands = 5;
  for (int i=0; i<numBands; i++) {
    if (anita1->bwslice_thresholds[i]>anita1->maxthreshold && anita1->bwslice_allowed[i]==1){
      anita1->maxthreshold = anita1->bwslice_thresholds[i];
    }
    if (anita1->BANDING==1) {
      if ((anita1->bwslice_max[i] - anita1->bwslice_min[i]) < anita1->bwmin && anita1->bwslice_allowed[i] == 1){

	anita1->bwmin = anita1->bwslice_max[i] - anita1->bwslice_min[i];
      }
    }
  }


  anita1->NOTCH_MIN = 1e6*notchFilterLimitsMHz.at(0);
  anita1->NOTCH_MAX = 1e6*notchFilterLimitsMHz.at(1);
  
  if (anita1->NOTCH_MIN>anita1->NOTCH_MAX) {
    std::cout << "Min of notch filter is greater than max. Try again." << std::endl;
  }
  if (anita1->NOTCH_MIN!=0 || anita1->NOTCH_MAX!=0){
    std::cout << "Applying a notch filter from " << anita1->NOTCH_MIN << " Hz to "
	      << anita1->NOTCH_MAX << " Hz" << std::endl;
  }	      
  anita1->trigRequirements[0] = trigRequirements[0];
  anita1->trigRequirements[1] = trigRequirements[1];	
  anita1->trigRequirements[2] = trigRequirements[2];	
  anita1->REQUIRE_CENTRE = REQUIRE_CENTRE;		
 for(unsigned int i=0; i < channelRequirePol.size(); i++){
    anita1->pol_required[i] = channelRequirePol.at(i); 
  }
 for(unsigned int i=0; i < channelAllowedPol.size(); i++){
    anita1->pol_allowed[i] = channelAllowedPol.at(i);
  }



  sig1->SetLPM(useLPM);
  if (sig1->GetLPM()!=1){
    std::cout << "Non-default setting:  LPM= " << sig1->GetLPM() << std::endl;
  }
  sig1->SetParameterization(askaryanParameterization);
  sig1->SetJaime_Factor(jamieFactor);
  sig1->SetMedium(medium);


  sec1->SECONDARIES=SECONDARIES;
  sec1->TAUDECAY=TAUDECAY;

  anita1->trigEffScanPhi = trigEffScanPhi;
  for (unsigned int i=0; i < efficiencyScanOffAxisAttenuations.size(); i++){
    anita1->trigEffScanAtt[i] = efficiencyScanOffAxisAttenuations.at(i);
  }

  for (unsigned int i=0; i < efficiencyScanPhiSectorDelay.size(); i++){
    //convert ns to s
    anita1->trigEffScanPhiDelay[i] = efficiencyScanPhiSectorDelay.at(i)*1e-9;
  }

  for (unsigned int i=0; i < efficiencyScanRingDelay.size(); i++){
    //convert ns to s
    anita1->trigEffScanRingDelay[i] = efficiencyScanRingDelay.at(i)*1e-9;
  }

  for (unsigned int i=0; i < efficiencyScanApplyRingDelay.size(); i++){
    anita1->trigEffScanApplyRingDelay[i] = efficiencyScanApplyRingDelay.at(i);
  }

  for (unsigned int i=0; i < efficiencyScanRingsUsed.size(); i++){
    anita1->trigEffScanRingsUsed[i] = efficiencyScanRingsUsed.at(i);
  }

  
  
  if (TRIGGEREFFSCAN){
    std::cout << "Let's do a trigger efficiency scan!" << std::endl;
    std::cout << "Apply pulse at AMPA (0) or SURF (1) : " << TRIGGEREFFSCAPULSE << std::endl;
    std::cout << "Central phi sector is " << anita1->trigEffScanPhi << std::endl;
    std::cout << "Attenuations are ";
    for (int i=0;i<5;i++) std::cout << anita1->trigEffScanAtt[i] << " ";
    std::cout << std::endl;
    std::cout << "Phi sector delays are ";
    for (int i=0;i<5;i++) std::cout << anita1->trigEffScanPhiDelay[i] << " ";
    std::cout << std::endl;
    std::cout << "The rings used in this scan are ";
    for (int i=0;i<3;i++) std::cout << anita1->trigEffScanRingsUsed[i] << " ";
    std::cout << std::endl;
    std::cout << "Ring delays are applie to : ";
    for (int i=0;i<5;i++) std::cout << anita1->trigEffScanApplyRingDelay[i] << " ";
    std::cout << std::endl;
    std::cout << "Ring delays are for T-M, M-B, T-B : ";
    for (int i=0;i<3;i++) std::cout << anita1->trigEffScanRingDelay[i] << " ";
    std::cout << std::endl;
  }
  


  for (unsigned int i=0; i< whichTUFFsON.size(); i++){
    anita1->TUFFstatus[i] = whichTUFFsON.at(i);
  }

  if (TUFFSON){
    std::cout << "The TUFFs are ON for the whole flight!" << std::endl;
    std::cout << "Notch 0 status " << anita1->TUFFstatus[0] << std::endl;
    std::cout << "Notch 1 status " << anita1->TUFFstatus[1] << std::endl;
    std::cout << "Notch 2 status " << anita1->TUFFstatus[2] << std::endl;
  }


}










void Settings::complainAboutNotFindingKey(const TString& key){
  std::cerr << "Warning in " << ANSI_COLOR_BLUE << __FILE__ << ANSI_COLOR_RESET
	    << ", unable to find setting " << ANSI_COLOR_RED << key << ANSI_COLOR_RESET << std::endl;
}


void Settings::getSetting(const char* key, int& value){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    value = atoi(it->second.Data());
  }
}

void Settings::getSetting(const char* key, float& value){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    value = atof(it->second.Data());
  }
}

void Settings::getSetting(const char* key, double& value){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    value = atof(it->second.Data());
  }
}

void Settings::getSetting(const char* key, std::vector<int>& valueArray){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    parseValueArray(it->second.Data(), valueArray);
  }
}

void Settings::getSetting(const char* key, std::vector<float>& valueArray){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    parseValueArray(it->second.Data(), valueArray);
  }
}

void Settings::getSetting(const char* key, std::vector<double>& valueArray){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    parseValueArray(it->second.Data(), valueArray);
  }
}

void Settings::parseValueArray(const char* valueString, std::vector<int>& values){
  TString theValueString(valueString);

  TObjArray* theValues = theValueString.Tokenize(",");
  for(int i=0; i < theValues->GetEntries(); ++i){

    TObjString* token = (TObjString*) theValues->At(i);
    int value = atoi(token->GetString().Data());
    values.push_back(value);
  }
}

void Settings::parseValueArray(const char* valueString, std::vector<float>& values){
  TString theValueString(valueString);

  TObjArray* theValues = theValueString.Tokenize(",");
  for(int i=0; i < theValues->GetEntries(); ++i){

    TObjString* token = (TObjString*) theValues->At(i);
    float value = atof(token->GetString().Data());
    values.push_back(value);
  }
}

void Settings::parseValueArray(const char* valueString, std::vector<double>& values){
  TString theValueString(valueString);

  TObjArray* theValues = theValueString.Tokenize(",");
  for(int i=0; i < theValues->GetEntries(); ++i){

    TObjString* token = (TObjString*) theValues->At(i);
    double value = atof(token->GetString().Data());
    values.push_back(value);
  }
}
