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
#include "trigger.hh"
#include "Spectra.h"
#include "signal.hh"
#include "secondaries.hh"
#include "ray.hh"
#include "counting.hh"
#include "Primaries.h"


#include "String.h"

#include "TString.h"
#include "TRegexp.h"
#include "TObjString.h"


/**
 * Default constructor
 *
 */
Settings::Settings() {
  Initialize();
}





/**
 * Constructor taking filename
 *
 * @param fileName
 */
Settings::Settings(const char* fileName){
  Initialize();
  readSettingsFile(fileName);
}





/**
 * Copies the contents of a yaml settings file into internal memory as strings
 * Later on these strings get turned into ints, floats, doubles...
 *
 * @param fileName the name of the settings file to read
 */
void Settings::readSettingsFile(const char* fileName){

  std::ifstream settingsFile(fileName);

  // Print error message if I can't read the file
  if(!settingsFile.is_open()){
    std::cerr << "Warning in " << __FILE__ << ", could not open file " << fileName << std::endl;
  }
  else{
    int lineNum = 1;

    // Read every line in the file...
    while(!settingsFile.eof()){

      std::string thisLine;
      std::getline(settingsFile, thisLine);

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
	  std::cerr << "Warning in " << __FILE__ << ". I couldn't parse line " << lineNum
		    << " in " << fileName << ". It said: " << std::endl;
	  std::cerr << thisLine << std::endl;
	}
      }
      delete tokens;

      lineNum++;
    }
  }
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
    std::cerr << "Warning in " << __FILE__ << ", " << fileName
	      << " has a variable with no name at line " << lineNum << "." << std::endl;
    isGood = false;
  }

  else if(value.Length()==0){
    std::cerr << "Warning in " << __FILE__ << ", " << fileName
	      << " has a variable with no value at line " << lineNum << "." << std::endl;
    isGood = false;
  }
  else{
    kvpMap::iterator it = keyValuePairStrings.find(key);

    if(it!=keyValuePairStrings.end()){
      std::cerr << "Warning in " << __FILE__ << ", " << fileName
		<< " already has a variable named " << key.Data() << std::endl;
      isGood = false;
    }
  }

  if(!isGood){
    std::cerr << "Will ignore line " << lineNum << std::endl;
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
  outputdir="outputs"; // directory where outputs go
  FREQ_LOW_SEAVEYS=200.E6;
  FREQ_HIGH_SEAVEYS=1200.E6;
  BW_SEAVEYS=FREQ_HIGH_SEAVEYS-FREQ_LOW_SEAVEYS;
  SIGMAPARAM=1;  // Connolly et al. 2011 default cross section parametrization
  YPARAM=1;  // Connolly et al. 2011 default y parametrization
  UNBIASED_SELECTION=1.; // (0) pick neutrino interaction in the ice and neutrino from any direction or (1) choose neutrino interaction point in the horizon on the balloon in the ice and neutrino direction on the cerenkov cone
  SIGMA_FACTOR=1;


  // Bunch of variables which were global in icemc.cc but are settings:
  FILLRAYTREES=1; // fill tree for each ray in roughness simulation
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


}


int lineCounter = 0;

void thisGetLine(ifstream& inFile, std::string& junk){
  getline(inFile, junk);
  lineCounter++;
}

void getNextNumberAsString(ifstream& fin, ofstream& fout, string& number) {
    string temp;
    thisGetLine(fin,temp); // get next line of the input file

    fout << temp << "\n"; // output this line to the summary file

    int place=0;
    place=temp.find_first_of(" \t"); // find where the first space or tab is

    number=temp.substr(0,place); // everything up until the first space is what we're looking for

    std::cout << lineCounter << "\t" << temp << std::endl;
} //GetNextNumberAsString







void Settings::ReadInputs(ifstream &inputsfile, ofstream &foutput, Anita* anita1, Secondaries* sec1, Signal* sig1, Balloon* bn1, Ray* ray1, int& NNU, double& RANDOMISEPOL) {
  // read inputs to the code.
  // See comments in input file

  // extern int NNU;
  // extern double RANDOMISEPOL;

  string number;
  string junk;

  thisGetLine(inputsfile,junk);

  foutput << "\n\n";

  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  foutput << "Current date and time are: " << asctime (timeinfo) << "\n";


  thisGetLine(inputsfile,junk);
  foutput << junk << "\n";
  // The following inputs have to do with event input/output

  //Fenfang's livetime
  //The LIVETIME I used is 34.78days*0.75

  getNextNumberAsString(inputsfile,foutput,number);
  NNU=(int)atoi(number.c_str());
  getSetting("Number of neutrinos", NNU);

  getNextNumberAsString(inputsfile,foutput,number);
  EXPONENT=(double)atof(number.c_str());
  getSetting("Energy exponent", EXPONENT);


  getNextNumberAsString(inputsfile,foutput,number);
  UNBIASED_SELECTION=(double)atof(number.c_str());
  getSetting("Neutrino position", UNBIASED_SELECTION);

  getNextNumberAsString(inputsfile,foutput,number);
  HIST=(int)atoi(number.c_str());
  getSetting("Write hists and trees", HIST);


  getNextNumberAsString(inputsfile,foutput,number);
  FILLRAYTREES=(int)atoi(number.c_str());
  getSetting("Write ray", FILLRAYTREES);

// ################################################################################################
// Number of neutrinos: 2000000 # How many neutrinos to generate
// Energy exponent: 20 # Select energy (just enter the exponent) or (30) for baseline ES&S (1) for E^-1 (2) for E^-2 (3) for E^-3 (4) for E^-4 (5) for ES&S flux with cosmological constant (6) for neutrino GZK flux from Iron nucleii (16-22)not using spectrum but just for a single energy (101-114)use all kinds of theoretical flux models
// Neutrino position: 0 # (0) Pick neutrino position in the ice within the horizon of the balloon and neutrino direction such that the specular ray hitting the balloon sits on the Cerenkov cone or (1) Pick neutrino interaction point from anywhere in the ice and neutrino direction from all directions equally
// Write hists and trees: 1 # 1=write to histograms and trees, 0=don't write to histograms & trees
// Write ray: 0      # 1=write to ray tree (for roughness study), 0 don't write
// Only final tree: 0     # 1=only write to final tree (for events that pass)
// Max histogram entries: 10000  # Maximum number of entries to store in histograms (10000=default).
// Random seed: 65546  # Random number seed
// Write neutrino position: 0      # Write neutrino position info to file
// Events map: 0 # (1) EVENTSMAP draw the events distribution map or (0) not



  getNextNumberAsString(inputsfile,foutput,number);
  ONLYFINAL=(int)atoi(number.c_str());
  getSetting("Only final tree", ONLYFINAL);


  getNextNumberAsString(inputsfile,foutput,number);
  HIST_MAX_ENTRIES=(int)atoi(number.c_str());

  getSetting("Max histogram entries", HIST_MAX_ENTRIES);
  getNextNumberAsString(inputsfile,foutput,number);

  SEED=(int)atoi(number.c_str());
  getSetting("Random seed", SEED);
  std::cout << "SEED is " << SEED << std::endl;
  gRandom->SetSeed(SEED);

  getNextNumberAsString(inputsfile,foutput,number);
  WRITEPOSFILE=(int)atof(number.c_str());
  getSetting("Write neutrino position", WRITEPOSFILE);

  if (WRITEPOSFILE==1){
    std::cout << "Non-default setting: WRITEPOSFILE= " << WRITEPOSFILE << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  EVENTSMAP=atoi(number.c_str());//draw the events map or not
  getSetting("Events map", EVENTSMAP);

  thisGetLine(inputsfile,junk);
  foutput << junk << "\n";
  // The following inputs have to do with the payload and balloon



// ################################################################################################
// # Balloon and payload
// ################################################################################################
// Which payload: 9 # Which payload:  ANITA-lite (0), Ross Payload (1), 2006-2007 flight with simple settings, Make Your Own (3), Ant Hill (4), 2006-2007 flight with Kurt measurements(6), EeVEX (one layer of antennas) (7), ANITA-II (8), ANITA-III (9), (SWORD) Satellite (10)
// Antenna layers: 4 # How many physical layers of antennas to include.  If you indicate N, it will consider the top N layers and ignore the rest.  Useful for turning off the nadirs.  4=default
// Inclination angle top three layers: 10     # Inclination angle in degrees for top three antenna layers 10=default
// Inclination angle fourth layer: 10     # Inclination angle in degrees for fourth layer of antennas 10=default
// Flight path: 8 # (0) fixed balloon position (1) randomized at 80deg S (2) anita-lite flight. (3) banana plot (note: banana plot forces energy, turns off noise and turns off slopyness.  "How many neutrinos to generate" is not used, nor are any inputs relating to a detection system. Also, constant crust density, ice thickness, and elevation are set.) (4) for comparison with Peter (5) Pick your own altitude, phi and theta 2=default (6) Anita 1 flight (7) Anita 2 flight (8) Anita 3 flight
// Fixed balloon latitude: 999  # (degrees south) (999 is default position)
// Fixed balloon longitude: 999 # (degrees east) (999 is default position) longitude is between -180 and 180
// Balloon orientation : 0  # Fixed (0) or randomized (1) balloon orientation.  0=default
// Balloon altitude: 120000. # Pick balloon altitude in ft. (0) 120000 ft. or whatever you choose.
// Gain setting: 1 # Constant gains (0) or Ped's measured gains (1)




  getNextNumberAsString(inputsfile,foutput,number);
  WHICH=(int)atoi(number.c_str());
  getSetting("Which payload", WHICH);

  getNextNumberAsString(inputsfile,foutput,number);
  NLAYERS=(int)atoi(number.c_str()); // this is number of layers, counting the upper 16 as 2 layers
  getSetting("Antenna layers", NLAYERS);

  if(((WHICH==1 || WHICH==6) && NLAYERS!=4) || (WHICH==0 && NLAYERS!=1) || (WHICH==7 && NLAYERS!=1)){
    std::cout << "Non-default setting: WHICH = " << WHICH << " and NLAYERS= " << NLAYERS << std::endl;
  }


  //When you look at the Anita payload there are 4 layers, with 8,8,16 and 8 antennas each.  But in the trigger, the top two become one layer of 16 antennas.  So that means for Anita 1 and Anita 2, there are one fewer trigger layers than physical layers.
  // anything anita like
  // anita 1 simple, anita 1 kurt, anita 2 kurt, anita 3, satellite
  if (WHICH==2 || WHICH==6 || WHICH==8 || WHICH==9 || WHICH==10){
    anita1->NTRIGGERLAYERS = NLAYERS - 1;
  }
  else{
    anita1->NTRIGGERLAYERS=NLAYERS;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  anita1->INCLINE_TOPTHREE=(double)atof(number.c_str());
  getSetting("Inclination top three layers", anita1->INCLINE_TOPTHREE);

  if(anita1->INCLINE_TOPTHREE!=10){
    std::cout << "Non-default setting: INCLINE_TOPTHREE= " << anita1->INCLINE_TOPTHREE << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  anita1->INCLINE_NADIR=(double)atof(number.c_str());
  getSetting("Inclination fourth layer", anita1->INCLINE_NADIR);

  if(anita1->INCLINE_NADIR!=10){
    std::cout << "Non-default setting: INCLINE_NADIR= " << anita1->INCLINE_NADIR << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  bn1->WHICHPATH=(int)atoi(number.c_str());
  getSetting("Flight path", bn1->WHICHPATH);

  if((WHICH==0 && bn1->WHICHPATH!=2) || (WHICH==2 && bn1->WHICHPATH!=6)){
    std::cout << "Non-default setting: bn1->WHICHPATH = " << bn1->WHICHPATH << " and WHICH = "
	      << WHICH << std::endl;
  }

  if(bn1->WHICHPATH==2){
    anita1->LIVETIME=45.*24.*3600.*0.75; // 45 days for anita-lite
  }
  else if (bn1->WHICHPATH==0){
    anita1->LIVETIME=6.02*24.*3600.; // anita-lite
  }
  else if (bn1->WHICHPATH==6){
    // kim's livetime for anita
    anita1->LIVETIME=17.*24.*3600.; // for anita, take 34.78 days * 0.75 efficiency
  }
  else{
    anita1->LIVETIME=14.*24.*3600.; // otherwise use 2 weeks by default
  }

  if (WHICH==7){
    // EeVEX
    anita1->LIVETIME=100.*24.*3600.; // ultra-long duration balloon flight of 100 days
  }

  getNextNumberAsString(inputsfile,foutput,number);
  bn1->BN_LATITUDE=(double)atof(number.c_str());
  getSetting("Balloon latitude", bn1->BN_LATITUDE);

  getNextNumberAsString(inputsfile,foutput,number);
  bn1->BN_LONGITUDE=(double)atof(number.c_str());
  getSetting("Balloon longitude", bn1->BN_LATITUDE);

  if((bn1->BN_LONGITUDE!=999 || bn1->BN_LATITUDE!=999) && bn1->WHICHPATH==0){
    std::cout << "BN_LATITUDE: "<< bn1->BN_LATITUDE << ", BN_LONGITUDE: " << bn1->BN_LONGITUDE << std::endl;
  }

  if (bn1->BN_LONGITUDE>180. && bn1->BN_LONGITUDE!=999){
    std::cout << "Entered balloon longitude wrong!  Should be between -180 and 180 degrees." << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  bn1->RANDOMIZE_BN_ORIENTATION=(int)atoi(number.c_str());
  getSetting("Balloon orientation", bn1->RANDOMIZE_BN_ORIENTATION);


  if (bn1->RANDOMIZE_BN_ORIENTATION==1 && (bn1->WHICHPATH==2 || bn1->WHICHPATH==6 ||
					   bn1->WHICHPATH==7 || bn1->WHICHPATH==8)){
    std::cout << "Warning:: Strangely you asked for a real flight path but a randomized balloon orientation.  WILL BE OVERRIDDEN." << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  bn1->BN_ALTITUDE=(double)atof(number.c_str());
  getSetting("Balloon altitude", bn1->BN_ALTITUDE);


  // whether to use constant gains as entered in GetBeamWidths (0) or to use Ped's measurements as entered in ReadGains (1)
  // GAINS is actually an int, not a double...
  getNextNumberAsString(inputsfile,foutput,number);
  anita1->GAINS=(int)atof(number.c_str());
  getSetting("Gain setting", anita1->GAINS);






  thisGetLine(inputsfile,junk);
  foutput << junk << "\n";
  // The following inputs have to do with event input/output





  getNextNumberAsString(inputsfile,foutput,number);
  TRIGGERSCHEME=atoi(number.c_str()); // whether it's frequency domain (0) or time domain
  getSetting("Trigger scheme", TRIGGERSCHEME);

  // Ben S, leaving this here... should it get its own config file entry?
  TRIGTYPE=1; // ANITA, not anita-lite.  But the anita-lite code back in later



  vector<string> vnumber;
  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,5);

  for (int n=0;n<5;n++) {
    anita1->bwslice_thresholds[n]=(double)atof(vnumber[n].c_str());
  }

  std::vector<double> tempThresholds;
  getSetting("Band thresholds", tempThresholds);
  for (unsigned int i=0;i<tempThresholds.size();i++) {
    anita1->bwslice_thresholds[i] = tempThresholds.at(i);
    std::cout << i << "\t" << anita1->bwslice_thresholds[i] << std::endl;
  }


  getNextNumberAsString(inputsfile,foutput,number);
  anita1->BANDING=atoi(number.c_str()); // whether you use anita-1 banding (0) or choose-your-own
  getSetting("Banding", anita1->BANDING);


  if(anita1->BANDING !=0 && anita1->BANDING!= 1 && anita1->BANDING!=2 && anita1->BANDING!=4) {
    std::cout << "Banding should be set to 0 (Anita 1), 1 (custum), 2 (Anita 2), 3 (Satellite) or 4 (Anita 3)."
	      << std::endl;
    exit(1);
  }

  if ((TRIGGERSCHEME==0 || TRIGGERSCHEME==1) && anita1->BANDING!=1) {
    std::cout << "Frequency domain trigger schemes can only be used with user-set sub-bands." << std::endl;
    exit(1);
  }

  if (TRIGGERSCHEME==2 && anita1->BANDING==1) {
    std::cout << "Time domain trigger scheme only works with Anita 1, Anita 2 or Anita 3 banding data, you can't set your own bands." << std::endl;
    exit(1);
  }


  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,5);
  vector<string> vnumber2;
  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber2,5);

  std::vector<double> bandLowEdgesMHz;
  getSetting("Lower band edges", bandLowEdgesMHz);
  std::vector<double> bandHighEdgesMHz;
  getSetting("Upper band edges", bandHighEdgesMHz);

  for (unsigned int i=0; i < bandLowEdgesMHz.size(); i++) {
    anita1->bwslice_min[i] = 1e6*bandLowEdgesMHz.at(i);
    anita1->bwslice_max[i] = 1e6*bandHighEdgesMHz.at(i);
    anita1->bwslice_center[i] = 0.5*(anita1->bwslice_min[i] + anita1->bwslice_max[i]);
    anita1->bwslice_width[i] = anita1->bwslice_max[i] - anita1->bwslice_min[i];
  }

  // for (int n=0;n<5;n++) {
  //   anita1->bwslice_min[n]=(double)atof(vnumber[n].c_str())*1.E6;
  //   anita1->bwslice_max[n]=(double)atof(vnumber2[n].c_str())*1.E6;
  //   anita1->bwslice_center[n]=(anita1->bwslice_min[n]+anita1->bwslice_max[n])/2.;
  //   anita1->bwslice_width[n]=(anita1->bwslice_max[n]-anita1->bwslice_min[n]);
  //   //cout << "center, width are " << anita1->bwslice_center[n] << " " << anita1->bwslice_width[n] << "\n";
  // }

  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,5);

  // for (int n=0;n<5;n++) {
  //   //anita1->bwslice is actaully an int
  //   anita1->bwslice_required[n]=(int)atof(vnumber[n].c_str());
  // }

  std::vector<int> requiredBands;
  getSetting("Required bands", requiredBands);
  for(unsigned int i=0; i < requiredBands.size(); i++){
    anita1->bwslice_required[i] = requiredBands.at(i);
  }

  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,5);
  for (int n=0;n<5;n++) {
    //anita1->bwslice_allowed[n] is still an int, not a double
    anita1->bwslice_allowed[n]=(int)atof(vnumber[n].c_str());
  }

  std::vector<int> allowedBands;
  getSetting("Allowed bands", allowedBands);
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


  // now get number of bands for banding option 3 (satellite)
  getNextNumberAsString(inputsfile,foutput,number);
  anita1->NBANDS=atoi(number.c_str());
  getSetting("Number of bands", anita1->NBANDS);


  getNextNumberAsString(inputsfile,foutput,number);
  anita1->PERCENTBW=atoi(number.c_str());
  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,2);

  anita1->NOTCH_MIN=(double)atof(vnumber[0].c_str())*1.E6; // min and max frequencies for notch filter
  anita1->NOTCH_MAX=(double)atof(vnumber[1].c_str())*1.E6;

  std::vector<double> notchFilterLimitsMHz;
  getSetting("Notch filter limits", notchFilterLimitsMHz);
  anita1->NOTCH_MIN = 1e6*notchFilterLimitsMHz.at(0);
  anita1->NOTCH_MAX = 1e6*notchFilterLimitsMHz.at(1);


  if (anita1->NOTCH_MIN>anita1->NOTCH_MAX) {
    std::cout << "Min of notch filter is greater than max. Try again." << std::endl;
  }
  if (anita1->NOTCH_MIN!=0 || anita1->NOTCH_MAX!=0){
    std::cout << "Applying a notch filter from " << anita1->NOTCH_MIN << " Hz to "
	      << anita1->NOTCH_MAX << " Hz" << std::endl;
  }


  getNextNumberAsString(inputsfile,foutput,number);
  anita1->trigRequirements[0]=atoi(number.c_str());
  getSetting("Num antenna channels for L1 trigger", anita1->trigRequirements[0]);

  getNextNumberAsString(inputsfile,foutput,number);
  anita1->trigRequirements[1]=atoi(number.c_str());
  getSetting("Num L1 hits to pass L2", anita1->trigRequirements[1]);

  getNextNumberAsString(inputsfile,foutput,number);
  antennaclump=atoi(number.c_str());
  getSetting("Num antenna for L2 trigger", antennaclump);

  getNextNumberAsString(inputsfile,foutput,number);
  anita1->REQUIRE_CENTRE=atoi(number.c_str()); // require centre antenna in clump is one of those his
  getSetting("Require centre antenna", anita1->REQUIRE_CENTRE);

  getNextNumberAsString(inputsfile,foutput,number);
  anita1->trigRequirements[2]=atoi(number.c_str());
  getSetting("L3 trigger requirement", anita1->trigRequirements[2]);

  getNextNumberAsString(inputsfile,foutput,number);
  LCPRCP=atoi(number.c_str());
  getSetting("LCP/RCP or V/H", LCPRCP);


  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,2);

  for (int n=0;n<2;n++) {
    //anita1->bwslice is actaully an int
    anita1->pol_required[n]=(int)atof(vnumber[n].c_str());
  }

  std::vector<int> channelRequirePol;
  getSetting("Channels required polarization", channelRequirePol);
  for(unsigned int i=0; i < channelRequirePol.size(); i++){
    anita1->pol_required[i] = channelRequirePol.at(i);
  }

  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,2);

  for (int n=0;n<2;n++) {
    //anita1->bwslice_allowed[n] is still an int, not a double
    anita1->pol_allowed[n]=(int)atof(vnumber[n].c_str());
  }

  std::vector<int> channelAllowedPol;
  getSetting("Channels allowed polarization", channelAllowedPol);
  for(unsigned int i=0; i < channelAllowedPol.size(); i++){
    anita1->pol_allowed[i] = channelAllowedPol.at(i);
  }




  getNextNumberAsString(inputsfile,foutput,number);
  DISCONES=(int)atof(number.c_str());
  getSetting("Exta antennas in trigger", DISCONES);

  getNextNumberAsString(inputsfile,foutput,number);
  anita1->INCLUDE_NADIRONLY=(double)atof(number.c_str());
  getSetting("Nadir only trigger", anita1->INCLUDE_NADIRONLY);

  if (anita1->INCLUDE_NADIRONLY!=0){
    std::cout << "Non-default setting:  INCLUDE_NADIRONLY = " << anita1->INCLUDE_NADIRONLY << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  CHMASKING=atoi(number.c_str());
  getSetting("ANITA-1 channel masking", CHMASKING);

  if (bn1->WHICHPATH!=6 && CHMASKING==1) {
    std::cout << "Cannot include masking for flights other than the ANITA-1 flight." << std::endl;
    std::cout << "For the ANITA-3 channel masking, it is implemented together with the phi masking and it's turned on whenever the PHIMASKING is ON." << std::endl;
    std::cout << "CHMASKING set to 0." << std::endl;
    CHMASKING=0;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  PHIMASKING=atoi(number.c_str());
  getSetting("ANITA-2 channel masking", PHIMASKING);



  thisGetLine(inputsfile,junk);
  foutput << junk << "\n";
  // The following are variables that have been used for modeling anita-lite


  getNextNumberAsString(inputsfile,foutput,number);
  SCALEDOWNLCPRX1=(int)atoi(number.c_str());
  getSetting("Scale down LCP voltage 1st ant", SCALEDOWNLCPRX1);

  getNextNumberAsString(inputsfile,foutput,number);
  SCALEDOWNEPOLRX1=(int)atoi(number.c_str());
  getSetting("Scale down E pol 1st ant", SCALEDOWNEPOLRX1);

  getNextNumberAsString(inputsfile,foutput,number);
  SCALEDOWNHPOLRX1=(int)atoi(number.c_str());
  getSetting("Scale down H pol 1st ant", SCALEDOWNHPOLRX1);

  getNextNumberAsString(inputsfile,foutput,number);
  SCALEDOWNEPOLRX2=(int)atoi(number.c_str());
  getSetting("Scale down E pol 2nd ant", SCALEDOWNEPOLRX2);

  getNextNumberAsString(inputsfile,foutput,number);
  SCALEFACTOREPOLRX2=(double)atof(number.c_str());
  getSetting("E pol scale down factor 2nd ant", SCALEFACTOREPOLRX2);

  getNextNumberAsString(inputsfile,foutput,number);
  SCALEDOWNHPOLRX2=(int)atoi(number.c_str());
  getSetting("H pol scale down factor 2nd ant", SCALEDOWNHPOLRX2);

  getNextNumberAsString(inputsfile,foutput,number);
  EPOLRX2ZERO=(int)atoi(number.c_str());
  getSetting("E pol 2nd ant dead", EPOLRX2ZERO);

  getNextNumberAsString(inputsfile,foutput,number);
  HPOLRX2ZERO=(int)atoi(number.c_str());
  getSetting("H pol 2nd ant dead", HPOLRX2ZERO);

  getNextNumberAsString(inputsfile,foutput,number);
  RCPRX2ZERO=(int)atoi(number.c_str());
  getSetting("RCP 2nd ant dead", RCPRX2ZERO);

  getNextNumberAsString(inputsfile,foutput,number);
  LCPRX2ZERO=(int)atoi(number.c_str());
  getSetting("LCP 2nd ant dead", LCPRX2ZERO);

  if (WHICH==0 && !(SCALEDOWNEPOLRX1==1 && RCPRX2ZERO==1)){
    std::cout << "Non-default setting:  WHICH= " << WHICH << " and EPOLRX2ZERO= " << EPOLRX2ZERO << std::endl;
  }













  thisGetLine(inputsfile,junk);
  foutput << junk << "\n";
  // Modify the following settings to make changes to the signal and noise
  // for testing


  getNextNumberAsString(inputsfile,foutput,number);
  SIGNAL_FLUCT=(int)atoi(number.c_str());
  getSetting("Add noise to signal", SIGNAL_FLUCT);

  if (SIGNAL_FLUCT!=1){
    std::cout << "Non-default setting:  SIGNAL_FLUCT= " << SIGNAL_FLUCT << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  ZEROSIGNAL=atoi(number.c_str()); // whether it's frequency domain (0) or time domain
  getSetting("Zero signal", ZEROSIGNAL);

  getNextNumberAsString(inputsfile,foutput,number);
  RANDOMISEPOL=atoi(number.c_str());
  getSetting("Random rotation polarization", RANDOMISEPOL);

  getNextNumberAsString(inputsfile,foutput,number);
  sig1->SetLPM((int)atoi(number.c_str()));
  int useLPM;
  getSetting("LPM effect", useLPM);
  sig1->SetLPM(useLPM);

  if (sig1->GetLPM()!=1)
    cout << "Non-default setting:  LPM= " << sig1->GetLPM() << "\n";

  getNextNumberAsString(inputsfile,foutput,number);
  double jamieFactor = 0;
  getSetting("E-field factor", jamieFactor);
  sig1->SetJaime_Factor(jamieFactor);


  getNextNumberAsString(inputsfile,foutput,number);
  THERMALNOISE_FACTOR=(double)atof(number.c_str());
  getSetting("Thermal noise factor", THERMALNOISE_FACTOR);

  if (THERMALNOISE_FACTOR!=1){
    std::cout << "Non-default setting:  THERMALNOISE_FACTOR= " << THERMALNOISE_FACTOR << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);

  REMOVEPOLARIZATION=(int)atof(number.c_str());
  getSetting("Disable polarization vectors", REMOVEPOLARIZATION);

  if (REMOVEPOLARIZATION==1){
    std::cout << "Non-default setting:  Polarizations turned off!" << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  anita1->PULSER=atoi(number.c_str());
  getSetting("Use pulser spectrum", anita1->PULSER);


  if (anita1->PULSER!=0){
    std::cout << "Warning! Injecting a pulser spectrum- not simulating neutrinos!  PULSER = "
	      << anita1->PULSER << std::endl;
  }



  getNextNumberAsString(inputsfile,foutput,number);
  bn1->CENTER=atoi(number.c_str());
  getSetting("Centre one phi-sector", bn1->CENTER);

  if (bn1->CENTER!=0){
    std::cout << "WARNING!!  Rotating payload to center one phi sector on the incoming signal for each event."
	      << std::endl;
  }


  getNextNumberAsString(inputsfile,foutput,number);
  ray1->MAKEVERTICAL=atoi(number.c_str());
  getSetting("Force vertical polarization", ray1->MAKEVERTICAL);

  if (ray1->MAKEVERTICAL!=0){
    std::cout << "WARNING!!  Rotating polarization so it is always vertical approaching the payload" << std::endl;
  }









// ################################################################################################
// # Ice properties
// ################################################################################################
// Slopeyness: 0.012  # This determines size of the slopeyness (0.10=5.4, 0.20=7.4 deg mean, 0.012=default)
// Enable slopeyness: 1      # Turn on (1) and off (0) surface slopeyness.  1=default
// Depth dependent refractive index: 1 # Turn on (1) and off (0) depth-dependent index of refraction.  1=default
// Variable attenuation length: 0 # Use constant (1) or variable (0) attenuation length in the ice. 0=default
// Constant ice thickness: 0 # Constant ice thickness (3 km)
// Antarctic ice model: 0 # Which model of the Antarctic ice to use.(0) Crust 2.0  (1) BEDMAP   (BEDMAP is much more finely binned, but will run more slowly.)
// Flat surface: 0 # Set the normal to the surface to always be straight up. (1)  Setting to (0) is default.
// Fixed ice elevation: 0 # Fixed ice elevation
// Medium: 0 # ice (0) or salt(1)
// Enable surface roughness: 0 # Include effects of surface roughness (1) or not (0)  For this, should also have the option "Skip making cuts on neutrinos we couldn't see anyway" set to 1
// Surface roughness: 0.025      # level of surface roughness: for ROUGHNESS==2, this setting is \sigma in note #375.  For ROUGHNESS==1, flatglass (0), 500 grit (1), 1000 grit (2), 1500 grit (3)
// FIRN: 1	# FIRN (1) yes or (0) no
// Which attenuation length: 0 #MOOREBAY (1)use Moore's Bay measured field attenuation length or (0)use South Pole measured data

  thisGetLine(inputsfile,junk);
  foutput << junk << "\n";
  // The following are variables that have been used for modeling anita-lite

  getNextNumberAsString(inputsfile,foutput,number);
  SLOPEYSIZE=(double)atof(number.c_str());
  getSetting("Slopeyness", SLOPEYSIZE);

  if (SLOPEYSIZE!=0.012){
    std::cout << "Non-default setting:  SLOPEYSIZE= " << SLOPEYSIZE << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  SLOPEY=(int)atoi(number.c_str());
  getSetting("Enable slopeyness", SLOPEY);

  if (SLOPEY!=1){
    std::cout << "Non-default setting:  SLOPEY= " << SLOPEY << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  NOFZ=(int)atoi(number.c_str());
  getSetting("Depth dependent refractive index", NOFZ);

  if (NOFZ!=1){
    std::cout << "Non-default setting:  NOFZ= " << NOFZ << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  VARIABLE_ATTEN=(int)atoi(number.c_str());
  getSetting("Variable attenuation length", VARIABLE_ATTEN);

  if (VARIABLE_ATTEN!=0){
    std::cout << "Non-default setting:  VARIABLE_ATTEN= " << VARIABLE_ATTEN << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  CONSTANTICETHICKNESS=(int)atof(number.c_str());
  getSetting("Constant ice thickness", CONSTANTICETHICKNESS);

  if (CONSTANTICETHICKNESS==1){
    std::cout << "Non-default setting:  CONSTANTICETHICKNESS= " << CONSTANTICETHICKNESS << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  ICE_MODEL=(int)atof(number.c_str());
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

  getNextNumberAsString(inputsfile,foutput,number);
  FLATSURFACE=(int)atof(number.c_str());
  getSetting("Flat surface", FLATSURFACE);

  if (FLATSURFACE==1){
    std::cout << "Non-default setting: all surface segments are flat." << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  FIXEDELEVATION=(int)atof(number.c_str());
  getSetting("Fixed ice elevation", FIXEDELEVATION);

  if (FIXEDELEVATION==1){
    std::cout << "Non-default setting:  FIXEDELEVATION= " << FIXEDELEVATION << std::endl;
  }

  //ice or salt
  getNextNumberAsString(inputsfile,foutput,number);
  sig1->SetMedium(atoi(number.c_str()));

  int medium = 0;
  getSetting("Medium", medium);
  sig1->SetMedium(medium);

  getNextNumberAsString(inputsfile,foutput,number);
  ROUGHNESS=(int)atoi(number.c_str());
  getSetting("Enable surface roughness", ROUGHNESS);

  getNextNumberAsString(inputsfile,foutput,number);
  ROUGHSIZE=(double)atof(number.c_str());
  getSetting("Surface roughness", ROUGHSIZE);

  getNextNumberAsString(inputsfile,foutput,number);
  FIRN=atoi(number.c_str());
  getSetting("FIRN", FIRN);

  if (FIRN==0){
    std::cout << "Warning!  Non-standard parameter setting.  FIRN = " << FIRN << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  MOOREBAY=atoi(number.c_str());
  getSetting("Which attenuation length", MOOREBAY);








  thisGetLine(inputsfile,junk);
  foutput << junk << std::endl;
  // The following are variables that have been used for modeling anita-lite

  getNextNumberAsString(inputsfile,foutput,number);
  SIGMA_FACTOR=(double)atof(number.c_str());
  getSetting("Cross-section factor", SIGMA_FACTOR);

  if (SIGMA_FACTOR!=1){
    std::cout << "Non-default setting:  settings->SIGMA_FACTOR= " << SIGMA_FACTOR << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  THETA_TH_FACTOR=(double)atof(number.c_str());
  getSetting("Theta_th factor", THETA_TH_FACTOR);

  if (THETA_TH_FACTOR!=1){
    std::cout << "Non-default setting:  THETA_TH_FACTOR= " << THETA_TH_FACTOR << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  CHANCEINHELL_FACTOR=(double)atof(number.c_str());
  getSetting("Chance in hell factor", CHANCEINHELL_FACTOR);

  if (CHANCEINHELL_FACTOR!=1){
    std::cout << "Non-default setting:  CHANCEINHELL_FACTOR= " << CHANCEINHELL_FACTOR << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  SKIPCUTS=(int)atof(number.c_str());
  getSetting("Skip neutrinos", SKIPCUTS);

  if (SKIPCUTS==1){
    std::cout << "Non-default setting:  Skipping all cuts!" << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  USEDIRECTIONWEIGHTS=(int)atof(number.c_str());
  getSetting("Restrict neutrino directions", USEDIRECTIONWEIGHTS);

  getNextNumberAsString(inputsfile,foutput,number);
  USEPOSITIONWEIGHTS=(int)atof(number.c_str());
  getSetting("Restrict neutrino positions", USEPOSITIONWEIGHTS);

  if (USEPOSITIONWEIGHTS==0){
    std::cout << "Non-default setting:  Not selecting events within the horizon." << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  WEIGHTABSORPTION=(int)atoi(number.c_str());
  getSetting("Weight on absorption", WEIGHTABSORPTION);

  getNextNumberAsString(inputsfile,foutput,number);
  horizontal_banana_points=(int)atof(number.c_str());
  getSetting("Phi points banana", horizontal_banana_points);

  getNextNumberAsString(inputsfile,foutput,number);
  vertical_banana_points=(int)atof(number.c_str());
  getSetting("Theta points banana", vertical_banana_points);

  getNextNumberAsString(inputsfile,foutput,number);
  FORSECKEL=(int)atoi(number.c_str());
  getSetting("Signal across frequencies", FORSECKEL);

  getNextNumberAsString(inputsfile,foutput,number);
  SHOWERTYPE=(int)atoi(number.c_str());
  getSetting("Shower type", SHOWERTYPE);

  getNextNumberAsString(inputsfile,foutput,number);
  BORESIGHTS=atoi(number.c_str());
  getSetting("Loop over boresights", BORESIGHTS);

  if (BORESIGHTS==0){
    std::cout << "Warning!  Non-standard parameter setting.  BORESIGHTS = " << BORESIGHTS << std::endl;
  }







  thisGetLine(inputsfile,junk);
  foutput << junk << std::endl;
  // Interactions

  getNextNumberAsString(inputsfile,foutput,number);
  sig1->SetParameterization(atoi(number.c_str()));
  int askaryanParameterization = 0;
  getSetting("Askaryan parameterization", askaryanParameterization);

  getNextNumberAsString(inputsfile,foutput,number);
  SIGMAPARAM=(int)atoi(number.c_str());
  getSetting("Cross-section parameterization", SIGMAPARAM);

  getNextNumberAsString(inputsfile,foutput,number);
  YPARAM=(int)atoi(number.c_str());
  getSetting("Inelasticity parameterization", YPARAM);

  getNextNumberAsString(inputsfile,foutput,number);
  sec1->SECONDARIES=(int)atoi(number.c_str());
  getSetting("Secondary interactions", sec1->SECONDARIES);

  if (sec1->SECONDARIES!=1){
    std::cout << "Non-default setting:  SECONDARIES= " << sec1->SECONDARIES << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  sec1->TAUDECAY=(int)atoi(number.c_str());
  getSetting("Tau decay as secondary interaction", sec1->TAUDECAY);

  if (sec1->TAUDECAY!=1){
    std::cout << "Non-default setting:  TAUDECAY= " << sec1->TAUDECAY << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  ATMOSPHERE=(int)atoi(number.c_str());
  getSetting("Include atmosphere", ATMOSPHERE);

  if (ATMOSPHERE!=1){
    std::cout << "Non-default setting:  ATMOSPHERE= " << ATMOSPHERE << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  CONSTANTCRUST=(int)atof(number.c_str());
  getSetting("Constant crust density", CONSTANTCRUST);

  if (CONSTANTCRUST==1){
    std::cout << "Non-default setting:  CONSTANTCRUST= " << CONSTANTCRUST << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  CONSTANTY=(int)atoi(number.c_str()); // whether to use contant of 0.2 for y (1) yes or (0) no
  getSetting("Constant y", CONSTANTY);

  getNextNumberAsString(inputsfile,foutput,number);
  bn1->MAXHORIZON=(int)atof(number.c_str()); // max distance from interaction point to horizon
  getSetting("Max interaction distance", bn1->MAXHORIZON);

  getNextNumberAsString(inputsfile,foutput,number);
  taumodes=(int)atoi(number.c_str()); // tau modes (1 for flat distribution for y)
  getSetting("Set tau modes", taumodes);

















  thisGetLine(inputsfile,junk);
  foutput << junk << std::endl;
  // General settings

  getNextNumberAsString(inputsfile,foutput,number);
  WHICHRAYS=(int)atoi(number.c_str());
  getSetting("Which rays", WHICHRAYS);

  if (WHICHRAYS!=1){
    std::cout << "Non-default setting:  WHICHRAYS= " << WHICHRAYS << std::endl;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  WRITE_FILE=(int)atof(number.c_str());
  getSetting("CreateHorizons file", WRITE_FILE);

  if (WRITE_FILE) std::cout<<"Writing CreateHorizons input file." << std::endl;

  //  std::cout << "WRITE_FILE is " << WRITE_FILE << std::endl;

  getNextNumberAsString(inputsfile,foutput,number);
  anita1->SIGMA_THETA=(double)atof(number.c_str());
  getSetting("Theta resolution", anita1->SIGMA_THETA);

  if (anita1->SIGMA_THETA==1){
    std::cout << "Non-default setting:  SIGMA_THETA = 1" << std::endl;
  }
  anita1->SIGMA_THETA*=RADDEG; // immediately convert degrees to radians


  getNextNumberAsString(inputsfile,foutput,number);
  anita1->FREQ_LOW=(double)atof(number.c_str());
  getSetting("Low frequency", anita1->FREQ_LOW);

  if (FREQ_LOW_SEAVEYS>anita1->FREQ_LOW){
    FREQ_LOW_SEAVEYS=anita1->FREQ_LOW;
  }

  getNextNumberAsString(inputsfile,foutput,number);
  anita1->FREQ_HIGH=(double)atof(number.c_str());
  getSetting("High frequency", anita1->FREQ_HIGH);

  BW = anita1->FREQ_HIGH - anita1->FREQ_LOW; // total bandwidth of simulation







  thisGetLine(inputsfile,junk);
  foutput << junk << std::endl;
  // Slac

  getNextNumberAsString(inputsfile,foutput,number);
  SLAC=atoi(number.c_str());
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
    bn1->BN_ALTITUDE=4.22/0.3; // balloon altitude in ft.!!
    foutput << "BN_ALTITUDE=4.22/0.3" << std::endl;
    bn1->RANDOMIZE_BN_ORIENTATION=0; // don't randomize the balloon orientation
    foutput << "RANDOMIZE_BN_ORIENTATION=0" << std::endl;
    SKIPCUTS=1; // don't make chance in hell cuts
    foutput << "SKIPCUTS=1" << std::endl;
    SLACSLOPE=5.8; // slope of the ice in degrees
    foutput << "SLACSLOPE=5.8" << std::endl;
    SLACICELENGTH=5.02; // length of the block of ice
    foutput << "SLACICELENGTH=5.02" << std::endl;
  }





  getNextNumberAsString(inputsfile,foutput,number);
  SLAC_HORIZDIST=(double)atof(number.c_str());
  getSetting("SLAC horizontal distance", SLAC_HORIZDIST);
  // horizontal distance from interaction point to payload


  getNextNumberAsString(inputsfile,foutput,number);
  SLACSLOPE=(double)atof(number.c_str());
  getSetting("SLAC ice slope", SLACSLOPE);
  // slope of the ice surface


  getNextNumberAsString(inputsfile,foutput,number);
  SLACICELENGTH=(double)atof(number.c_str());
  getSetting("SLAC block length", SLACICELENGTH);
  // length of the block of ice


  getNextNumberAsString(inputsfile,foutput,number);
  SLAC_HORIZ_DEPTH=(double)atof(number.c_str());
  getSetting("SLAC interaction depth", SLAC_HORIZ_DEPTH);
  // horizontal distance from interaction point to surface


  SLAC_DEPTH=tan(SLACSLOPE*RADDEG)*(SLACICELENGTH-SLAC_HORIZ_DEPTH) // height from lowest point of ice
    +21.375*CMINCH/100.; // height from beam to lowest point of ice









// ################################################################################################
// # Coherent Sum Trigger settings
// ################################################################################################
// Coherent power threshold: 1300	# power threshold (arbitrary units) for coherent sum trigger and summed power trigger

// ################################################################################################
// # Settings for digitiser path
// ################################################################################################
// Apply inpulse response: 1 # apply impulse response (default for Anita-3 is 1)

// ################################################################################################
// # Settings for time-dependent thresholds
// ################################################################################################
// Time dependent thresholds: 1 # use time-dependent thresholds (only available for Anita-3)






  thisGetLine(inputsfile,junk);
  foutput << junk << std::endl;


  getNextNumberAsString(inputsfile,foutput,number);
  COHERENT_THRESHOLD = double (atof(number.c_str()));
  getSetting("Coherent power threshold", COHERENT_THRESHOLD );


  // default values are 0
  APPLYIMPULSERESPONSEDIGITIZER=0;
  APPLYIMPULSERESPONSETRIGGER=0;
  USETIMEDEPENDENTTHRESHOLDS=0;

  thisGetLine(inputsfile,junk);
  foutput << junk << std::endl;
  getNextNumberAsString(inputsfile,foutput,number);
  APPLYIMPULSERESPONSEDIGITIZER=atoi(number.c_str());
  getSetting("Digitizer path impulse response", APPLYIMPULSERESPONSEDIGITIZER);
  std::cout << "Apply impulse response to digitizer path: " << APPLYIMPULSERESPONSEDIGITIZER << std::endl;



  getNextNumberAsString(inputsfile,foutput,number);
  APPLYIMPULSERESPONSETRIGGER=atoi(number.c_str());
  getSetting("Trigger path impulse response", APPLYIMPULSERESPONSETRIGGER);
  std::cout << "Apply impulse response to trigger path: " << APPLYIMPULSERESPONSETRIGGER << std::endl;

#ifdef ANITA_UTIL_EXISTS
  if ( (APPLYIMPULSERESPONSEDIGITIZER || APPLYIMPULSERESPONSETRIGGER) && WHICH!=8 && WHICH!=9) {
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


  thisGetLine(inputsfile,junk);
  foutput << junk << std::endl;

  getNextNumberAsString(inputsfile,foutput,number);
  USETIMEDEPENDENTTHRESHOLDS=atoi(number.c_str());
  getSetting("Time dependent thresholds", USETIMEDEPENDENTTHRESHOLDS);
  std::cout << "Use time-dependent thresholds: " << USETIMEDEPENDENTTHRESHOLDS << std::endl;

  if ( USETIMEDEPENDENTTHRESHOLDS && WHICH!=9) {
    std::cout << "Time-dependent thresholds are only available for anita-3." << std::endl;
    exit(1);
  }



  thisGetLine(inputsfile,junk);
  foutput << junk << std::endl;


  getNextNumberAsString(inputsfile,foutput,number);
  NOISEFROMFLIGHTDIGITIZER=atoi(number.c_str());
  getSetting("Digitizer noise from flight", NOISEFROMFLIGHTDIGITIZER);
  std::cout << "Use noise from flight for digitizer path: " << NOISEFROMFLIGHTDIGITIZER << std::endl;


  getNextNumberAsString(inputsfile,foutput,number);
  NOISEFROMFLIGHTTRIGGER=atoi(number.c_str());
  getSetting("Trigger noise from flight", NOISEFROMFLIGHTTRIGGER);
  std::cout << "Use noise from flight for trigger path: " << NOISEFROMFLIGHTTRIGGER << std::endl;

#ifdef ANITA_UTIL_EXISTS
  if ( (NOISEFROMFLIGHTDIGITIZER || NOISEFROMFLIGHTTRIGGER) && WHICH!=9) {
    std::cout << "Noise from flight only available for anita-3." << std::endl;
    exit(1);
  }
  if (!APPLYIMPULSERESPONSETRIGGER && NOISEFROMFLIGHTTRIGGER ){
    std::cout << "Noise from flight can only be applied to trigger path if impulse reponse is also used " << std::endl;
    exit(1);
  }
  if (NOISEFROMFLIGHTDIGITIZER || NOISEFROMFLIGHTTRIGGER){
    std::cout << "Noise from flight can only be applied when the Anita tools are sourced." << std::endl;
    exit(1);
  }
#endif


  getNextNumberAsString(inputsfile,foutput,number);
  MINBIAS=atoi(number.c_str());
  getSetting("Min bias", MINBIAS);
  if (MINBIAS){
    std::cout << "Generate Minimum Bias sample: " << MINBIAS << std::endl;
  }




} //method ReadInputs
















void Settings::getSetting(const char* key, int& value){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to find setting " << key << std::endl;
  }
  else{
    // found a match for the key
    value = atoi(it->second.Data());
  }
}

void Settings::getSetting(const char* key, float& value){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to find setting " << key << std::endl;
  }
  else{
    // found a match for the key
    value = atof(it->second.Data());
  }
}

void Settings::getSetting(const char* key, double& value){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to find setting " << key << std::endl;
  }
  else{
    // found a match for the key
    value = atof(it->second.Data());
  }
}

void Settings::getSetting(const char* key, std::vector<int>& valueArray){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to find setting " << key << std::endl;
  }
  else{
    // found a match for the key
    parseValueArray(it->second.Data(), valueArray);
  }
}

void Settings::getSetting(const char* key, std::vector<float>& valueArray){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to find setting " << key << std::endl;
  }
  else{
    // found a match for the key
    parseValueArray(it->second.Data(), valueArray);
  }
}

void Settings::getSetting(const char* key, std::vector<double>& valueArray){

  kvpMap::iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to find setting " << key << std::endl;
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
