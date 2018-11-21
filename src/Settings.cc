#include "Settings.h"

#include <string.h>
#include <iostream>
#include "EnvironmentVariable.h"
#include "TFile.h"

#include "TString.h"
#include "TRegexp.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TRandom3.h"

// Prettify warnings because, why not?
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

ClassImp(icemc::Settings);




/**
 * Default constructor
 *
 */
icemc::Settings::Settings() : jamieFactor(0), medium(0), askaryanParameterization(0)
{
  Initialize();
}


/**
 * Default destructor
 *
 */
icemc::Settings::~Settings() {

}





/**
 * Copies the contents of a yaml settings file into internal memory as strings
 * Later on these strings get turned into ints, floats, doubles...
 *
 * @param fileName the name of the settings file to read
 */
void icemc::Settings::parseSettingsFile(const char* fileName, std::ofstream& outputFile) {


  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  outputFile << "Current date and time are: " << asctime(timeinfo) << std::endl;

  int lineNum = 1;
  wholeSettingsFile = "";


  bool foundFile = false;

  const std::string thisFile = __FILE__;
  const std::string::size_type n = thisFile.rfind("/");  

  const TString prettyFileName = TString(ANSI_COLOR_BLUE) + fileName + ANSI_COLOR_RESET;
  const TString prettySourceFile = TString(ANSI_COLOR_RED) + thisFile.substr(n+1) + ANSI_COLOR_RESET;


  // here we loop through possible prefixes for a text settings file or a root output file containing settings
  TString icemc_src_dir = EnvironmentVariable::ICEMC_SRC_DIR();
  std::vector<TString> prefixes {".", "./config", icemc_src_dir + "/config"};
  for(const auto& prefix : prefixes){

    std::cout << prettySourceFile << " is searching for " << prettyFileName << " in " << prefix << "/ \n";
    
    TString fileName2 = prefix + "/" + fileName;
    std::ifstream testExists(fileName2.Data());
    
    TObjArray* tkns = NULL;
  
    if(fileName2.Contains(".root")){
      TFile* f = TFile::Open(fileName);
      TNamed* n = (TNamed*) f->Get("Settings");

      if(n){
	std::cout << "Found " << prettyFileName << " in " << prefix << "/!"  << std::endl;
	foundFile=true;
	
	TString oldSettings = n->GetTitle();
	tkns = oldSettings.Tokenize("\n");

	for(int i=0; i < tkns->GetEntries(); i++){
	  TObjString* tkn = (TObjString*) tkns->At(i);
	  std::string thisLine(tkn->String().Data());
	  processLine(thisLine, outputFile, fileName, lineNum);
	  lineNum++;
	}
      }

      if(n){
	delete n;
	n = NULL;
      }
      if(f){
	f->Close();
	delete f;
	f = NULL;
      }
    }

    if(!tkns){ // if we never got anything useful from the ROOT file attempt, try to parse as a text file
      std::ifstream settingsFile(fileName2);
      // Print error message if I can't read the file
      if(settingsFile.is_open()){
	std::cout << "Found " << prettyFileName << " in " << prefix << "/!"  << std::endl;
	foundFile=true;
	while(!settingsFile.eof()){
	  std::string thisLine;
	  std::getline(settingsFile, thisLine);
	  processLine(thisLine, outputFile, fileName, lineNum);
	  lineNum++;
	}
      }
    }

    if(tkns){
      delete tkns;
    }

    if(foundFile){
      outputFile << std::endl << std::endl;
      outputFile << __FILE__ << " has finished parsing " << fileName << std::endl;
      outputFile << std::endl << std::endl;
      break;
    }
  }
}



void icemc::Settings::processLine(const std::string& thisLine, std::ofstream& outputFile, const char* fileName, int lineNum){

  // here we keep track of the whole file so we can store it as a TNamed later
  wholeSettingsFile += thisLine + "\n";

  // Copy to silly text output file
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
Bool_t icemc::Settings::newKvpPassesSanityChecks(const TString& key, const TString& value, const char* fileName, int lineNum) const{

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
    kvpMap::const_iterator it = keyValuePairStrings.find(key);

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
void icemc::Settings::printAllKeyValuePairStrings() const {

  kvpMap::const_iterator it;
  for(it = keyValuePairStrings.begin(); it!=keyValuePairStrings.end(); ++it){
    std::cout << it->first << "\t" << it->second << std::endl;
  }
}



/**
 * Set member variables to default values
 *
 */
void icemc::Settings::Initialize() {
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
  PAYLOAD_USE_SPECIFIC_TIME = 0; 
  PAYLOAD_USE_SPECIFIC_TIME_DELTA = 3600; 

  SPECIFIC_NU_POSITION = 0; 
  SPECIFIC_NU_POSITION_LATITUDE = -77.7; 
  SPECIFIC_NU_POSITION_LONGITUDE = 166.7; 
  SPECIFIC_NU_POSITION_ALTITUDE = 0; 
  SPECIFIC_NU_POSITION_DISTANCE = 100e3; 


}









void icemc::Settings::ReadInputs(const char* inputFileName, std::ofstream &foutput){

  parseSettingsFile(inputFileName, foutput);


//################################################################################################
//# Event input output
//################################################################################################

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

  // getSetting("Centre one phi-sector", CENTER);

  // if (CENTER!=0){
  //   std::cout << "WARNING!!  Rotating payload to center one phi sector on the incoming signal for each event."
  // 	      << std::endl;
  // }

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
  
  getSetting("Fixed ice elevation", FIXEDELEVATION);
  if (FIXEDELEVATION==1){
    std::cout << "Non-default setting:  FIXEDELEVATION= " << FIXEDELEVATION << std::endl;
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

  getSetting("Medium", medium);

  getSetting("Enable surface roughness", ROUGHNESS);
  getSetting("Surface roughness", ROUGHSIZE);
  getSetting("Screen edge length [meters]", SCREENEDGELENGTH);

  getSetting("Screen step size [meters]", SCREENSTEPSIZE);

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

  getSetting("Shower type", SHOWERTYPE);







  getSetting("Askaryan parameterization", askaryanParameterization);
  getSetting("Max interaction distance", MAX_HORIZON_DISTANCE);
  if(MAX_HORIZON_DISTANCE != 800e3){
    std::cout << "Non-default settings: MAX_HORIZON_DISTANCE=" << MAX_HORIZON_DISTANCE << std::endl;
  }

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
  getSetting("Set tau modes", taumodes);














  getSetting("Which rays", WHICHRAYS);
  if (WHICHRAYS!=1){
    std::cout << "Non-default setting:  WHICHRAYS= " << WHICHRAYS << std::endl;
  }

  // Set minray/maxray from whichrays, moved from icemc main
  if (WHICHRAYS==1) {
    MINRAY=0;
    MAXRAY=0;
  }
  if (WHICHRAYS==2) {
    MINRAY=0;
    MAXRAY=1;
  } 
  if (WHICHRAYS==3) {
    MINRAY=1;
    MAXRAY=1;
  }

  
  getSetting("CreateHorizons file", WRITE_FILE);
  if (WRITE_FILE){
    std::cout<<"Writing CreateHorizons input file." << std::endl;
  }

  
  getSetting("Min bias", MINBIAS);
  if (MINBIAS){
    std::cout << "Generate Minimum Bias sample: " << MINBIAS << std::endl;
  }



  getSetting("Specific Payload Unix Time", PAYLOAD_USE_SPECIFIC_TIME); 
  getSetting("Specific Payload Delta Time", PAYLOAD_USE_SPECIFIC_TIME_DELTA); 
  
  getSetting("Use Specific Interaction Location", SPECIFIC_NU_POSITION); 
  std::vector<double> specific_place; 
  getSetting("Specific Interaction Location", specific_place); 
  if (specific_place.size() == 3)
  {
    SPECIFIC_NU_POSITION_LATITUDE = specific_place[0];
    SPECIFIC_NU_POSITION_LONGITUDE = specific_place[1]; 
    SPECIFIC_NU_POSITION_ALTITUDE = specific_place[2]; 
  }
     
  getSetting("Specific Interaction Location Maximum Distance", SPECIFIC_NU_POSITION_DISTANCE); 

} //method ReadInputs



  












void icemc::Settings::complainAboutNotFindingKey(const TString& key) const {
  std::cerr << "Warning in " << ANSI_COLOR_BLUE << __FILE__ << ANSI_COLOR_RESET
	    << ", unable to find setting " << ANSI_COLOR_RED << key << ANSI_COLOR_RESET << std::endl;
}


void icemc::Settings::getSetting(const char* key, int& value) const {

  kvpMap::const_iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    value = atoi(it->second.Data());
  }
}

void icemc::Settings::getSetting(const char* key, float& value) const {

  kvpMap::const_iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    value = atof(it->second.Data());
  }
}

void icemc::Settings::getSetting(const char* key, double& value) const { 

  kvpMap::const_iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    value = atof(it->second.Data());
  }
}

void icemc::Settings::processStrings(const std::string& raw, std::vector<std::string >& processed) const {

  bool stringAccumulation = false;
  std::string::const_iterator it = raw.begin();

  // empty the processed vector
  processed.clear();
  processed.push_back("");

  while (it != raw.end()){
    char c = *it++;
    bool escapedDoubleQuote = false;
    // we got an escape character...
    if (c == '\\' && it != raw.end()){
      // c is now the next char, decide what to do based on that
      switch (*it++) {
      case '\\':
	c = '\\'; // this is just a \, but it needs escaping in this source code
	break;
      case 'n' :
	c = '\n';
	break;
      case 't' :
	c = '\t';
	break;
      case '"' :
	c = '"';
	escapedDoubleQuote = true;
	break;

	// Add any other escapes here?
      default: 
	break;
      }
    }

    // if it's a double quote (that's not escaped)
    if(c=='"' && !escapedDoubleQuote){
      // toggle the string accumulation
      stringAccumulation = !stringAccumulation;
    }
    else { // it's a normal character in the string

      // if we're not storing this char, and it's an array delimiter, make a new string
      if(!stringAccumulation && c==','){
	processed.push_back("");
      }
      // if we're storing, just append it to the string...
      else if(stringAccumulation){
	processed.back().append(1, c);
      }
    }
  }
}


void icemc::Settings::getSetting(const char* key, std::string& value) const {
  std::vector<std::string> tempVec;
  getSetting(key, tempVec);
  value = tempVec[0];
}


void icemc::Settings::getSetting(const char* key, std::vector<std::string>& value) const {

  kvpMap::const_iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    std::string raw = it->second.Data();
    processStrings(raw, value);
  }
}



void icemc::Settings::getSetting(const char* key, std::vector<int>& valueArray) const {

  kvpMap::const_iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    parseValueArray(it->second.Data(), valueArray);
  }
}

void icemc::Settings::getSetting(const char* key, std::vector<float>& valueArray) const {

  kvpMap::const_iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    parseValueArray(it->second.Data(), valueArray);
  }
}

void icemc::Settings::getSetting(const char* key, std::vector<double>& valueArray) const {

  kvpMap::const_iterator it = keyValuePairStrings.find(key);
  if(it == keyValuePairStrings.end()){
    complainAboutNotFindingKey(key);
  }
  else{
    // found a match for the key
    parseValueArray(it->second.Data(), valueArray);
  }
}

void icemc::Settings::parseValueArray(const char* valueString, std::vector<int>& values) const{
  TString theValueString(valueString);

  TObjArray* theValues = theValueString.Tokenize(",");
  for(int i=0; i < theValues->GetEntries(); ++i){

    TObjString* token = (TObjString*) theValues->At(i);
    int value = atoi(token->GetString().Data());
    values.push_back(value);
  }
}

void icemc::Settings::parseValueArray(const char* valueString, std::vector<float>& values) const{
  TString theValueString(valueString);

  TObjArray* theValues = theValueString.Tokenize(",");
  for(int i=0; i < theValues->GetEntries(); ++i){

    TObjString* token = (TObjString*) theValues->At(i);
    float value = atof(token->GetString().Data());
    values.push_back(value);
  }
}

void icemc::Settings::parseValueArray(const char* valueString, std::vector<double>& values) const{
  TString theValueString(valueString);

  TObjArray* theValues = theValueString.Tokenize(",");
  for(int i=0; i < theValues->GetEntries(); ++i){

    TObjString* token = (TObjString*) theValues->At(i);
    double value = atof(token->GetString().Data());
    values.push_back(value);
  }
}


TNamed* icemc::Settings::makeRootSaveableSettings() const {
  return new TNamed("Settings", wholeSettingsFile.Data());
}
