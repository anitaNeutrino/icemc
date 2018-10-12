#include "TChain.h"
#include "TFile.h"
#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "Crust2.h"
#include "Antarctica.h"
#include "TVector3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMath.h"

#include "Geoid.h"
#include "anita.hh"
#include "RayTracer.h"
#include "FancyTTreeInterpolator.h"

#include "balloon.hh"
#include "Settings.h"
#include "ConnollyEtAl2011.h"
#include "EnvironmentVariable.h"
#include "Report.h"

std::ostream& operator<<(std::ostream& os, const icemc::FlightPath& fp){
  switch (fp){
  case icemc::FlightPath::FixedPosition:
    return os << "FlightPath::FixedPosition";
  case icemc::FlightPath::Circle80DegreesSouth:
    return os << "FlightPath::Circle80DegreesSouth";
  case icemc::FlightPath::AnitaLite:
    return os << "FlightPath::AnitaLite";
  case icemc::FlightPath::Custom:
    return os << "FlightPath::Custom";
  case icemc::FlightPath::Anita1:
    return os << "FlightPath::Anita1";
  case icemc::FlightPath::Anita2:
    return os << "FlightPath::Anita2";
  case icemc::FlightPath::Anita3:
    return os << "FlightPath::Anita3";
  case icemc::FlightPath::Anita4:
    return os << "FlightPath::Anita4";
  default:
    return os << "Unknown FlightPath";
  }
}



// ANITA-2 (WHICHPATH 7) analysis fitted fixed const pitch (-0.29) and const roll (0.89)
constexpr double fixedAnita2Pitch = -0.29; // degrees
constexpr double fixedAnita2Roll  =  0.89; // degrees
// ANITA-3 (WHICHPATH 8) analysis forced pitch and roll to be zero
constexpr double fixedAnita3Pitch =  0.00; // degrees
constexpr double fixedAnita3Roll  =  0.00; // degrees
// @todo have we got ANITA-4 pitch and roll numbers?


double icemc::Balloon::getPitch() const {
  switch(WHICHPATH){
  case FlightPath::Anita2:
    return fixedAnita2Pitch;
  case FlightPath::Anita3:
    return fixedAnita3Pitch;
  default:
    return pitch;
  }
}

double icemc::Balloon::getRoll() const {
  switch(WHICHPATH){
  case FlightPath::Anita2:
    return fixedAnita2Roll;
  case FlightPath::Anita3:
    return fixedAnita3Roll;
  default:
    return roll;
  }
}


icemc::Balloon::Balloon(icemc::FlightPath path, const Settings* settings) : fSettings(settings), WHICHPATH(path)
{
  InitializeBalloon(nullptr);
}

icemc::Balloon::Balloon(const Settings* settings)
  : WHICHPATH(settings ? static_cast<icemc::FlightPath>(settings->WHICHPATH) : FlightPath::FixedPosition){

  if(!settings){
    icemc::report() << severity::warning << __PRETTY_FUNCTION__ << " was given nullptr for const Settings*. "
	       << "Assuming " << FlightPath::FixedPosition << std::endl;
  }
  InitializeBalloon(settings);
}

icemc::Balloon::~Balloon(){
  if(fInterp) delete fInterp;
  
  // icemc::report() << severity::info << __PRETTY_FUNCTION__ << " for " << WHICHPATH << std::endl;
}







// void icemc::Balloon::SetDefaultBalloonPosition() { // position of surface of earth under balloon
    
//   // // set the default balloon position
//   // // if you are using real Anita-lite path, these get overwritten for each event
    
//   if(BN_LATITUDE==999){
//     fPosition.SetLonLatAlt(0, -90, 40e3);
//   }
//   else{
//     if (BN_ALTITUDE==0){ // if the altitude isn't set in the input file
//       altitude_bn=120000*12.*constants::CMINCH/100.; // 120000 ft.=36.6 m
//     }
//     else{
//       altitude_bn=BN_ALTITUDE*12.*constants::CMINCH/100.; // converts the altitude in the input file to meters
//     }
//     fPosition.SetLonLatAlt(BN_LONGITUDE, BN_LATITUDE, altitude_bn);
//   }  
// } // set default balloon position


void icemc::Balloon::ReadAnitaliteFlight() {

  const std::string ICEMC_SRC_DIR = icemc::EnvironmentVariable::ICEMC_SRC_DIR();
  const std::string ICEMC_DATA_DIR = ICEMC_SRC_DIR+"/data/";  
  const std::string anitaliteflight = ICEMC_DATA_DIR+"/BalloonGPS.txt"; // the gps path of the anita-lite flight
  std::ifstream flightfile(anitaliteflight.c_str()); //set file to the right one
  if (!flightfile) {
    std::cout << "Flight file not found.\n";
    exit(1);
  }
    
    


  for(int line=0; line < 4; line++){
    std::string junk;
    getline(flightfile,junk); // first few lines are junk
  }
    
  int NPOINTS=0; // number of points we have for the flight path

  const int plentyLong = 100000;
  latitude_bn_anitalite.reserve(plentyLong);
  longitude_bn_anitalite.reserve(plentyLong);
  altitude_bn_anitalite.reserve(plentyLong);
  heading_bn_anitalite.reserve(plentyLong);
  realtime_bn_anitalite.reserve(plentyLong);
  
  while (!flightfile.eof()) {
    std::string slatitude; // latitude in string format
    std::string slongitude; //longitude in string format
    std::string saltitude; // altitude in string format
    std::string sheading; // heading in string format
    std::string srealtime; // real life unix time in string format    
    std::string junk;

    double latitude; // same quantities as above but in double format
    double longitude;
    double altitude;
    double heading;
    double realtime;
    
    flightfile >> srealtime >> junk >> slatitude >> slongitude >> saltitude >> sheading >> junk >> junk >> junk;
    latitude  = (double)atof(slatitude.c_str());
    longitude = (double)atof(slongitude.c_str());
    altitude  = (double)atof(saltitude.c_str());
    heading   = (double)atof(sheading.c_str());
    realtime  = (double)atof(srealtime.c_str());

    // latitude and longitude of each balloon position
    latitude_bn_anitalite.push_back(latitude);
    longitude_bn_anitalite.push_back(longitude);
    altitude_bn_anitalite.push_back(altitude);
    heading_bn_anitalite.push_back(heading);
    realtime_bn_anitalite.push_back(realtime);
    NPOINTS++;
		
    getline(flightfile,junk);
		
  }//while

  // plenty long is probably too long

  int NPOINTS_MAX=NPOINTS-140; // exclude the fall
  // int NPOINTS_MIN=0;

  while(latitude_bn_anitalite.size() > NPOINTS_MAX){
    latitude_bn_anitalite.pop_back();
    longitude_bn_anitalite.pop_back();
    altitude_bn_anitalite.pop_back();
    heading_bn_anitalite.pop_back();
    realtime_bn_anitalite.pop_back();
  }

  latitude_bn_anitalite.shrink_to_fit();
  longitude_bn_anitalite.shrink_to_fit();
  altitude_bn_anitalite.shrink_to_fit();
  heading_bn_anitalite.shrink_to_fit();
  realtime_bn_anitalite.shrink_to_fit();
  
  fFirstRealTime = realtime_bn_anitalite.front();
  fLastRealTime = realtime_bn_anitalite.back();

}//ReadAnitaliteFlight



void icemc::Balloon::InitializeBalloon(const Settings* settings) {

  
  MAXHORIZON=800000.; // pick the interaction within this distance from the balloon so that it is within the horizon
  ibnposition=0;
  igps=0;
  BN_LONGITUDE=999;   // balloon longitude for fixed balloon location
  BN_LATITUDE=999;    // balloon latitude for fixed balloon location
  if(settings){
    BN_LATITUDE              = settings->BN_LATITUDE;
    BN_LONGITUDE             = settings->BN_LONGITUDE;
    BN_ALTITUDE              = settings->BN_ALTITUDE;
    RANDOMIZE_BN_ORIENTATION = settings->RANDOMIZE_BN_ORIENTATION;
    MAXHORIZON               = settings->MAXHORIZON;    
  }

  
  for (int i=0;i<9;i++) {
    for (int j=0;j<2;j++) {
      surfTrigBandMask[i][j]=0;
    }
    for (int j=0;j<32;j++) {
      powerthresh[i][j]=3.E-12;
      meanp[i][j]=1.E-12;
    }
  }

  
  const std::string ICEMC_SRC_DIR = icemc::EnvironmentVariable::ICEMC_SRC_DIR();
  const std::string ICEMC_DATA_DIR = ICEMC_SRC_DIR+"/data/";


  // GPS positions of Anita or Anita-lite balloon flight
  if (WHICHPATH==FlightPath::AnitaLite){
    ReadAnitaliteFlight();
  }


  MINALTITUDE=30e3; // balloon has to be 30 km altitude at least for us to read the event from the flight data file

  // initialisation of igps_previous
  if (WHICHPATH==FlightPath::AnitaLite || WHICHPATH==FlightPath::Anita1 || WHICHPATH==FlightPath::Anita2 || WHICHPATH==FlightPath::Anita3 || WHICHPATH==FlightPath::Anita4){
    igps_previous=0; // which entry from the flight data file the previous event was
  }

  
  // This for Anita 1 flight
  if (WHICHPATH==FlightPath::Anita1) {

    fChain = new TChain("foricemc");
    fChain->Add((ICEMC_DATA_DIR+"/anita1flightdata.root").c_str());

    fChain->SetBranchAddress("surfTrigBandMask",surfTrigBandMask);
    fChain->SetBranchAddress("powerthresh",powerthresh);
    fChain->SetBranchAddress("meanp",meanp);
    fChain->SetBranchAddress("longitude",&flongitude);
    fChain->SetBranchAddress("latitude",&flatitude);
    fChain->SetBranchAddress("altitude",&faltitude);
    const char* whichRealTime = "realTime_surfhk"; // there are several realTime variables in this tree
    fChain->SetBranchAddress(whichRealTime,&realTime);
    fChain->SetBranchAddress("heading",&fheading);

    fChain->BuildIndex(whichRealTime);

    Long64_t firstGoodEntry = -1;
    faltitude = 0;
    while(faltitude < MINALTITUDE || heading < 0){
      firstGoodEntry++;
      fChain->GetEntry(firstGoodEntry);
      fFirstRealTime = realTime;
    }

    faltitude = 0;
    Long64_t lastGoodEntry = fChain->GetEntries();
    while(faltitude < MINALTITUDE || heading < 0){
      lastGoodEntry--;
      fChain->GetEntry(lastGoodEntry);
      fLastRealTime = realTime;
    }
    // fInterp = new FancyTTreeInterpolator(fChain,  whichRealTime);
    // TString cut = TString::Format("Entry$ >= %lld && Entry$ < %lld", firstGoodEntry, lastGoodEntry);
    // fInterp->add("heading", cut, 360);
    // fInterp->add("longitude", cut);
    // fInterp->add("latitude", cut);
    // fInterp->add("altitude", cut);

    std::cout << "Loaded chain " << whichPath() << " with first good entry " << firstGoodEntry << " at " << fFirstRealTime << " and last good entry " << lastGoodEntry << " at " << fLastRealTime << std::endl;
  }

  else if (WHICHPATH==FlightPath::Anita2 || WHICHPATH==FlightPath::Anita3 || WHICHPATH==FlightPath::Anita4) { // for anita-3 and 4 flights

    TString balloonFile = ICEMC_DATA_DIR;
    switch(WHICHPATH){
    case FlightPath::Anita2:
      balloonFile += "/anita2gps_pitchandroll.root";
      break;
    case FlightPath::Anita3:
      balloonFile += "/anita3gps_pitchroll.root";
      break;
    case FlightPath::Anita4:
    default:
      balloonFile += "/anita4gps_pitchroll.root";
      break;
    }

    
    fChain = new TChain("adu5PatTree");
    
    // fChain->SetMakeClass(1);
    fChain->Add(balloonFile);//created to include pitch and roll.
    fChain->SetBranchAddress("longitude",&flongitude);
    fChain->SetBranchAddress("latitude",&flatitude);
    fChain->SetBranchAddress("altitude",&faltitude);
    fChain->SetBranchAddress("heading",&fheading);
    fChain->SetBranchAddress("realTime",&realTime);

    if(WHICHPATH==FlightPath::Anita2){// someone was really stupid
      fChain->SetBranchAddress("pitch",&pitch);
      fChain->SetBranchAddress("roll",&roll);
    }
    else{
      fChain->SetBranchAddress("pitch",&fpitch);
      fChain->SetBranchAddress("roll",&froll);
    }
    
    fChain->BuildIndex("realTime");

    Long64_t firstGoodEntry = -1;
    faltitude = 0;
    while(faltitude < MINALTITUDE || heading < 0){
      firstGoodEntry++;
      fChain->GetEntry(firstGoodEntry);
      fFirstRealTime = realTime;      
    }

    faltitude = 0;
    Long64_t lastGoodEntry = fChain->GetEntries();
    while(faltitude < MINALTITUDE || heading < 0){
      lastGoodEntry--;
      fChain->GetEntry(lastGoodEntry);
      fLastRealTime = realTime;
    }

    // fInterp = new FancyTTreeInterpolator(fChain,"realTime");
    // TString cut = TString::Format("Entry$ >= %lld && Entry$ < %lld", firstGoodEntry, lastGoodEntry);
    // fInterp->add("heading", cut, 360);
    // fInterp->add("longitude", cut, 180, -180);
    // fInterp->add("latitude", cut);
    // fInterp->add("altitude", cut);

    std::cout << "Loaded chain " << whichPath() << " with first good entry " << firstGoodEntry << " at " << fFirstRealTime << " and last good entry " << lastGoodEntry << " at " << fLastRealTime << std::endl;
  }
  
  // REDUCEBALLOONPOSITIONS=100;
}





// for tuffs for anita-4
int getTuffIndex(int Curr_time) {
  // all the Tuffconfig stuff is in the icemc::constants namespace
  // but the following gets a bit verbose if you actually specify it
  using namespace icemc::constants;
  
  if((TUFFconfig_B_end_3 < Curr_time) && (Curr_time <= TUFFconfig_A_end_1)) {// config A trigconfigA.imp
    return 0;
  }
  else if(((0 < Curr_time) && (Curr_time <= TUFFconfig_B_end_1)) || ((TUFFconfig_P_end_3 < Curr_time) && (Curr_time <= TUFFconfig_B_end_2)) || ((TUFFconfig_P_end_4 < Curr_time) && (Curr_time <= TUFFconfig_B_end_3)) || ((TUFFconfig_A_end_1 < Curr_time) && (Curr_time <= TUFFconfig_B_end_4)) || ((TUFFconfig_P_end_5 < Curr_time) && (Curr_time <= TUFFconfig_B_end_5)) || ((TUFFconfig_P_end_6 < Curr_time) && (Curr_time <= TUFFconfig_B_end_6)) || (TUFFconfig_P_end_7 < Curr_time) ) { // config B trigconfigB.imp
    return 1;
  }
  else if((TUFFconfig_P_end_1 < Curr_time) && (Curr_time <= TUFFconfig_C_end_1)) { // config C trigconfigC.imp
    return 2;
  }
  else if( ((TUFFconfig_P_end_2 < Curr_time) && (Curr_time <= TUFFconfig_G_end_1)) || ((TUFFconfig_O_end_1 < Curr_time) && (Curr_time <= TUFFconfig_G_end_2)) ) { // config G trigconfigG.imp
    return 3;
  }
  else if( ((TUFFconfig_G_end_1 < Curr_time) && (Curr_time <= TUFFconfig_O_end_1)) || ((TUFFconfig_G_end_2 < Curr_time) && (Curr_time <= TUFFconfig_O_end_2)) ) { // config O trigconfigO.imp
    return 4;
  }
  else if( ((TUFFconfig_B_end_1 < Curr_time) && (Curr_time <= TUFFconfig_P_end_1)) || ((TUFFconfig_C_end_1 < Curr_time) && (Curr_time <= TUFFconfig_P_end_2)) || ((TUFFconfig_O_end_2 < Curr_time) && (Curr_time <= TUFFconfig_P_end_3)) || ((TUFFconfig_B_end_2 < Curr_time) && (Curr_time <= TUFFconfig_P_end_4)) || ((TUFFconfig_B_end_4 < Curr_time) && (Curr_time <= TUFFconfig_P_end_5)) || ((TUFFconfig_B_end_5 < Curr_time) && (Curr_time <= TUFFconfig_P_end_6)) || ((TUFFconfig_B_end_6 < Curr_time) && (Curr_time <= TUFFconfig_P_end_7)) ) { // config P trigconfigP.imp
    return 5;
  }
  icemc::report() << icemc::severity::warning << __PRETTY_FUNCTION__
	     << " could not get TUFF index from current time "
	     << Curr_time << ", returning -1." << std::endl;
  return -1;
}


// this is called for each neutrino
// void icemc::Balloon::PickBalloonPosition(const Antarctica *antarctica1, const Settings *settings1, int inu, Anita *anita1, double randomNumber) {
// void icemc::Balloon::PickBalloonPosition(const Settings *settings1, int inu, Anita *anita1, double randomNumber) {
// void icemc::Balloon::PickBalloonPosition(double eventTime, const Settings* settings1 ,Anita *anita1) {
void icemc::Balloon::getBalloonPosition(double eventTime, Anita *anita1) {  


  pitch=0.;
  roll=0.;
  Long64_t entry = 0;
  if (WHICHPATH==FlightPath::AnitaLite ||
      WHICHPATH==FlightPath::Anita1 ||
      WHICHPATH==FlightPath::Anita2 ||
      WHICHPATH==FlightPath::Anita3 ||
      WHICHPATH==FlightPath::Anita4) {

    if (WHICHPATH==FlightPath::AnitaLite) {
      
      auto it = std::upper_bound(realtime_bn_anitalite.begin(),  realtime_bn_anitalite.end(), eventTime);
      
      int igps = std::distance(realtime_bn_anitalite.begin(), it);
      flatitude=(float)latitude_bn_anitalite.at(igps);
      flongitude=(float)longitude_bn_anitalite.at(igps);
      faltitude=(float)altitude_bn_anitalite.at(igps);
      fheading=(float)heading_bn_anitalite.at(igps);
    }
    else if (WHICHPATH==FlightPath::Anita1 ||
	     WHICHPATH==FlightPath::Anita2 ||
	     WHICHPATH==FlightPath::Anita3 ||
	     WHICHPATH==FlightPath::Anita4) {

      entry = fChain->GetEntryNumberWithBestIndex((UInt_t)eventTime);
      fChain->GetEntry(entry);
      
      // For Anita 1 and Anita 2 and Anita 3:
      // igps = (igps_previous+1)%fChain->GetEntries(); // pick which event in the tree we want
      // static int start_igps = 0; 
      // static int ngps = fChain->GetEntries(); 
      // static int init_best = 0;


      ///@todo restore payload use specific time
      // if (fSettings->PAYLOAD_USE_SPECIFIC_TIME && !init_best) 
      // {
      //    int N = fChain->Draw("realTime","","goff"); 
      //    double * times = fChain->GetV1(); 

      //    int best_igps =  TMath::BinarySearch(N, times, (double) fSettings->PAYLOAD_USE_SPECIFIC_TIME); 
      //    start_igps = best_igps;
      //    int end_igps = best_igps;

      //    while (times[start_igps] > fSettings->PAYLOAD_USE_SPECIFIC_TIME - fSettings->PAYLOAD_USE_SPECIFIC_TIME_DELTA)
      //    {
      //      start_igps--; 
      //    }

      //    while (times[end_igps] < fSettings->PAYLOAD_USE_SPECIFIC_TIME + fSettings->PAYLOAD_USE_SPECIFIC_TIME_DELTA)
      //    {
      //      end_igps++; 
      //    }

      //    ngps = end_igps - start_igps + 1; 
      //    init_best = 1; 
      // }

      // igps = start_igps + int(randomNumber*ngps); // use random position 

      //////////////////////////// TEMPORARY HACKS FOR ANITA4 !!!!!!      
      // if (WHICHPATH==FlightPath::Anita4 && ((igps>870 && igps<880) || (igps>7730 && igps<7740) || (igps>23810 && igps<23820) || (igps>31630 && igps<31660)) || (igps==17862) ){
      // 	igps = igps+30;
      // }
      
      // fChain->GetEvent(igps); // this grabs the balloon position data for this event
      if(fSettings && anita1 && fSettings->TUFFSON){
	anita1->tuffIndex = getTuffIndex(realTime);
      }// end if tuffson 


      ///@todo figure out if I need to switch this back on!
      // while (faltitude<MINALTITUDE || fheading<0) { // if the altitude is too low, pick another event.

      // 	igps++; // increment by 1
      // 	igps=igps%fChain->GetEntries(); // make sure it's not beyond the maximum entry number

      // 	fChain->GetEvent(igps);	  // get new event
      // }
      
      // set phi Masking for Anita 2 or Anita 3
      // the deadtime is read from the same tree
      if ((WHICHPATH==FlightPath::Anita2 ||
	   WHICHPATH==FlightPath::Anita3 ||
	   WHICHPATH==FlightPath::Anita4) &&
	  fSettings && anita1 && 
	  (fSettings->PHIMASKING==1 || fSettings->USEDEADTIME)){
	anita1->setphiTrigMask(realTime);
      }
      if ((WHICHPATH==FlightPath::Anita3 || WHICHPATH==FlightPath::Anita4) && fSettings && anita1 && fSettings->USETIMEDEPENDENTTHRESHOLDS==1){ // set time-dependent thresholds
	anita1->setTimeDependentThresholds(realTime);
      }
    }
    igps_previous = igps;
    latitude  = fInterp ? fInterp->interp("latitude", eventTime) : (double) flatitude;    
    longitude = fInterp ? fInterp->interp("longitude", eventTime) : (double) flongitude;
    altitude  = fInterp ? fInterp->interp("altitude", eventTime) : (double) faltitude;
    heading   = fInterp ? fInterp->interp("heading", eventTime) : (double) fheading;
    // latitude  = (double) flatitude;
    // longitude = (double) flongitude;
    // altitude  = (double) faltitude;
    // heading   = (double) fheading;


    if(WHICHPATH!=FlightPath::Anita2){
      // someone set these as doubles in the A2 tree
      // but everything else is a float... (I have some
      // opinions about that, putting those to one side)
      // for A2 we read the doubles directly so no need
      // to convert
      roll      = (double) froll;
      pitch     = (double) fpitch;
    }
    

    static double lastLat =  -999;
    static double lastLon =  -999;
    static UInt_t lastTime = 0;
    static Long64_t lastEntry = -1;
    if(lastLat != -999 && lastLon != -999){
      // if(fabs(latitude - lastLat) > 0.2 || (fabs(longitude - lastLon) > 1 && fabs(longitude - lastLon) < 359)){ 
      // 	std::cout << lastEntry << "\t" << lastTime << "\t" << lastLon << "\t" << lastLat << std::endl;
      // 	std::cout << entry << "\t" << realTime << "\t" << longitude << "\t" << latitude << std::endl;	
      // 	std::cout << entry - lastEntry << "\t" << realTime - lastTime << "\t" << longitude - lastLon << "\t" << latitude - lastLat << std::endl;
      // 	std::cout << std::endl << std::endl;;	
      // 	// std::cout << entry << "\t" << dt << "\t" << realTime << "\t" << eventTime << std::endl;
      // }
    }

    lastLat = latitude;
    lastLon = longitude;
    lastTime = realTime;
    lastEntry = entry;
    

    if (WHICHPATH==FlightPath::AnitaLite){
      altitude_bn=altitude*12.*constants::CMINCH/100.;
    }
    else{
      altitude_bn=altitude; // get the altitude of the balloon in the right units
    }

    fPosition.SetLonLatAlt(longitude, latitude, altitude_bn);

    
  }

  else if (WHICHPATH==FlightPath::Circle80DegreesSouth){ // pick random phi at 80 deg S
    longitude = gRandom->Rndm()*360;
    latitude = -80;
    fPosition.SetLonLatAlt(longitude, latitude, 40e3);
    igps=0;
  } //else if(random position at 80 deg S)
  else if (WHICHPATH==FlightPath::FixedPosition){
		
    igps=0; // set gps position to be zero - we're at a fixed balloon position
		
  }
  else if (WHICHPATH==FlightPath::Custom) {
    if (BN_ALTITUDE!=0){
      altitude_bn=BN_ALTITUDE*12.*constants::CMINCH/100.; // set the altitude of the balloon to be what you pick.  This isn't in time for CreateHorizons though!
    }

    icemc::report() << severity::error << FlightPath::Custom << " is currently broken!" << std::endl;
    fPosition.SetLonLatAlt(0, -90, altitude_bn);
    
    // surface_under_balloon = antarctica1->Surface(r_bn);
    // r_bn_shadow = surface_under_balloon * r_bn.Unit();
    //    r_bn = (antarctica->Geoid(r_bn)+altitude_bn) * r_bn.Unit();
		
    // r_bn = (antarctica1->Surface(r_bn)+altitude_bn) * r_bn.Unit();
  } // you pick it
  else{
    // icemc::report() << severity::error << "Can't get position for " << static_cast<int>(whichPath()) << std::endl;
  }
  // if (!fSettings->UNBIASED_SELECTION && dtryingposition!=-999){
  //   dtryingposition=antarctica1->GetBalloonPositionWeight(ibnposition);
  // }
  // else{
  //   dtryingposition=1.;
  // }
    
    
} 








#ifdef ANITA_UTIL_EXISTS
Adu5Pat icemc::Balloon::pat() const{
  Adu5Pat pat;
  pat.latitude  = getLatitude();
  pat.longitude = getLongitude();
  pat.altitude  = getAltitude();
  pat.realTime  = getRealTime();
  pat.heading   = getHeading();
  pat.pitch     = getPitch();
  pat.roll      = getRoll();
  pat.run       = -1;
  return pat;
}
#endif
