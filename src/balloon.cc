#include "TChain.h"
#include "TFile.h"
#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "Earth.h"
#include "Antarctica.h"
#include "TVector3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMath.h"

#include "GeoidModel.h"
#include "anita.hh"
#include "RayTracer.h"

#include "balloon.hh"
#include "Settings.h"
#include "Primaries.h"
#include "EnvironmentVariable.h"
#include "IcemcLog.h"

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


icemc::Balloon::Balloon(const Settings* settings)
  : WHICHPATH(settings ? static_cast<icemc::FlightPath>(settings->WHICHPATH) : FlightPath::FixedPosition){

  
  MAXHORIZON=800000.; // pick the interaction within this distance from the balloon so that it is within the horizon
  ibnposition=0;
  igps=0;
  BN_LONGITUDE=999;   // balloon longitude for fixed balloon location
  BN_LATITUDE=999;    // balloon latitude for fixed balloon location
  if(!settings){
    icemcLog() << icemc::warning << __PRETTY_FUNCTION__ << " was given nullptr for const Settings*. "
	       << "Assuming " << FlightPath::FixedPosition << std::endl;
  }
  else{
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

  InitializeBalloon();
}





void icemc::Balloon::SetDefaultBalloonPosition(const Antarctica *antarctica1) { // position of surface of earth under balloon
    
  // set the default balloon position
  // if you are using real Anita-lite path, these get overwritten for each event
  //std::cout << "BN_LATITUDE is " << BN_LATITUDE << "\n";
    
  int BN_LATITUDE_SETTING = BN_LATITUDE;
  int BN_LONGITUDE_SETTING = BN_LONGITUDE;
    
  if(BN_LATITUDE_SETTING==999){
    theta_bn=10*constants::RADDEG; // wrt south pole
  }
  else{
    theta_bn=(90-BN_LATITUDE_SETTING)*constants::RADDEG;
  }
  
  if(BN_LONGITUDE_SETTING==999){
    phi_bn = constants::PI/4; //wrt 90E longitude
  }
  else {
    phi_bn = Earth::LongtoPhi_0isPrimeMeridian(BN_LONGITUDE_SETTING); //remember input of LongtoPhi is between -180 and 180
  }    
    
  r_bn = GeoidModel::Position(0, 0, 1);
  r_bn.SetTheta(theta_bn);
  r_bn.SetPhi(phi_bn); // direction of balloon- right now this is a unit vector

  if (BN_ALTITUDE==0){ // if the altitude isn't set in the input file
    altitude_bn=120000*12.*constants::CMINCH/100.; // 120000 ft.=36.6 m
  }
  else{
    altitude_bn=BN_ALTITUDE*12.*constants::CMINCH/100.; // converts the altitude in the input file to meters
  }    
  surface_under_balloon = antarctica1->Surface(r_bn); // distance between center of earth and surface under balloon

  r_bn_shadow = surface_under_balloon * r_bn.Unit(); // position of surface under balloon

  r_bn = (antarctica1->Geoid(r_bn)+altitude_bn) * r_bn; // position of balloon

} // set default balloon position


void icemc::Balloon::ReadAnitaliteFlight() {

  const std::string ICEMC_SRC_DIR = icemc::EnvironmentVariable::ICEMC_SRC_DIR();
  const std::string ICEMC_DATA_DIR = ICEMC_SRC_DIR+"/data/";  
  const std::string anitaliteflight = ICEMC_DATA_DIR+"/BalloonGPS.txt"; // the gps path of the anita-lite flight
  std::ifstream flightfile(anitaliteflight.c_str()); //set file to the right one
  if (!flightfile) {
    std::cout << "Flight file not found.\n";
    exit(1);
  }
    
  std::string slatitude; // latitude in string format
  std::string slongitude; //longitude in string format
  std::string saltitude; // altitude in string format
  std::string sheading; // heading in string format
  std::string srealtime; // real life unix time in string format
  std::string junk; // string not used for anything
    
  std::string line; // for reading a whole line
  double latitude; // same quantities as above but in double format
  double longitude;
  double altitude;
  double heading;
  double realtime;
    
  getline(flightfile,junk); // first few lines are junk
  getline(flightfile,junk);
  getline(flightfile,junk);
  getline(flightfile,junk);
    
  NPOINTS=0; // number of points we have for the flight path

  const int plentyLong = 100000;
  latitude_bn_anitalite.reserve(plentyLong);
  longitude_bn_anitalite.reserve(plentyLong);
  altitude_bn_anitalite.reserve(plentyLong);
  heading_bn_anitalite.reserve(plentyLong);
  realtime_bn_anitalite.reserve(plentyLong);
  
  while (!flightfile.eof()) {
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
    // latitude_bn_anitalite[NPOINTS]=latitude;
    // longitude_bn_anitalite[NPOINTS]=longitude;
    // altitude_bn_anitalite[NPOINTS]=altitude;
    // heading_bn_anitalite[NPOINTS]=heading;
    // realtime_bn_anitalite[NPOINTS]=realtime;
		
    NPOINTS++;
		
    //    std::cout << "NPOINTS, altitude are " << NPOINTS << " " << altitude << "\n";
		
    getline(flightfile,junk);
		
  }//while

  // plenty long is probably too long
  latitude_bn_anitalite.shrink_to_fit();
  longitude_bn_anitalite.shrink_to_fit();
  altitude_bn_anitalite.shrink_to_fit();
  heading_bn_anitalite.shrink_to_fit();
  realtime_bn_anitalite.shrink_to_fit();
    
  NPOINTS_MAX=NPOINTS-140; // exclude the fall
  NPOINTS_MIN=0;

}//ReadAnitaliteFlight



void icemc::Balloon::InitializeBalloon() {
  const std::string ICEMC_SRC_DIR = icemc::EnvironmentVariable::ICEMC_SRC_DIR();
  const std::string ICEMC_DATA_DIR = ICEMC_SRC_DIR+"/data/";
  const std::string anitaflight = ICEMC_DATA_DIR+"/anitagps.txt";// gps path of anita flight    

  // GPS positions of Anita or Anita-lite balloon flight
  if (WHICHPATH==FlightPath::AnitaLite){
    ReadAnitaliteFlight();
  }

  MINALTITUDE=30000; // balloon has to be 30 km altitude at least for us to read the event from the flight data file

  // initialisation of igps_previous
  if (WHICHPATH==FlightPath::Anita1 || WHICHPATH==FlightPath::Anita2 || WHICHPATH==FlightPath::Anita3 || WHICHPATH==FlightPath::Anita4){
    igps_previous=0; // which entry from the flight data file the previous event was
  }
  if (WHICHPATH==FlightPath::AnitaLite){
    igps_previous=NPOINTS_MIN; // initialise here to avoid times during launch
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
    fChain->SetBranchAddress("realTime_surfhk",&realTime_flightdata);
    fChain->SetBranchAddress("heading",&fheading);
  }

  else if (WHICHPATH==FlightPath::Anita2 || WHICHPATH==FlightPath::Anita3 || WHICHPATH==FlightPath::Anita4) { // for anita-3 and 4 flights

    std::string balloonFile = ICEMC_DATA_DIR;
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
    fChain->SetMakeClass(1);
    fChain->Add(balloonFile.c_str());//created to include pitch and roll.
    fChain->SetBranchAddress("longitude",&flongitude);
    fChain->SetBranchAddress("latitude",&flatitude);
    fChain->SetBranchAddress("altitude",&faltitude);
    fChain->SetBranchAddress("heading",&fheading);
    fChain->SetBranchAddress("realTime",&realTime_flightdata_temp);
    fChain->SetBranchAddress("pitch",&fpitch);
    fChain->SetBranchAddress("roll",&froll);
  }
  
  NPOINTS=0;
  REDUCEBALLOONPOSITIONS=100;

}

// double icemc::Balloon::GetBalloonSpin(double heading) const { // get the azimuth of the balloon
      
//   double phi_spin;
//   if (WHICHPATH==FlightPath::AnitaLite ||
//       WHICHPATH==FlightPath::Anita1    ||
//       WHICHPATH==FlightPath::Anita2    ||
//       WHICHPATH==FlightPath::Anita3    ||
//       WHICHPATH==FlightPath::Anita4){

//     phi_spin=heading*constants::RADDEG;
//   }
//   else {
//     if (RANDOMIZE_BN_ORIENTATION==1){
//       phi_spin=gRandom->Rndm()*2*constants::PI;
//     }
//     else{
//       phi_spin=0.;
//     }
//   }
  
//   return phi_spin;
// }


int icemc::Balloon::Getibnposition() {
  int ibnposition_tmp;
  if (WHICHPATH==FlightPath::Circle80DegreesSouth){
    double lon = r_bn.Longitude();
    if(lon < 0 ){lon += 360;}
    ibnposition_tmp = (int)(lon/2);
  }
  else if (WHICHPATH==FlightPath::AnitaLite ||
	   WHICHPATH==FlightPath::Anita1 ||
	   WHICHPATH==FlightPath::Anita2 ||
	   WHICHPATH==FlightPath::Anita3 ||
	   WHICHPATH==FlightPath::Anita4){
    ibnposition_tmp=(int)((double)igps/(double)REDUCEBALLOONPOSITIONS);
    // std::cout << __FUNCTION__ << " " << igps << " " << REDUCEBALLOONPOSITIONS << " " << ibnposition_tmp << std::endl;
  }
  else{
    ibnposition_tmp=0;
  }
    
  return ibnposition_tmp;
    
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
  icemcLog() << icemc::warning << __PRETTY_FUNCTION__
	     << " could not get TUFF index from current time "
	     << Curr_time << ", returning -1." << std::endl;
  return -1;
}


// this is called for each neutrino
void icemc::Balloon::PickBalloonPosition(const Antarctica *antarctica1, const Settings *settings1, int inu, Anita *anita1, double randomNumber) {

  pitch=0.;
  roll=0.;

  if (WHICHPATH==FlightPath::AnitaLite ||
      WHICHPATH==FlightPath::Anita1 ||
      WHICHPATH==FlightPath::Anita2 ||
      WHICHPATH==FlightPath::Anita3 ||
      WHICHPATH==FlightPath::Anita4) {
        
    if (WHICHPATH==FlightPath::AnitaLite) {
      igps = NPOINTS_MIN+(igps_previous+1-NPOINTS_MIN)%(NPOINTS_MAX-NPOINTS_MIN); //Note: ignore last 140 points, where balloon is falling - Stephen
      flatitude=(float)latitude_bn_anitalite.at(igps);
      flongitude=(float)longitude_bn_anitalite.at(igps);
      faltitude=(float)altitude_bn_anitalite.at(igps);
      fheading=(float)heading_bn_anitalite.at(igps);
    }
    else if (WHICHPATH==FlightPath::Anita1 ||
	     WHICHPATH==FlightPath::Anita2 ||
	     WHICHPATH==FlightPath::Anita3 ||
	     WHICHPATH==FlightPath::Anita4) {
      // For Anita 1 and Anita 2 and Anita 3:
      // igps = (igps_previous+1)%fChain->GetEntries(); // pick which event in the tree we want
      static int start_igps = 0; 
      static int ngps = fChain->GetEntries(); 
      static int init_best = 0;
      
      if (settings1->PAYLOAD_USE_SPECIFIC_TIME && !init_best) 
      {
         int N = fChain->Draw("realTime","","goff"); 
         double * times = fChain->GetV1(); 

         int best_igps =  TMath::BinarySearch(N, times, (double) settings1->PAYLOAD_USE_SPECIFIC_TIME); 
         start_igps = best_igps;
         int end_igps = best_igps;

         while (times[start_igps] > settings1->PAYLOAD_USE_SPECIFIC_TIME - settings1->PAYLOAD_USE_SPECIFIC_TIME_DELTA)
         {
           start_igps--; 
         }

         while (times[end_igps] < settings1->PAYLOAD_USE_SPECIFIC_TIME + settings1->PAYLOAD_USE_SPECIFIC_TIME_DELTA)
         {
           end_igps++; 
         }

         ngps = end_igps - start_igps + 1; 
         init_best = 1; 
      }

      igps = start_igps + int(randomNumber*ngps); // use random position 

      //////////////////////////// TEMPORARY HACKS FOR ANITA4 !!!!!!      
      if (WHICHPATH==FlightPath::Anita4 && ((igps>870 && igps<880) || (igps>7730 && igps<7740) || (igps>23810 && igps<23820) || (igps>31630 && igps<31660)) || (igps==17862) ){
	igps = igps+30;
      }
      
      fChain->GetEvent(igps); // this grabs the balloon position data for this event
      realTime_flightdata = realTime_flightdata_temp;
      if(settings1->TUFFSON){
       anita1->tuffIndex = getTuffIndex(realTime_flightdata);
      }// end if tuffson 
      
      while (faltitude<MINALTITUDE || fheading<0) { // if the altitude is too low, pick another event.
		    
	igps++; // increment by 1
	igps=igps%fChain->GetEntries(); // make sure it's not beyond the maximum entry number
		    
	fChain->GetEvent(igps);	  // get new event
      }
      // set phi Masking for Anita 2 or Anita 3
      // the deadtime is read from the same tree
      if ((WHICHPATH==FlightPath::Anita2 ||
	   WHICHPATH==FlightPath::Anita3 ||
	   WHICHPATH==FlightPath::Anita4) &&
	  (settings1->PHIMASKING==1 || settings1->USEDEADTIME)){
	anita1->setphiTrigMask(realTime_flightdata);
      }
      if ((WHICHPATH==FlightPath::Anita3 || WHICHPATH==FlightPath::Anita4) && settings1->USETIMEDEPENDENTTHRESHOLDS==1){ // set time-dependent thresholds
	anita1->setTimeDependentThresholds(realTime_flightdata);
      }
    }
    igps_previous=igps;
    latitude  = (double) flatitude;
    longitude = (double) flongitude;
    altitude  = (double) faltitude;
    heading   = (double) fheading;
    roll      = (double) froll;
    pitch     = (double) fpitch;

    setr_bn(latitude,longitude); // sets theta_bn, phi_bn and r_bn.  r_bn is a unit vector pointing in the right direction

    if (WHICHPATH==FlightPath::AnitaLite){
      altitude_bn=altitude*12.*constants::CMINCH/100.;
    }
    else if (WHICHPATH==FlightPath::Anita1 ||
	     WHICHPATH==FlightPath::Anita2 ||
	     WHICHPATH==FlightPath::Anita3 ||
	     WHICHPATH==FlightPath::Anita4){
      altitude_bn=altitude; // get the altitude of the balloon in the right units
    }
    surface_under_balloon = antarctica1->Surface(r_bn); // get altitude of the surface under the balloon

    r_bn_shadow = surface_under_balloon * r_bn.Unit(); // this is a vector pointing to spot just under the balloon on the surface (its shadow at high noon)
    r_bn = (antarctica1->Geoid(r_bn)+altitude_bn) * r_bn.Unit();
    //r_bn = (antarctica->Surface(r_bn)+altitude_bn) * r_bn.Unit(); //this points to balloon position (not a unit vector)
  } //if (ANITA-lite path) or anita 1 or anita 2
    
    
    
  else if (WHICHPATH==FlightPath::Circle80DegreesSouth){ // pick random phi at 80 deg S
    phi_bn=gRandom->Rndm()*constants::TWOPI;
		
    // r_bn = Position(theta_bn,phi_bn);
    r_bn = GeoidModel::Position(0, 0, 1);
    r_bn.SetTheta(theta_bn);
    r_bn.SetPhi(phi_bn);    
    surface_under_balloon = antarctica1->Surface(r_bn);
		
    r_bn_shadow = surface_under_balloon * r_bn.Unit();
    //    r_bn = (antarctica->Geoid(r_bn)+altitude_bn) * r_bn.Unit();
    r_bn = (antarctica1->Surface(r_bn)+altitude_bn) * r_bn.Unit();
    igps=0;
  } //else if(random position at 80 deg S)
  else if (WHICHPATH==FlightPath::FixedPosition){
		
    igps=0; // set gps position to be zero - we're at a fixed balloon position
		
  }
  else if (WHICHPATH==FlightPath::Custom) {
    igps=0;
    theta_bn=1.*constants::RADDEG; // 1deg
    phi_bn=1.*constants::RADDEG; // 1deg
    // r_bn=Position(theta_bn,phi_bn); // sets r_bn
    r_bn = GeoidModel::Position(0, 0, 1);
    r_bn.SetTheta(theta_bn);
    r_bn.SetPhi(phi_bn);    
    if (BN_ALTITUDE!=0){
      altitude_bn=BN_ALTITUDE*12.*constants::CMINCH/100.; // set the altitude of the balloon to be what you pick.  This isn't in time for CreateHorizons though!
    }
    surface_under_balloon = antarctica1->Surface(r_bn);
    r_bn_shadow = surface_under_balloon * r_bn.Unit();
    //    r_bn = (antarctica->Geoid(r_bn)+altitude_bn) * r_bn.Unit();
		
    r_bn = (antarctica1->Surface(r_bn)+altitude_bn) * r_bn.Unit();
  } // you pick it
    
  ibnposition = Getibnposition();

  if (!settings1->UNBIASED_SELECTION && dtryingposition!=-999){
    dtryingposition=antarctica1->GetBalloonPositionWeight(ibnposition);
  }
  else{
    dtryingposition=1.;
  }
    
  // // normalized balloon position
  TVector3 n_bn = r_bn.Unit();
    
  if (settings1->SLAC){
    icemcLog() << icemc::error << "SLAC settings are currently disabled!" << std::endl;
    // AdjustSlacBalloonPosition(inu); // move payload around like we did at slac
  }
   
} // end PickBalloonPosition




void icemc::Balloon::setr_bn(double latitude,double longitude) {

  r_bn.SetLonLatAlt(longitude, latitude, 0);
  r_bn = r_bn.Unit();
    
  // // latitude is between -90 and 0.
  // // theta_bn measured from the SP and is between 0 and constants::PI/2.
  // theta_bn = (90+latitude)*constants::RADDEG;

  // // this is the payload's longitude, not the azimuth of the balloon like it sounds.
  // // longitude is between -180 to 180 with 0 at prime meridian
  // // phi is from 0 to 360 with 0 at +90 longitude
  // phi_bn = (-1*longitude+90.);
    
  // if (phi_bn<0){
  //   phi_bn += 360.;
  // }    
  // phi_bn *= constants::RADDEG;

  // r_bn = Position(theta_bn,phi_bn);  //r_bn is a unit vector pointing in the right direction
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
