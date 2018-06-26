#include "TChain.h"
#include "TFile.h"
#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "earthmodel.hh"
#include "icemodel.hh"
#include "vector.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMath.h"

#include "position.hh"
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
  case icemc::FlightPath::PeterEvent:
    return os << "FlightPath::PeterEvent";
  case icemc::FlightPath::SLAC:
    return os << "FlightPath::SLAC";
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


icemc::Balloon::Balloon(const Settings* settings)
  : WHICHPATH(settings ? static_cast<icemc::FlightPath>(settings->WHICHPATH) : FlightPath::FixedPosition){

  
  MAXHORIZON=800000.; // pick the interaction within this distance from the balloon so that it is within the horizon
  ibnposition=0;
  igps=0;
  horizcoord_bn=0;    // x component of balloon position
  vertcoord_bn=0;     // y component of balloon position
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

  InitializeBalloon();
}



void  icemc::Balloon::setObservationLocation(Interaction *interaction1,int inu, const IceModel *antarctica, const Settings *settings1) {
  interaction1->banana_volts = 0; //Zero the variable
  interaction1->banana_obs = Vector(0,0,Interaction::banana_observation_distance);
    
  // First pick theta of the obs vector from nu vector

  interaction1->banana_theta_obs = -0.5*constants::PI/settings1->vertical_banana_points * ((int)(inu/settings1->horizontal_banana_points));
  interaction1->banana_phi_obs = 2*constants::PI/settings1->horizontal_banana_points * (inu%settings1->horizontal_banana_points);
  interaction1->banana_weight = (double)fabs(sin(interaction1->banana_theta_obs)); //Set "weight" - purely a phase factor
  interaction1->banana_obs = interaction1->banana_obs.RotateY(interaction1->banana_theta_obs);
    
  // Now the phi of the obs vector relative to the nu vector
  interaction1->banana_obs = interaction1->banana_obs.RotateZ(interaction1->banana_phi_obs);
    
  // Finally rotate the obs vector to the nu vector's orientation (change coordinate systems, in effect)
  interaction1->banana_obs = interaction1->banana_obs.RotateY(constants::PI-Interaction::theta_nu_banana);
  interaction1->banana_obs = interaction1->banana_obs.RotateZ(Interaction::phi_nu_banana);
    
  //Finally, set the balloon position to the location which we would like to observe
  r_bn = interaction1->nu_banana_surface + interaction1->banana_obs;
    
  theta_bn = r_bn.Theta();
  phi_bn = r_bn.Phi();
  surface_under_balloon = antarctica->Surface(r_bn);
  // position of balloon at earth's surface (the shadow it would cast with the sun overhead)
  r_bn_shadow = surface_under_balloon * r_bn.Unit();
  altitude_bn = r_bn.Mag() - surface_under_balloon;
    
  //Finished setting observation location
}


void icemc::Balloon::SetDefaultBalloonPosition(const IceModel *antarctica1) { // position of surface of earth under balloon
    
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
    phi_bn=constants::PI/4; //wrt 90E longitude
  }
  else {
    phi_bn=EarthModel::LongtoPhi_0isPrimeMeridian(BN_LONGITUDE_SETTING); //remember input of LongtoPhi is between -180 and 180
  }    
    
  r_bn = Position(theta_bn,phi_bn); // direction of balloon- right now this is a unit vector
    
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

  const string ICEMC_SRC_DIR = icemc::EnvironmentVariable::ICEMC_SRC_DIR();
  const string ICEMC_DATA_DIR = ICEMC_SRC_DIR+"/data/";  
  const string anitaliteflight = ICEMC_DATA_DIR+"/BalloonGPS.txt"; // the gps path of the anita-lite flight
  std::ifstream flightfile(anitaliteflight.c_str()); //set file to the right one
  if (!flightfile) {
    std::cout << "Flight file not found.\n";
    exit(1);
  }
    
  string slatitude; // latitude in string format
  string slongitude; //longitude in string format
  string saltitude; // altitude in string format
  string sheading; // heading in string format
  string srealtime; // real life unix time in string format
  string junk; // string not used for anything
    
  string line; // for reading a whole line
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

  const int oldArrayLength = 100000;
  latitude_bn_anitalite.reserve(oldArrayLength);
  longitude_bn_anitalite.reserve(oldArrayLength);
  altitude_bn_anitalite.reserve(oldArrayLength);
  heading_bn_anitalite.reserve(oldArrayLength);
  realtime_bn_anitalite.reserve(oldArrayLength);
  
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
    
  latitude_bn_anitalite.shrink_to_fit();
  longitude_bn_anitalite.shrink_to_fit();
  altitude_bn_anitalite.shrink_to_fit();
  heading_bn_anitalite.shrink_to_fit();
  realtime_bn_anitalite.shrink_to_fit();
    
  NPOINTS_MAX=NPOINTS-140; // exclude the fall
  NPOINTS_MIN=0;

}//ReadAnitaliteFlight



void icemc::Balloon::InitializeBalloon() {
  const string ICEMC_SRC_DIR = icemc::EnvironmentVariable::ICEMC_SRC_DIR();
  const string ICEMC_DATA_DIR = ICEMC_SRC_DIR+"/data/";
  const string anitaflight = ICEMC_DATA_DIR+"/anitagps.txt";// gps path of anita flight    
    
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
  flightdatachain = nullptr; 

  // This for Anita 1 flight
  if (WHICHPATH==FlightPath::Anita1) {

    flightdatachain = new TChain("foricemc");
    flightdatachain->Add((ICEMC_DATA_DIR+"/anita1flightdata.root").c_str());
		
    flightdatachain->SetBranchAddress("surfTrigBandMask",surfTrigBandMask);
    flightdatachain->SetBranchAddress("powerthresh",powerthresh);
    flightdatachain->SetBranchAddress("meanp",meanp);
    flightdatachain->SetBranchAddress("longitude",&flongitude);
    flightdatachain->SetBranchAddress("latitude",&flatitude);
    flightdatachain->SetBranchAddress("altitude",&faltitude);
    flightdatachain->SetBranchAddress("realTime_surfhk",&realTime_flightdata);
    flightdatachain->SetBranchAddress("heading",&fheading);
  }
    
  // This for Anita 2 flight
  if (WHICHPATH==FlightPath::Anita2) {
		
    flightdatachain = new TChain("adu5PatTree");
    flightdatachain->SetMakeClass(1);
    flightdatachain->Add((ICEMC_DATA_DIR+"/anita2gps_pitchandroll.root").c_str());//created to include pitch and roll.
    flightdatachain->SetBranchAddress("longitude",&flongitude);
    flightdatachain->SetBranchAddress("latitude",&flatitude);
    flightdatachain->SetBranchAddress("altitude",&faltitude);
    flightdatachain->SetBranchAddress("heading",&fheading);
    flightdatachain->SetBranchAddress("realTime",&realTime_flightdata_temp);
    flightdatachain->SetBranchAddress("pitch",&fpitch);
    flightdatachain->SetBranchAddress("roll",&froll);
    //std::cout << "Loading file.  n events is " << flightdatachain->GetEntries() << "\n";
		
  }
  else if (WHICHPATH==FlightPath::Anita3 || WHICHPATH==FlightPath::Anita4) { // for anita-3 and 4 flights

    std::string balloonFile = ICEMC_DATA_DIR;
    switch(WHICHPATH){
    case FlightPath::Anita3:
      balloonFile+="/anita3gps_pitchroll.root";
      break;
    case FlightPath::Anita4:
    default:
      balloonFile+="/anita4gps_pitchroll.root";
      break;
    }
    
    flightdatachain = new TChain("adu5PatTree");
    flightdatachain->SetMakeClass(1);
    flightdatachain->Add(balloonFile.c_str());//created to include pitch and roll.
    flightdatachain->SetBranchAddress("longitude",&flongitude);
    flightdatachain->SetBranchAddress("latitude",&flatitude);
    flightdatachain->SetBranchAddress("altitude",&faltitude);
    flightdatachain->SetBranchAddress("heading",&fheading);
    flightdatachain->SetBranchAddress("realTime",&realTime_flightdata_temp);
    flightdatachain->SetBranchAddress("pitch",&fpitch);
    flightdatachain->SetBranchAddress("roll",&froll);
  }

  // for (int i=0;i<10000;i++) {
  //   latitude_bn_anitalite[i]=0;
  //   longitude_bn_anitalite[i]=0;
  //   altitude_bn_anitalite[i]=0;
  // }
  NPOINTS=0;
  REDUCEBALLOONPOSITIONS=100;

  if (WHICHPATH!=FlightPath::Anita1) { // @todo is this if statement correct?
    // if it's not the anita path, you won't be able to read in thresholds and masks
    for (int i=0;i<9;i++) {
      for (int j=0;j<2;j++) {
	surfTrigBandMask[i][j]=0;
      }
      for (int j=0;j<32;j++) {
	powerthresh[i][j]=3.E-12;
	meanp[i][j]=1.E-12;
      }			
    }		
  }
}

double icemc::Balloon::GetBalloonSpin(double heading) { // get the azimuth of the balloon
      
  double phi_spin;
  if (WHICHPATH==FlightPath::AnitaLite     || WHICHPATH==FlightPath::Anita1        || WHICHPATH==FlightPath::Anita2        || WHICHPATH==FlightPath::Anita3        || WHICHPATH==FlightPath::Anita4       )
    phi_spin=heading*constants::RADDEG;
  else {
    if (RANDOMIZE_BN_ORIENTATION==1){
      phi_spin=gRandom->Rndm()*2*constants::PI;
    }
    else{
      phi_spin=0.;
    }
  }
  
  return phi_spin;
}


int icemc::Balloon::Getibnposition() {
  int ibnposition_tmp;
  if (WHICHPATH==FlightPath::Circle80DegreesSouth){
    ibnposition_tmp = (int)(r_bn.Lon() / 2);
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

void icemc::Balloon::PickBalloonPosition(Vector straightup,const IceModel *antarctica1,const Settings *settings1,Anita *anita1) {
  // takes a 3d vector pointing along the z axis
  Vector thetazero(0.,0.,1.);
  Vector phizero(1.,0.,0.);
  igps=0;
  theta_bn=acos(straightup.Dot(thetazero)); // 1deg
 
  if (straightup.GetX()==0 && straightup.GetY()==0){
    phi_bn=0.;
  }
  else{
    phi_bn=acos((straightup.Cross(thetazero)).Dot(phizero)); // 1deg
  }

  r_bn=Position(theta_bn,phi_bn); // sets r_bn
  if (BN_ALTITUDE!=0){
    altitude_bn=BN_ALTITUDE*12.*constants::CMINCH/100.; // set the altitude of the balloon to be what you pick.  This isn't in time for CreateHorizons though!
  }
  surface_under_balloon = antarctica1->Surface(r_bn);
  r_bn_shadow = surface_under_balloon * r_bn.Unit();
  r_bn = (antarctica1->Surface(r_bn)+altitude_bn) * r_bn.Unit();


  ibnposition=Getibnposition();
  
  dtryingposition=1.;
  
  phi_spin=0.; // get the azimuth of the balloon.
    
  // normalized balloon position
  n_bn = r_bn.Unit();
  
  // finding which direction is east under the balloon
  n_east = Vector(sin(phi_bn), -1*cos(phi_bn), 0.);
  

  // find position of each antenna boresight
  
  // now finding north
  n_north = n_bn.Cross(n_east);
  
  // these coordinates are for filling ntuples.
  horizcoord_bn=r_bn[0]/1000.; //m to km
  vertcoord_bn=r_bn[1]/1000.;
  
  
  calculate_antenna_positions(settings1,anita1);
  //  std::cout << "calculated antenna positions.\n";
  if (settings1->BORESIGHTS){
    GetBoresights(settings1,anita1);
  }    
}


// for tuffs for anita-4
int getTuffIndex(int Curr_time) {
  if((icemc::constants::TUFFconfig_B_end_3 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_A_end_1)) {// config A trigconfigA.imp
    return 0;
  }
  else if(((0 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_B_end_1)) || ((icemc::constants::TUFFconfig_P_end_3 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_B_end_2)) || ((icemc::constants::TUFFconfig_P_end_4 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_B_end_3)) || ((icemc::constants::TUFFconfig_A_end_1 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_B_end_4)) || ((icemc::constants::TUFFconfig_P_end_5 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_B_end_5)) || ((icemc::constants::TUFFconfig_P_end_6 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_B_end_6)) || (icemc::constants::TUFFconfig_P_end_7 < Curr_time) ) { // config B trigconfigB.imp
    return 1;
  }
  else if((icemc::constants::TUFFconfig_P_end_1 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_C_end_1)) { // config C trigconfigC.imp
    return 2;
  }
  else if( ((icemc::constants::TUFFconfig_P_end_2 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_G_end_1)) || ((icemc::constants::TUFFconfig_O_end_1 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_G_end_2)) ) { // config G trigconfigG.imp
    return 3;
  }
  else if( ((icemc::constants::TUFFconfig_G_end_1 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_O_end_1)) || ((icemc::constants::TUFFconfig_G_end_2 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_O_end_2)) ) { // config O trigconfigO.imp
    return 4;
  }
  else if( ((icemc::constants::TUFFconfig_B_end_1 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_P_end_1)) || ((icemc::constants::TUFFconfig_C_end_1 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_P_end_2)) || ((icemc::constants::TUFFconfig_O_end_2 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_P_end_3)) || ((icemc::constants::TUFFconfig_B_end_2 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_P_end_4)) || ((icemc::constants::TUFFconfig_B_end_4 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_P_end_5)) || ((icemc::constants::TUFFconfig_B_end_5 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_P_end_6)) || ((icemc::constants::TUFFconfig_B_end_6 < Curr_time) && (Curr_time <= icemc::constants::TUFFconfig_P_end_7)) ) { // config P trigconfigP.imp
    return 5;
  }
  std::cerr << "Error in" << __PRETTY_FUNCTION__ << ", could not get TUFF index from current time... returning -1" << std::endl;
  return -1;
}


// this is called for each neutrino
void icemc::Balloon::PickBalloonPosition(const IceModel *antarctica1, const Settings *settings1, int inu, Anita *anita1, double randomNumber) {

  // r_bn_shadow=position of spot under the balloon on earth's surface

  //std::cout << "calling pickballoonposition.\n";
  pitch=0.;
  roll=0.;
  phi_spin=0.;

  //flightdatatree->SetBranchAddress("surfTrigBandMask",&surfTrigBandMask);
  //unsigned short test[9][2];
    
  //double latitude,longitude,altitude;
    
  //  std::cout << "I'm in pickballoonposition. whichpath is " << WHICHPATH << "\n";
  //Pick balloon position
  if (WHICHPATH==FlightPath::AnitaLite ||
      WHICHPATH==FlightPath::Anita1 ||
      WHICHPATH==FlightPath::Anita2 ||
      WHICHPATH==FlightPath::Anita3 ||
      WHICHPATH==FlightPath::Anita4) { // anita-lite or anita-I or -II path
        
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
      // igps = (igps_previous+1)%flightdatachain->GetEntries(); // pick which event in the tree we want
      static int start_igps = 0; 
      static int ngps = flightdatachain->GetEntries(); 
      static int init_best = 0;
      
      if (settings1->PAYLOAD_USE_SPECIFIC_TIME && !init_best) 
      {
         int N = flightdatachain->Draw("realTime","","goff"); 
         double * times = flightdatachain->GetV1(); 

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

      
      flightdatachain->GetEvent(igps); // this grabs the balloon position data for this event
      realTime_flightdata = realTime_flightdata_temp;
      if(settings1->TUFFSON){
       anita1->tuffIndex = getTuffIndex(realTime_flightdata);
      }// end if tuffson 
      
      while (faltitude<MINALTITUDE || fheading<0) { // if the altitude is too low, pick another event.
		    
	igps++; // increment by 1
	igps=igps%flightdatachain->GetEntries(); // make sure it's not beyond the maximum entry number
		    
	flightdatachain->GetEvent(igps);	  // get new event
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
		
    r_bn = Position(theta_bn,phi_bn);
    surface_under_balloon = antarctica1->Surface(r_bn);
		
    r_bn_shadow = surface_under_balloon * r_bn.Unit();
    //    r_bn = (antarctica->Geoid(r_bn)+altitude_bn) * r_bn.Unit();
    r_bn = (antarctica1->Surface(r_bn)+altitude_bn) * r_bn.Unit();
    igps=0;
  } //else if(random position at 80 deg S)
  else if (WHICHPATH==FlightPath::FixedPosition){
		
    igps=0; // set gps position to be zero - we're at a fixed balloon position
		
  }
  else if (WHICHPATH==FlightPath::SLAC) { // for slac?
    igps=0;
    theta_bn=1.*constants::RADDEG; // 1deg
    phi_bn=1.*constants::RADDEG; // 1deg
    r_bn=Position(theta_bn,phi_bn); // sets r_bn
    if (BN_ALTITUDE!=0){
      altitude_bn=BN_ALTITUDE*12.*constants::CMINCH/100.; // set the altitude of the balloon to be what you pick.  This isn't in time for CreateHorizons though!
    }
    surface_under_balloon = antarctica1->Surface(r_bn);
    r_bn_shadow = surface_under_balloon * r_bn.Unit();
    //    r_bn = (antarctica->Geoid(r_bn)+altitude_bn) * r_bn.Unit();
		
    r_bn = (antarctica1->Surface(r_bn)+altitude_bn) * r_bn.Unit();
  } // you pick it
  else if (WHICHPATH==FlightPath::PeterEvent){ // for making comparison with Peter's event
    //phi_bn=0.767945;
    // phi_bn=0.95993;
    //166.73 longitude
    //    phi_bn=-1.339;// this is McM phi
    // 156.73
    phi_bn=-1.1646;// this is McM phi
    //phi_bn=-0.523;
    // theta_bn has been set to 10 degrees in setdefaultballoonposition
    // we want theta_bn to be on a path along same longitude as
    // McM but up to 10 deg further south in latitude
    // 12.12 deg. is McM
    double theta_max=0.2115; // this is McM theta
    double theta_min=0.2115-0.110; // this approximately one horizon north of McM
    //This is 701.59 km along the surface of the earth
    int npositions=1000; // number of possible positions between those two bounds
    theta_bn=theta_min+(theta_max-theta_min)*(double)(inu%npositions)/(double)npositions;
        
    r_bn = Position(theta_bn,phi_bn);
    surface_under_balloon = antarctica1->Surface(r_bn);
		
    r_bn_shadow = surface_under_balloon * r_bn.Unit();
    //    r_bn = (antarctica->Geoid(r_bn)+altitude_bn) * r_bn.Unit();
		
    r_bn = (antarctica1->Surface(r_bn)+altitude_bn) * r_bn.Unit();
    igps=0;
  } // comparing with Peter's event
    
  ibnposition=Getibnposition();
    
  if (!settings1->UNBIASED_SELECTION && dtryingposition!=-999){
    dtryingposition=antarctica1->GetBalloonPositionWeight(ibnposition);
  }
  else{
    dtryingposition=1.;
  }
  phi_spin=GetBalloonSpin(heading); // get the azimuth of the balloon.
    
  // normalized balloon position
  n_bn = r_bn.Unit();
    
  // finding which direction is east under the balloon
  n_east = Vector(sin(phi_bn), -1*cos(phi_bn), 0.);
    
  if (settings1->SLAC){
    icemcLog() << icemc::error << "SLAC settings are currently disabled!" << std::endl;
    // AdjustSlacBalloonPosition(inu); // move payload around like we did at slac
  }
  // find position of each antenna boresight
   
  // now finding north
  n_north = n_bn.Cross(n_east);
    
  // these coordinates are for filling ntuples.
  horizcoord_bn=r_bn[0]/1000.; //m to km
  vertcoord_bn=r_bn[1]/1000.;

  calculate_antenna_positions(settings1,anita1);
  if (settings1->BORESIGHTS){
    GetBoresights(settings1,anita1);
  }



} // end PickBalloonPosition




// void icemc::Balloon::AdjustSlacBalloonPosition(int inu) { // move payload around like we did at slac
    
//   if (inu<=MAX_POSITIONS)  // use the event number for the position if we c an
//     islacposition=inu;
//   else
//     islacposition=0; // otherwise just use the first positions
//   // adjust the payload position
//   r_bn+=slacpositions[islacposition].GetX()*n_east
//     +slacpositions[islacposition].GetY()*n_bn;
    
// }


void icemc::Balloon::CenterPayload(double& hitangle_e) {
    
  // allow an option to rotate the payload so the signal is
  // always on the boresight of one phi sector
    
  phi_spin-=hitangle_e;

}


void icemc::Balloon::GetAntennaOrientation(const Settings *settings1, Anita *anita1, int ilayer, int ifold, Vector& n_eplane, Vector& n_hplane, Vector& n_normal) const {

  // rotate for antenna's orientation on payload
  // face of antenna starts out relative to +x because phi is relative to +x
  // n_normal points northwards for antenna 0 (0-31)
  // const vectors const_z (n_eplane), const_y (-n_hplane), const_x (n_normal) defined under Constants.h   -- oindree

  if(settings1->WHICH==6 || settings1->WHICH==8 || settings1->WHICH==9) {
    n_eplane = constants::const_z.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
    n_hplane = (-constants::const_y).RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
    n_normal = constants::const_x.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
  }
  else {
    n_eplane = constants::const_z.RotateY(anita1->THETA_ZENITH[ilayer] - constants::PI/2);
    n_hplane = (-constants::const_y).RotateY(anita1->THETA_ZENITH[ilayer] - constants::PI/2);
    n_normal = constants::const_x.RotateY(anita1->THETA_ZENITH[ilayer] - constants::PI/2);
  }

  double phi;
  // rotate about z axis for phi
  if (settings1->CYLINDRICALSYMMETRY==1) {
    phi=(double)ifold/(double)anita1->NRX_PHI[ilayer]*2*constants::PI + anita1->PHI_OFFSET[ilayer] + phi_spin;
  }
  else{
    //phi=anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer] +phi_spin;
    phi=anita1->PHI_EACHLAYER[ilayer][ifold];
  }

  n_eplane = n_eplane.RotateZ(phi);
  n_hplane = n_hplane.RotateZ(phi);
  n_normal = n_normal.RotateZ(phi);

  n_eplane = RotatePayload(n_eplane); 
  n_hplane = RotatePayload(n_hplane);
  n_normal = RotatePayload(n_normal);

} //end void GetAntennaOrientation




void icemc::Balloon::setr_bn(double latitude,double longitude) {
    
  // latitude is between -90 and 0.
  // theta_bn measured from the SP and is between 0 and constants::PI/2.
  theta_bn = (90+latitude)*constants::RADDEG;

  // this is the payload's longitude, not the azimuth of the balloon like it sounds.
  // longitude is between -180 to 180 with 0 at prime meridian
  // phi is from 0 to 360 with 0 at +90 longitude
  phi_bn = (-1*longitude+90.);
    
  if (phi_bn<0){
    phi_bn += 360.;
  }    
  phi_bn *= constants::RADDEG;

  r_bn = Position(theta_bn,phi_bn);  //r_bn is a unit vector pointing in the right direction
}


void icemc::Balloon::PickDownwardInteractionPoint(Interaction *interaction1, Anita *anita1, const Settings *settings1, const IceModel *antarctica1, RayTracer *ray1, int &beyondhorizon) {
    
  // double distance=1.E7;
  double phi=0,theta=0;
  double lon=0;
  double latfromSP=0;
  
  
  if (settings1->UNBIASED_SELECTION==1) {

    if (antarctica1->PickUnbiased(interaction1,antarctica1)) { // pick neutrino direction and interaction point
      interaction1->dtryingdirection=1.;
      interaction1->iceinteraction=1;
    }
    else{
      interaction1->iceinteraction=0;
    }
  }
  else {
    interaction1->iceinteraction=1;
    if (WHICHPATH==FlightPath::Custom) { //Force interaction point if we want to make a banana plot
      interaction1->posnu = interaction1->nu_banana;
    } //if (making banana plot)
    else if (WHICHPATH==FlightPath::PeterEvent) {// Force interaction point for comparison with Peter.
      
      // want interaction location at Taylor Dome
      // According to lab book it's 77 deg, 52.818' S=77.8803
      // 158 deg, 27.555' east=158.45925
      latfromSP=12.1197; // this is degrees latitude from the south pole
      //lon=180.+166.73;
      lon=180.+158.45925;
      //lon=180.+120.
      phi=EarthModel::LongtoPhi_0is180thMeridian(lon);
      
      theta = latfromSP*constants::RADDEG;
      
      double elevation=antarctica1->SurfaceAboveGeoid(lon,latfromSP)-500.;
      
      interaction1->posnu = Vector((elevation+antarctica1->Geoid(latfromSP))*sin(theta)*cos(phi),(elevation+antarctica1->Geoid(latfromSP))*sin(theta)*sin(phi),(elevation+antarctica1->Geoid(latfromSP))*cos(theta));
			
      
      
    } //if (WHICHPATH==FlightPath::PeterEvent)
    
    else if (settings1->SLAC) {
      
      Vector zaxis(0.,0.,1.); // start with vector pointing in the +z direction
      
      interaction1->posnu=zaxis.RotateY(r_bn.Theta()-settings1->SLAC_HORIZDIST/EarthModel::EarthRadiusMeters); // rotate to theta of balloon, less the distance from the interaction to the balloon
      
      interaction1->posnu=interaction1->posnu.RotateZ(r_bn.Phi()); // rotate to phi of the balloon
      
      interaction1->posnu=(antarctica1->Surface(interaction1->posnu)-settings1->SLAC_DEPTH)*interaction1->posnu; // put the interaction position at depth settings1->SLAC_DEPTH in the ice
      
    }
    else{

      // If we require neutrinos from a particular position
      // we generate that cartesian position here

      static Vector specific_position; 

      if (settings1->SPECIFIC_NU_POSITION) 
      {
        double R = settings1->SPECIFIC_NU_POSITION_ALTITUDE + antarctica1->Geoid(settings1->SPECIFIC_NU_POSITION_LATITUDE); 
        double theta = settings1->SPECIFIC_NU_POSITION_LATITUDE * constants::RADDEG; 
        double phi = EarthModel::LongtoPhi_0isPrimeMeridian(settings1->SPECIFIC_NU_POSITION_LONGITUDE); 
        specific_position.SetXYZ(R * sin(theta) * cos(phi), R * sin(theta) * sin(phi), R * cos(theta)); 
      }
        
      do
      {
        interaction1->posnu = antarctica1->PickInteractionLocation(ibnposition, settings1, r_bn, interaction1);
      } while(settings1->SPECIFIC_NU_POSITION &&  (interaction1->posnu - specific_position).Mag() > settings1->SPECIFIC_NU_POSITION_DISTANCE); 

    }
  }
  
  ray1->rfexit[0] = antarctica1->Surface(interaction1->posnu) * interaction1->posnu.Unit();  
  
  // unit vector pointing to antenna from exit point.
  ray1->n_exit2bn[0] = (r_bn - ray1->rfexit[0]).Unit();
  
  if (settings1->BORESIGHTS) {
    for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
      for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
	ray1->rfexit_eachboresight[0][ilayer][ifold] = antarctica1->Surface(interaction1->posnu) * interaction1->posnu.Unit();// this first guess rfexit is the same for all antennas too
	ray1->n_exit2bn_eachboresight[0][ilayer][ifold] = (r_boresights[ilayer][ifold]- ray1->rfexit_eachboresight[0][ilayer][ifold]).Unit();
	// std::cout << "ilayer, ifold, n_exit2bn are " << ilayer << "\t" << ifold << " ";
      }
    }
  }
  
  // first pass at direction of vector from interaction to exit point
  // just make the ray go radially outward away from center of earth.
  // still for first guess
  // should incorporate this into PickInteractionPoint
  ray1->nrf_iceside[0] = interaction1->posnu.Unit();
  
  if (settings1->BORESIGHTS) {
    // this is the same for all of the antennas too
    for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) { // loop over layers on the payload
      for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
	ray1->nrf_iceside_eachboresight[0][ilayer][ifold]=interaction1->posnu.Unit();
      } // end loop over fold
    } // end loop over payload layers
  } // end if boresights
  
  
  
  
  double r_down = 2*(antarctica1->Surface(interaction1->posnu)-antarctica1->IceThickness(interaction1->posnu))-interaction1->posnu.Mag();
  interaction1->posnu_down = r_down * interaction1->posnu.Unit();
  //position of the mirror point of interaction
  
  //interaction1->posnu is downward interaction1->posnu.
  // distance=interaction1->posnu.Distance(r_bn);
  
  
  // depth of interaction
  // gets distance between interaction and exit point, this time it's same as depth
  // because our first guess at exit point is radially outward from interaction.
  // negative means below surface
  interaction1->altitude_int=-1*ray1->rfexit[0].Distance(interaction1->posnu);
  interaction1->altitude_int_mirror=-1*ray1->rfexit[0].Distance(interaction1->posnu_down);//get depth of mirror point
  
  
  interaction1->r_fromballoon[0]=r_bn.Distance(interaction1->posnu);
  
  //distance from the mirror point to the balloon because it is equal to the path that signals pass
  interaction1->r_fromballoon[1]=r_bn.Distance(interaction1->posnu_down);
  
  
  // std::cout << (interaction1->posnu) << std::endl;
  if (ray1->n_exit2bn[0].Angle(interaction1->posnu) > constants::PI/2 &&
      !(WHICHPATH==FlightPath::Custom || WHICHPATH==FlightPath::PeterEvent)){
    beyondhorizon = 1;
  }
  
  return;
}//PickDownwardInteractionPoint



void icemc::Balloon::GetBoresights(const Settings *settings1, Anita *anita1, Position r_bn, double phi_spin, Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX]){//,Vector n_north,Vector n_east) {
    
  // this fills r_boresights with the position of the antenna boresights.
  for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
    for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
      
      Vector n_boresight(1,0,0); // this will be a unit vector that points from the center axis of the payload to the boresight on each level of the payload (it will stay in the x-y plane)
      
      double phi;
      // rotate about z axis for phi
      if (settings1->CYLINDRICALSYMMETRY==1) {
	phi=(double)ifold/(double)anita1->NRX_PHI[ilayer]*2*constants::PI + anita1->PHI_OFFSET[ilayer]+phi_spin;
      }
      else{
	phi=anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer] +phi_spin;
      }
      
      
      n_boresight = n_boresight.RotateZ(phi);
      
      // start with a vector that points in the +z direction but with the same length as r_bn
      r_boresights[ilayer][ifold].SetX(0.); // this will be a unit vector that points from the center axis of the payload to the boresight on each level of the payload (it will stay in the x-y plane)
      r_boresights[ilayer][ifold].SetY(0.);
      r_boresights[ilayer][ifold].SetZ(r_bn.Mag());
      
      Vector zaxis(0.,0.,1.);
      Vector xaxis(1.,0.,0.);
      Vector yaxis(0.,1.,0.);
      
      // Add the positions of the antennas relative to the center of the payload
      r_boresights[ilayer][ifold] += (anita1->RRX[ilayer]*n_boresight
				      + anita1->LAYER_VPOSITION[ilayer]*zaxis
				      + anita1->LAYER_HPOSITION[ilayer]*cos(anita1->LAYER_PHIPOSITION[ilayer])*xaxis
				      + anita1->LAYER_HPOSITION[ilayer]*sin(anita1->LAYER_PHIPOSITION[ilayer])*yaxis);

      //   now rotate to balloon's position on the continent
      r_boresights[ilayer][ifold] = r_boresights[ilayer][ifold].ChangeCoord(n_north,-1.*n_east);
    }
  }
}


void icemc::Balloon::GetBoresights(const Settings *settings1,Anita *anita1) {
  Vector ant_pos;
  for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
    for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
      ant_pos = anita1->ANTENNA_POSITION_START[0][ilayer][ifold];
      ant_pos = RotatePayload(ant_pos);
      r_boresights[ilayer][ifold] = ant_pos + r_bn;
    }
  }
}



void icemc::Balloon::calculate_antenna_positions(const Settings *settings1, Anita *anita1) const {
  int number_all_antennas = 0;
  Vector antenna_position;

  for (int ipol=0; ipol<2; ipol++){
    number_all_antennas=0;
    for (int ilayer = 0; ilayer < settings1->NLAYERS; ilayer++){
      for (int ifold = 0; ifold < anita1->NRX_PHI[ilayer]; ifold++){
	double phi = 0;
	if (settings1->WHICH==6 || settings1->WHICH==8 || settings1->WHICH == 9 || settings1->WHICH == 10){ //If payload is either
	  antenna_position = anita1->ANTENNA_POSITION_START[ipol][ilayer][ifold];
	}
	else {
	  if (settings1->CYLINDRICALSYMMETRY==1){ // for timing code
	    // phi is 0 for antenna 0 (0-31) and antenna 16 (0-31)
	    // antenna 1 (1-32) and antenna 18 (1-32)
	    phi = (double) ifold / (double) anita1->NRX_PHI[ilayer] * 2 * constants::PI + anita1->PHI_OFFSET[ilayer];
	  }
	  else{
	    phi = anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer];
	  }
	  antenna_position = Vector(anita1->RRX[ilayer]*cos(phi) + anita1->LAYER_HPOSITION[ilayer]*cos(anita1->LAYER_PHIPOSITION[ilayer]), anita1->RRX[ilayer]*sin(phi)+anita1->LAYER_HPOSITION[ilayer]*sin(anita1->LAYER_PHIPOSITION[ilayer]), anita1->LAYER_VPOSITION[ilayer]);

	}//else

	//std::cout<<"antenna_position start for "<<number_all_antennas<<" is  "<<antenna_position<<"\n";
	antenna_position = RotatePayload(antenna_position);

	anita1->antenna_positions[ipol][number_all_antennas] = antenna_position;
	//std::cout<<"antenna_position for "<<number_all_antennas<<" is  "<<antenna_position<<"\n";
	number_all_antennas++;
      }
    }
  }
  return;
}

icemc::Vector icemc::Balloon::RotatePayload(Vector ant_pos_pre) const {

  const double thisPitch = WHICHPATH == FlightPath::Anita2 ? fixedAnita2Pitch : WHICHPATH == FlightPath::Anita3 ? fixedAnita3Pitch : pitch;
  const double thisRoll  = WHICHPATH == FlightPath::Anita2 ? fixedAnita2Roll  : WHICHPATH == FlightPath::Anita3 ? fixedAnita3Roll  : roll;

  const Vector BalloonPos = r_bn;
  const double BalloonTheta = BalloonPos.Theta();
  double BalloonPhi = BalloonPos.Phi();
  
  if(BalloonPhi > constants::PI){
    BalloonPhi = BalloonPhi-constants::TWOPI;
  }
  
  Vector ant_pos = ant_pos_pre;

  const Vector zaxis(0.,0.,-1.);
  Vector xaxis(1.,0.,0.);  // roll axis
  Vector yaxis(0.,-1.,0.); // pitch axis for positive rotation to the clockwise of roll

  // do heading...
  ant_pos = ant_pos.Rotate(heading*constants::RADDEG,zaxis);
  xaxis = xaxis.Rotate(heading*constants::RADDEG,zaxis);
  yaxis = yaxis.Rotate(heading*constants::RADDEG,zaxis);

  // do pitch...
  ant_pos = ant_pos.Rotate(thisPitch*constants::RADDEG,yaxis);
  xaxis = xaxis.Rotate(thisPitch*constants::RADDEG,yaxis);

  // do roll...
  ant_pos = ant_pos.Rotate(thisRoll*constants::RADDEG,xaxis);//roll and pitch

  ////now place balloon at proper lat and lon
  // BalloonPhi =latitude*constants::RADDEG;
  ant_pos = ant_pos.RotateY(BalloonTheta);
  ant_pos =  ant_pos.RotateZ(BalloonPhi);
  
  // this->x_axis_rot = xaxis;
  // this->y_axis_rot = yaxis;
  // this->z_axis_rot = zaxis;
  return ant_pos;
}  

icemc::Vector icemc::Balloon::unRotatePayload(Vector ant_pos_pre) const {//rotate back to Payload Centric coordinates
  ///@todo does we have static pitch and roll offsets for A4?
  const double thisPitch = WHICHPATH == FlightPath::Anita2 ? fixedAnita2Pitch : WHICHPATH == FlightPath::Anita3 ? fixedAnita3Pitch : pitch;
  const double thisRoll  = WHICHPATH == FlightPath::Anita2 ? fixedAnita2Roll  : WHICHPATH == FlightPath::Anita3 ? fixedAnita3Roll  : roll;
  
  Vector BalloonPos;

  BalloonPos = r_bn;
  const double BalloonTheta = BalloonPos.Theta();
  double BalloonPhi = BalloonPos.Phi();
  
  //std::cout <<"BalloonTheta is "<<BalloonTheta<<" phi is "<<BalloonPhi<<"\n";
  if(BalloonPhi > constants::PI){
    BalloonPhi =BalloonPhi-constants::TWOPI;
  }
  
 
  Vector ant_pos = ant_pos_pre;

  // regenerate the pitch/roll axes here to allow constification...
  const Vector zaxis(0.,0.,-1.);
  Vector xaxis(1.,0.,0.);  // roll axis
  Vector yaxis(0.,-1.,0.); // pitch axis for positive rotation to the clockwise of roll
  // do heading...
  xaxis = xaxis.Rotate(heading*constants::RADDEG,zaxis);
  yaxis = yaxis.Rotate(heading*constants::RADDEG,zaxis);

  // do pitch...
  xaxis = xaxis.Rotate(thisPitch*constants::RADDEG,yaxis);  
  
  //rotate to correct heading, roll and pitch

  
  ant_pos = ant_pos.RotateZ(-1*BalloonPhi);

  ant_pos = ant_pos.RotateY(-1*BalloonTheta);

  ant_pos = ant_pos.Rotate(-1*thisRoll*constants::RADDEG,xaxis);//roll and pitch

  ant_pos = ant_pos.Rotate(-1*thisPitch*constants::RADDEG,yaxis);

  ant_pos = ant_pos.Rotate(-1*heading*constants::RADDEG,zaxis);

  return ant_pos;
}  


