#include "Earth.h"
#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include <cmath>
#include "Tools.h"
#include "TVector3.h"
#include "GeoidModel.h"
#include <iostream>
#include <fstream>

#include "AskaryanFreqsGenerator.h"
#include "Primaries.h"
#include "secondaries.hh"
#include "EnvironmentVariable.h"
#include "Constants.h"

#include "TFile.h" ///@todo remove after debugging

#include "RampdemReader.h"



const std::string ICEMC_SRC_DIR=icemc::EnvironmentVariable::ICEMC_SRC_DIR();
const std::string ICEMC_DATA_DIR=ICEMC_SRC_DIR+"/data/";
// input files for Crust 2.0
const std::string crust20_in=ICEMC_DATA_DIR+"/outcr"; // Crust 2.0 data

const double icemc::Earth::COASTLINE(-60); //30.);
// const double icemc::Earth::MAXTHETA(180.);
// const int icemc::Earth::ILAT_COASTLINE((int)((COASTLINE/MAXTHETA)*(double)NLAT+0.00001)); // corresponding latitude bin to "coastline"
// const double icemc::Earth::GEOID_MAX(6.378137E6); // parameters of geoid model
// const double icemc::Earth::GEOID_MIN(6.356752E6); // from Geodetic Reference System 1980, Bulletin Geodesique, Vol 54:395,1980. // The previous reference gave issue number 3 instead of page number 395


icemc::Earth::Earth(int model,int WEIGHTABSORPTION_SETTING) {

  fLayerNames.emplace(CrustLayer::Water        , std::string("water"));
  fLayerNames.emplace(CrustLayer::Ice          , std::string("ice"));
  fLayerNames.emplace(CrustLayer::SoftSediment , std::string("hard sed"));
  fLayerNames.emplace(CrustLayer::HardSediment , std::string("soft sed"));
  fLayerNames.emplace(CrustLayer::UpperCrust   , std::string("upper crust"));
  fLayerNames.emplace(CrustLayer::MiddleCrust  , std::string("middle crust"));
  fLayerNames.emplace(CrustLayer::LowerCrust   , std::string("lower crust"));

  fPropertyNames.emplace(CrustProperty::ThicknessKm, std::string("thickness (km)"));
  fPropertyNames.emplace(CrustProperty::ElevationLowerBoundary, std::string("Elevation at the lower boundary (m)"));
  fPropertyNames.emplace(CrustProperty::Density,     std::string("density"));

  radii[0]=1.2e13;
  radii[1]=(Earth::BulgeRadius-4.0E4)*(Earth::BulgeRadius-4.0E4);
  radii[2]=(Earth::BulgeRadius*Earth::BulgeRadius); // average radii of boundaries between earth layers

  //  cout << "In Earth, model is " << model << "\n";
  weightabsorption= WEIGHTABSORPTION_SETTING;

  CONSTANTICETHICKNESS = (int) (model / 1000);
  model -= CONSTANTICETHICKNESS * 1000;

  CONSTANTCRUST = (int) (model / 100);
  model -= CONSTANTCRUST * 100;

  FIXEDELEVATION = (int) (model / 10);
  model -= FIXEDELEVATION * 10;

  EARTH_MODEL = model;
  //cout<<"CONSTICETHK = "<<CONSTANTICETHICKNESS<<", CNSTCRST = "<<CONSTANTCRUST<<", FIXDELV = "<<FIXEDELEVATION<<", EARTH_MODEL = "<<EARTH_MODEL<<endl;
  for (int i=0;i<NLON;i++) {
    // Tools::Zero(elevationarray[i],NLAT);
    // Tools::Zero(waterthkarray[i],NLAT);
    // Tools::Zero(icethkarray[i],NLAT);
    // Tools::Zero(softsedthkarray[i],NLAT);
    // Tools::Zero(hardsedthkarray[i],NLAT);
    // Tools::Zero(uppercrustthkarray[i],NLAT);
    // Tools::Zero(middlecrustthkarray[i],NLAT);
    // Tools::Zero(lowercrustthkarray[i],NLAT);
    // Tools::Zero(crustthkarray[i],NLAT);

    // Tools::Zero(surfacer[i],NLAT);
    // Tools::Zero(icer[i],NLAT);
    // Tools::Zero(waterr[i],NLAT);
    // Tools::Zero(softsedr[i],NLAT);
    // Tools::Zero(hardsedr[i],NLAT);
    // Tools::Zero(uppercrustr[i],NLAT);
    // Tools::Zero(middlecrustr[i],NLAT);
    // Tools::Zero(lowercrustr[i],NLAT);

    // Tools::Zero(waterdensityarray[i],NLAT);
    // Tools::Zero(icedensityarray[i],NLAT);
    // Tools::Zero(softseddensityarray[i],NLAT);
    // Tools::Zero(hardseddensityarray[i],NLAT);
    // Tools::Zero(uppercrustdensityarray[i],NLAT);
    // Tools::Zero(middlecrustdensityarray[i],NLAT);
    // Tools::Zero(lowercrustdensityarray[i],NLAT);
  } //Zero Earth model arrays
    
    // see monte carlo note #17
  // for (int i=0;i<NLAT;i++) {
  //   geoid[i]=GEOID_MIN*GEOID_MAX/sqrt(pow(GEOID_MIN,2.)-(pow(GEOID_MIN,2.)-pow(GEOID_MAX,2.))*pow(cos(dGetTheta(i)),2.));
  // } //for

  // Crust 2.0 is binned in 2deg x 2deg bins, area of bin depends on latitude.
  // calculating surface area of bins
  // phistep=2*constants::PI/(double)NLON;
  // thetastep=(MAXTHETA*constants::RADDEG)/NLAT;
  // for (int i=0;i<NLAT;i++) {
  //   area[i]=phistep*(cos(dGetTheta(i))-cos(dGetTheta(i+1)))*pow(geoid[i],2.);
  // } //for
    
  if (EARTH_MODEL == 0){
    ReadCrust(crust20_in);
  }
  else {
    std::cout<<"Error!  Unknown Earth model requested!  Defaulting to Crust 2.0 model.\n";
    ReadCrust(crust20_in);
  } //else
} //Earth constructor (int mode)


icemc::Earth::~Earth() {
  std::map<CrustKey, TH2D*>& m = fRawCrustData;
  for(auto it = m.begin(); it!=m.end(); ++it){
    if(it->second){
      delete it->second;
    }
  }

  if(fCrustElevation){
    delete fCrustElevation;
  }
  if(fSurfaceAboveGeoid){
    delete fSurfaceAboveGeoid;
  }
} 


const std::string& icemc::Earth::getLayerName(CrustLayer layer) const {
  const std::string& name = findThingInMapToString(fLayerNames, layer);
  if(name=="unknown"){
    icemcLog() << icemc::error << "Requested name of unknown layer!" << std::endl;
  }
  return name;
}

const std::string& icemc::Earth::getPropertyName(CrustProperty property) const {
  const std::string& name = findThingInMapToString(fPropertyNames, property);
  if(name=="unknown"){
    icemcLog() << icemc::error << "Requested name of unknown propery!" << std::endl;
  }
  return name;
}


// const std::string& icemc::Earth::getLayerName(CrustLayer layer) const {
//   auto it = fLayerNames.find(layer);
//   if(it!=fLayerNames.end()){
//     return it->second;
//   }
//   else{
//     icemcLog() << icemc::error << "Requested name of unknown layer!" << std::endl;
//     static const std::string unknownLayer("unknown layer");
//     return unknownLayer;
//   }
// }


icemc::Earth::CrustLayer icemc::Earth::getLayerFromString(const std::string& layerType) const {  
  for(auto& it : fLayerNames){
    const std::string& name = it.second;
    if(layerType.substr(0, name.length())==name){      
      return it.first;
    }    
  }
  icemcLog() << icemc::error << "Could not deduce CrustLayer from string: " << layerType << std::endl;
  return CrustLayer::LowerCrust;  
}







constexpr double minHistLon = -180;
constexpr double maxHistLon = 182;
constexpr int nBinsLon = (maxHistLon - minHistLon)/2;
constexpr double minHistLat = -90;
constexpr double maxHistLat = 90;
constexpr int nBinsLat = (maxHistLat - minHistLat)/2;

inline TH2D* makeRawDataHist(const TString& name){
  TH2D* h = new TH2D(name, name, nBinsLon, minHistLon, maxHistLon, nBinsLat, minHistLat, maxHistLat);
  h->GetXaxis()->SetTitle("Longitude (Degrees)");
  h->GetYaxis()->SetTitle("Latitude (Degrees)");
  return h;
}

inline void fillWrapped(TH2D* h, double lon, double lat, double val){
  h->Fill(lon, lat, val);
  if(lon + 360 < maxHistLon){
    h->Fill(lon+360, lat, val);	
  }
}



TH2D* icemc::Earth::makeRawHist(CrustLayer layer, CrustProperty property){
  CrustKey key(layer, property);
  TString name("h_");
  name += getLayerName(layer);
  name += "_";
  name += getPropertyName(property);
  name.ReplaceAll(" ",  "_");
  name.ReplaceAll("(",  "");
  name.ReplaceAll(")",  "");

  TH2D* h = makeRawDataHist(name);
  fRawCrustData[key] = h;
  return h;
}



TH2D* icemc::Earth::getRawHist(CrustLayer layer, CrustProperty property) const {
  CrustKey key(layer, property);
  auto it = fRawCrustData.find(key);
  if(it!=fRawCrustData.end()){
    return it->second;
  }
  else{
    return nullptr;
  }
}





double icemc::Earth::IceThickness(double lon,double lat) const {
  return getRawHist(CrustLayer::Ice, CrustProperty::ThicknessKm)->Interpolate(lon, lat);
}

double icemc::Earth::IceThickness(const GeoidModel::Position&pos) const {
  return IceThickness(pos.Longitude(), pos.Latitude());
}


int icemc::Earth::InFirn(const GeoidModel::Position&pos) const {
  if (pos.Mag()-Surface(pos)<constants::FIRNDEPTH){
    return 0;
  }
  return 1;
}


double icemc::Earth::SurfaceDeepIce(const GeoidModel::Position&pos) const { // surface of the deep ice (where you reach the firn)
  return  Surface(pos) + constants::FIRNDEPTH;
}

double icemc::Earth::Surface(double lon,double lat) const {
  GeoidModel::Position p;
  p.SetLonLatAlt(lon, lat, 0);
  return Surface(p);
} //Surface(lon,lat)

double icemc::Earth::Surface(const GeoidModel::Position&pos) const {
  return pos.Mag() + SurfaceAboveGeoid(pos);
  // return surfacer[(int)((pos.Longitude() < 0 ? pos.Longitude() + 360 : pos.Longitude())/2)][(int)(-pos.Latitude()/2)] + geoid[(int)(-pos.Latitude()/2)];
} //Surface(Position)




double icemc::Earth::RockSurface(double lon,double lat) const {
  return fCrustElevation->Interpolate(lon, lat); //(Surface(lon,lat) - IceThickness(lon,lat) - WaterDepth(lon,lat));
} //RockSurface(lon,lat)

double icemc::Earth::RockSurface(const GeoidModel::Position&pos) const {
  return RockSurface(pos.Longitude(), pos.Latitude());
} //RockSurface(lon,lat)




double icemc::Earth::SurfaceAboveGeoid(double lon,double lat) const {
  return fSurfaceAboveGeoid->Interpolate(lon, lat);
  // return surfacer[(int)(lon/2)][(int)(lat/2)];
} //SurfaceAboveGeoid(lon,lat)



double icemc::Earth::SurfaceAboveGeoid(const GeoidModel::Position&pos) const {
  return fSurfaceAboveGeoid->Interpolate(pos.Longitude(), pos.Latitude());
} //SurfaceAboveGeoid(Position)

double icemc::Earth::WaterDepth(double lon,double lat) const {
  return getRawHist(CrustLayer::Water, CrustProperty::ThicknessKm)->Interpolate(lon, lat);  
  // return waterthkarray[(int)(lon/2)][(int)(lat/2)]*1000;
} //WaterDepth(lon,lat)

double icemc::Earth::WaterDepth(const GeoidModel::Position&pos) const {
  return WaterDepth(pos.Longitude(), pos.Latitude());
} //WaterDepth(Position)


double icemc::Earth::GetDensity(const GeoidModel::Position& pos,
				int& crust_entered /* 1 or 0 */) const{
  
  double lon = pos.Longitude();
  double lat = pos.Latitude();
  double alt = pos.Altitude();
  double ddensity =0; //initilize ddensity

  if(alt > SurfaceAboveGeoid(pos)){
    ddensity = 1.25; /// @todo better model than this?
  }
  else {
    // this is actually pretty inefficient, but reads nicely and will do for now...
    for(auto layer : CrustLayers() ){
      double elevationLow = getRawHist(layer,  CrustProperty::ElevationLowerBoundary)->Interpolate(lon, lat);
      if(alt >= elevationLow){ //  
	ddensity = getRawHist(layer,  CrustProperty::Density)->Interpolate(lon, lat);
      }
      else{ // if you're below the lower boundary then the last density was correct
	break;
      }
    }
  }
  if(ddensity == 0){
    ddensity = densities[1]; // I guess this a mantle thing for below the moho?
  }
  return ddensity;

}//Get Density


GeoidModel::Position icemc::Earth::PickInteractionLocation(const GeoidModel::Position &balloon, Interaction *interaction1) const {

  GeoidModel::Position nuPos;

  const double MAX_HORIZON_DIST = 800e3; ///@todo Get this from settings or something?
  bool iceThickEnough = false;

  // This is a pretty important statement about icemc, we sample the ICE uniformly.
  // To do that we first sample the x/y plane uniformly within the horizon radius,
  // then we weight randomly chosen positions by ice thickness such that places where
  // the ice is twice as thick are twice as likely to be chosen.
  while(!iceThickEnough){
    double dx, dy;
    do{ // first pick any point within the horizon radius
      dx = MAX_HORIZON_DIST*gRandom->Rndm();
      dy = MAX_HORIZON_DIST*gRandom->Rndm();
    } while(dx*dx* + dy*dy > MAX_HORIZON_DIST*MAX_HORIZON_DIST);

    TVector3 deltaPos(dx, dy, 0);
    nuPos = balloon + deltaPos;

    // then roll a dice to see if it's deep enough
    double iceThicknessHere = this->IceThickness(nuPos);
    double randomHeight = gRandom->Uniform()*GetMaxIceThickness();
    if(randomHeight < iceThicknessHere){
      iceThickEnough = true;
      double elevationIceBottom = get(CrustLayer::Ice, CrustProperty::ElevationLowerBoundary, nuPos);
      const int kmToMetres = 1000;
      nuPos.SetAltitude(elevationIceBottom + kmToMetres*randomHeight);
    }
  }
  return nuPos;
}



int icemc::Earth::Getchord(const Settings *settings1,
				double len_int_kgm2,
				const GeoidModel::Position &earth_in, // place where neutrino entered the earth
				const GeoidModel::Position &r_enterice,
				const GeoidModel::Position &nuexitice,
				const GeoidModel::Position &posnu, // position of the interaction
				int inu,
				// ChordInfo& chordInfo){
				double& chord, // chord length
				double& probability_tmp, // weight
				double& weight1_tmp,
				double& nearthlayers, // core, mantle, crust
				double myair,
				double& total_kgm2, // length in kg m^2
				int& crust_entered, // 1 or 0
				int& mantle_entered, // 1 or 0
				int& core_entered)  {
    
  TVector3 chord3;
  TVector3 nchord;
  double x=0;
  double lat,lon;
  // int ilon,ilat;
    
  total_kgm2 = 0; //Initialize column density
  nearthlayers=0; // this counts the number of earth layers the neutrino traverses.
  // Want to find probability that the neutrino survives its trip
  // through the earth.
    
  //Find the chord, its length and its unit vector.
  chord3 = posnu - earth_in;
  chord=chord3.Mag();
  nchord = chord3*(1./chord);
    
  if (chord<=1) {
    std::cout << "short chord " << chord << "\n";
    return 0;
  }
  // if (chord>2.*R_EARTH+1000) {
  if (chord>2.*BulgeRadius+1000) {      
    std::cout << "bad chord" << " " << chord << ".  Event is " << inu << "\n";
  }
    
  GeoidModel::Position where=earth_in;
  //cout <<"where(1) is "<<where;
  // the sin of the angle between the neutrino path and the 
  // radial vector to its earth entrance point determines
  // if it will get to the next layer down.
  double costh=where.Dot(nchord)/where.Mag();
  double sinth=sqrt(1-costh*costh);
  double distance=0;
  double halfchord=0;

  if (getchord_method<1 || getchord_method>3){
    std::cout << "Bogus method!\n";
  }
    
  // we are really focusing on method 2 - method 1 has not been maintenanced in a long time!!
  // use at your own risk.
  if (getchord_method==1) {
    double L=0;
    weight1_tmp=0;
	
    if (sinth>sqrt(radii[1]/radii[2])) {
      nearthlayers++;
	    
      // these only skim the first layer.
      L=len_int_kgm2/densities[2];
	    
      weight1_tmp=exp(-posnu.Distance(where)/L);
    }
    else {
      nearthlayers++;
	    
      // these get to the second layer down.
      L=len_int_kgm2/densities[2];
      // compute distance before it gets to the next layer.
      halfchord=sqrt(radii[1]-radii[2]*sinth*sinth);
      distance=sqrt(radii[2])*costh-halfchord;
	    
      weight1_tmp=exp(-distance/L);
	    
      // get position where it enters the second layer.
      where = earth_in + distance*nchord;
	    
      // determine if it enters the core or not.
      costh=(where.Dot(nchord))/where.Mag();
      sinth=sqrt(1-costh*costh);
	    
      if (sinth>sqrt(radii[0]/radii[1])) {
		
		
	halfchord=sqrt(radii[1])*costh;
	nearthlayers++;
		
	// these do not enter the core.
	L=len_int_kgm2/densities[1];
		
		
	weight1_tmp *= exp(-2*halfchord/L); 
		
	L=len_int_kgm2/densities[2];
	// this is where it exits the second layer and enters the crust again.
	where = where + 2*halfchord*nchord;
	weight1_tmp*=exp(-where.Distance(posnu)/L);
      }
      else {
	nearthlayers++;
	// these enter the core.
	L=len_int_kgm2/densities[1];
		
	// compute distance before entering the core.
	halfchord=sqrt(radii[0]-radii[1]*sinth*sinth);
	distance=sqrt(radii[1])*costh-halfchord;
	weight1_tmp*=exp(-distance/L);
		
	// go through the core.
	L=len_int_kgm2/densities[0];
	weight1_tmp*=exp(-2*halfchord/L);
		
	// go through 2nd layer again.
	L=len_int_kgm2/densities[1];
	weight1_tmp*=exp(-distance/L);
		
	// through the crust and end up at posnu.
	L=len_int_kgm2/densities[2];
		
	where = where + (2*distance+2*halfchord)*nchord;
	weight1_tmp*=exp(-where.Distance(posnu)/L);
      } //else
    } //else
  } //if getchord_method==1
  if (getchord_method==2) {

    x=0; // x is the distance you move as you step through the earth.

    lon = where.Longitude();
    lat = where.Latitude();
    // ilon = (int)(lon/2);
    // ilat = (int)(lat/2);
	
    double surface_elevation = this->SurfaceAboveGeoid(where);  //lon,lat); // altitude of surface relative to geoid at earth entrance point
	
    // double local_icethickness = this->IceThickness(lon,lat);
    // double local_waterdepth = WaterDepth(lon,lat);
    // double altitude=0;
    weight1_tmp=1;
    probability_tmp=1;
    double step=TMath::Min(len_int_kgm2/densities[1]/10,500.); //how big is the step size
    // either 1/10 of an interaction length in the mantle or 500 m, whichever is smaller.
    // 500 m is approximately the feature size in Crust 2.0.
    //------------------added on Dec 8------------------------
    weight1_tmp*=exp(-myair/len_int_kgm2);//add atmosphere attenuation // fenfang's atten. due to atmosphere
    //------------------added on Dec 8------------------------
    total_kgm2+=myair;
	
	
	
    double L_ice=len_int_kgm2/AskaryanFreqsGenerator::RHOICE;
	
    if (settings1->UNBIASED_SELECTION){
      probability_tmp*=1.-exp(-1.*(r_enterice.Distance(nuexitice)/L_ice)); // probability it interacts in ice along its path
    }	
    double L=0;
	
    double ddensity=AskaryanFreqsGenerator::RHOAIR;
    nearthlayers=1;

    if (where.Dot(nchord)>0.)  { // look at direction of neutrino where it enters the earth.
      std::cout << "This one's trouble.  Neutrino exit point looks more like an entrance point.  Event is " << inu << "\n";
      std::cout << "where is " << where[0] << " " << where[1] << " " << where[2] << "\n";
      std::cout << "nchord is " << nchord[0] << " " << nchord[1] << " " << nchord[2] << "\n";
      std::cout << "dot product is " << where.Dot(nchord)/sqrt(where.Dot(where)) << "\n";
      std::cout << "posnu is " << posnu[0] << " " << posnu[1] << " " << posnu[2] << "\n";
      std::cout << "Length of chord is : "<<chord<<std::endl;
    } //end if
	
    double altitude=where.Altitude(); //Mag()-Geoid(lat); // what is the altitude of the entrance point
	
    if(altitude>surface_elevation+0.1){ // if it is above the surface, it's messed up
      std::cout << "neutrino entrance point is above the surface.  Event is " << inu << "\n";
    }
    //cout <<"altitude is "<<altitude<<"\n";

    while(altitude>MIN_ALTITUDE_CRUST && x<posnu.Distance(earth_in)) { // starting at earth entrance point, step toward interaction position until you've reached the interaction or you are below the crust.
      //    while(altitude>MIN_ALTITUDE_CRUST && x<dDistance(enterice,earth_in)) {
	  
      ddensity = this->GetDensity(where,crust_entered);

      L=len_int_kgm2/ddensity; // get the interaction length for that density
      weight1_tmp*=exp(-step/L);  // adjust the weight accordingly

      total_kgm2+=ddensity*step; //increase column density accordingly
      if (exp(-step/L) > 1){
	std::cout<<"Oops! len_int_kgm2, ddensity, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<exp(-step/L)<<std::endl;
      }
      x+=step; // distance you have stepped through the earth so far.

      where += step*nchord;// find where you are now along the neutrino's path 
	    
      lon = where.Longitude();
      lat = where.Latitude();
      // ilon = (int)(lon/2);
      // ilat = (int)(lat/2);
      altitude=where.Altitude(); //Mag()-Geoid(lat); //what is the altitude
      surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
      // local_icethickness = this->IceThickness(lon,lat);
      // local_waterdepth = WaterDepth(lon,lat);
	    
    } //end while
	
    if (x>posnu.Distance(earth_in) && weightabsorption) { // if you left the loop because you have already stepped the whole distance from the entrance point to the neutrino interaction position
      probability_tmp*=weight1_tmp;
      return 1;
    }
    // if you left the loop because you're not in the crust anymore
    if (altitude<=MIN_ALTITUDE_CRUST) {
	    
      mantle_entered=1; //Switch that lets us know we're into the mantle
	    
      // determine if it enters the core or not.
      sinth=sin(where.Angle(nchord));
      costh=sqrt(1-sinth*sinth);
	    
      if (sinth>sqrt(radii[0]/radii[1])) { // it does not enter the core, just the mantle
		
	nearthlayers++;  // count the mantle as a layer traversed.
		
	L=len_int_kgm2/densities[1]; // interaction length in the mantle
	halfchord=sqrt(radii[1])*costh; // 1/2 chord the neutrino goes through in the mantle
		
	weight1_tmp *= exp(-2*halfchord/L); // adjust the weight for path through mantle
	total_kgm2+= 2*halfchord*densities[1];  //add column density for path through mantle 
	where += (2*halfchord)*nchord; // neutrino's new position once it reaches the crust again
		
      } //end if (not entering core)
      // these enter the core
      else {
	core_entered=1; //Switch that lets us know we've entered the core
		
	nearthlayers+=2; // count the mantle and core as a layer traversed.
		
	L=len_int_kgm2/densities[1]; // interaction length in mantle
		
	// compute distance before entering the core.
	halfchord=sqrt(radii[0]-radii[1]*sinth*sinth); // find distance it travels in the mantle
	distance=sqrt(radii[1])*costh-halfchord;
	weight1_tmp*=exp(-distance/L); // adjust the weight
	total_kgm2 += 2*distance*densities[1]; //Add column density for trip through mantle
	// go through the core.
	L=len_int_kgm2/densities[0]; // interaction length in core
	weight1_tmp*=exp(-2*halfchord/L); // adjust the weight
	total_kgm2 += 2*halfchord*densities[0]; //Add column density for trip through core
	// go through 2nd layer again.
	L=len_int_kgm2/densities[1];
	weight1_tmp*=exp(-distance/L);
	if (exp(-distance/L) > 1){
	  std::cout<<"Oops2! len_int_kgm2, ddensity, distance, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<distance<<" , "<<exp(-distance/L)<<std::endl;
	}
	where += (2*distance+2*halfchord)*nchord;  // neutrino's new position once it reaches the crust again
		
      } //end else(enter core)
    } //end if(left crust)
	
    lon = where.Longitude();
    lat = where.Latitude();
    // ilon = (int)(lon/2);
    // ilat = (int)(lat/2);
    altitude=where.Altitude(); //.Mag()-Geoid(lat); //what is the altitude
    surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
    // local_icethickness = this->IceThickness(lon,lat);
    // local_waterdepth = WaterDepth(lon,lat);
	
    double distance_remaining=where.Distance(posnu); // how much farther you need to travel before you reach the neutrino interaction point
	
    x=0; // this keeps track of how far you've stepped along the neutrino path, starting at the crust entrance.
    while(x<=distance_remaining) { // keep going until you have reached the interaction position
	    
      ddensity=this->GetDensity(where,crust_entered);
	  
      L=len_int_kgm2/ddensity; // get the interaction length for that density
      weight1_tmp*=exp(-step/L);  // adjust the weight accordingly
      total_kgm2 += step*ddensity;
	    
      if (exp(-step/L) > 1){
	std::cout<<"Oops3! len_int_kgm2, ddensity, step, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<step<<" , "<<exp(-step/L)<<std::endl;
      }
      x+=step; // increment how far you've stepped through crust

      // possible for a neutrino to go through the air but not likely because they aren't the most extreme skimmers (they went through the mantle)
      where += step*nchord; // where you are now along neutrino's path
	    
      lon = where.Longitude();
      lat = where.Latitude();
      // ilon = (int)(lon/2);
      // ilat = (int)(lat/2);
      altitude=where.Altitude(); //where.Mag()-Geoid(lat); //what is the altitude
      surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
      // local_icethickness = this->IceThickness(lon,lat);
      // local_waterdepth = WaterDepth(lon,lat);
    } //while
	
  } //if (getchord_method == 2)
    
  probability_tmp*=weight1_tmp;
  //cout <<"probability_tmp(non-tau) is "<<probability_tmp<<".\n";
    
  if (weightabsorption==0) {
    if (gRandom->Rndm()>weight1_tmp) { 
      weight1_tmp=0.;
      return 0;
    }
    else {
      weight1_tmp=1.;
      return 1;
    }
  }
  else {
    return 1;
  }
  std::cout << "made it this far.\n";    
  return 1;
} //end Getchord

TVector3 icemc::Earth::GetSurfaceNormal(const GeoidModel::Position &r_out) const
{
  ///@todo rewrite this method!

  GeoidModel::Position justAbove = r_out;
  justAbove.SetAltitude(r_out.Altitude()+1);
  
  TVector3 n_surf = justAbove - r_out;
  
  // if (FLATSURFACE){
  //   return n_surf;
  // }    
    
  return n_surf;
    
} //method GetSurfaceNormal

// double icemc::Earth::SmearPhi(int ilon, double rand) const {
    
    
//   double phi=((double)(360.*3./4.-((double)ilon+rand)*360/180))*constants::RADDEG;
//   if (phi<0 && phi>-1*constants::PI/2)
//     phi+=2*constants::PI;
    
    
//   return phi;
// } //SmearPhi

// double icemc::Earth::SmearTheta(int ilat, double rand) const {
    
//   // remember that we should smear it evenly in cos(theta).
//   // first get the cos(theta)'s at the boundaries.
    
//   double theta1=dGetTheta(ilat)-constants::PI/(double)NLAT/2.;
//   double theta2=dGetTheta(ilat+1)-constants::PI/(double)NLAT/2.;
    
//   double costheta1=cos(theta1);
//   double costheta2=cos(theta2);
    
    
    
//   double costheta=rand*(costheta2-costheta1)+costheta1;
    
//   double theta=acos(costheta);
    
//   return theta;
// } //SmearTheta








void icemc::Earth::ReadCrust(const std::string& fName) {

  // reads in altitudes of 7 layers of crust, ice and water
  // puts data in arrays

  std::fstream infile(fName.c_str(),std::ios::in);

  std::string thisline; // for reading in file
  std::string slon; //longitude as a string
  std::string slat; // latitude as a string
  std::string selev; // elevation (km relative to geoid)
  std::string sdepth; // depth (km)
  std::string sdensity; // density (g/cm^3)
  int endindex; // index along thisline for parsing
  int beginindex; // same

  std::string layertype; // water, ice, etc.

  int indexlon = 0;
  int indexlat = 0;
  double dlon, dlat;

  while(!infile.eof()) {
    getline(infile,thisline,'\n');

    int loc=thisline.find("type, latitude, longitude,"); 
	
    if (loc!=(int)(std::string::npos)) {      
 
      beginindex=thisline.find_first_not_of(" ",57);

      endindex=thisline.find_first_of(" ",61);

      slat=thisline.substr(beginindex,endindex-beginindex);
      dlat=(double)atof(slat.c_str());

      beginindex=thisline.find_first_not_of(" ",68);
      endindex=thisline.find_first_of(" ",72);

      slon=thisline.substr(beginindex,endindex-beginindex);
      dlon=(double)atof(slon.c_str());

      indexlon=(int)((dlon+180)/2);
      indexlat=(int)((90+dlat)/2);

      beginindex=thisline.find_first_not_of(" ",78);
      endindex=thisline.find_first_of(" ",83);

      selev=thisline.substr(beginindex,endindex-beginindex);
      double elevation = (double)atof(selev.c_str());

      if(!fCrustElevation) {
	fCrustElevation = makeRawDataHist("fCrustElevation");
      }
      fillWrapped(fCrustElevation, dlon, dlat, elevation);

    } //if

    for (int i=0;i<4;i++) {
      getline(infile,thisline,'\n');
    } //for

    for (int i=0;i<7;i++) {
      getline(infile,thisline,'\n');

      endindex = thisline.length()-1;
      beginindex = thisline.find_last_of("0123456789",1000);
      layertype = thisline.substr(beginindex+3,endindex-beginindex);

      beginindex = thisline.find_first_not_of(" ",0);
      endindex = thisline.find_first_of(" ",beginindex);

      sdepth = thisline.substr(beginindex, endindex-beginindex-1);
      
      beginindex=thisline.find_first_not_of(" ",endindex);
      endindex=thisline.find_first_of(" ",beginindex);

      beginindex=thisline.find_first_not_of(" ",endindex);
      endindex=thisline.find_first_of(" ",beginindex);

      beginindex=thisline.find_first_not_of(" ",endindex);
      endindex=thisline.find_first_of(" ",beginindex);

      sdensity=thisline.substr(beginindex,endindex-beginindex);

      double density=(double)atof(sdensity.c_str());
      double depth = (double)atof(sdepth.c_str());;

      CrustLayer layer = getLayerFromString(layertype);
      // region where Ross Ice Shelf was not accounted for in Crust 2.0
      // add it in by hand
      if (layer==CrustLayer::Ice && indexlat==5 && (indexlon<=5 || indexlon>=176)){ // Ross Ice Shelf
	depth = 0.5;
      }

      TH2D* hThickness = getRawHist(layer, CrustProperty::ThicknessKm);
      if(!hThickness){
	hThickness = makeRawHist(layer, CrustProperty::ThicknessKm);
      }
      fillWrapped(hThickness, dlon, dlat, depth);

      if(layer == CrustLayer::Ice && depth > fMaxIceThickness){
	fMaxIceThickness = depth;
      }

      double centerLon = hThickness->GetXaxis()->GetBinCenter(hThickness->GetXaxis()->FindBin(dlon));
      double centerLat = hThickness->GetYaxis()->GetBinCenter(hThickness->GetYaxis()->FindBin(dlat));
      if(centerLon - dlon != 0 || centerLat != dlat){
	icemcLog() << icemc::warning << "Binning error in " << hThickness->GetName()
		   << ", raw lon = " << dlon << " binCenter = " << centerLon
		   << ", raw lat = " << dlat << " binCenter = " << centerLat
		   << std::endl;
      }

      TH2D* hDensity = getRawHist(layer, CrustProperty::Density);
      if(!hDensity){
	hDensity = makeRawHist(layer, CrustProperty::Density);
      }
      fillWrapped(hDensity,  dlon, dlat, density);
            
    } //for (reading all lines for one location given in Crust 2.0 input file)


    ///@todo restore constant crust values
    // if (CONSTANTCRUST) {
    //   softsedthkarray[indexlon][indexlat]=40.;
    //   hardsedthkarray[indexlon][indexlat]=0;
    //   uppercrustthkarray[indexlon][indexlat]=0;
    //   middlecrustthkarray[indexlon][indexlat]=0;
    //   lowercrustthkarray[indexlon][indexlat]=0;
    //   crustthkarray[indexlon][indexlat]=0;
    //   softseddensityarray[indexlon][indexlat]=2.9;
    // } //if (set crust thickness to constant everywhere)
    // if (CONSTANTICETHICKNESS) {
    //   icethkarray[indexlon][indexlat]=3.;
    //   waterthkarray[indexlon][indexlat]=0.;
    // } //if (set ice thickness to constant everywhere)

    // adds up total thickness of crust
    
    
    // crustthkarray[indexlon][indexlat] = (softsedthkarray[indexlon][indexlat]+
    // 					 hardsedthkarray[indexlon][indexlat]+
    // 					 uppercrustthkarray[indexlon][indexlat]+
    // 					 middlecrustthkarray[indexlon][indexlat]+
    // 					 lowercrustthkarray[indexlon][indexlat]);

    if (indexlon==179 && indexlat==0){
      break;
    }
  }  // done reading file

  TFile* fTest = new TFile("earthTest.root", "recreate");
  for(auto it : fRawCrustData){
    TH2D* h = it.second;
    h->SetDirectory(gDirectory);
  }
  fCrustElevation->SetDirectory(gDirectory);

  fTest->Write();
  fTest->Close();
  exit(1);



  
  if(!fSurfaceAboveGeoid){
    fSurfaceAboveGeoid = makeRawDataHist("fSurfaceAboveGeoid");
  }
  else{
    fSurfaceAboveGeoid->Clear();
  }
  const double kmToMetres = 1000;
  for(int by=1; by < fCrustElevation->GetNbinsY(); by++){
    for(int bx=1; bx < fCrustElevation->GetNbinsX(); bx++){
      double surfaceAboveGeoid = fCrustElevation->GetBinContent(bx,  by);
      for(auto layer : {CrustLayer::Ice, CrustLayer::Water}){
	surfaceAboveGeoid += kmToMetres*(getRawHist(layer, CrustProperty::ThicknessKm)->GetBinContent(bx, by));
      }
      fSurfaceAboveGeoid->SetBinContent(bx, by, surfaceAboveGeoid);
    }
  }

  // find the place where the crust is the deepest.
  // for finding where to start stepping in Getchord  
  MIN_ALTITUDE_CRUST = DBL_MAX;
  
  // now we have the surface which has everything underneath it, we can go back and deduce the elevation of each layer
  // by the subtracting the cumulative thickness of the layers from the surfaceAboveGeoid 
  // so the CurstPropery::ElevationLowerBoundary is, descriptively, the elevation at the *BOTTOM* of each layer.
  for(int by=1; by < fCrustElevation->GetNbinsY(); by++){
    for(int bx=1; bx < fCrustElevation->GetNbinsX(); bx++){
      double layerElevation = fSurfaceAboveGeoid->GetBinContent(bx, by);
      for(auto layer : CrustLayers() ){	
	double thicknessMetres = kmToMetres*(getRawHist(layer, CrustProperty::ThicknessKm)->GetBinContent(bx, by));
	layerElevation -= thicknessMetres;
	TH2D* hLayerElev = getRawHist(layer, CrustProperty::ElevationLowerBoundary);
	if(!hLayerElev){
	  hLayerElev = makeRawHist(layer, CrustProperty::ElevationLowerBoundary);
	}
	hLayerElev->SetBinContent(bx, by, layerElevation);

	if(layerElevation < MIN_ALTITUDE_CRUST){
	  MIN_ALTITUDE_CRUST = layerElevation;
	}
      }
    }
  }
 
  // calculate ice volume for comparison with expectation
  
  

  // volume=0;
  // ice_area=0.;
  // double sumarea=0; // sum of surface area of ice
  // average_iceth = 0;
  // max_icevol_perbin=0.;
  // max_icethk_perbin=0.;
  // for (int j=0;j<ILAT_COASTLINE;j++) {
  //   for (int i=0;i<NLON;i++) {
  //     volume+=icethkarray[i][j]*1000.*area[j];
  //     if (icethkarray[i][j]>0){
  // 	ice_area+=area[j];
  //     }
  //     if (icethkarray[i][j]*area[j]>max_icevol_perbin){
  // 	max_icevol_perbin=icethkarray[i][j]*area[j];
  //     }
  //     if (icethkarray[i][j]>max_icethk_perbin){
  // 	max_icethk_perbin=icethkarray[i][j];
  //     }	    
  //     /*
  //     // j=6 corresponds to 80deg S
  //     if (j==6) {
  //     // fill output file, just for plotting
  //     outfile << surfacer[i][j] << "\t" << waterr[i][j] << "\t" << icer[i][j] << "\t" << icethkarray[i][j] << " " << waterr[i][j]-icer[i][j] << " " << (waterr[i][j]-icer[i][j])*area[j] << " " << area[j] << " " << volume << "\n";
  //     }//ifx
  //     */
	    
  //     // find average ice thickness
  //     average_iceth+=(surfacer[i][j]-icer[i][j])*area[j];
  //     sumarea+=area[j];
  //   } //for
  // } //for
  // average_iceth=average_iceth/sumarea;
  
  //record depth of crust-mantle interface
  // radii[1]=(Geoid(0.)+MIN_ALTITUDE_CRUST)*(Geoid(0.)+MIN_ALTITUDE_CRUST);
  radii[1]=(GeoidModel::GEOID_MIN+MIN_ALTITUDE_CRUST)*(GeoidModel::GEOID_MIN+MIN_ALTITUDE_CRUST);

}//ReadCrust


// GeoidModel::Position icemc::Earth::PickPosnuForaLonLat(double lon,double lat) const {
    
    
//   double surface_above_geoid = this->SurfaceAboveGeoid(lon,lat);
//   double local_ice_thickness = this->IceThickness(lon,lat);
    
//   double rnd3=gRandom->Rndm();

//   double elevation = surface_above_geoid - rnd3*local_ice_thickness; // elevation of interaction
//   //cout << "Inside PickInteractionLocation, lon, lat, local_ice_thickness are " << lon << " " << lat << " " << local_ice_thickness << "\n";

//   if (elevation > surface_above_geoid || elevation < (surface_above_geoid - local_ice_thickness)){
//     std::cout<<"elevation > surface_above_geoid || elevation < (surface_above_geoid - local_ice_thickness)\n";
//   }
  
//   GeoidModel::Position posnu;
//   posnu.SetLonLatAlt(lon, lat, elevation);

//   return posnu;
// }





GeoidModel::Position icemc::Earth::WhereDoesItEnter(const GeoidModel::Position &posnu,const TVector3 &nnu) const {
  // now get neutrino entry point...
  double p = posnu.Mag(); // radius of interaction
  double costheta = (nnu.Dot(posnu)) / p; // theta of neutrino at interaction position
  double sintheta = sqrt(1-costheta*costheta);
    
  double lon = posnu.Longitude();
  double lat = posnu.Latitude();

  double a=0; // length of chord

  double R = Surface(lon,lat);
  double delta = R - p; // depth of the interaction
  // if interaction occurs below surface, as it should
    
  if (delta>-0.001) {
    a=p*costheta+sqrt(R*R*costheta*costheta+2*delta*R*sintheta*sintheta); // chord length
    if (a<0) {
      std::cout << "Negative chord length: " << a << "\n";
    } //end if
  } //end if (interaction below surface)  
  else if (delta<=-0.001) {

    std::cout << "lon, lat from WhereDoesItEnter is " << lon << " " << lat << "\n";
    std::cout << "geoid, surface, p, surface-p are " << GeoidModel::getGeoidRadiusAtLatitude(lat) << " " << Surface(lon,lat) << " " << p << " , "<<(Surface(lon,lat)-p)<<"\n";
	
  } //else if: error: interaction takes place above the surface


    // first approx
    // find where nnu intersects spherical earth surface
  GeoidModel::Position r_in = posnu - a*nnu;

  int iter = 0;
  // now do correction 3 times
  //for (iter=0; iter<3; iter++) {
  //    delta = r_in.Mag() - antarctica1->Surface( r_in );
  //    r_in = r_in + (delta * nnu);
  //}
    
  // how far away r_in is from surface, along line to center of earth
  delta = r_in.Mag() - Surface( r_in );
    
  while ( fabs(delta) >= 0.1 ) {  //if it's greater than 10 cm
    r_in = r_in + (delta * nnu); // move distance delta along nnu  for the next iteration
    delta = r_in.Mag() - Surface( r_in ); // find new delta
    iter++;
    if ( iter > 10 ) { // stop after 10 iterations and just revert to simple answer
      //cout<<"\n r_in iteration more than 10 times!!! delta : "<<delta<<". now set r_in as before."<<endl;
      r_in = Surface( r_in ) * r_in.Unit();   // the way before
      delta = r_in.Mag() - Surface( r_in );
    }
  }

    
  //lon = r_in.Lon();
  //lat = r_in.Lat();
    
  //r_in = antarctica1->Surface(lon,lat) * r_in.Unit();
    
    
  return r_in;
} //method WhereDoesItEnter
