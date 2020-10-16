#include "Crust2.h"
#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include <cmath>
#include "Tools.h"
#include "TVector3.h"
#include "Geoid.h"
#include <iostream>
#include <fstream>

#include "AskaryanRadiationModel.h"
#include "ConnollyEtAl2011.h"
#include "ShowerModel.h"
#include "EnvironmentVariable.h"
#include "Constants.h"

#include "Math/Delaunay2D.h"

#include "RampdemReader.h"
#include "LocalCoordinateSystem.h"


const double icemc::Crust2::COASTLINE(-60); //30.);

icemc::Crust2::Crust2(int model,int WEIGHTABSORPTION_SETTING) :
  fSurfaceAboveGeoid("fSurfaceAboveGeoid", "Surface above geoid")
{

  fLayerNames.emplace(Layer::Air,          std::string("air"));
  fLayerNames.emplace(Layer::Water,        std::string("water"));
  fLayerNames.emplace(Layer::Ice,          std::string("ice"));
  fLayerNames.emplace(Layer::SoftSediment, std::string("soft sed"));
  fLayerNames.emplace(Layer::HardSediment, std::string("hard sed"));
  fLayerNames.emplace(Layer::UpperCrust,   std::string("upper crust"));
  fLayerNames.emplace(Layer::MiddleCrust,  std::string("middle crust"));
  fLayerNames.emplace(Layer::LowerCrust,   std::string("lower crust"));
  fLayerNames.emplace(Layer::Mantle,       std::string("mantle"));
  fLayerNames.emplace(Layer::OuterCore,    std::string("outer core"));
  fLayerNames.emplace(Layer::InnerCore,    std::string("inner core"));
  
  namespace G=Geoid;
  radii[0]=1.2e13;
  radii[1]=(G::R_EARTH-4.0E4)*(G::R_EARTH-4.0E4);
  radii[2]=(G::R_EARTH*G::R_EARTH); // average radii of boundaries between earth layers

  //  cout << "In Earth, model is " << model << "\n";
  weightabsorption= WEIGHTABSORPTION_SETTING;

  CONSTANTICETHICKNESS = (int) (model / 1000);
  model -= CONSTANTICETHICKNESS * 1000;

  CONSTANTCRUST = (int) (model / 100);
  model -= CONSTANTCRUST * 100;

  FIXEDELEVATION = (int) (model / 10);
  model -= FIXEDELEVATION * 10;

  EARTH_MODEL = model;

  totalIceVolume = 0;
  totalIceArea = 0;
  
  const std::string ICEMC_SRC_DIR=icemc::EnvironmentVariable::ICEMC_SRC_DIR();
  const std::string crust20_in=ICEMC_SRC_DIR+"/data/outcr"; // Crust 2.0 data file

  if (EARTH_MODEL == 0){
    ReadCrust(crust20_in);
  }
  else {
    std::cout<<"Error!  Unknown Earth model requested!  Defaulting to Crust 2.0 model.\n";
    ReadCrust(crust20_in);
  } //else
} //Earth constructor (int mode)



icemc::Crust2::Layer icemc::Crust2::getLayerFromString(const std::string& layerType) const {  
  for(auto& it : fLayerNames){
    const std::string& name = it.second;
    // std::cout << "trying... " << name << " for " <<  layerType << std::endl;
    if(layerType.substr(0, name.length())==name){
      // std::cout << "Matched!" << std::endl;
      return it.first;
    }
  }
  icemc::report() << severity::error << "Could not deduce Layer from string: " << layerType << std::endl;
  return Layer::LowerCrust;
}


///@todo handle out of layer bounds somehow
icemc::Crust2::Layer icemc::Crust2::layerAbove(Layer layer) {
  switch (layer){
  case Layer::Air:          return Layer::Air;
  case Layer::Water:        return Layer::Air;
  case Layer::Ice:          return Layer::Water;
  case Layer::SoftSediment: return Layer::Ice;
  case Layer::HardSediment: return Layer::SoftSediment;
  case Layer::UpperCrust:   return Layer::HardSediment;
  case Layer::MiddleCrust:  return Layer::UpperCrust;
  case Layer::LowerCrust:   return Layer::MiddleCrust;
  case Layer::Mantle:       return Layer::LowerCrust;
  case Layer::OuterCore:    return Layer::Mantle;
  case Layer::InnerCore:    return Layer::OuterCore;
  }
  return layer;
}

///@todo handle out of layer bounds somehow
icemc::Crust2::Layer icemc::Crust2::layerBelow(Layer layer) {
  switch (layer){
  case Layer::Air:          return Layer::Water;    
  case Layer::Water:        return Layer::Ice;
  case Layer::Ice:          return Layer::SoftSediment;
  case Layer::SoftSediment: return Layer::HardSediment;
  case Layer::HardSediment: return Layer::UpperCrust;
  case Layer::UpperCrust:   return Layer::MiddleCrust;
  case Layer::MiddleCrust:  return Layer::LowerCrust;
  case Layer::LowerCrust:   return Layer::Mantle;
  case Layer::Mantle:       return Layer::OuterCore;
  case Layer::OuterCore:    return Layer::InnerCore;
  case Layer::InnerCore:    return Layer::InnerCore;
  }
  return layer;
}

double icemc::Crust2::IceVolumeWithinHorizon(const Geoid::Position& p, double horizonDistanceMeters) const {

  const double maxDeltaRSq = horizonDistanceMeters*horizonDistanceMeters;
  const double delta = 5e3;
  const LocalCoordinateSystem lc(p); // local coordate system tangential to geoid centered at p
  double v = 0;
  for(double dy = -horizonDistanceMeters; dy <= horizonDistanceMeters; dy += delta){
    for(double dx = -horizonDistanceMeters; dx <= horizonDistanceMeters; dx += delta){
      double drSq = dy*dy+dx*dx;
      if(drSq <= maxDeltaRSq){
	Geoid::Position p2 = lc.localPositionToGlobal(TVector3(dx, dy, 0));
	double t = IceThickness(p2);
	v += delta*delta*t;
      }
    }
  }
  
  return v;
}

double icemc::Crust2::IceThickness(const Geoid::Position& p) const {
  return fThicknesses.at(Layer::Ice).eval(p);
}

double icemc::Crust2::maxIceThicknessWithinDistance(const Geoid::Position& p, double distanceMeters) const {
  return fThicknesses.at(Layer::Ice).findMaxWithinDistance(p,  distanceMeters);
}


int icemc::Crust2::InFirn(const Geoid::Position&pos) const {
  if (pos.Mag()-Surface(pos)<constants::FIRNDEPTH){
    return 0;
  }
  return 1;
}


double icemc::Crust2::SurfaceDeepIce(const Geoid::Position&pos) const { // surface of the deep ice (where you reach the firn)
  return  Surface(pos) + constants::FIRNDEPTH;
}


double icemc::Crust2::SurfaceAboveGeoid(const Geoid::Position&pos) const {
  return fSurfaceAboveGeoid.eval(pos);
} //SurfaceAboveGeoid(Position)

double icemc::Crust2::WaterDepth(const Geoid::Position&pos) const {
  return fThicknesses.at(Layer::Water).eval(pos);  
} //WaterDepth(Position)




std::pair<Geoid::Position, double> icemc::Crust2::integratePath(const Geoid::Position& start, const TVector3& direction) const {

  double cumulativeColumnDepth = 0;
  Geoid::Position posNow = start;
  Geoid::Position posEnd = posNow;
  double fractionalStep = 1.0;
  Layer layerNow = getLayer(posNow);
  // int nSteps = 0, nSteps2 = 0;
  // std::cout << fLayerNames.at(layerNow) << "\n";

  while(layerNow > Layer::Air){
    layerNow = getLayer(posNow);
    // const std::string& ln = fLayerNames.find(l)->second;

    double density = Density(posNow);

    double cosTheta = posNow.Dot(direction)/posNow.Mag();
    const double minStepSize = 1.0;    
    const double maxStepSize = 1000e3;
    double distGuess = maxStepSize;
    double dr = 1;
    // const char* dir = "horizontal";    
    
    if(cosTheta > 0){
      // neutrino path is upgoing at this position
      // where up is radially outwards at p
      // but we are integrating downwards, so find distance to layer below
      dr = posNow.Mag() - getLayerSurfaceAtPosition(posNow, layerBelow(layerNow));
      distGuess = layerNow >= Layer::Mantle ? dr/cosTheta : distGuess;
      // dir = "(nu up) integrating down";
    }
    else if(cosTheta < 0){
      // neutrino path is downgoing at this position
      // where down is radially indwards at p
      // but we are integrating upwards, so find the distance to the surface of this layer.
      dr = getLayerSurfaceAtPosition(posNow, layerNow) - posNow.Mag();
      distGuess = layerNow >= Layer::Mantle ?  dr/-cosTheta : dr;
      // dir = "(nu down) integrating up";
    }

    double stepSize = fractionalStep*distGuess;


    stepSize = TMath::Max(stepSize, minStepSize);
    // std::cout << dir << "\t" << nSteps2 << "\t" << cosTheta << "\t" << TMath::ACos(cosTheta)*TMath::RadToDeg() << "\t" << ln << " dr = " << dr << "\t distGuess = " << distGuess << ", stepSize = " << stepSize << "\t...";
    posEnd = posNow + stepSize*direction;

    Layer l = getLayer(posEnd);
    if(l == layerNow || stepSize <= minStepSize){
      // if you're in the same layer
      // or you crossed with a step size of 1
      // update posNow,layernow,  reset step size fraction
      // and update columnDepth

      // if(l!=layerNow){
      // 	std::cout << fLayerNames.at(l) << "\n";
      // }
      
      posNow = posEnd;
      layerNow = l; // this allows us to break out of the loop on this iteration
      fractionalStep = 1;
      cumulativeColumnDepth += density*stepSize;
      // nSteps++;
    }
    else{
      // try again with a smaller step
      fractionalStep *= 0.5;
    }
    // nSteps2++;
  }

  // std::cout << nSteps << "\t" << nSteps2 << "\t" << cumulativeColumnDepth << "\n" << std::endl;
  
  return {posEnd, cumulativeColumnDepth};
}




double icemc::Crust2::getLayerSurfaceAtPosition(const Geoid::Position& pos, Layer layer) const {  

  switch(layer){
  case Layer::InnerCore: return 1220e3;
  case Layer::OuterCore: return 3400e3;
  default:               return fSurfaceMag.at(layer).eval(pos);
  }  
}



icemc::Crust2::Layer icemc::Crust2::getLayer(const Geoid::Position& pos, Layer startLayer) const {

  Layer inLayer = startLayer;
  double posMag = pos.Mag();
  for(auto layer : Layers()){
    if(layer <= startLayer) continue;

    double layerSurface = getLayerSurfaceAtPosition(pos, layer);

    // std::cout << "layer name = " << fLayerNames.at(layer) << "\n";
    // std::cout << "layer surface = " << layerSurface << "\n";
    // std::cout << "pos = " << posMag << "\n";    
    // std::cout << "delta (km) " << 1e-3*(posMag - layerSurface)  << "\n";
    // std::cout << "surface radius " << pos.EllipsoidSurface() + fSurfaceAboveGeoid.eval(pos) <<"\n";
    // std::cout << std::endl;
    
    if(posMag > layerSurface){
      // then we're in the layer from the previous loop!
      return inLayer;
    }
    inLayer = layer;
  }
  return inLayer;
}


bool icemc::Crust2::InsideIce(const Geoid::Position& pos) const {
  Layer l = getLayer(pos);
  // std::cout << fLayerNames.at(l) << std::endl;
  return l == Crust2::Layer::Ice;
}

double icemc::Crust2::Density(const Geoid::Position& pos,
				int& crust_entered /* 1 or 0 */) const{

  Layer l = getLayer(pos);

  switch(l){
  case Layer::Air:
    return 0; ///@todo Air model...
  case Layer::Mantle:
    return 3400; //densities[1]; ///@todo check this!

    // these values were just copied from hyperphysics...
    // it looks like the mantle density was too...
    // http://hyperphysics.phy-astr.gsu.edu/hbase/Geophys/earthstruct.html
  case Layer::OuterCore:
    return 9900;
  case Layer::InnerCore:
    return 1280;
  default:
    return fDensities.at(l).eval(pos);
  }
  
  // if(posMag > surfaceHere){
  //   // we're above the surface...
  //   /// @todo handle water!!!
  //   densityHere = 1.25;
  //   return densityHere;
  // }
  // else{
    
  //   for(auto layer : Layers()){
  //     double depth = fLayers.find(layer)->second.eval(pos);

  //     if(posMag < surfaceHere - depth){
  // 	// then we're in this layer, get the density!
  // 	densityHere = fDensities.find(layer)->second.eval(pos);
  // 	return densityHere;
  //     }
  //   }
  // }
  // // mantle below moho?
  // return densities[1];  
  
  // double lon = pos.Longitude();
  // lon = lon < 0 ? lon + 360 : lon;  
  // double lat = pos.Latitude();
  // double alt = pos.Altitude();
  // double ddensity =0; //initilize ddensity

  // if(alt > SurfaceAboveGeoid(pos)){
  //   ddensity = 1.25; /// @todo better model than this?
  // }
  // else {
  //   // this is actually pretty inefficient, but reads nicely and will do for now...
  //   for(auto layer : Layers() ){
  //     double elevationLow = getRawHist(layer,  CrustProperty::ElevationLowerBoundary)->Interpolate(lon, lat);
  //     if(alt >= elevationLow){ //  
  // 	ddensity = getRawHist(layer,  CrustProperty::Density)->Interpolate(lon, lat);
  //     }
  //     else{ // if you're below the lower boundary then the last density was correct
  // 	break;
  //     }
  //   }
  // }
  // if(ddensity == 0){
  //   ddensity = densities[1]; // I guess this a mantle thing for below the moho?
  // }
  // return ddensity;

}//Get Density







int icemc::Crust2::Getchord(const Settings *settings1,
			    double len_int_kgm2,
			    const Geoid::Position &earth_in, // place where neutrino entered the earth
			    const Geoid::Position &r_enterice,
			    const Geoid::Position &nuexitice,
			    const Geoid::Position &posnu, // position of the interaction
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
  if (chord>2.*Geoid::R_EARTH+1000) {      
    std::cout << "bad chord" << " " << chord << ".  Event is " << inu << "\n";
  }
    
  Geoid::Position where=earth_in;
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

    double L_ice=len_int_kgm2/AskaryanRadiationModel::RHOICE;

    if (settings1->UNBIASED_SELECTION){
      probability_tmp *= 1.-exp(-1.*(r_enterice.Distance(nuexitice)/L_ice)); // probability it interacts in ice along its path
    }
    double L = 0;

    double ddensity = AskaryanRadiationModel::RHOAIR;
    nearthlayers=1;

    if (where.Dot(nchord)>0.)  { // look at direction of neutrino where it enters the earth.
      std::cout << "This one's trouble.  Neutrino exit point looks more like an entrance point.  Event is " << inu << "\n";
      std::cout << "where is " << where[0] << " " << where[1] << " " << where[2] << "\n";
      std::cout << "nchord is " << nchord[0] << " " << nchord[1] << " " << nchord[2] << "\n";
      std::cout << "dot product is " << where.Dot(nchord)/sqrt(where.Dot(where)) << "\n";
      std::cout << "posnu is " << posnu[0] << " " << posnu[1] << " " << posnu[2] << "\n";
      std::cout << "Length of chord is : "<<chord<<std::endl;
    } //end if

    double altitude = where.Altitude(); //Mag()-Geoid(lat); // what is the altitude of the entrance point

    if(altitude>surface_elevation+0.1){ // if it is above the surface, it's messed up
      std::cout << "neutrino entrance point is above the surface.  Event is " << inu << "\n";
    }
    //cout <<"altitude is "<<altitude<<"\n";

    while(altitude>MIN_ALTITUDE_CRUST && x<posnu.Distance(earth_in)) {
      // starting at earth entrance point, step toward interaction position until you've reached the interaction or you are below the crust.

      ddensity = this->Density(where,crust_entered);
      
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
      altitude=where.Altitude(); //Mag()-Geoid(lat); //what is the altitude
      surface_elevation = this->SurfaceAboveGeoid(where); // altitude of surface relative to geoid at earth entrance point
      // local_icethickness = this->IceThickness(lon,lat);
      // local_waterdepth = WaterDepth(lon,lat);

    } //end while

    if (x>posnu.Distance(earth_in) && weightabsorption) {
      // if you left the loop because you have already stepped the whole distance from the entrance point to the neutrino interaction position
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
    surface_elevation = this->SurfaceAboveGeoid(where); // altitude of surface relative to geoid at earth entrance point
    // local_icethickness = this->IceThickness(lon,lat);
    // local_waterdepth = WaterDepth(lon,lat);

    double distance_remaining=where.Distance(posnu); // how much farther you need to travel before you reach the neutrino interaction point
	
    x=0; // this keeps track of how far you've stepped along the neutrino path, starting at the crust entrance.
    while(x<=distance_remaining) { // keep going until you have reached the interaction position

      ddensity=this->Density(where,crust_entered);

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
      surface_elevation = this->SurfaceAboveGeoid(where); // altitude of surface relative to geoid at earth entrance point
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





void icemc::Crust2::ReadCrust(const std::string& fName) {

  // reads in altitudes of 7 layers of crust, ice and water

  std::fstream infile(fName.c_str(),std::ios::in);
  std::string thisline; // for reading in file
  while(!infile.eof()) {
    getline(infile,thisline,'\n');

    Geoid::Position ellipsoidPos;
    Geoid::Position surfacePos;    
    int indexlon = 0;
    int indexlat = 0;
    double elevation = 0;
    
    if (thisline.find("type, latitude, longitude,")!=(int)(std::string::npos)){
      // a new point in the model  
      int beginindex=thisline.find_first_not_of(" ",57);
      int endindex=thisline.find_first_of(" ",61);
      
      std::string slat=thisline.substr(beginindex,endindex-beginindex);
      double latitude = (double)atof(slat.c_str());

      beginindex=thisline.find_first_not_of(" ",68);
      endindex=thisline.find_first_of(" ",72);

      std::string slon=thisline.substr(beginindex,endindex-beginindex);
      double longitude =(double)atof(slon.c_str());

      indexlon=(int)((longitude+180)/2);
      indexlat=(int)((90+latitude)/2);

      beginindex=thisline.find_first_not_of(" ",78);
      endindex=thisline.find_first_of(" ",83);

      std::string selev=thisline.substr(beginindex,endindex-beginindex);
      elevation = (double)atof(selev.c_str());

      ellipsoidPos.SetLonLatAlt(longitude, latitude, 0);
      fSurfaceAboveGeoid.addPoint(ellipsoidPos, elevation);
      surfacePos.SetLonLatAlt(longitude, latitude, elevation);
    } // new data point

    // skip next 4 lines of file
    for (int i=0;i<4;i++) {
      getline(infile,thisline,'\n');
    }

    // read crust values at this point
    double depthBelowSurfaceMeters = 0;
    
    for (int i=0;i<7;i++) {
      getline(infile,thisline,'\n');

      int endindex = thisline.length()-1;
      int beginindex = thisline.find_last_of("0123456789",1000);
      std::string layertype = thisline.substr(beginindex+3,endindex-beginindex);
      Layer layer = getLayerFromString(layertype);

      beginindex = thisline.find_first_not_of(" ",0);
      endindex = thisline.find_first_of(" ",beginindex);

      std::string sdepth = thisline.substr(beginindex, endindex-beginindex-1);
      
      for(int i=0; i < 3;  i++){ // skip the next three data numbers
	beginindex = thisline.find_first_not_of(" ",endindex);
	endindex = thisline.find_first_of(" ",beginindex);
      }

      std::string sdensity=thisline.substr(beginindex,endindex-beginindex);
      double density=(double)atof(sdensity.c_str());
      const double kilometersToMeters = 1e3;      
      double depthMeters = kilometersToMeters*(double)atof(sdepth.c_str());;
      
      // region where Ross Ice Shelf was not accounted for in Crust 2.0 add it in by hand
      if (layer==Layer::Ice && indexlat==5 && (indexlon<=5 || indexlon>=176)){ // Ross Ice Shelf
	depthMeters = kilometersToMeters*0.5;
      }

      // Sum up ice volume
      if(layer==Layer::Ice && ellipsoidPos.Latitude() <= COASTLINE){
	//Ice within the continent

	double areaOnGeoid = Geoid::getAreaOnGeoid(ellipsoidPos.Longitude(), ellipsoidPos.Latitude(), DLAT, DLON);
	totalIceVolume += depthMeters*areaOnGeoid;
	
	if (depthMeters > 0){
	  totalIceArea += areaOnGeoid;
	}

      }
      
      auto find_or_make_mesh = [&](std::map<Crust2::Layer, Mesh>& m, Crust2::Layer l, const std::string& n) -> Mesh& {
				 auto it = m.find(l);
				 if(it==m.end()){
				   auto& newMesh = m[l];
				   const std::string& name = fLayerNames.find(l)->second;
				   const std::string title = name + n;
				   newMesh.SetName(name.c_str());
				   newMesh.SetTitle(title.c_str());
				   return newMesh;
				 }
				 return it->second;
			       };

      icemc::Mesh& thicknessMesh = find_or_make_mesh(fThicknesses, layer,  "thickness"); //fThicknesses[layer];
      thicknessMesh.addPoint(ellipsoidPos, depthMeters);
      icemc::Mesh& surfaceMesh = find_or_make_mesh(fSurfaceMag, layer, "surface magnitude");
      
      // don't add water to cumulative depth, as it is ABOVE surface!
      if(layer==Layer::Water){
	surfaceMesh.addPoint(ellipsoidPos, surfacePos.Mag() + depthMeters);
      }
      else {
	surfaceMesh.addPoint(ellipsoidPos, surfacePos.Mag() - depthBelowSurfaceMeters);
	depthBelowSurfaceMeters += depthMeters;
	if(layer==Layer::LowerCrust){
	  icemc::Mesh& mantleMesh = find_or_make_mesh(fSurfaceMag, Layer::Mantle, "surface magnitude");
	  mantleMesh.addPoint(ellipsoidPos, surfacePos.Mag() - depthBelowSurfaceMeters);
	}
      }
      
      icemc::Mesh& densityMesh = find_or_make_mesh(fDensities, layer,  "density"); ////fDensities[layer];
      densityMesh.addPoint(ellipsoidPos, density);      
    } // for each data point in the crust file


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

  // find the place where the crust is the deepest.
  // for finding where to start stepping in Getchord
  MIN_ALTITUDE_CRUST = DBL_MAX;
  
  radii[1]=(Geoid::GEOID_MIN+MIN_ALTITUDE_CRUST)*(Geoid::GEOID_MIN+MIN_ALTITUDE_CRUST);
  
  auto build_all_meshes = [&](std::map<Crust2::Layer, Mesh>& m){
			    for(auto it = m.begin(); it != m.end(); ++it){
			      it->second.build();
			    }
			  };
  build_all_meshes(fThicknesses);
  build_all_meshes(fSurfaceMag);
  build_all_meshes(fDensities);
  fSurfaceAboveGeoid.build();
}//ReadCrust



Geoid::Position icemc::Crust2::WhereDoesItEnter(const Geoid::Position &posnu,const TVector3 &nnu) const {
  // now get neutrino entry point...
  double p = posnu.Mag(); // radius of interaction
  double costheta = (nnu.Dot(posnu)) / p; // theta of neutrino at interaction position
  double sintheta = sqrt(1-costheta*costheta);
    
  double lon = posnu.Longitude();
  double lat = posnu.Latitude();

  double a=0; // length of chord

  double R = Surface(posnu);
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
    // std::cout << "geoid, surface, p, surface-p are " << Geoid::getGeoidRadiusAtLatitude(lat) << " " << Surface(lon,lat) << " " << p << " , "<<(Surface(lon,lat)-p)<<"\n";
	
  } //else if: error: interaction takes place above the surface


    // first approx
    // find where nnu intersects spherical earth surface
  Geoid::Position r_in = posnu - a*nnu;

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
