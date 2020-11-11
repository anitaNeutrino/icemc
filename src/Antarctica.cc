#include "TVector3.h"
#include "TChain.h"
#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "Crust.h"
#include "Antarctica.h"

#include "AskaryanRadiationModel.h"
#include "Geoid.h"
#include "ConnollyEtAl2011.h"
#include "RayTracer.h"
#include "Constants.h"
#include "EnvironmentVariable.h"
#include "InteractionGenerator.h"


#include <fstream>
#include <iostream>

const std::string ICEMC_SRC_DIR = icemc::EnvironmentVariable::ICEMC_SRC_DIR();
const std::string ICEMC_DATA_DIR = ICEMC_SRC_DIR+"/data/";


icemc::Antarctica::Antarctica(int model,
			      int earth_model,
			      int WEIGHTABSORPTION_SETTING)
  : Crust(WEIGHTABSORPTION_SETTING, model)
{
  // double bedmap_R = scale_factor*bedmap_c_0 * pow(( (1 + eccentricity*sin(71*constants::RADDEG)) / (1 - eccentricity*sin(71*constants::RADDEG)) ),eccentricity/2) * tan((constants::PI/4) - (71*constants::RADDEG)/2); //varies with latitude, defined here for 71 deg S latitude
  
  // bedmap_nu = bedmap_R / cos(71*constants::RADDEG);
    
  //Parameters of the BEDMAP ice model. (See http://www.antarctica.ac.uk/aedc/bedmap/download/)
  // nCols_ice         = 1200; //number of columns in data, set by header file (should be 1200)
  // nRows_ice         = 1000; //number of rows in data, set by header file (should be 1000)
  // cellSize          = 5000; //in meters, set by header file (should be 5000) - same for both files
  // xLowerLeft_ice    = -3000000;
  // yLowerLeft_ice    = -2500000;
  // nCols_ground      = 1068;
  // nRows_ground      = 869;
  // xLowerLeft_ground = -2661600;
  // yLowerLeft_ground = -2149967;
  // nCols_water       = 1200;
  // nRows_water       = 1000;
  // xLowerLeft_water  = -3000000;
  // yLowerLeft_water  = -2500000;
  // NODATA            = -9999;

  ice_model=model;
  DEPTH_DEPENDENT_N = (int) (model / 10);
  ice_model -= DEPTH_DEPENDENT_N * 10;

  if (ice_model != 0 && ice_model != 1 && ice_model !=2) {
    icemc::report() << severity::error << "Unknown ice model requested!  Defaulting to Crust 2.0." << std::endl;
    ice_model = 2;
  }
  else if (ice_model==0) {
    // Set up layer names
    fLayerNames.emplace(Layer::Air,          std::string("air"));
    fLayerNames.emplace(Layer::Ice,          std::string("ice"));
    fLayerNames.emplace(Layer::Water,        std::string("water")); // Water only below ice in Bedmap
    fLayerNames.emplace(Layer::Bed,          std::string("bed"));
    ReadBedmap();
  }

  //read in attenuation length data for direct signals
  int i=0;
  std::ifstream sheetup((ICEMC_DATA_DIR+"/icesheet_attenlength_up.txt").c_str());
  if(sheetup.fail())
    {
      std::cerr << "Failed to open icesheet_attenlength_up.txt" << std::endl;
      exit(1);
    }
    
  i=0;
  while(sheetup>>d_sheetup[i]>>l_sheetup[i]){
    i++;
  }
  sheetup.close();
    
  std::ifstream shelfup((ICEMC_DATA_DIR+"/iceshelf_attenlength_up.txt").c_str());
  if(shelfup.fail()){
    std::cerr << "Failed to open iceshelf_attenlength_up.txt" << std::endl;
    exit(1);
  }
    
  i=0;
  while(shelfup>>d_shelfup[i]>>l_shelfup[i]){
    i++;
  }
  shelfup.close();
    
  std::ifstream westlandup((ICEMC_DATA_DIR+"/westland_attenlength_up.txt").c_str());
    
  if(westlandup.fail())
    {std::cerr << "Failed to open westland_attenlength_up.txt";
      exit(1);
    }
  i=0;
  while(westlandup>>d_westlandup[i]>>l_westlandup[i])
    {
      i++;
    }
  westlandup.close();

    
  //read in attenuation length for downgoing signals
  std::ifstream sheetdown((ICEMC_DATA_DIR+"/icesheet_attenlength_down.txt").c_str());
  if(sheetdown.fail())
    {
      std::cerr << "Failed to open icesheet_attenlength_down.txt" << std::endl;
      exit(1);
    }
    
  i=0;
  while(sheetdown>>d_sheetdown[i]>>l_sheetdown[i])
    {
      i++;
    }
  sheetdown.close();
    
  std::ifstream shelfdown((ICEMC_DATA_DIR+"/iceshelf_attenlength_down.txt").c_str());
  if(shelfdown.fail())
    {
      std::cerr << "Failed to open iceshelf_attenlength_down.txt" << std::endl;
      exit(1);
    }
    
  i=0;
  while(shelfdown>>d_shelfdown[i]>>l_shelfdown[i])
    {
      i++;
    }
  shelfdown.close();
    
  std::ifstream westlanddown((ICEMC_DATA_DIR+"/westland_attenlength_down.txt").c_str());
  if(westlanddown.fail())
    {std::cerr << "Failed to open westland_attenlength_down.txt";
      exit(1);
    }
  i=0;
  while(westlanddown>>d_westlanddown[i]>>l_westlanddown[i])
    {
      i++;
    }
  westlanddown.close();

}
//constructor Antarctica(int model)





int icemc::Antarctica::WhereDoesItEnterIce(const Geoid::Position &posnu,
					 const TVector3 &nnu,
					 double stepsize,
					 Geoid::Position &r_enterice) const {
  // now get exit point...
  //   see my geometry notes.
  // parameterize the neutrino trajectory and just see where it
  // crosses the earth radius.
    
  //  Geoid::Position r_enterice;
  double distance=0;
  int left_edge=0;
  Geoid::Position x = posnu;
  double x2;
    
  Geoid::Position x_previous = posnu;
    
  double x_previous2= x_previous.Dot(x_previous);
  x2=x_previous2;
    
  double lon = x.Longitude(),lat = x.Latitude();
  double lon_old = lon,lat_old = lat;
  double local_surface = Surface(lon,lat);
  double rock_previous2= pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
  double surface_previous2=pow(local_surface,2);
    
  double rock2=rock_previous2;
  double surface2=surface_previous2;
  int foundit=0;  // keeps track of whether you found an ice entrance point
    
  //  std::cout << "lon, lat are " << posnu.Longitude() << " " << posnu.Latitude() << "\n";
  //std::cout << "x2 at start is " << x2 << "\n";
  while (distance<2*local_surface+1000) {
	
    distance+=stepsize;
	
    x -= stepsize*nnu;
    x2=x.Dot(x);
    //std::cout << "x2 is " << x2 << "\n";
    lon = x.Longitude();
    lat = x.Latitude();
	
    double ice_thickness=IceThickness(lon,lat);
    if (lon!=lon_old || lat!=lat_old) {
      local_surface = Surface(lon,lat);
	    
      //if (lat>COASTLINE) 
      //left_edge=1;
	    
      rock2=pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
      surface2=pow(local_surface,2);    
	    
      if (ice_model==2) {
	if ((int)(lat)==COASTLINE && rock_previous2 < x2 && surface2 > x2)
	  left_edge=1;
      } //if (Crust 2.0)
    } //if (neutrino has stepped into new lon/lat bin)
	
    if ((((x_previous2>rock_previous2 && x2<rock2) // crosses rock boundary from above
	  || (x_previous2<surface_previous2 && x2>surface2)) && ice_thickness>0 && lat<COASTLINE) // crosses surface boundary from below
	|| left_edge) {
      //  std::cout << "lat, COASTLINE, left_edge is " << lat << " " << COASTLINE<< " " << left_edge << "\n";
      //std::cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
      r_enterice = x;
      // this gets you out of the loop.
      //continue;
      distance=3*Geoid::getGeoidRadiusAtLatitude(lat);
      foundit=1;
      //std::cout << "foundit is " << foundit << "\n";
      //std::cout << "r_enterice is ";r_enterice.Print();
      //continue;
    } //if
	
    x_previous = x;
    x_previous2 = x2;
    //std::cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";
	
    if (lon!=lon_old || lat!=lat_old) {
      rock_previous2 = rock2;
      surface_previous2 = surface2;
      lat_old = lat;
      lon_old = lon;
    } //if
	
  } //while
    
  return foundit;
}//WhereDoesItEnterIce

int icemc::Antarctica::WhereDoesItExitIce(const Geoid::Position &posnu,
					const TVector3 &nnu,
					double stepsize,
					Geoid::Position &r_enterice) const {
  // now get exit point...
  //   see my geometry notes.
  // parameterize the neutrino trajectory and just see where it
  // crosses the earth radius.
    
  //  Geoid::Position r_enterice;
  double distance=0;
  int left_edge=0;
  Geoid::Position x = posnu;
  double x2;   
    
  Geoid::Position x_previous = posnu;
    
  double x_previous2= x_previous.Dot(x_previous);
  x2=x_previous2;
    
  double lon = x.Longitude(),lat = x.Latitude();
  double lon_old = lon,lat_old = lat;
  double local_surface = Surface(lon,lat);
  double rock_previous2= pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
  double surface_previous2=pow(local_surface,2);
    
  double rock2=rock_previous2;
  double surface2=surface_previous2;
  int foundit=0;  // keeps track of whether you found an ice entrance point
    
    
    
  //  std::cout << "lon, lat are " << posnu.Longitude() << " " << posnu.Latitude() << "\n";
  //std::cout << "x2 at start is " << x2 << "\n";
  int nsteps=0;
  while (distance<2*local_surface+1000) {
    //std::cout << "another step.\n";
    distance+=stepsize;
    nsteps++;
    //    std::cout << "inu, nsteps is " << inu << " " << nsteps << "\n";
    x -= stepsize*nnu;
    x2=x.Dot(x);
    //std::cout << "x2 is " << x2 << "\n";
    lon = x.Longitude();
    lat = x.Latitude();
	
    double ice_thickness=IceThickness(lon,lat);
    if (lon!=lon_old || lat!=lat_old) {
      local_surface = Surface(lon,lat);
	    
      //if (lat>COASTLINE) 
      //left_edge=1;
	    
      rock2=pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
      surface2=pow(local_surface,2);    
	    
      if (ice_model==2) {
	if ((int)(lat)==COASTLINE && rock_previous2 < x2 && surface2 > x2)
	  left_edge=1;
      } //if (Crust 2.0)
    } //if (neutrino has stepped into new lon/lat bin)
	
    if ((((x_previous2<rock_previous2 && x2>rock2) // crosses rock boundary from above
	  || (x_previous2>surface_previous2 && x2<surface2)) && ice_thickness>0 && lat<COASTLINE) // crosses surface boundary from above
	|| left_edge) {
      //  std::cout << "lat, COASTLINE, left_edge is " << lat << " " << COASTLINE<< " " << left_edge << "\n";
      //std::cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
      r_enterice = x;
      // this gets you out of the loop.
      //continue;
      distance=3*Geoid::getGeoidRadiusAtLatitude(lat);//Geoid(lat);
      foundit=1;
      //std::cout << "foundit is " << foundit << "\n";
      //continue;
    } //if
	
    x_previous = x;
    x_previous2 = x2;
    //std::cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";
	
    if (lon!=lon_old || lat!=lat_old) {
      rock_previous2 = rock2;
      surface_previous2 = surface2;
      lat_old = lat;
      lon_old = lon;
    } //if
	
  } //while
    
  return foundit;
}//WhereDoesItExitIce

int icemc::Antarctica::WhereDoesItExitIceForward(const Geoid::Position &posnu,
					       const TVector3 &nnu,
					       double stepsize,
					       Geoid::Position &r_enterice) const {
  // now get exit point...
  //   see my geometry notes.
  // parameterize the neutrino trajectory and just see where it
  // crosses the earth radius.
    
  //  Geoid::Position r_enterice;
  double distance=0;
  int left_edge=0;
  Geoid::Position x = posnu;
  double x2;
    
  Geoid::Position x_previous = posnu;
    
  double x_previous2= x_previous.Dot(x_previous);
  x2=x_previous2;
    
  double lon = x.Longitude(),lat = x.Latitude();
  double lon_old = lon,lat_old = lat;
  double local_surface = Surface(lon,lat);
  double rock_previous2= pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
  double surface_previous2=pow(local_surface,2);
    
  double rock2=rock_previous2;
  double surface2=surface_previous2;
  int foundit=0;  // keeps track of whether you found an ice entrance point
    
  //  std::cout << "lon, lat are " << posnu.Longitude() << " " << posnu.Latitude() << "\n";
  //std::cout << "x2 at start is " << x2 << "\n";
  while (distance<2*local_surface+1000) {
	
    distance+=stepsize;
	
    x += stepsize*nnu;
    x2=x.Dot(x);
    //std::cout << "x2 is " << x2 << "\n";
    lon = x.Longitude();
    lat = x.Latitude();
	
    double ice_thickness=IceThickness(lon,lat);
    if (lon!=lon_old || lat!=lat_old) {
      local_surface = Surface(lon,lat);
	    
      //if (lat>COASTLINE) 
      //left_edge=1;
	    
      rock2=pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
      surface2=pow(local_surface,2);    
	    
      if (ice_model==2) {
	if ((int)(lat)==COASTLINE && rock_previous2 < x2 && surface2 > x2)
	  left_edge=1;
      } //if (Crust 2.0)
    } //if (neutrino has stepped into new lon/lat bin)
	
    if ((((x_previous2<rock_previous2 && x2>rock2) // enters rock boundary from above
	  || (x_previous2>surface_previous2 && x2<surface2)) && ice_thickness>0 && lat<COASTLINE) // crosses surface boundary from below
	|| left_edge) {
      //  std::cout << "lat, COASTLINE, left_edge is " << lat << " " << COASTLINE<< " " << left_edge << "\n";
      //std::cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
      r_enterice = x;
      // this gets you out of the loop.
      //continue;
      distance=3*Geoid::getGeoidRadiusAtLatitude(lat);//Geoid(lat);
      foundit=1;
      //std::cout << "foundit is " << foundit << "\n";
      //continue;
    } //if
	
    x_previous = x;
    x_previous2 = x2;
    //std::cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";
	
    if (lon!=lon_old || lat!=lat_old) {
      rock_previous2 = rock2;
      surface_previous2 = surface2;
      lat_old = lat;
      lon_old = lon;
    } //if
	
  } //while
    
  return foundit;
}//WhereDoesItExitIceForward




double icemc::Antarctica::IceThickness(double lon, double lat) const {
  //This method returns the thickness of the ice in meters at a location specified by a latitude and longitude (in degrees).  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  double ice_thickness=0;
  Geoid::Position p;
  p.SetLonLatAlt(lon, lat, 0);
  
  return Crust::IceThickness(p);
} //method IceThickness

double icemc::Antarctica::IceThickness(const Geoid::Position &pos) const {
  //This method returns the thickness of the ice in meters at a location under a given position vector.  Code by Stephen Hoover.
  return IceThickness(pos.Longitude(),pos.Latitude());
} //method IceThickness(position)

double icemc::Antarctica::SurfaceAboveGeoid(double lon, double lat) const {
  //This method returns the elevation above the geoid of the surface of the ice (or bare ground, if no ice is present) in meters, at a location specified by a latitude and longitude (in degrees).  In areas covered by water where no ice present, the method returns 0.  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  // lon must be 0 to 360
  double surface=0;

  Geoid::Position p;
  p.SetLonLatAlt(lon, lat, 0);

  return Crust::SurfaceAboveGeoid(p);
} //method SurfaceAboveGeoid

double icemc::Antarctica::SurfaceAboveGeoid(const Geoid::Position &pos) const {
  //This method returns the elevation above the geoid of the surface of the ice (or bare ground, if no ice is present) in meters, at a location specified by a position vector.  Code by Stephen Hoover.
    
  return SurfaceAboveGeoid(pos.Longitude(),pos.Latitude());
} //method SurfaceAboveGeoid(position)

double icemc::Antarctica::Surface(double lon,double lat) const {
  Geoid::Position p(lon, lat, SurfaceAboveGeoid(lon,lat));
  return p.Mag();
} //Surface

// double icemc::Antarctica::Surface(const Geoid::Position&pos) const {
//   return Surface(pos.Longitude(),pos.Latitude());
// } //Surface

double icemc::Antarctica::WaterDepth(double lon, double lat) const {
  //This method returns the depth of water beneath ice shelves in meters, at a location specified by a latitude and longitude (in degrees).  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  double water_depth_value=0;

  Geoid::Position p;
  p.SetLonLatAlt(lon, lat, 0);
 
  return Crust::WaterDepth(p);
} //method WaterDepth(longitude, latitude)
double icemc::Antarctica::WaterDepth(const Geoid::Position &pos) const {
  //This method returns the depth of water beneath ice shelves in meters, at a location specified by a position vector.  Code by Stephen Hoover.
    
  return WaterDepth(pos.Longitude(),pos.Latitude());
} //method WaterDepth(position)

int icemc::Antarctica::IceOnWater(const Geoid::Position &pos) const{
  if(IceThickness(pos)>0.&&WaterDepth(pos)>0.)
    return 1;
  else return 0;
    
}
int icemc::Antarctica::RossIceShelf(const Geoid::Position &pos) const {

  ///@todo RESTORE ME!
  return 0;
  // int ilon,ilat;
    
  // GetILonILat(pos,ilon,ilat);
    
  // if ((ilat==2 && ilon>=5 && ilon<=14) ||
  //     (ilat==3 && (ilon>=168 || ilon<=14)) ||
  //     (ilat==4 && (ilon>=168 || ilon<=13)) ||
  //     (ilat==5 && (ilon>=168 || ilon<=14))){
  //   return 1;
  // }
  // else{
  //   return 0;
  // }
}//RossIceShelf

int icemc::Antarctica::RossExcept(const Geoid::Position &pos) const{
  ///@todo RESTORE ME!
  // int ilon,ilat;
  // GetILonILat(pos,ilon,ilat);
  // if(ilon<=178&&ilon>=174&&ilat>=4&&ilat<=5){
  //   return 1;
  // }
  // else {
    return 0;
  // }
}


int icemc::Antarctica::RonneIceShelf(const Geoid::Position &pos) const {
  ///@todo RESTORE ME!
  return 0;
  // int ilon,ilat;
  
  // GetILonILat(pos,ilon,ilat);
    
  // if ((ilat==4 && ilon>=52 && ilon<=74) ||
  //     (ilat==5 && ilon>=50 && ilon<=71) ||
  //     (ilat==6 && ilon>=55 && ilon<=64)){
  //   return 1;
  // }
  // else{
  //   return 0;
  // }    
}//RonneIceShelf

int icemc::Antarctica::WestLand(const Geoid::Position &pos) const {
  double lon = pos.Longitude() , lat = pos.Latitude();
    
  if((lat>=4&&lat<=26)&&((lon>=0&&lon<=180)||lon>=336))
    return 1;
  else return 0;
    
}//WestLand


double icemc::Antarctica::GetBalloonPositionWeight(int ibnpos) const {
  //  std::cout << "ibnpos, volume_inhorizon, volume are " << ibnpos << " " << volume_inhorizon[ibnpos] << " " << volume << "\n";
  ///@todo The functionality here needs to be restored somehow...
  
  // if (volume_inhorizon[ibnpos]==0) {
  //   icemc::report() << severity::error << "Volume in horizon is zero!\n";
  //   exit(1);
  // }
  return 1;
  // return volume/volume_inhorizon[ibnpos];
} //GetBalloonPositionWeight


int icemc::Antarctica::OutsideAntarctica(const Geoid::Position &pos) const {
  return (pos.Latitude() >= COASTLINE);
} //OutsideAntarctica(Position)


int icemc::Antarctica::OutsideAntarctica(double lat) const {
  return (lat >= COASTLINE);
} //OutsideAntarctica(double lat)


int icemc::Antarctica::AcceptableRfexit(const TVector3 &nsurf_rfexit,const Geoid::Position &rfexit,const TVector3 &n_exit2rx) {

  //Make sure there's actually ice where the ray leaves
  if (rfexit.Latitude()>COASTLINE || IceThickness(rfexit)<0.0001) {
    std::cout << "latitude is " << rfexit.Latitude() << " compared to COASTLINE at " << COASTLINE << "\n";
    std::cout << "ice thickness is " << IceThickness(rfexit) << "\n";
    return 0;
  } //if
    
  if (nsurf_rfexit.Dot(n_exit2rx)<0) {
    //std::cout << "dot product is " << nsurf_rfexit*n_exit2rx << "\n";
    return 0;
  } //if
    
  return 1;
} //AcceptableRfexit

double icemc::Antarctica::GetN(double altitude) const {
  // these are Peter's fit parameters
  double a1=0.463251;
  double b1=0.0140157;
  double n=0;

  if (altitude < constants::FIRNDEPTH) {
    n=AskaryanRadiationModel::NICE;
  }
  else if (altitude >= constants::FIRNDEPTH && altitude <=0 && DEPTH_DEPENDENT_N) {
    //    N_DEPTH=NFIRN-(4.6198+13.62*(altitude_int/1000.))*
    //(altitude_int/1000.);   // Besson's equation for n(z)
    n=constants::NFIRN+a1*(1.0-exp(b1*altitude));   // Peter's equation for n(z)
  }
  else if (altitude > 0){
    std::cout<<"Error!  N requested for position in air!\n";
  }
  else if (!DEPTH_DEPENDENT_N){
    n = constants::NFIRN;
  }    
  return n;
} //GetN(altitude)

double icemc::Antarctica::GetN(const Geoid::Position &pos) const {
  return GetN(pos.Mag() - Surface(pos.Longitude(),pos.Latitude()));
} //GetN(Position)



double icemc::Antarctica::EffectiveAttenuationLength(const Settings *settings1,const Geoid::Position &pos,const int &whichray) const {
  double localmaxdepth = IceThickness(pos);
  double depth = WorldModel::Surface(pos) - pos.Mag();
  //std::cout << "depth is " << depth << "\n";
  int depth_index=0;
  double attenuation_length=0.0;
    
  if(WestLand(pos) && !CONSTANTICETHICKNESS) {
    depth_index=int(depth*419.9/localmaxdepth);//use 420 m ice shelf attenuation length data as the standard, squeeze or stretch if localmaxdepth is longer or shorter than 420m.
    if(RossIceShelf(pos) || RonneIceShelf(pos)) {	  
      if(whichray==0){
	attenuation_length=l_shelfup[depth_index];
      }
      else if(whichray==1){
	attenuation_length=l_shelfdown[depth_index];
      }
      else{
	std::cerr << " wrong attenuation length " <<std::endl;
      }	    
      //for sanity check
      if((depth_index+0.5)!=d_shelfup[depth_index])
	{
	  std::cerr << "the index of the array l_iceshelfup is wrong!" << std::endl;
	  exit(1);
	}
    }
    else{ //in ice sheet of westland
      if(whichray==0){
	attenuation_length=l_westlandup[depth_index]; 
      }
      else if(whichray==1){
	attenuation_length=l_westlanddown[depth_index];
      }
      else{
	std::cerr << " wrong attenuation length " <<std::endl;
      }
    }
	
    if(settings1->MOOREBAY)//if use Moore's Bay measured data for the west land
      attenuation_length*=1.717557; //about 450 m (field attenuation length) for one whole way when assuming -3dB for the power loss at the bottom
  }
  else{ //in east antarctica or constant ice thickness

    depth_index =int(depth*(2809.9/localmaxdepth));

    if(whichray==0){
      attenuation_length =l_sheetup[depth_index];
    }
    else if(whichray==1){
      attenuation_length =l_sheetdown[depth_index];
    }
    else {
      std::cerr << " wrong attenuation length " <<std::endl;
    }
  }
    
  return attenuation_length;
} //EffectiveAttenuationLengthUp

double icemc::Antarctica::Area(double latitude) const {
  //Returns the area of one square of the BEDMAP data at a given latitude. 
  double lat_rad = (90 - latitude) * constants::RADDEG;
    
  return (pow(cellSize* ((1 + sin(71*constants::RADDEG)) / (1 + sin(lat_rad))),2));
} //method Area

void icemc::Antarctica::LonLattoEN(double lon, double lat, double xLowerLeft, double yLowerLeft, int& e_coord, int& n_coord) const {
  //takes as input a latitude and longitude (in degrees) and converts to indicies for BEDMAP matricies. Needs a location for the corner of the matrix, as not all the BEDMAP files cover the same area.  Code by Stephen Hoover.
    
  double easting=0;
  double northing=0;
    
  double lon_rad = lon * constants::RADDEG; //convert to radians, and shift origin to conventional spot
  double lat_rad = lat * constants::RADDEG;
  // double lon_rad = (lon - 180) * constants::RADDEG; //convert to radians, and shift origin to conventional spot
  // double lat_rad = (90 - lat) * constants::RADDEG;
    
  double bedmap_R = Geoid::scale_factor*Geoid::c_0 * pow(( (1 + Geoid::eccentricity*sin(lat_rad)) / (1 - Geoid::eccentricity*sin(lat_rad)) ),Geoid::eccentricity/2) * tan((constants::PI/4) - lat_rad/2);
    
  easting = bedmap_R * sin(lon_rad);
  northing = bedmap_R * cos(lon_rad);
    
  //  std::cout << "bedmap_R is " << bedmap_R << "\n";
  //std::cout << "easting, northing are " << easting << " " << northing << "\n";
    
  e_coord = (int)((easting - xLowerLeft) / cellSize);
  n_coord = (int)((-1*northing - yLowerLeft) / cellSize);
    
  return;
} //method LonLattoEN

// void icemc::Antarctica::IceLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
//   //Converts a latitude and longitude (in degrees) to indicies for BEDMAP ice thickness data.  Code by Stephen Hoover.
//   LonLattoEN(lon, lat, xLowerLeft_ice, yLowerLeft_ice, e_coord, n_coord);
// }//IceLonLattoEN
// void icemc::Antarctica::GroundLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
//   //Converts a latitude and longitude (in degrees) to indicies for BEDMAP ground elevation data.  Code by Stephen Hoover.
//   LonLattoEN(lon, lat, xLowerLeft_ground, yLowerLeft_ground, e_coord, n_coord);
// }//GroundLonLattoEN
// void icemc::Antarctica::WaterLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
//   //Converts a latitude and longitude (in degrees) to indicies for BEDMAP water depth data.  Code by Stephen Hoover.
//   LonLattoEN(lon, lat, xLowerLeft_water, yLowerLeft_water, e_coord, n_coord);
// }//WaterLonLattoEN

void icemc::Antarctica::ENtoLonLat(int e_coord, int n_coord, double xLowerLeft, double yLowerLeft, double& lon, double& lat) const {
  //Takes as input the indicies from a BEDMAP data set, and turns them into latitude and longitude coordinates.  Information on which data set (surface data, ice depth, water depth) is necessary, in the form of coordinates of a corner of the map.  Code by Stephen Hoover.
    
  double isometric_lat=0;
  double easting = xLowerLeft+(cellSize*(e_coord+0.5)); //Add offset of 0.5 to get coordinates of middle of cell instead of edges.
  double northing = -1*(yLowerLeft+(cellSize*(n_coord+0.5)));
    
  //  std::cout << "easting, northing are " << easting << " " << northing << "\n";
    
  //first set longitude
    
  if (northing!=0){
    lon = atan2(easting,northing);
  }
  else{
    lon = 90*constants::RADDEG;
  }    
  // this puts lon between -pi and pi
  if (easting > 0 && lon < 0){ //adjust sign of longitude
    lon += constants::PI;
  }
  else if (easting < 0 && lon > 0){
    lon -= constants::PI;
  }
  else if (easting == 0 && northing < 0){
    lon += constants::PI;
  }
    
  //  now find latitude
  double bedmap_R = Geoid::scale_factor*Geoid::c_0 * pow(( (1 + Geoid::eccentricity*sin(71*constants::RADDEG)) / (1 - Geoid::eccentricity*sin(71*constants::RADDEG)) ),Geoid::eccentricity/2) * tan((constants::PI/4) - (71*constants::RADDEG)/2); //varies with latitude, defined here for 71 deg S latitude
  if (easting != 0){
    bedmap_R = fabs(easting/sin(lon));
  }
  else if (easting == 0 && northing != 0){
    bedmap_R = fabs(northing);
  }
  else {
    lat = 0; //at the pole, set lat=0 degrees
    lon = lon*constants::DEGRAD; // now put lon between 180 and 180 (only at pol)
    return;
  } //else
    
  isometric_lat = (constants::PI/2) - 2*atan(bedmap_R/(Geoid::scale_factor*Geoid::c_0));
    
  lat = isometric_lat + Geoid::a_bar*sin(2*isometric_lat) + Geoid::b_bar*sin(4*isometric_lat) + Geoid::c_bar*sin(6*isometric_lat) + Geoid::d_bar*sin(8*isometric_lat);

  // Make both negative for coordinate flip -- south pole at +90, East longtitudes positive
  lon = lon * constants::DEGRAD;  //convert to degrees, shift 0 to line up with bin 0 of Crust 2.0
  lat = -lat * constants::DEGRAD; //convert to degrees, with 0 degrees at the south pole
    
  //  if (lon>160 && lon<165)
  //std::cout << "e_coord, n_coord, easting, northing, lon are " << e_coord << " " << n_coord << " " << easting << " " << northing << " " << lon << "\n";
  return;
    
} //method ENtoLonLat

// void icemc::Antarctica::IceENtoLonLat(int e, int n, double& lon, double& lat) const {
//   //Converts indicies of the BEDMAP ice thickness matrix into longitude and latitude.  Code by Stephen Hoover.
//   // std::cout << "I'm inside IceENtoLonLat.\n";
//   ENtoLonLat(e,n,xLowerLeft_ice,yLowerLeft_ice,lon,lat);
// }//IceENtoLonLat
// void icemc::Antarctica::GroundENtoLonLat(int e, int n, double& lon, double& lat) const {
//   //Converts indicies of the BEDMAP ground elevation matrix into longitude and latitude.  Code by Stephen Hoover.
//   ENtoLonLat(e,n,xLowerLeft_ground,yLowerLeft_ground,lon,lat);
// }//GroundENtoLonLat
// void icemc::Antarctica::WaterENtoLonLat(int e, int n, double& lon, double& lat) const {
//   //Converts indicies of the BEDMAP water depth matrix into longitude and latitude.  Code by Stephen Hoover.
//   ENtoLonLat(e,n,xLowerLeft_water,yLowerLeft_water,lon,lat);
// }//WaterENtoLonLat

void icemc::Antarctica::ReadBedmap() {

  Geoid::Position ellipsoidPos;
  Geoid::Position surfacePos;
  double lon = 0;
  double lat = 0;
  double elevation = 0;
  double readValue;

  auto find_or_make_mesh = [&](std::map<Crust::Layer, Mesh>& m, Crust::Layer l, const std::string& n) -> Mesh& {
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
  auto build_mesh = [&](std::map<Crust::Layer, Mesh>& m, Crust::Layer l){
    auto it = m.find(l);
    if(it==m.end()){
      std::cout << "Couldn't find mesh to build" << std::endl;
    }
    else{
      it->second.build();
    }
  };
  auto build_all_meshes = [&](std::map<Crust::Layer, Mesh>& m){
    for(auto it = m.begin(); it != m.end(); ++it){
      it->second.build();
    }
  };

  //Fill the earth not covered by Bedmap with zeros on 2 degree grid
  for (int lat=-60.5; lat<90; lat++){
    for (int lon=-179.5; lon<180; lon++){
      ellipsoidPos.SetLonLatAlt(lon, lat, 0);

      icemc::Mesh& bedSurfaceMesh = find_or_make_mesh(fSurfaceMag, Layer::Bed, "surface magnitude");
      bedSurfaceMesh.addPoint(ellipsoidPos, 0);
      
      icemc::Mesh& waterThicknessMesh = find_or_make_mesh(fThicknesses, Layer::Water,  "thickness");
      waterThicknessMesh.addPoint(ellipsoidPos, 0);
      icemc::Mesh& waterSurfaceMesh = find_or_make_mesh(fSurfaceMag, Layer::Water, "surface magnitude");
      waterSurfaceMesh.addPoint(ellipsoidPos, 0);

      icemc::Mesh& iceThicknessMesh = find_or_make_mesh(fThicknesses, Layer::Ice,  "thickness"); //fThicknesses[layer];
      iceThicknessMesh.addPoint(ellipsoidPos, 0);
      icemc::Mesh& iceSurfaceMesh = find_or_make_mesh(fSurfaceMag, Layer::Ice, "surface magnitude");
      iceSurfaceMesh.addPoint(ellipsoidPos, 0);
      fSurfaceAboveGeoid.addPoint(ellipsoidPos, 0);
    }
  }
  
  //===================== Groundbed ==================//
   //Reads the BEDMAP data on the elevation of the ground beneath the ice.  If there is water beneath the ice, the ground elevation is given the value 0.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
  std::ifstream GroundBedFile((ICEMC_DATA_DIR+"/groundbed.asc").c_str());
  if(!GroundBedFile) {
    std::cerr << "Couldn't open: ICEMC_DATA_DIR/groundbed.asc" << std::endl;
    exit(1);
  }

  std::string thisline;
  int nCols, nRows, xLowerLeft, yLowerLeft, NODATA;
  Layer layer = Layer::Bed;
  std::cout<<"Reading in BEDMAP data on elevation of ground" << std::flush;
  //Read in header lines
  //@todo make function to read this header?
  for (int i=0; i<6; i++){
    getline(GroundBedFile, thisline, '\n');
    if(thisline.find("ncols")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string snCols = thisline.substr(beginindex, endindex-beginindex);
      nCols = (int)atoi(snCols.c_str());
    }
    else if(thisline.find("nrows")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string snRows = thisline.substr(beginindex, endindex-beginindex);
      nRows = (int)atoi(snRows.c_str());
    }
    else if(thisline.find("xllcorner")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string sxll = thisline.substr(beginindex, endindex-beginindex);
      xLowerLeft = (int)atoi(sxll.c_str());
    }
    else if(thisline.find("yllcorner")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string syll = thisline.substr(beginindex, endindex-beginindex);
      yLowerLeft = (int)atoi(syll.c_str());
    }
    else if(thisline.find("cellsize")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string ssize = thisline.substr(beginindex, endindex-beginindex);
      cellSize = (int)atoi(ssize.c_str());
    }
    else if(thisline.find("NODATA_value")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string sNoData = thisline.substr(beginindex, endindex-beginindex);
      NODATA = (int)atoi(sNoData.c_str());
    }    
  }
  //std::cout << "-----Ground-----" << std::endl <<"  nCols=" << nCols << " nRows= " << nRows << std::endl;
  //std::cout << "  xLowerLeft=" << xLowerLeft << " yLowerLeft=" << yLowerLeft << " cellSize=" << cellSize << " NODATA=" << NODATA << std::endl;
    
  for(int rowNum=0;rowNum<nRows;rowNum++) {
    if( rowNum%(nRows/10)==0)
      std::cout << "...." << (int)100*rowNum/nRows << "%" << std::flush;
    for(int colNum=0;colNum<nCols;colNum++) {
      ENtoLonLat(colNum, rowNum, xLowerLeft, yLowerLeft, lon,lat);
      ellipsoidPos.SetLonLatAlt(lon, lat, 0);

      GroundBedFile >> readValue;	    
      if(readValue==NODATA)
	readValue=0; //Set elevation to 0 where we have no data.

      icemc::Mesh& surfaceMesh = find_or_make_mesh(fSurfaceMag, layer, "surface magnitude");
      surfaceMesh.addPoint(ellipsoidPos, readValue);
    }
  }
  std::cout << std::endl;
  GroundBedFile.close();
  
  //Build groundbed surface mesh so it can be called
  build_mesh(fSurfaceMag, layer);
  
  //===================== Water Depth ==================//
  //Reads BEDMAP data on the depth of water beneath the ice.  Where no water is present, the value 0 is entered.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
  std::ifstream WaterDepthFile((ICEMC_DATA_DIR+"/water.asc").c_str());
  if(!WaterDepthFile) {
    std::cerr << "Couldn't open: ICEMC_DATA_DIR/water.asc" << std::endl;
    exit(1);
  }

  layer = Layer::Water;
  std::cout<<"Reading in BEDMAP data on water depth" << std::flush;
  for (int i=0; i<6; i++){
    getline(WaterDepthFile, thisline, '\n');

    if(thisline.find("ncols")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string snCols = thisline.substr(beginindex, endindex-beginindex);
      nCols = (int)atoi(snCols.c_str());
    }
    else if(thisline.find("nrows")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string snRows = thisline.substr(beginindex, endindex-beginindex);
      nRows = (int)atoi(snRows.c_str());
    }
    else if(thisline.find("xllcorner")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string sxll = thisline.substr(beginindex, endindex-beginindex);
      xLowerLeft = (int)atoi(sxll.c_str());
    }
    else if(thisline.find("yllcorner")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string syll = thisline.substr(beginindex, endindex-beginindex);
      yLowerLeft = (int)atoi(syll.c_str());
    }
    else if(thisline.find("cellsize")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string ssize = thisline.substr(beginindex, endindex-beginindex);
      cellSize = (int)atoi(ssize.c_str());
    }
    else if(thisline.find("NODATA_value")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string sNoData = thisline.substr(beginindex, endindex-beginindex);
      NODATA = (int)atoi(sNoData.c_str());
    } 
  }
  //std::cout << "-----Water-----" << std::endl <<"  nCols=" << nCols << " nRows= " << nRows << std::endl;
  //std::cout << "  xLowerLeft=" << xLowerLeft << " yLowerLeft=" << yLowerLeft << " cellSize=" << cellSize << " NODATA=" << NODATA << std::endl;
  
  for(int rowNum=0;rowNum<nRows;rowNum++) {
    if( rowNum%(int)(nRows/10)==0)
      std::cout << "...." << (int)100*rowNum/nRows << "%" << std::flush;
    for(int colNum=0;colNum<nCols;colNum++) {
      ENtoLonLat(colNum, rowNum, xLowerLeft, yLowerLeft, lon, lat);
      ellipsoidPos.SetLonLatAlt(lon, lat, 0);

      WaterDepthFile >> readValue;
      if(readValue==NODATA)
	readValue=0; //Set depth to 0 where we have no data.

      icemc::Mesh& thicknessMesh = find_or_make_mesh(fThicknesses, layer,  "thickness"); //fThicknesses[layer];
      thicknessMesh.addPoint(ellipsoidPos, readValue);
      icemc::Mesh& surfaceMesh = find_or_make_mesh(fSurfaceMag, layer, "surface magnitude");
      surfaceMesh.addPoint(ellipsoidPos, fSurfaceMag.at(Layer::Bed).eval(ellipsoidPos)+readValue);
    }
  }
  std::cout << std::endl;
  WaterDepthFile.close();

  build_mesh(fThicknesses, layer);
  build_mesh(fSurfaceMag, layer);
  


  //===================== Ice Thickness ==================//
  //Reads the BEDMAP ice thickness data.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
    
  std::ifstream IceThicknessFile((ICEMC_DATA_DIR+"/icethic.asc").c_str());
  if(!IceThicknessFile) {
    std::cerr << "Couldn't open: $ICEMC_DATA_DIR/icethic.asc" << std::endl;
    exit(1);
  }

  layer = Layer::Ice;
  std::cout<<"Reading in BEDMAP data on ice thickness" << std::flush;
  for (int i=0; i<6; i++){
    getline(IceThicknessFile, thisline, '\n');

    if(thisline.find("ncols")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string snCols = thisline.substr(beginindex, endindex-beginindex);
      nCols = (int)atoi(snCols.c_str());
    }
    else if(thisline.find("nRows")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string snRows = thisline.substr(beginindex, endindex-beginindex);
      nRows = (int)atoi(snRows.c_str());
    }
    else if(thisline.find("xllcorner")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string sxll = thisline.substr(beginindex, endindex-beginindex);
      xLowerLeft = (int)atoi(sxll.c_str());
    }
    else if(thisline.find("yllcorner")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string syll = thisline.substr(beginindex, endindex-beginindex);
      yLowerLeft = (int)atoi(syll.c_str());
    }
    else if(thisline.find("cellsize")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string ssize = thisline.substr(beginindex, endindex-beginindex);
      cellSize = (int)atoi(ssize.c_str());
    }
    else if(thisline.find("NODATA_value")!=(int)(std::string::npos)){
      int beginindex = thisline.find_first_not_of(" ", 13);
      int endindex = thisline.find_first_of(" ", beginindex+1);
      std::string sNoData = thisline.substr(beginindex, endindex-beginindex);
      NODATA = (int)atoi(sNoData.c_str());
    }
  }
  //std::cout << "-----Ice-----" << std::endl <<"  nCols=" << nCols << " nRows= " << nRows << std::endl;
  //std::cout << "  xLowerLeft=" << xLowerLeft << " yLowerLeft=" << yLowerLeft << " cellSize=" << cellSize << " NODATA=" << NODATA << std::endl;
  
  layer = Layer::Ice;
  volume=0.;
  ice_area=0.;
  // max_icevol_perbin=0.;
  // max_icethk_perbin=0.;
  for(int rowNum=0;rowNum<nRows;rowNum++) {
    if( rowNum%(nRows/10)==0)
      std::cout << "...." << (int)100*rowNum/nRows << "%" << std::flush;
    for(int colNum=0;colNum<nCols;colNum++) {
      ENtoLonLat(colNum, rowNum, xLowerLeft, yLowerLeft, lon, lat);
      ellipsoidPos.SetLonLatAlt(lon, lat, 0);

     
      IceThicknessFile >> readValue;
      if(readValue==NODATA){
	readValue=0; //Set ice depth to 0 where we have no data.
      }

      icemc::Mesh& thicknessMesh = find_or_make_mesh(fThicknesses, layer,  "thickness"); //fThicknesses[layer];
      thicknessMesh.addPoint(ellipsoidPos, readValue);
      icemc::Mesh& surfaceMesh = find_or_make_mesh(fSurfaceMag, layer, "surface magnitude");
      surfaceMesh.addPoint(ellipsoidPos, fSurfaceMag.at(Layer::Water).eval(ellipsoidPos)+readValue);
      fSurfaceAboveGeoid.addPoint(ellipsoidPos, fSurfaceMag.at(Layer::Water).eval(ellipsoidPos)+readValue);
      totalIceVolume += readValue*cellSize*cellSize; // 5 km^2 cells
      if(readValue>0)
	totalIceArea += cellSize*cellSize;
    }
  }
  std::cout << std::endl;
  IceThicknessFile.close();

  build_mesh(fThicknesses, layer);
  build_mesh(fSurfaceMag, layer);
  fSurfaceAboveGeoid.build();
  

}// ReadBedmap
