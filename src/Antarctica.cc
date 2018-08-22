#include "TVector3.h"
#include "TChain.h"
#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "Crust2.h"
#include "Antarctica.h"

#include "AskaryanFreqsGenerator.h"
#include "Geoid.h"
#include "Primaries.h"
#include "anita.hh"
#include "RayTracer.h"
#include "balloon.hh"
#include "Constants.h"
#include "EnvironmentVariable.h"

#include <fstream>
#include <iostream>

const std::string ICEMC_SRC_DIR = icemc::EnvironmentVariable::ICEMC_SRC_DIR();
const std::string ICEMC_DATA_DIR = ICEMC_SRC_DIR+"/data/";


icemc::Antarctica::Antarctica(int model,
			      int earth_model,
			      int WEIGHTABSORPTION_SETTING)
  : Crust2(earth_model,WEIGHTABSORPTION_SETTING)
{

  // double bedmap_R = scale_factor*bedmap_c_0 * pow(( (1 + eccentricity*sin(71*constants::RADDEG)) / (1 - eccentricity*sin(71*constants::RADDEG)) ),eccentricity/2) * tan((constants::PI/4) - (71*constants::RADDEG)/2); //varies with latitude, defined here for 71 deg S latitude
  
  // bedmap_nu = bedmap_R / cos(71*constants::RADDEG);
    
  //Parameters of the BEDMAP ice model. (See http://www.antarctica.ac.uk/aedc/bedmap/download/)
  nCols_ice         = 1200; //number of columns in data, set by header file (should be 1200)
  nRows_ice         = 1000; //number of rows in data, set by header file (should be 1000)
  cellSize          = 5000; //in meters, set by header file (should be 5000) - same for both files
  xLowerLeft_ice    = -3000000;
  yLowerLeft_ice    = -2500000;
  nCols_ground      = 1068;
  nRows_ground      = 869;
  xLowerLeft_ground = -2661600;
  yLowerLeft_ground = -2149967;
  nCols_water       = 1200;
  nRows_water       = 1000;
  xLowerLeft_water  = -3000000;
  yLowerLeft_water  = -2500000;
  NODATA            = -9999;

  ice_model=model;
  DEPTH_DEPENDENT_N = (int) (model / 10);
  ice_model -= DEPTH_DEPENDENT_N * 10;

  if (ice_model != 0 && ice_model != 1) {
    icemcLog() << icemc::error << "Unknown ice model requested!  Defaulting to Crust 2.0." << std::endl;
    ice_model = 0;
  } //if
  else if (ice_model==1) {
    ReadIceThickness();
    ReadWaterDepth();
    ReadGroundBed();
  } //else if (BEDMAP)
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



Geoid::Position icemc::Antarctica::PickInteractionLocation(const Geoid::Position &balloon) const {

  
  // // random numbers used
  // double rnd1=0;
  // double rnd2=2.;
  // double rnd3=1.;
  
  // double vol_bybin=0.; // volume of ice in each bin
  // int whichbin_forcrust20=0; // choice of a bin of ice for the interaction to occur in, for Crust 2.0
  // int which_coord=0;// choice of coordinates for the interaction, for Bedmap
  // double phi=0,theta=0; // phi and theta for interaction location
  // double lon=0,lat=0; // longitude and latitude corresponding to those phi and theta
  // // int ilat=0,ilon=0; // ilat and ilon bins where interaction occurs
  // int e_coord=0; //east coordinate of interaction, from Bedmap
  // int n_coord=0; // north coordinate of interaction, from Bedmap

  // //std::cout << "in pickinteractionlocation, size of ilat_inhorizon is " << ilat_inhorizon[ibnposition].size() << "\n";
  // if (ice_model == 0) { // this is Crust 2.0
  //   //std::cout << "Inside Crust 2.0 if statement.\n";
  //   // vol_bybin is initialized to 0
  //   while(vol_bybin/maxvol_inhorizon[ibnposition]<rnd3) {
  //     rnd3=gRandom->Rndm(); // pick random numbers between 0 and 1
  //     rnd1=gRandom->Rndm();
  //     rnd2=gRandom->Rndm();

  //     whichbin_forcrust20=(int)(rnd1*(double)ilat_inhorizon[ibnposition].size());
  //     ilat=ilat_inhorizon[ibnposition][whichbin_forcrust20];

  //     ilon=ilon_inhorizon[ibnposition][whichbin_forcrust20];

  //     // vol_bybin=icethkarray[ilon][ilat]*1000.*area[ilat];
  //   } //while
  //   phi=SmearPhi(ilon, gRandom->Rndm());
  //   theta=SmearTheta(ilat, gRandom->Rndm());
  //   lon = GetLon(phi);
  //   lat = GetLat(theta);
  //   //        std::cout << "ibnposition, phi, theta, lon, lat are " << ibnposition << " " << phi << " " << theta << " " << lon << " " << lat << "\n";
  // } //end if(Crust 2.0)
  // else if (ice_model==1) { // this is Bedmap
  //   //std::cout << "Inside Bedmap if statement.\n";
  //   while(vol_bybin/maxvol_inhorizon[ibnposition]<rnd3) {
  //     rnd3=gRandom->Rndm();
  //     rnd1=gRandom->Rndm();
  //     rnd2=gRandom->Rndm();

  //     which_coord=(int)(rnd1*(double)easting_inhorizon[ibnposition].size());

  //     e_coord=easting_inhorizon[ibnposition][which_coord];
  //     n_coord=northing_inhorizon[ibnposition][which_coord];

  //     //      std::cout << "e_coord, n_coord are " << e_coord << " " << n_coord << "\n";

  //     GroundENtoLonLat(e_coord,n_coord,lon,lat); //Recall that the e / n coordinates in horizon were picked from the ground bed array.
  //     //lon+=180.;
  //     //std::cout << "lon, lat are " << lon << " " << lat << "\n";

  //     if (e_coord > 1068 || e_coord < 0 || n_coord < 0 || n_coord > 869){
  // 	std::cout<<"Problem in PickDownward: e_coord, n_coord : "<<e_coord<<" , "<<n_coord<<std::endl;
  //     }
  //     vol_bybin=IceThickness(lon,lat)*Area(lat);
  //   } //while
  //   theta = lat*constants::RADDEG;
  //   phi=LongtoPhi_0is180thMeridian(lon); // convert longitude to phi
  // } //end if(BEDMAP)
  //   //}
  //   //all routines (roughness or no) \it{should} deliver a lon, lat, theta, phi
  //   //roughness may not sometimes, should add checks or something
  // TVector3 posnu=PickPosnuForaLonLat(lon,lat,theta,phi);

  std::cerr << __PRETTY_FUNCTION__ << " Disabled during refactor!" << std::endl;

  ///@todo do this again!
  return Crust2::PickInteractionLocation(balloon);
  // Geoid::Position posNu = balloon;
  // return posNu;
} //PickInteractionLocation


int icemc::Antarctica::PickUnbiased(Interaction *interaction1) const {

  interaction1->PickAnyDirection(); // first pick the neutrino direction
  return 0;
  // double mincos=cos(COASTLINE*constants::RADDEG);
  // double maxcos=cos(0.);
  // double minphi=0.;
  // double maxphi=2.*constants::PI;
  // double thisphi,thiscos,thissin;
        
  // thisphi = gRandom->Rndm()*(maxphi-minphi) + minphi;
  // thiscos = gRandom->Rndm()*(maxcos-mincos) + mincos;
  // thissin = sqrt(1.-thiscos*thiscos);
  // Geoid::Position thisr_in;// entrance point
  // Geoid::Position thisr_enterice;
  // Geoid::Position thisr_enterice_tmp;
  // Geoid::Position thisnuexitearth;
  // Geoid::Position thisnuexitice;
  // Geoid::Position thisr_exitice;
  // interaction1->noway=0;
  // interaction1->wheredoesitleave_err=0;
  // interaction1->neverseesice=0;
  // interaction1->wheredoesitenterice_err=0;
  // interaction1->toohigh=0;
  // interaction1->toolow=0;
    
  // thisr_in.SetXYZ(G::R_EARTH*thissin*thiscos,G::R_EARTH*thissin*thissin,G::R_EARTH*thiscos);
  // // interaction1->r_in = thisr_in;

  // if (thisr_in.Dot(interaction1->nnu)>0){
  //   interaction1->nnu=-1.*interaction1->nnu;
  // }
  // // does this intersect any ice
  // //std::cout << "lat, coastline, cos are " << thisr_in.Latitude() << " " << COASTLINE << " " << cos(interaction1->nnu.Theta()) << "\n";
  // if (-thisr_in.Latitude()>COASTLINE && cos(interaction1->nnu.Theta())<0) {
  //   interaction1->noway=1;

  //   return 0; // there is no way it's going through the ice
  // }
    
  // int count1=0;
  // int count2=0;

  // if (RayTracer::WhereDoesItLeave(thisr_in, interaction1->nnu, this, thisnuexitearth)) { // where does it leave Earth
  //   // really want to find where it leaves ice
  //   // Does it leave in an ice bin
  //   if (IceThickness(thisnuexitearth) && -thisnuexitearth.Latitude()<COASTLINE) { // if this is an ice bin in the Antarctic
  //     //std::cout << "inu is " << inu << " it's in ice.\n";
  //     //std::cout << "this is an ice bin.\n";
  //     thisnuexitice = thisnuexitearth;
  //     thisr_exitice = thisnuexitearth;
  //     if (thisnuexitice.Mag()>Surface(thisnuexitice)) { // if the exit point is above the surface
  // 	if ((thisnuexitice.Mag()-Surface(thisnuexitice))/cos(interaction1->nnu.Theta())>5.E3) { 
  // 	  WhereDoesItExitIce(thisnuexitearth,interaction1->nnu,5.E3, // then back up and find it more precisely
  // 			     thisr_exitice);
  // 	  thisnuexitice=(5000.)*interaction1->nnu;
  // 	  thisnuexitice+=thisr_exitice;
  // 	  count1++;
  // 	}
  // 	if ((thisnuexitice.Mag()-Surface(thisnuexitice))/cos(interaction1->nnu.Theta())>5.E2) {
		    
  // 	  WhereDoesItExitIce(thisnuexitice,interaction1->nnu,5.E2, // then back up and find it more precisely
  // 			     thisr_exitice);
  // 	  thisnuexitice=5.E2*interaction1->nnu;
  // 	  thisnuexitice+=thisr_exitice;
  // 	  count1++;
  // 	}
  // 	if ((thisnuexitice.Mag()-Surface(thisnuexitice))/cos(interaction1->nnu.Theta())>50.) {
		    
  // 	  WhereDoesItExitIce(thisnuexitice,interaction1->nnu,50., // then back up and find it more precisely
  // 			     thisr_exitice);
  // 	  count1++;
  // 	} // end third wheredoesitexit
  // 	thisnuexitice = thisr_exitice;
  //     } // if the exit point overshoots
  //     else{
  // 	thisnuexitice = thisnuexitearth;
  //     }	    
  //     // should also correct for undershooting
  //     if (count1>10){
  // 	std::cout << "count1 is " << count1 << "\n";
  //     }
  //   } // if it's an Antarctic ice bin
  //   else { // it leaves a rock bin so back up and find where it leaves ice
  //     //std::cout << "inu is " << inu << " it's in rock.\n";
  //     if (thisr_in.Distance(thisnuexitearth)>5.E4) {
  // 	count2++;
  // 	if (WhereDoesItExitIce(thisnuexitearth,interaction1->nnu,5.E4, // then back up and find it more precisely
  // 			       thisr_exitice)) {
		    
  // 	  thisnuexitice=(5.E4)*interaction1->nnu;
  // 	  thisnuexitice+=thisr_exitice;
  // 	  //std::cout << "inu is " << inu << " I'm here 1.\n";
		    
  // 	}
  // 	else {
  // 	  interaction1->neverseesice=1;
  // 	  return 0;
  // 	}
  //     }
  //     else
  // 	thisnuexitice=thisnuexitearth;
  //     //   WhereDoesItExitIce(thisnuexit,interaction1->nnu,5.E4, // then back up and find it more precisely
  //     // 			     thisr_exitice);
  //     // 	  thisnuexit=5.E4*interaction1->nnu;
  //     // 	  thisnuexit+=thisr_exitice;
  //     if (thisr_in.Distance(thisnuexitice)>5.E3) {
		
		
  // 	if (WhereDoesItExitIce(thisnuexitice,interaction1->nnu,5.E3, // then back up and find it more precisely
  // 			       thisr_exitice)) {
  // 	  count2++;
  // 	  //interaction1->neverseesice=1;
  // 	  thisnuexitice=5.E3*interaction1->nnu;
  // 	  thisnuexitice+=thisr_exitice;
  // 	  //std::cout << "inu is " << inu << " I'm here 2\n";
  // 	  //return 0;
		    
  // 	}
  //     }
  //     if (thisr_in.Distance(thisnuexitice)>5.E2) {
		
		
  // 	if (WhereDoesItExitIce(thisnuexitice,interaction1->nnu,5.E2, // then back up and find it more precisely
  // 			       thisr_exitice)) {
  // 	  count2++;
  // 	  //interaction1->neverseesice=1;
		    
  // 	  thisnuexitice=5.E2*interaction1->nnu;
  // 	  thisnuexitice+=thisr_exitice;
  // 	  //std::cout << "inu is " << inu << " I'm here 3\n";
  // 	  //return 0;
  // 	}
		
  //     }
  //     if (thisr_in.Distance(thisnuexitice)>50.) {
		
		
  // 	if (WhereDoesItExitIce(thisnuexitice,interaction1->nnu,50., // then back up and find it more precisely
  // 			       thisr_exitice)) {
  // 	  //interaction1->neverseesice=1;
  // 	  count2++;
  // 	  //std::cout << "inu is " << inu << " I'm here 4\n";
  // 	  //return 0;
  // 	}
  //     }
  //     thisnuexitice=thisr_exitice;
  //     if (count2>10){
  // 	std::cout << "count1 is " << count2 << "\n";
  //     }
  //     //	else return 0;  // never reaches any ice or is it because our step is too big
  //   } // if the nu leaves a rock bin
  // } // end wheredoesitleave
  // else {
  //   interaction1->wheredoesitleave_err=1;
  //   return 0;
  // }
  // // end finding where it leaves ice
    
  // // 	if (thisnuexit.Mag()<Surface(thisnuexit)) { // if the exit point is below the surface
  // // 	  WhereDoesItExitIceForward(thisnuexit,interaction1->nnu,20., // then find it more finely
  // // 			     thisr_exitice);
  // // 	  thisnuexit=thisr_enterice;
  // // 	  // then back up and find it more precisely
  // // 	}
    
  // if (WhereDoesItEnterIce(thisnuexitearth,interaction1->nnu,5.E3, // first pass with sort of course binning
  // 			  thisr_enterice)) {
  //   thisr_enterice_tmp = thisr_enterice+5.E3*interaction1->nnu;

  //   if (WhereDoesItEnterIce(thisr_enterice_tmp,interaction1->nnu,20., // second pass with finer binning
  // 			    thisr_enterice)) {
  //     interaction1->pathlength_inice = thisr_enterice.Distance(thisnuexitice);
  //     interaction1->posnu = interaction1->pathlength_inice*gRandom->Rndm()*interaction1->nnu;
  //     interaction1->posnu = interaction1->posnu+thisr_enterice;
  //   }
  // }
  // else {
  //   thisr_enterice = thisr_in;
  //   interaction1->wheredoesitenterice_err = 1;
  //   return 0;
  // }
  // interaction1->nuexitice = thisnuexitice;
  // interaction1->r_enterice = thisr_enterice;
    
  // if (interaction1->posnu.Mag()-Surface(interaction1->posnu)>0) {
  //   interaction1->toohigh=1;
  //   //std::cout << "inu, toohigh is " << inu << " " << interaction1->toohigh << "\n";
  //   return 0;
  // }
  // if (interaction1->posnu.Mag()-Surface(interaction1->posnu)+IceThickness(interaction1->posnu)<0) {
  //   interaction1->toolow=1;
  //   //std::cout << "inu, toolow is " << inu << " " << interaction1->toolow << "\n";
  //   return 0;
  // }    
  // return 1;
    
}


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
	    
      if (ice_model==0) {
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
	    
      if (ice_model==0) {
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
	    
      if (ice_model==0) {
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

  if (ice_model==1) {
    int e_coord=0;
    int n_coord=0;
    IceLonLattoEN(lon,lat,e_coord,n_coord);
    std::cout << lon << "\t" << e_coord << "\t" << n_coord << std::endl;
    if (e_coord <= 1200 && e_coord >= 0 && n_coord <= 1000 && n_coord > 0){      
      ice_thickness = ice_thickness_array[e_coord][n_coord]; //if this region has BEDMAP data, use it.
    }
    else{
      ice_thickness = Crust2::IceThickness(lon, lat);
    }
  } //BEDMAP ice thickness
  else if (ice_model==0) {
    ice_thickness = Crust2::IceThickness(lon, lat);
    //std::cout << "ilon, ilat are " << (int)(lon/2) << " " << (int)(lat/2) << "\n";
  } //Crust 2.0 ice thickness
    
  return ice_thickness;
} //method IceThickness

double icemc::Antarctica::IceThickness(const Geoid::Position &pos) const {
  //This method returns the thickness of the ice in meters at a location under a given position vector.  Code by Stephen Hoover.
    
  // return IceThickness(pos.Longitude(),pos.Latitude());
  return Crust2::IceThickness(pos);
} //method IceThickness(position)

double icemc::Antarctica::SurfaceAboveGeoid(double lon, double lat) const {
  //This method returns the elevation above the geoid of the surface of the ice (or bare ground, if no ice is present) in meters, at a location specified by a latitude and longitude (in degrees).  In areas covered by water where no ice present, the method returns 0.  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  // lon must be 0 to 360
  double surface=0;

  if (ice_model==1) {
    int e_coord_ice=0;
    int n_coord_ice=0;
    int e_coord_ground=0;
    int n_coord_ground=0;
    IceLonLattoEN(lon,lat,e_coord_ice,n_coord_ice);
    GroundLonLattoEN(lon,lat,e_coord_ground,n_coord_ground);
    if (e_coord_ground <= 1068 && e_coord_ground >= 0 && n_coord_ground <= 869 && n_coord_ground >= 0 && e_coord_ice <= 1200 && e_coord_ice >= 0 && n_coord_ice <= 1000 && n_coord_ice >= 0){
      surface = ground_elevation[e_coord_ground][n_coord_ground] + ice_thickness_array[e_coord_ice][n_coord_ice] + water_depth[e_coord_ice][n_coord_ice];
    }
    else{
      surface = Crust2::SurfaceAboveGeoid(lon, lat);
    }
  } //Elevation of surface above geoid according to BEDMAP
  else if (ice_model==0) {
      surface = Crust2::SurfaceAboveGeoid(lon, lat);
  } //Elevation of surface above geoid according to Crust 2.0
    
  return surface;
} //method SurfaceAboveGeoid

double icemc::Antarctica::SurfaceAboveGeoid(const Geoid::Position &pos) const {
  //This method returns the elevation above the geoid of the surface of the ice (or bare ground, if no ice is present) in meters, at a location specified by a position vector.  Code by Stephen Hoover.
    
  return SurfaceAboveGeoid(pos.Longitude(),pos.Latitude());
} //method SurfaceAboveGeoid(position)

double icemc::Antarctica::Surface(double lon,double lat) const {
  Geoid::Position p(lon, lat, SurfaceAboveGeoid(lon,lat));
  return p.Mag();
} //Surface

double icemc::Antarctica::Surface(const Geoid::Position&pos) const {
  return Surface(pos.Longitude(),pos.Latitude());
} //Surface

double icemc::Antarctica::WaterDepth(double lon, double lat) const {
  //This method returns the depth of water beneath ice shelves in meters, at a location specified by a latitude and longitude (in degrees).  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  double water_depth_value=0;
    
  if (ice_model==0) {
    water_depth_value = Crust2::WaterDepth(lon, lat);
  } //if(Crust 2.0)
  else if (ice_model==1) {
    int e_coord=0;
    int n_coord=0;
    WaterLonLattoEN(lon,lat,e_coord,n_coord);
    if (e_coord <= 1200 && e_coord >= 0 && n_coord <= 1000 && n_coord >= 0)
      water_depth_value = water_depth[e_coord][n_coord];
    else
      water_depth_value = Crust2::WaterDepth(lon, lat);
  } //else if(BEDMAP)
    
  return water_depth_value;
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
  //   icemcLog() << icemc::error << "Volume in horizon is zero!\n";
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
    n=AskaryanFreqsGenerator::NICE;
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
  double depth = Surface(pos) - pos.Mag();
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

void icemc::Antarctica::IceLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
  //Converts a latitude and longitude (in degrees) to indicies for BEDMAP ice thickness data.  Code by Stephen Hoover.
  LonLattoEN(lon, lat, xLowerLeft_ice, yLowerLeft_ice, e_coord, n_coord);
}//IceLonLattoEN
void icemc::Antarctica::GroundLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
  //Converts a latitude and longitude (in degrees) to indicies for BEDMAP ground elevation data.  Code by Stephen Hoover.
  LonLattoEN(lon, lat, xLowerLeft_ground, yLowerLeft_ground, e_coord, n_coord);
}//GroundLonLattoEN
void icemc::Antarctica::WaterLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
  //Converts a latitude and longitude (in degrees) to indicies for BEDMAP water depth data.  Code by Stephen Hoover.
  LonLattoEN(lon, lat, xLowerLeft_water, yLowerLeft_water, e_coord, n_coord);
}//WaterLonLattoEN

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
    
  lon = lon * constants::DEGRAD;  //convert to degrees, shift 0 to line up with bin 0 of Crust 2.0
  lat = lat * constants::DEGRAD; //convert to degrees, with 0 degrees at the south pole
    
  //  if (lon>160 && lon<165)
  //std::cout << "e_coord, n_coord, easting, northing, lon are " << e_coord << " " << n_coord << " " << easting << " " << northing << " " << lon << "\n";
  return;
    
} //method ENtoLonLat

void icemc::Antarctica::IceENtoLonLat(int e, int n, double& lon, double& lat) const {
  //Converts indicies of the BEDMAP ice thickness matrix into longitude and latitude.  Code by Stephen Hoover.
  // std::cout << "I'm inside IceENtoLonLat.\n";
  ENtoLonLat(e,n,xLowerLeft_ice,yLowerLeft_ice,lon,lat);
}//IceENtoLonLat
void icemc::Antarctica::GroundENtoLonLat(int e, int n, double& lon, double& lat) const {
  //Converts indicies of the BEDMAP ground elevation matrix into longitude and latitude.  Code by Stephen Hoover.
  ENtoLonLat(e,n,xLowerLeft_ground,yLowerLeft_ground,lon,lat);
}//GroundENtoLonLat
void icemc::Antarctica::WaterENtoLonLat(int e, int n, double& lon, double& lat) const {
  //Converts indicies of the BEDMAP water depth matrix into longitude and latitude.  Code by Stephen Hoover.
  ENtoLonLat(e,n,xLowerLeft_water,yLowerLeft_water,lon,lat);
}//WaterENtoLonLat

void icemc::Antarctica::GetMAXHORIZON(Balloon *bn1) const {

  double altitude_inmeters=bn1->BN_ALTITUDE*12.*constants::CMINCH/100.;
  if (bn1->BN_ALTITUDE==0.){
    bn1->MAXHORIZON=8.E5; // if it is a standard flight then use a horizon of 800 km
  }
  else{
    bn1->MAXHORIZON=(sqrt((Geoid::R_EARTH+altitude_inmeters)*(Geoid::R_EARTH+altitude_inmeters)-Geoid::R_EARTH*Geoid::R_EARTH))*1.1;
    // find distance from hrizon to balloon, increase it by 10% to be conservative.
  }
  icemcLog() << icemc::info << "MAXHORIZON is " << bn1->MAXHORIZON << std::endl;

}


void icemc::Antarctica::CreateHorizons(const Settings *settings1, const Balloon *bn1) {

  // add up volume of ice within horizon of payload
  // goes a little beyond horizon.

  // when we select a path in a circle at 80deg S latitude,
  // vectors are binned in phi (longitude).

  // when we select the Anita-lite path,
  // vectors are binned in time.

  //for (int i=0; i<60;i++)
  //cout<<"area at lat "<<(double)i/2<<" is "<<Area((double)i/2)<<std::endl;

  volume = 0.; // initialize volume to zero


  icemcLog() << icemc::warning << "Temporarily disabled " << __PRETTY_FUNCTION__ << " during refactor! Volume = 2.3e8!" << std::endl;
  volume = 2.3e8;

  // double total_area=0; // initialize total area to zero
  // int NBALLOONPOSITIONS; // number of balloon positions considered here
  // if (bn1->whichPath()==FlightPath::AnitaLite){ // if anita-lite
  //   NBALLOONPOSITIONS=(int)((double)bn1->NPOINTS/(double)bn1->REDUCEBALLOONPOSITIONS); //only take 1/100 of the total balloon positions that we have because otherwise it's overkill.
  // }
  // else if (bn1->whichPath()==FlightPath::Anita1 ||
  // 	   bn1->whichPath()==FlightPath::Anita2 ||
  // 	   bn1->whichPath()==FlightPath::Anita3 ||
  // 	   bn1->whichPath()==FlightPath::Anita4) {

  //   NBALLOONPOSITIONS=(int)((double)bn1->fChain->GetEntries()/(double)bn1->REDUCEBALLOONPOSITIONS)+1;
  // }

  // else if (bn1->whichPath()==FlightPath::Circle80DegreesSouth){ // for picking random point along 80 deg south

  //   NBALLOONPOSITIONS=NPHI; // NPHI is the number of bins in phi for the visible ice in the horizon.  This is not the same as NLON, the number of bins in longitude for crust 2.0
  // }
  // else{ // includes fixed position (bn1->whichPath()=0)
  //   NBALLOONPOSITIONS = 1;
  // }

  // if (NBALLOONPOSITIONS>NBNPOSITIONS_MAX) {
  //   icemcLog() << icemc::error << "Number of balloon positions is greater than max allowed.\n";
  //   exit(1);
  // }

  // double phi_bn_temp=0; //phi of balloon, temporary variable
  // Geoid::Position r_bn_temp; //position of balloon
  // Geoid::Position r_bin; // position of each bin
    
  // double surface_elevation=0;
  // double local_ice_thickness=0;
  // double lat=0;
  // double lon=0;
  // //double r=0; // r,theta and phi are temporary variables
  // double theta=0;
  // double phi=0;
  // int volume_found=0;
  // char horizon_file[80];
  // FILE *bedmap_horizons = new FILE();
  // char line[200];
  // int e_coord = 0;
  // int n_coord = 0;
    
  // sprintf(horizon_file,
  // 	  "bedmap_horizons_whichpath%i_weights%i.dat",
  // 	  bn1->whichPath(),
  // 	  settings1->USEPOSITIONWEIGHTS);

  // bool writeFile = settings1->WRITE_FILE > 0;

  // if (ice_model==1 && !writeFile) { // for bedmap model, need to be able to read file
  //   if(!(bedmap_horizons = fopen(horizon_file, "r"))) {
  //     printf("Error: unable to open %s.  Creating new file.\n", horizon_file);
  //     writeFile=1;
  //   }//if
  // } //if
  // if (ice_model==1 && writeFile) { // for bedmap model, need to be able to write to file
  //   if(!(bedmap_horizons = fopen(horizon_file, "w"))) {
  //     printf("Error: unable to open %s\n", horizon_file);
  //     exit(1);
  //   }//if
  // } //if

  // //@todo does this if statement make sense?
  // if (bn1->whichPath()!=FlightPath::AnitaLite && bn1->whichPath()!=FlightPath::Anita1){ // not anita and not anita-lite
  //   lat = GetLat(theta_bn); //get index of latitude, same for all balloon positions
  // }    

  // for (int i=0;i<NBALLOONPOSITIONS;i++) { // loop over balloon positions
	
  //   maxvol_inhorizon[i]=-1.; // volume of bin with the most ice in the horizon
	
  //   if (bn1->whichPath()==FlightPath::AnitaLite) { // anita or anita-lite path
  //     theta_bn=(90+bn1->latitude_bn_anitalite.at(i*100))*constants::RADDEG; //theta of the balloon wrt north pole
  //     lat = GetLat(theta_bn); // latitude  
  //     phi_bn_temp =  (-1*bn1->longitude_bn_anitalite.at(i*100)+90.); //phi of the balloon, with 0 at +x and going counter clockwise looking down from the south pole
  //     if (phi_bn_temp<0){ //correct phi_bn if it's negative
  // 	phi_bn_temp+=360.;
  //     }
  //     phi_bn_temp*=constants::RADDEG;// turn it into radians

  //     altitude_bn=bn1->altitude_bn_anitalite.at(i*100)*12.*constants::CMINCH/100.; // get the altitude for this balloon posistion
  //     //altitude_bn=altitude_bn_anitalite[i*100]*12.*CMINCH/100.; // for anita, altitude in is meters
	    
  //   } // end if anita-lite

  //   else if (bn1->whichPath()==FlightPath::Anita1 || 
  // 	     bn1->whichPath()==FlightPath::Anita2 || 
  // 	     bn1->whichPath()==FlightPath::Anita3 || 
  // 	     bn1->whichPath()==FlightPath::Anita4) {

  //     bn1->fChain->GetEvent(i*100);

  //     theta_bn = (90+(double)bn1->flatitude)*constants::RADDEG; //theta of the balloon wrt north pole

  //     lat=GetLat(theta_bn); // latitude  
  //     phi_bn_temp=(-1*(double)bn1->flongitude+90.); //phi of the balloon, with 0 at +x and going counter clockwise looking down from the south pole
  //     if (phi_bn_temp<0){ //correct phi_bn if it's negative
  // 	phi_bn_temp+=360.;
  //     }
  //     phi_bn_temp*=constants::RADDEG;// turn it into radians

  //     altitude_bn = (double) bn1->faltitude; // for anita, altitude in is meters

  //   } // end if anita or anita-lite
  //   else if (bn1->whichPath()==FlightPath::Circle80DegreesSouth){ // for picking random position along 80 deg south
  //     phi_bn_temp=dGetPhi(i); //get phi of the balloon just based on the balloon positions.  Remember nballoonpositions runs from 0 to 179 here.
  //   // the output of dGetPhi is in radians and runs from -pi/2 to 3*pi/2 relative to 
  //   }
  //   else{ // includes bn1->whichPath()=0 
  //     phi_bn_temp=phi_bn;
  //   }
  //   // altitude_bn has already been set in SetDefaultBalloonPosition
  //   // same for theta_bn
  //   // lat set above
  //   lon = GetLon(phi_bn_temp); //get longitude for this phi position
  //   // lon goes from 0 (at -180 longitude) to 360 (at 180 longitude)

  //   // position of balloon
  //   surface_elevation = this->Surface(lon,lat); // distance from center of the earth to surface at this lon and lat
    
  //   r_bn_temp = TVector3(sin(theta_bn)*cos(phi_bn_temp)*(surface_elevation+altitude_bn),
  // 			 sin(theta_bn)*sin(phi_bn_temp)*(surface_elevation+altitude_bn),
  // 			 cos(theta_bn)*(surface_elevation+altitude_bn)); // vector from center of the earth to balloon
	
  //   if (ice_model==0) { // Crust 2.0
  //     for (int j=0;j<NLON;j++) { // loop over bins in longitude
  // 	for (int k=0;k<ILAT_COASTLINE;k++) { // loop over bins in latitude

  // 	  // get position of each bin
  // 	  r_bin.SetLonLatAlt();
  // 	  r_bin = TVector3(sin(dGetTheta(k))*cos(dGetPhi(j))*(geoid[k]+surfacer[j][k]),
  // 			   sin(dGetTheta(k))*sin(dGetPhi(j))*(geoid[k]+surfacer[j][k]),
  // 			   cos(dGetTheta(k))*(geoid[k]+surfacer[j][k])); // vector from center of the earth to the surface of this bin

  // 	  if (!volume_found){ 
  // 	    volume += icethkarray[j][k]*1000*area[k]; // add this to the total volume of ice in antarctica
  // 	  }
  // 	  if (!volume_found && icethkarray[j][k] > 0){
  // 	    total_area += area[k]; // add this to the total area of ice in antarctica
  // 	  }
  // 	  // if the bin is within the maximum horizon of the balloon or if we don't care
  // 	  //r=(geoid[k]+surfacer[j][k]);
  // 	  //phi=dGetPhi(j);
  // 	  //theta=dGetTheta(k);
  // 	  //	  cout << "USEWEIGHTS is " << USEWEIGHTS << "\n";
  // 	  if ((settings1->USEPOSITIONWEIGHTS && r_bin.Distance(r_bn_temp) < bn1->MAXHORIZON)
  // 	      || !settings1->USEPOSITIONWEIGHTS) {
  // 	    // then put this latitude and longitude in vector

  // 	    ilat_inhorizon[i].push_back(k); 
  // 	    ilon_inhorizon[i].push_back(j);
  // 	    // add up volume in horizon
			
  // 	    volume_inhorizon[i]+=icethkarray[j][k]*1000*area[k];

  // 	    // finding which bin in horizon has maximum volume
  // 	    if (icethkarray[j][k]*1000*area[k]>maxvol_inhorizon[i]) {
  // 	      maxvol_inhorizon[i]=icethkarray[j][k]*1000.*area[k];
  // 	    }
  // 	  } //end if (distance < 800 km)
  // 	} //end for (k loop)
  //     } //end for (j loop)
  //     //      cout << "i, volume_inhorizon are " << i << " " << volume_inhorizon[i] << "\n";

  //     // ifi the balloon is too close to the ice, it will think there aren't any
  //     // bins in the horizon, so force it the include the ones just below the balloon
  //     int ilat_bn,ilon_bn;
  //     GetILonILat(r_bn_temp,ilon_bn,ilat_bn); // find which longitude and latitude the balloon is above

  //     if (ilat_inhorizon[i].size()==0 || ilon_inhorizon[i].size()==0) {
  // 	maxvol_inhorizon[i]=icethkarray[ilon_bn][ilat_bn]*1000.*area[ilat_bn]; // need to give it a maximum volume for a bin in horizon 
  // 	volume_inhorizon[i]=icethkarray[ilon_bn][ilat_bn]*1000.*area[ilat_bn]; // and a total volume for the horizon
  //     }

  //     if (ilat_inhorizon[i].size()==0){ // for the ith balloon position, if it didn't find a latitude bin in horizon
  // 	ilat_inhorizon[i].push_back(ilat_bn); // force it to be the one below the balloon       
  //     }	    
  //     if (ilon_inhorizon[i].size()==0){ // for the ith balloon position, if it didn't find a longitude bin in horizon
  // 	ilon_inhorizon[i].push_back(ilon_bn); // force it to be the one below the balloon
  //     }
  //   } //end if (ice_model==0) Crust 2.0

  //   else if (ice_model==1 && !writeFile) { // for bedmap model
  //     fgets(line,200,bedmap_horizons);
  //     while(line[0] != 'X') {
  // 	e_coord = atoi(strtok(line,","));
  // 	n_coord = atoi(strtok(NULL,","));
  // 	easting_inhorizon[i].push_back(e_coord);
  // 	northing_inhorizon[i].push_back(n_coord);
  // 	fgets(line,200,bedmap_horizons);
  //     } //while there are more eastings and northings
  //     strtok(line,",");
  //     volume_inhorizon[i] = atof(strtok(NULL,",")); // volume in the horizon
  //     maxvol_inhorizon[i] = atof(strtok(NULL,","));

  //     if (!volume_found) {
  // 	total_area = atof(fgets(line,200,bedmap_horizons)); // total area on the continent
  // 	volume = atof(fgets(line,200,bedmap_horizons)); // total volume on the continent
  //     } //if
  //   } //end if (ice_model==1) && !writeFile

  //   else if (ice_model==1 && writeFile) { //for BEDMAP model, look through all easting and northing coordinates in the groundbed map (our smallest).  Output what we find to a file for later use.
  //     for (int n_coord=0;n_coord<nRows_ground;n_coord++) {
  // 	for (int e_coord=0;e_coord<nCols_ground;e_coord++) {

  // 	  GroundENtoLonLat(e_coord,n_coord,lon,lat);

  // 	  theta = lat * constants::RADDEG;
  // 	  phi=LongtoPhi_0is180thMeridian(lon);

  // 	  surface_elevation = this->Surface(lon,lat);
  // 	  local_ice_thickness = this->IceThickness(lon,lat);
		    
  // 	  // get position of each bin
  // 	  r_bin = TVector3(sin(theta)*cos(phi)*surface_elevation,
  // 			   sin(theta)*sin(phi)*surface_elevation,
  // 			   cos(theta)*surface_elevation);
		    
  // 	  if (!volume_found){
  // 	    volume += local_ice_thickness*Area(lat);
  // 	  }
  // 	  if (!volume_found && local_ice_thickness > 5){
  // 	    total_area += Area(lat);
  // 	  }
  // 	  if ((settings1->USEPOSITIONWEIGHTS && r_bn_temp.Distance(r_bin)<bn1->MAXHORIZON) || !settings1->USEPOSITIONWEIGHTS) {
  // 	    fprintf(bedmap_horizons,"%i,%i,\n",e_coord,n_coord);
  // 	    // then put this latitude and longitude in vector
  // 	    easting_inhorizon[i].push_back(e_coord);
  // 	    northing_inhorizon[i].push_back(n_coord);
  // 	    // add up volume in horizon

  // 	    GroundENtoLonLat(e_coord,n_coord,lon,lat);

  // 	    volume_inhorizon[i]+=local_ice_thickness*Area(lat);

  // 	    // finding which bin in horizon has maximum volumey
  // 	    if (local_ice_thickness*Area(lat)>maxvol_inhorizon[i]) {
  // 	      maxvol_inhorizon[i]=local_ice_thickness*Area(lat);
  // 	    }
  // 	  } //end if (distance < 800 km & ice present)
  // 	} //end for (e_coord loop)
  //     } //end for (n_coord loop)

  //     fprintf(bedmap_horizons,"X,%f,%f,\n",volume_inhorizon[i],maxvol_inhorizon[i]);
  //     if (!volume_found) {
  // 	fprintf(bedmap_horizons,"%f\n",total_area);
  // 	fprintf(bedmap_horizons,"%f\n",volume);
  //     } //if
  //   } //end if (ice_model==1) && writeFile

  //   if (!volume_found) {
  //     std::cout << "Total surface area covered with ice (in m^2) is : "<<total_area<<std::endl;
  //     volume_found=1;
  //   } //if
  // } //end loop over balloon positions
    
  //   // finding average volume in horizon over all balloon positions.
  // volume_inhorizon_average=0;
    
  // for (int i=0;i<NPHI;i++) {
  //   volume_inhorizon_average+=volume_inhorizon[i];
	
  // } //for
  // volume_inhorizon_average/=(double)NBALLOONPOSITIONS;
    
  // //cout << "Total volume of ice in Antarctica with this ice model (m^3): " << volume << "\n";
  // //cout << "Average volume of ice within a horizon is " << volume_inhorizon_average << "\n";
    
  // icemcLog() << "Average volume of ice within a horizon is " << volume_inhorizon_average << "\n";
    
  // icemcLog() << "Average thickness of ice within horizon, averages over balloon positions " << volume_inhorizon_average/constants::PI/pow(bn1->MAXHORIZON,2) << "\n";
} //method CreateHorizons 


void icemc::Antarctica::ReadIceThickness() {
  //Reads the BEDMAP ice thickness data.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
    
  std::ifstream IceThicknessFile((ICEMC_DATA_DIR+"/icethic.asc").c_str());
  if(!IceThicknessFile) {
    std::cerr << "Couldn't open: $ICEMC_DATA_DIR/icethic.asc" << std::endl;
    exit(1);
  }
    
  std::cout<<"Reading in BEDMAP data on ice thickness.\n";
  std::string tempBuf1;
  std::string tempBuf2;
  std::string tempBuf3;
  std::string tempBuf4;
  std::string tempBuf5;
  std::string tempBuf6;
  int temp1,temp2,temp3,temp4,temp5,temp6;
    
  IceThicknessFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
		   >> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
		   >> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
    
  if(tempBuf1 == std::string("ncols")) {
    nCols_ice=temp1;
  }
  if(tempBuf2 == std::string("nrows")) {
    nRows_ice=temp2;
  }
  if(tempBuf3 == std::string("xllcorner")) {
    xLowerLeft_ice=temp3;
  }
  if(tempBuf4 == std::string("yllcorner")) {
    yLowerLeft_ice=temp4;
  }
  if(tempBuf5 == std::string("cellsize")) {
    cellSize=temp5;
  }
  if(tempBuf6 == std::string("NODATA_value")) {
    NODATA=temp6;
  }
  //cout<<"nCols_ice, nRows_ice "<<nCols_ice<<" , "<<nRows_ice<<std::endl;
  //cout<<"xLL_ice, yLL_ice, cellsize "<<xLowerLeft_ice<<" , "<<yLowerLeft_ice<<" , "<<cellSize<<endl<<std::endl;
    
  double theValue;
  volume=0.;
  ice_area=0.;
  // max_icevol_perbin=0.;
  // max_icethk_perbin=0.;
  double lon_tmp,lat_tmp;
  for(int rowNum=0;rowNum<nRows_ice;rowNum++) {
    for(int colNum=0;colNum<nCols_ice;colNum++) {
      IceENtoLonLat(colNum,rowNum,lon_tmp,lat_tmp);
      IceThicknessFile >> theValue;
      if(theValue==NODATA){
	theValue=0; //Set ice depth to 0 where we have no data.
      }
      ice_thickness_array[colNum][rowNum] = double(theValue); //This stores data as ice_thickness_array[easting][northing]
      volume+=ice_thickness_array[colNum][rowNum]*Area(lat_tmp);
      if (ice_thickness_array[colNum][rowNum]>0){
	ice_area+=Area(lat_tmp);
      }
      // if (ice_thickness_array[colNum][rowNum]*Area(lat_tmp)>max_icevol_perbin){
      // 	max_icevol_perbin=ice_thickness_array[colNum][rowNum]*Area(lat_tmp);
      // }
      // if (ice_thickness_array[colNum][rowNum]>max_icethk_perbin){
      // 	max_icethk_perbin=ice_thickness_array[colNum][rowNum];
      // }
      if(theValue > fMaxIceThickness ){
	fMaxIceThickness = theValue;
      }
    }//for
  }//for
    
  IceThicknessFile.close();
  return;
} //method ReadIceThickness

void icemc::Antarctica::ReadGroundBed() {
  //Reads the BEDMAP data on the elevation of the ground beneath the ice.  If there is water beneath the ice, the ground elevation is given the value 0.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
  std::ifstream GroundBedFile((ICEMC_DATA_DIR+"/groundbed.asc").c_str());
  if(!GroundBedFile) {
    std::cerr << "Couldn't open: ICEMC_DATA_DIR/groundbed.asc" << std::endl;
    exit(1);
  }
    
  std::cout<<"Reading in BEDMAP data on elevation of ground.\n";
    
  std::string tempBuf1;
  std::string tempBuf2;
  std::string tempBuf3;
  std::string tempBuf4;
  std::string tempBuf5;
  std::string tempBuf6;
  int temp1,temp2,temp3,temp4,temp5,temp6;
    
  GroundBedFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
		>> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
		>> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
    
  if(tempBuf1 == std::string("ncols")) {
    nCols_ground=temp1;
  }
  if(tempBuf2 == std::string("nrows")) {
    nRows_ground=temp2;
  }
  if(tempBuf3 == std::string("xllcorner")) {
    xLowerLeft_ground=temp3;
  }
  if(tempBuf4 == std::string("yllcorner")) {
    yLowerLeft_ground=temp4;
  }
  if(tempBuf5 == std::string("cellsize")) {
    cellSize=temp5;
  }
  if(tempBuf6 == std::string("NODATA_value")) {
    NODATA=temp6;
  }
    
  //cout<<"nCols_ground, nRows_ground "<<nCols_ground<<" , "<<nRows_ground<<std::endl;
  //cout<<"xLL_ground, yLL_ground, cellsize "<<xLowerLeft_ground<<" , "<<yLowerLeft_ground<<" , "<<cellSize<<endl<<std::endl;
    
  double theValue;
  for(int rowNum=0;rowNum<nRows_ground;rowNum++) {
    for(int colNum=0;colNum<nCols_ground;colNum++) {
      GroundBedFile >> theValue;
	    
      if(theValue==NODATA)
	theValue=0; //Set elevation to 0 where we have no data.
      ground_elevation[colNum][rowNum] = double(theValue);
      //if (theValue != -96 && theValue != 0)
      //cout<<"ground_elevation: "<<theValue<<std::endl;
    }//for
  }//for
    
  GroundBedFile.close();
  return;
} //method ReadGroundBed

void icemc::Antarctica::ReadWaterDepth() {
  //Reads BEDMAP data on the depth of water beneath the ice.  Where no water is present, the value 0 is entered.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
  std::ifstream WaterDepthFile((ICEMC_DATA_DIR+"/water.asc").c_str());
  if(!WaterDepthFile) {
    std::cerr << "Couldn't open: ICEMC_DATA_DIR/water.asc" << std::endl;
    exit(1);
  }
    
  std::cout<<"Reading in BEDMAP data on water depth.\n";
    
  std::string tempBuf1;
  std::string tempBuf2;
  std::string tempBuf3;
  std::string tempBuf4;
  std::string tempBuf5;
  std::string tempBuf6;
  int temp1,temp2,temp3,temp4,temp5,temp6;
    
  WaterDepthFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
		 >> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
		 >> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
    
  if(tempBuf1 == std::string("ncols")) {
    nCols_water=temp1;
  }
  if(tempBuf2 == std::string("nrows")) {
    nRows_water=temp2;
  }
  if(tempBuf3 == std::string("xllcorner")) {
    xLowerLeft_water=temp3;
  }
  if(tempBuf4 == std::string("yllcorner")) {
    yLowerLeft_water=temp4;
  }
  if(tempBuf5 == std::string("cellsize")) {
    cellSize=temp5;
  }
  if(tempBuf6 == std::string("NODATA_value")) {
    NODATA=temp6;
  }
    
  //cout<<"nCols_water, nRows_water "<<nCols_water<<" , "<<nRows_water<<std::endl;
  //cout<<"xLL_water, yLL_water, cellsize "<<xLowerLeft_water<<" , "<<yLowerLeft_water<<" , "<<cellSize<<endl<<std::endl;
    
  double theValue;
  for(int rowNum=0;rowNum<nRows_water;rowNum++) {
    for(int colNum=0;colNum<nCols_water;colNum++) {
      WaterDepthFile >> theValue;
	    
      if(theValue==NODATA)
	theValue=0; //Set depth to 0 where we have no data.
      water_depth[colNum][rowNum] = double(theValue);
    }//for
  }//for
    
  WaterDepthFile.close();
  return;
} //method ReadWaterDepth