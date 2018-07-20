#ifndef ICEMC_ANTARCTICA_H
#define ICEMC_ANTARCTICA_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "Earth.h"

class TRandom3;

namespace icemc {
  class Settings;
  class Interaction;
  class Ray;  
  class Balloon;
  class Earth;


  //Variables for conversion between BEDMAP polar stereographic coordinates and lat/lon.  Conversion equations from ftp://164.214.2.65/pub/gig/tm8358.2/TM8358_2.pdf
  const double scale_factor=0.97276901289;  //scale factor at pole corresponding to 71 deg S latitude of true scale (used in BEDMAP)
  const double ellipsoid_inv_f = 298.257223563; //of Earth
  // const double ellipsoid_b = Earth::R_EARTH*(1-(1/ellipsoid_inv_f));
  const double ellipsoid_b = Earth::BulgeRadius*(1-(1/ellipsoid_inv_f));  
  const double eccentricity = sqrt((1/ellipsoid_inv_f)*(2-(1/ellipsoid_inv_f)));
  const double bedmap_a_bar = pow(eccentricity,2)/2 + 5*pow(eccentricity,4)/24 + pow(eccentricity,6)/12 + 13*pow(eccentricity,8)/360;
  const double bedmap_b_bar = 7*pow(eccentricity,4)/48 + 29*pow(eccentricity,6)/240 + 811*pow(eccentricity,8)/11520;
  const double bedmap_c_bar = 7*pow(eccentricity,6)/120 + 81*pow(eccentricity,8)/1120;
  const double bedmap_d_bar = 4279*pow(eccentricity,8)/161280;
  const double bedmap_c_0 = (2*Earth::BulgeRadius / sqrt(1-pow(eccentricity,2))) * pow(( (1-eccentricity) / (1+eccentricity) ),eccentricity/2);
  // const double bedmap_c_0 = (2*Earth::R_EARTH / sqrt(1-pow(eccentricity,2))) * pow(( (1-eccentricity) / (1+eccentricity) ),eccentricity/2);  

  //! Ice thicknesses and water depth
  class Antarctica : public Earth {


  public:

    //BEDMAP data
    double ice_thickness_array[1200][1000];  //thickness of the ice
    double ground_elevation[1068][869]; //elevation above geoid at which ice starts
    double water_depth[1200][1000]; //depth of water under ice


    // double bedmap_R; //varies with latitude, defined here for 71 deg S latitude
    // double bedmap_nu;


    //Parameters of the BEDMAP ice model. (See http://www.antarctica.ac.uk/aedc/bedmap/download/)
    int nCols_ice; //number of columns in data, set by header file (should be 1200)
    int nRows_ice; //number of rows in data, set by header file (should be 1000)
    int cellSize; //in meters, set by header file (should be 5000) - same for both files
    int xLowerLeft_ice; 
    int yLowerLeft_ice;
    int nCols_ground;
    int nRows_ground;
    int xLowerLeft_ground;
    int yLowerLeft_ground;
    int nCols_water;
    int nRows_water;
    int xLowerLeft_water;
    int yLowerLeft_water;
    int NODATA;


    void IceENtoLonLat(int e, int n, double& lon, double& lat) const;  
    void GroundENtoLonLat(int e, int n,	double& lon, double& lat) const;

    const static int NBNPOSITIONS_MAX=26000;
    double volume_inhorizon[NBNPOSITIONS_MAX]; // volume of ice within horizon for each balloon phi position 
    // const static int NBNPOSITIONS_MAX=26000;
    // double volume_inhorizon[NBNPOSITIONS_MAX]; // volume of ice within horizon for each balloon phi position 
    Antarctica(int model=0,int earth_mode=0,int WEIGHTABSORPTION_SETTING=1);
    double IceThickness(double lon,double lat) const;
    double IceThickness(const GeoidModel::Position& pos) const;
    double Surface(double lon,double lat) const;
    double Surface(const GeoidModel::Position& pos) const;
    double SurfaceAboveGeoid(double lon,double lat) const;
    double SurfaceAboveGeoid(const GeoidModel::Position& pos) const;
    double WaterDepth(double lon,double lat) const;
    double WaterDepth(const GeoidModel::Position& pos) const;
    GeoidModel::Position PickInteractionLocation(const GeoidModel::Position &balloon) const;
    void GetMAXHORIZON(Balloon *bn1) const; // get upper limit on the horizon wrt the balloon.
    int RossIceShelf(const GeoidModel::Position &position) const; 
    int IceOnWater(const GeoidModel::Position &postition) const;
    int RossExcept(const GeoidModel::Position &position) const;
    int RonneIceShelf(const GeoidModel::Position &position) const;
    int WestLand(const GeoidModel::Position &pos) const; 
    int AcceptableRfexit(const TVector3 &nsurf_rfexit,const GeoidModel::Position &rfexit,const TVector3 &n_exit2rx); 
    double GetBalloonPositionWeight(int ibnpos) const;
    int OutsideAntarctica(const GeoidModel::Position &pos) const;
    int OutsideAntarctica(double lat) const;
    int WhereDoesItEnterIce(const GeoidModel::Position &posnu,
			    const TVector3 &nnu,
			    double stepsize,
			    GeoidModel::Position &r_enterice) const;

    int WhereDoesItExitIce(const GeoidModel::Position &posnu,
			   const TVector3 &nnu,
			   double stepsize,
			   GeoidModel::Position &r_enterice) const;
    int WhereDoesItExitIceForward(const GeoidModel::Position &posnu,
				  const TVector3 &nnu,
				  double stepsize,
				  GeoidModel::Position &r_enterice) const;
    void CreateHorizons(const Settings *settings1,Balloon *bn1,double theta_bn,double phi_bn,double altitude_bn);
    TVector3 GetSurfaceNormal(const GeoidModel::Position &r_out) const; //overloaded from Earth to include procedures for new ice models.
    double GetN(double depth) const;
    double GetN(const GeoidModel::Position &pos) const;
    double EffectiveAttenuationLength(const Settings *settings1,const GeoidModel::Position &pos, const int &whichray) const;
  
    void IceLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const;

    int PickUnbiased(Interaction *interaction1) const;


  protected:
    std::string fDataDir;
    int ice_model;
    int DEPTH_DEPENDENT_N;

    //Information on horizons - what ice the balloon can see at each position along its path.
 
    double volume_inhorizon_average; // average volume of ice seen by balloon
    std::vector<int> ilon_inhorizon[NBNPOSITIONS_MAX]; // indices in lon and lat for bins in horizon for NPHI balloon positions along 80 deg latitude line.
    std::vector<int> ilat_inhorizon[NBNPOSITIONS_MAX];
    std::vector<int> easting_inhorizon[NBNPOSITIONS_MAX]; //indicies in easting and northing for bins in horizon for NPHI balloon positions along 80 deg latitude line.
    std::vector<int> northing_inhorizon[NBNPOSITIONS_MAX];
    double maxvol_inhorizon[NBNPOSITIONS_MAX]; // maximum volume of ice for a bin 

    //BEDMAP utility methods
    double Area(double latitude) const;

    void ENtoLonLat(int e_coord, int n_coord, double xLowerLeft, double yLowerLeft, double& lon, double& lat) const;



    void WaterENtoLonLat(int e,
			 int n,
		       
			 double& lon,
			 double& lat) const;
    void LonLattoEN(double lon, 
		    double lat,
		    double xLowerLeft,
		    double yLowerLeft,
		
		    int& e_coord, 
		    int& n_coord) const;
 
    void GroundLonLattoEN(double lon, 
			  double lat,
			
			  int& e_coord, 
			  int& n_coord) const;
    void WaterLonLattoEN(double lon,
			 double lat,
		       
			 int& e_coord,
			 int& n_coord) const;

    //BEDMAP data input methods
    void ReadIceThickness();
    void ReadGroundBed();
    void ReadWaterDepth();

  private:

    const static int N_sheetup=2810;
    double d_sheetup[N_sheetup], l_sheetup[N_sheetup];
    const static int N_shelfup=420;
    double d_shelfup[N_shelfup], l_shelfup[N_shelfup];
    const static int N_westlandup=420;
    double d_westlandup[N_westlandup],l_westlandup[N_westlandup];

    const static int N_sheetdown=2810;
    double d_sheetdown[N_sheetup], l_sheetdown[N_sheetdown];
    const static int N_shelfdown=420;
    double d_shelfdown[N_shelfdown], l_shelfdown[N_shelfdown];
    const static int N_westlanddown=420;
    double d_westlanddown[N_westlanddown],l_westlanddown[N_westlanddown];


  };
}


#endif //ICEMC_ANTARCTICA_H
