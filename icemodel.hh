#ifndef ICEMODEL_HH_
#define ICEMODEL_HH_

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "TVector3.h" 

#include "earthmodel.hh"
#include "source.hh" 
#include "TH2.h" 

class Balloon;
class Position;
class Vector;
class EarthModel;
class Settings;
class Interaction;
class Ray;

//Constants relating to all ice models
const double FIRNDEPTH=-150.;                // depth of the firn, in meters: currently a constant over all ice
using std::ofstream;
using std::vector;
using std::string;



//Variables for conversion between BEDMAP polar stereographic coordinates and lat/lon.  Conversion equations from ftp://164.214.2.65/pub/gig/tm8358.2/TM8358_2.pdf
const double scale_factor=0.97276901289;  //scale factor at pole corresponding to 71 deg S latitude of true scale (used in BEDMAP)
const double ellipsoid_inv_f = 298.257223563; //of Earth
const double ellipsoid_b = EarthModel::R_EARTH*(1-(1/ellipsoid_inv_f));
const double eccentricity = sqrt((1/ellipsoid_inv_f)*(2-(1/ellipsoid_inv_f)));
const double bedmap_a_bar = pow(eccentricity,2)/2 + 5*pow(eccentricity,4)/24 + pow(eccentricity,6)/12 + 13*pow(eccentricity,8)/360;
const double bedmap_b_bar = 7*pow(eccentricity,4)/48 + 29*pow(eccentricity,6)/240 + 811*pow(eccentricity,8)/11520;
const double bedmap_c_bar = 7*pow(eccentricity,6)/120 + 81*pow(eccentricity,8)/1120;
const double bedmap_d_bar = 4279*pow(eccentricity,8)/161280;
const double bedmap_c_0 = (2*EarthModel::R_EARTH / sqrt(1-pow(eccentricity,2))) * pow(( (1-eccentricity) / (1+eccentricity) ),eccentricity/2);

struct icemodel_debug
{
  TVector3 balloon;
  TVector3 ref; 
  int nintersect; 
  TVector3 int1; 
  TVector3 int2; 
  TVector3 pint; 
  TVector3 nudir;
  bool rightdir1; 
  bool rightdir2; 
  bool noway; 
  bool leave_err; 
  TVector3 nupos; 
  TVector3 enterice; 
  TVector3 exitice; 
  TVector3 exitearth; 
  double x,y; 

  void reset() 
  {
    balloon.SetXYZ(0,0,0);
    ref.SetXYZ(0,0,0);
    int1.SetXYZ(0,0,0);
    int2.SetXYZ(0,0,0);
    pint.SetXYZ(0,0,0);
    nudir.SetXYZ(0,0,0);
    nupos.SetXYZ(0,0,0);
    enterice.SetXYZ(0,0,0);
    exitice.SetXYZ(0,0,0);
    exitearth.SetXYZ(0,0,0);
    nintersect = 0;
    rightdir1 = false;
    rightdir2 = false;
    noway = false;
    x = 0;
    y = 0; 
    leave_err = false; 
  }
}; 



//! Ice thicknesses and water depth
class IceModel : public EarthModel {


public:

  //BEDMAP data
  TH2D h_ice_thickness;
  TH2D h_ground_elevation;;
  TH2D h_water_depth;;


double bedmap_R; //varies with latitude, defined here for 71 deg S latitude
double bedmap_nu;


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


  void IceENtoLonLat(int e,
		     int n,
		     
		     double& lon,
		     double& lat);  

  void GroundENtoLonLat(int e,
			int n,
			double& lon,
			double& lat);
  void WaterENtoLonLat(int e,
			int n,
			double& lon,
			double& lat);



  vector<double> volume_inhorizon; // volume of ice within horizon for each balloon phi position 
  IceModel(int model=0,int earth_mode=0,int WEIGHTABSORPTION_SETTING=1);
  virtual ~IceModel(); 
  double IceThickness(double lon,double lat);
  double IceThickness(const Position& pos) ;
  double Surface(double lon,double lat) ;
  double Surface(const Position& pos) ;
  double SurfaceAboveGeoid(double lon,double lat);
  double SurfaceAboveGeoid(const Position& pos) ;
  double WaterDepth(double lon,double lat);
  double WaterDepth(const Position& pos);
  Position PickInteractionLocation(int ibnposition, Settings *settings1, const Position &rbn, Interaction *interaction1);
  Position PickBalloonPosition();
  void GetMAXHORIZON(Balloon *bn1); // get upper limit on the horizon wrt the balloon.
  int RossIceShelf(const Position &position); 
  int IceOnWater(const Position &postition);
  int RossExcept(const Position &position);
  int RonneIceShelf(const Position &position);
  int WestLand(const Position &pos); 
  int AcceptableRfexit(const Vector &nsurf_rfexit,const Position &rfexit,const Vector &n_exit2rx); 
  double GetBalloonPositionWeight(int ibnpos);
  int OutsideAntarctica(const Position &pos);
  int OutsideAntarctica(double lat);
  int WhereDoesItEnterIce(const Position &posnu,
			       const Vector &nnu,
			       double stepsize,
			       Position &r_enterice);

  int WhereDoesItExitIce(const Position &posnu,
			 const Vector &nnu,
			 double stepsize,
			 Position &r_enterice);
  int WhereDoesItExitIceForward(const Position &posnu,
			 const Vector &nnu,
			 double stepsize,
			 Position &r_enterice);


  //resolution is in cartesian meters 
  void CreateCartesianTopAndBottom(int resolution, bool force_new =false); 
  const TH2 * GetCartesianTop() const { return  cart_ice_top; }
  const TH2 * GetCartesianBottom() const { return  cart_ice_bot; }
  bool CartesianIsInIce(double x, double y, double z); 
  
  int GetIceIntersectionsCartesian(const Position &posnu,  const Vector &nnu, 
      std::vector<std::pair<double,double> > & intersections, double initial_step_size = 50, int map_resolution=1000); 


  void CreateHorizons(Settings *settings1,Balloon *bn1,double theta_bn,double phi_bn,double altitude_bn,ofstream &foutput);
  Vector GetSurfaceNormal(const Position &r_out); //overloaded from EarthModel to include procedures for new ice models.
  double GetN(double depth);
  double GetN(const Position &pos);
  double EffectiveAttenuationLength(Settings *settings1,const Position &pos, const int &whichray);
  
  void FillArraysforTree(double lon_ground[1068][869],double lat_ground[1068][869],double lon_ice[1200][1000],double lat_ice[1200][1000],double lon_water[1200][1000],double lat_water[1200][1000]);
  int PickUnbiased(Interaction *interaction1, double len_int_kgm2, double & position_weight, double chord_step, Vector * force_dir = 0);

  int PickUnbiasedPointSourceNearBalloon(Interaction *interaction1, 
      const Position * balloon_position, double max_ps_distance, double chord_step, 
      double len_int_kgm2, const Vector * force_dir = 0);

  double getSampleX() const { return sample_x; } 
  double getSampleY() const { return sample_y; } 

  void LonLattoEN(double lon, 
		  double lat,
		  double& E, 
		  double& N);


protected:
  int ice_model;
  int DEPTH_DEPENDENT_N;


  //Information on horizons - what ice the balloon can see at each position along its path.
 
 
  double volume_inhorizon_average; // average volume of ice seen by balloon
  vector< vector<int> > ilon_inhorizon; // indices in lon and lat for bins in horizon for NPHI balloon positions along 80 deg latitude line.
  vector< vector<int> >ilat_inhorizon;
  vector< vector<int> >easting_inhorizon; //indicies in easting and northing for bins in horizon for NPHI balloon positions along 80 deg latitude line.
  vector< vector<int> >northing_inhorizon;
  vector<double> maxvol_inhorizon; // maximum volume of ice for a bin 

  double cart_max_z; 
  double cart_min_z; 
  double sample_x, sample_y; 

  TH2 *cart_ice_top; 
  TH2 *cart_ice_bot; 
  TFile *cart_ice_file; 
  int cart_resolution; 

  //BEDMAP utility methods
  double Area(double latitude);

  void ENtoLonLat(int e_coord, 
		  int n_coord,
		  double xLowerLeft,
		  double yLowerLeft,
		  
		  double& lon, 
		  double& lat) ;




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


}; //class IceModel

#endif
