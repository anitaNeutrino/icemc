#ifndef EARTHMODEL_H
#define EARTHMODEL_H

#include <string>
#include <cstdlib>
#include <map>

#include "TRandom3.h"
#include "TH2D.h"


#include "GeoidModel.h"
#include "TVector3.h"

namespace icemc{

  class Primaries;
  class Interaction;
  class IceModel;
  class Settings;
  


  ////////////////////////////////////////////////////////////////////////////////////////////////
  //class EarthModel
  //This class contains a model of the physical structure of the Earth, with information on shape,
  //density at all points inside, and depth of water and ice on the surface.
  //
  //   methods:
  //
  //IceThickness : Returns the thickness of ice in meters at a given GeoidModel::Position or lon/lat.
  //WaterDepth  : Returns the depth of water in meters at a given GeoidModel::Position or lon/lat.
  //Surface   : Returns the distance in meters from the center of the Earth to the surface,
  //            i.e. to the start of air.
  //SurfaceAboveGeoid : Returns the distance in meters from the local geoid to the start of air.
  //Geoid  : Returns the height in meters of the geoid at a given GeoidModel::Position or lon/lat.
  //RockSurface  : Returns the distance in meters from the center of the Earth to the end of rock
  //               (the begninning of the ice/water/air)
  //GetSurfaceNormal  : Returns a unit vector pointing in the direction of the surface normal at
  //                    a given Position.
  //Getchord
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * @class ChordInfo
   * @brief A little class to hold information about a neutrino's
   * journey through a cord length of the Earth
   */

  class ChordInfo {
  public:
    ChordInfo() : chord(0), probability_tmp(0),  weight1_tmp(0),  nearthlayers(0), 
		  myair(0), total_kgm2(0), crust_entered(0), mantle_entered(0), 
		  core_entered(0)
    { }

    double chord; ///< chord length
    double probability_tmp; ///< weight
    double weight1_tmp;
    double nearthlayers; ///< core; mantle; crust
    double myair;
    double total_kgm2; ///< length in kg m^2
    int crust_entered; ///< 1 or 0
    int mantle_entered; ///< 1 or 0
    int core_entered; ///< 1 or 0
    ClassDef(ChordInfo, 1);
  };

  //! Shape of the earth, ice thicknesses, profiles of earth layers, densities, neutrino absorption
  class Earth {

  public:
    Earth(int model = 0,int WEIGHTABSORPTION_SETTING=1);
    virtual ~Earth();

    // Rename to avoid collision with AnitaGeomTool, although it's that lib's fault
    static constexpr double BulgeRadius=6.378140E6; ///< radius of Earth in m at bulge
  
    double radii[3];
    // = {1.2e13,(Earth::EarthRadiusMeters-4.0E4)*(Earth::EarthRadiusMeters-4.0E4),Earth::EarthRadiusMeters*Earth::EarthRadiusMeters}; // average radii of boundaries between earth layers

    double volume; // sums the volume of medium (ice or salt)
    double ice_area; // sums the area of the earth's surface that has antarctic ice underneath
    // double max_icevol_perbin; // maximum ice volume in any bin
    // double max_icethk_perbin;
    // virtual double Geoid(double latitude) const;
    // virtual double Geoid(const GeoidModel::Position &pos) const;

    virtual double GetMaxIceThickness() const {
      return fMaxIceThickness;
    }
    virtual double IceThickness(double lon,double lat) const;
    virtual double IceThickness(const GeoidModel::Position& pos) const;
    virtual double Surface(double lon,double lat) const;    
    virtual double Surface(const GeoidModel::Position& pos) const;    
    virtual int InFirn(const GeoidModel::Position& pos) const;
    virtual double SurfaceDeepIce(const GeoidModel::Position& pos) const;
    virtual double SurfaceAboveGeoid(double lon,double lat) const;
    virtual double SurfaceAboveGeoid(const GeoidModel::Position& pos) const;
    virtual double WaterDepth(double lon,double lat) const;
    virtual double WaterDepth(const GeoidModel::Position& pos) const;
    virtual double RockSurface(double lon,double lat) const;
    virtual double RockSurface(const GeoidModel::Position& pos) const;
    double GetDensity(const GeoidModel::Position& pos, int& crust_entered) const;

    /** 
     * Figures out whether a neutrino will make it through thea Earth along a chord
     * (Or gives the journey a weight)
     * 
     * @param settings1 The simulation settings
     * @param len_int_kgm2 
     * @param earth_in place where neutrino entered the earth
     * @param r_enterice 
     * @param nuexitice 
     * @param posnu position of the interaction
     * @param inu index of generated neutrino
     * @param chordInfo hold quantites calculated from the neutrino journey
     * 
     * @return 1 if it makes it, 0 otherwise
     */
    int Getchord(const Settings *settings1, double len_int_kgm2, const GeoidModel::Position &earth_in, const GeoidModel::Position &r_enterice,
		 const GeoidModel::Position &nuexitice, const GeoidModel::Position &posnu, int inu, ChordInfo& ci){
      return Getchord(settings1, len_int_kgm2, earth_in, r_enterice,
		      nuexitice, posnu, inu, ci.chord, ci.probability_tmp, ci.weight1_tmp,
		      ci.nearthlayers, ci.myair, ci.total_kgm2, ci.crust_entered, ci.mantle_entered, ci.core_entered);
    }

    int Getchord(const Settings *settings1, double len_int_kgm2, const GeoidModel::Position &earth_in, const GeoidModel::Position &r_enterice,
		 const GeoidModel::Position &nuexitice, const GeoidModel::Position &posnu, int inu,  double& chord, double& probability_tmp, double& weight1_tmp,
		 double& nearthlayers, double myair, double& total_kgm2, int& crust_entered, int& mantle_entered, int& core_entered);

    TVector3 GetSurfaceNormal(const GeoidModel::Position &r_out) const;
    // static double LongtoPhi_0isPrimeMeridian(double longitude); // convert longitude to phiwith 0 longitude being the prime meridian
    // static double LongtoPhi_0is180thMeridian(double longitude); // convert longitude to phi with 0 longitude being at the 180th meridian
    // void EarthCurvature(double *array,double depth_temp) const; // adjusts coordinates within the mine to account for the curvature of the earth.
    GeoidModel::Position WhereDoesItEnter(const GeoidModel::Position &posnu,const TVector3 &nnu) const;

 
  protected:
    int EARTH_MODEL;
    int CONSTANTICETHICKNESS;
    int CONSTANTCRUST;
    int FIXEDELEVATION;
    int FLATSURFACE;
    int weightabsorption;


    // pick method to step neutrino through the earth and get attenuation length
    static constexpr int getchord_method=2;	///< 1=core,mantle,crust; 2=crust 2.0 model

    static const double COASTLINE;		///< if the rf leaves from beyond this "coastline" (in degrees of latitude relative to south pole) it's not in antarctica.  Real coastline is always closer than this.

						///< parameters of the Crust 2.0 earth model
    static constexpr int NLON=180;		///< number of bins in longitude for crust 2.0
    static constexpr int NLAT=90;		///< number of bins in latitude
    static constexpr int NPHI=180;		///< bins in longitude for visible ice in horizon
    // static const double MAXTHETA;		///< maximum value of theta in degrees in polar coordinates
    // double thetastep;				///< how big do you step in theta-> always 2deg with Crust 2.0
    // double phistep;				///< how big do you step in phi->always 2deg in Crust 2.0
    // static const int ILAT_COASTLINE;
    // double surfacer[NLON][NLAT];		///< elevation at the surface (top of ice) compared to geoid (in meters)
    // double icer[NLON][NLAT];			///< elevation at the *bottom* of ice layer (in meters)
    // double waterr[NLON][NLAT];			///< elevation at the bottom of water layer (in meters)
    // double softsedr[NLON][NLAT];		///< elevation at the bottom of soft set layer (in meters)
    // double hardsedr[NLON][NLAT];		///< elev at bottom of hard sed layer (in meters)
    // double uppercrustr[NLON][NLAT];		///< elev at bottom of upper crust layer (in meters)
    // double middlecrustr[NLON][NLAT];		///< elev at bottom of middle crust layer (in meters)
    // double lowercrustr[NLON][NLAT];		///< elev at bottom of lower crust layer (in meters)
    // double geoid[NLAT];				///< realistic shape of earth-radius at each latitude (in meters)
    double MIN_ALTITUDE_CRUST;			///< maximum depth of crust- determines when to start stepping
						///<double MAX_VOL; ///< maximum volume of ice in a bin in Crust 2.0 - not used


    // double elevationarray[NLON][NLAT];		///< If no water, measures the elevation (relative to geoid, in meters) of the top of the ice or rock (i.e., air interface).  If there is water, measures elevation to bottom of water. (There may or may not be ice on top of the water.)
    // double waterthkarray[NLON][NLAT];		///< thickness of water layer (in km)
    // double icethkarray[NLON][NLAT];		///< thickness of ice layer (in km)
    // double softsedthkarray[NLON][NLAT];		///< thickness of soft sed layer (in km)
    // double hardsedthkarray[NLON][NLAT];		///< thickness of hard sed layer (in km)
    // double uppercrustthkarray[NLON][NLAT];	///< thickness of upper crust layer (in km)
    // double middlecrustthkarray[NLON][NLAT];	///< thickness of middle crust layer (in km)
    // double lowercrustthkarray[NLON][NLAT];	///< thickness of lower crust layer (in km)
    // double crustthkarray[NLON][NLAT];		///< total thickness of crust (in km)
    // double waterdensityarray[NLON][NLAT];	///< density of water layer bin by bin
    // double icedensityarray[NLON][NLAT];		///< density of ice layer bin by bin
    // double softseddensityarray[NLON][NLAT];	///< density of soft sed layer
    // double hardseddensityarray[NLON][NLAT];	///< density of hard sed layer
    // double uppercrustdensityarray[NLON][NLAT];	///< density of upper crust layer
    // double middlecrustdensityarray[NLON][NLAT]; ///< density of middle crust layer
    // double lowercrustdensityarray[NLON][NLAT];	///< density of lower crust layer
    // double area[NLAT];				///< area of a bin at a given latitude- calculated once
    double average_iceth;			///< average ice thickness over the continent-calculated once

    /////////////////////////////////////
    //methods
    void ReadCrust(const std::string&);
    // double SmearPhi(int ilon, double rand) const;
    // double SmearTheta(int ilat, double rand) const;
    // GeoidModel::Position PickPosnuForaLonLat(double lon,double lat) const; // given that an interaction occurs at a lon and lat, pick an interaction position in the ice
    GeoidModel::Position PickInteractionLocation(const GeoidModel::Position &balloon, Interaction *interaction1) const;


    static const int numCrustLayers = 7;
    enum class CrustLayer {Water,
			   Ice,
			   SoftSediment,
			   HardSediment,
			   UpperCrust,
			   MiddleCrust,
			   LowerCrust};

    /** 
     * Utility function for ease of for-range based for loops and iteration
     * If the CrustLayer enum class is expanded this array should too! 
     * Allows for nice code like:
     * for(auto layer : CrustLayers() ){
     *    // something with a layer
     * }
     * @return const reference to a static array.
     */
    inline const std::array<CrustLayer, numCrustLayers> CrustLayers() const {
      static std::array<CrustLayer, numCrustLayers> allLayers {CrustLayer::Water,
							       CrustLayer::Ice,
							       CrustLayer::SoftSediment,
							       CrustLayer::HardSediment,
							       CrustLayer::UpperCrust,
							       CrustLayer::MiddleCrust,
							       CrustLayer::LowerCrust};
      return allLayers;
    }

    enum class CrustProperty {ThicknessKm,
			      ElevationLowerBoundary, // this is the elevation at the BOTTOM of the layer!
			      Density};

  protected:
    double fMaxIceThickness = 0;    
  private:



    const std::string& getLayerName(CrustLayer layer) const;
    CrustLayer getLayerFromString(const std::string& layerType) const;

    const std::string& getPropertyName(CrustProperty property) const;

    TH2D* makeRawHist(CrustLayer layer, CrustProperty property);
    TH2D* getRawHist(CrustLayer layer, CrustProperty property) const;
    inline double get(CrustLayer layer, CrustProperty property, const GeoidModel::Position& p) const{
      return getRawHist(layer, property)->Interpolate(p.Longitude(), p.Latitude());
    }

    std::map<CrustLayer, std::string> fLayerNames;
    std::map<CrustProperty, std::string> fPropertyNames;
    typedef std::pair<CrustLayer, CrustProperty> CrustKey;
    std::map<CrustKey, TH2D*> fRawCrustData;

    inline TH2D* getRawHist(CrustKey k) const {
      return getRawHist(k.first,  k.second);
    }
    TH2D* fCrustElevation = nullptr;
    TH2D* fSurfaceAboveGeoid = nullptr; // the top of of everything





    

    template <class T>
    static const std::string& findThingInMapToString(const std::map<T, std::string>& m, T t){
      auto it = m.find(t);
      if(it!=m.end()){
	return it->second;
      }
      else{
	static const std::string unknown("unknown");
	return unknown;
      }      
    }



    

  }; //class Earth

  /**
   * @todo put this somewhere in a function so it gets accessed by a descriptive name
   * And maybe with a bounds check too...
   */
  constexpr double densities[3]={14000.,3400.,2900.}; // average density of each earth layer
}

#endif
