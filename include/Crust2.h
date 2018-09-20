#ifndef CRUST2_H
#define CRUST2_H

#include <string>
#include <cstdlib>
#include <map>

#include "Geoid.h"
#include "TVector3.h"

#include "WorldModel.h"

#include "TKDTree.h"

class TH2D;

namespace ROOT {
  namespace Math {
    class Delaunay2D;
  }
}


namespace icemc{

  class Primaries;
  class Interaction;
  class IceModel;
  class Settings;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //class Crust2.0
  //This class contains a model of the physical structure of the Earth, with information on shape,
  //density at all points inside, and depth of water and ice on the surface.
  //
  //   methods:
  //
  //IceThickness : Returns the thickness of ice in meters at a given Geoid::Position or lon/lat.
  //WaterDepth  : Returns the depth of water in meters at a given Geoid::Position or lon/lat.
  //Surface   : Returns the distance in meters from the center of the Earth to the surface,
  //            i.e. to the start of air.
  //SurfaceAboveGeoid : Returns the distance in meters from the local geoid to the start of air.
  //Geoid  : Returns the height in meters of the geoid at a given Geoid::Position or lon/lat.
  //RockSurface  : Returns the distance in meters from the center of the Earth to the end of rock
  //               (the begninning of the ice/water/air)
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
  class Crust2 : public WorldModel {

  public:
    Crust2(int model = 0,int WEIGHTABSORPTION_SETTING=1);
    virtual ~Crust2(){;}

    double radii[3];
    // = {1.2e3,(Earth::EarthRadiusMeters-4.0E4)*(Earth::EarthRadiusMeters-4.0E4),Earth::EarthRadiusMeters*Earth::EarthRadiusMeters}; // average radii of boundaries between earth layers

    double volume; // sums the volume of medium (ice or salt)
    double ice_area; // sums the area of the earth's surface that has antarctic ice underneath
    // double max_icevol_perbin; // maximum ice volume in any bin
    // double max_icethk_perbin;
    // virtual double Geoid(double latitude) const;
    // virtual double Geoid(const Geoid::Position &pos) const;

    virtual double GetMaxIceThickness() const {
      return fMaxIceThickness;
    }
    virtual double IceThickness(const Geoid::Position& pos) const;
    virtual int InFirn(const Geoid::Position& pos) const;
    virtual double SurfaceDeepIce(const Geoid::Position& pos) const;
    virtual double SurfaceAboveGeoid(const Geoid::Position& pos) const;
    virtual double WaterDepth(const Geoid::Position& pos) const;
    // virtual double RockSurface(const Geoid::Position& pos) const;
    virtual double Density(const Geoid::Position& pos) const {
      int ce;
      return Density(pos, ce);
    }
    virtual double Density(const Geoid::Position& pos, int& crust_entered) const;

    // virtual double fractionalIceVolumeWithinHorizon(const Geoid::Position& centeredOn, double horizonDistance) const;
    virtual double IceVolumeWithinHorizon(const Geoid::Position& p, double horizonDistanceMeters) const ;

    /** 
     * Figures out whether a neutrino will make it through the a Earth along a chord
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
    int Getchord(const Settings *settings1, double len_int_kgm2, const Geoid::Position &earth_in, const Geoid::Position &r_enterice,
		 const Geoid::Position &nuexitice, const Geoid::Position &posnu, int inu, ChordInfo& ci){
      return Getchord(settings1, len_int_kgm2, earth_in, r_enterice,
		      nuexitice, posnu, inu, ci.chord, ci.probability_tmp, ci.weight1_tmp,
		      ci.nearthlayers, ci.myair, ci.total_kgm2, ci.crust_entered, ci.mantle_entered, ci.core_entered);
    }

    int Getchord(const Settings *settings1, double len_int_kgm2, const Geoid::Position &earth_in, const Geoid::Position &r_enterice,
		 const Geoid::Position &nuexitice, const Geoid::Position &posnu, int inu,  double& chord, double& probability_tmp, double& weight1_tmp,
		 double& nearthlayers, double myair, double& total_kgm2, int& crust_entered, int& mantle_entered, int& core_entered);

    Geoid::Position WhereDoesItEnter(const Geoid::Position &posnu,const TVector3 &nnu) const;

    double integratePath(const Geoid::Position& interaction, const TVector3& neutrinoDir) const;
 
  protected:
    double fMaxIceThickness;
    
    int EARTH_MODEL;
    int CONSTANTICETHICKNESS;
    int CONSTANTCRUST;
    int FIXEDELEVATION;
    int FLATSURFACE;
    int weightabsorption;

    // pick method to step neutrino through the earth and get attenuation length
    static constexpr int getchord_method=2;	///< 1=core,mantle,crust; 2=crust 2.0 model

    static const double COASTLINE;		///< if the rf leaves from beyond this "coastline" (in degrees of latitude relative to south pole) it's not in antarctica.  Real coastline is always closer than this.						///< parameters of the Crust 2.0 earth model
    static constexpr int NLON=180;		///< number of bins in longitude for crust 2.0
    static constexpr int NLAT=90;		///< number of bins in latitude
    static constexpr int NPHI=180;		///< bins in longitude for visible ice in horizon
    double MIN_ALTITUDE_CRUST;			///< maximum depth of crust- determines when to start stepping
    double average_iceth;			///< average ice thickness over the continent-calculated once

    /////////////////////////////////////
    //methods
    void ReadCrust(const std::string&);
    Geoid::Position PickInteractionLocation(const Geoid::Position &detector) const;

    static const int numLayers = 11;
    enum class Layer {Air, ///@todo think about this one
		      Water,
		      Ice,
		      SoftSediment,
		      HardSediment,
		      UpperCrust,
		      MiddleCrust,
		      LowerCrust,
		      Mantle,
		      OuterCore,
		      InnerCore///@todo also special case
    };
    /** 
     * Utility function for ease of for-range based for loops and iteration
     * If the Layer enum class is expanded this array should too! 
     * Allows for nice code like:
     * for(auto layer : Layers() ){
     *    // something with a layer
     * }
     * @return const reference to a static array.
     */
    const std::array<Crust2::Layer, numLayers>& Layers() const {

      static std::array<Crust2::Layer, numLayers> allLayers {Crust2::Layer::Air,
							     Crust2::Layer::Water,
							     Crust2::Layer::Ice,
							     Crust2::Layer::SoftSediment,
							     Crust2::Layer::HardSediment,
							     Crust2::Layer::UpperCrust,
							     Crust2::Layer::MiddleCrust,
							     Crust2::Layer::LowerCrust,
							     Crust2::Layer::Mantle,
							     Crust2::Layer::OuterCore,
							     Crust2::Layer::InnerCore
      };
      return allLayers;
    }
    
    static Layer layerAbove(Layer);
    static Layer layerBelow(Layer);
    Layer getLayer(const Geoid::Position& pos, Layer startLayer=Layer::Air) const;
    double getLayerSurfaceAtPosition(const Geoid::Position& pos, Layer l) const;
    
  private:

    const std::string& getLayerName(Layer layer) const;
    Layer getLayerFromString(const std::string& layerType) const;
    std::map<Layer, std::string> fLayerNames;

    std::map<Layer, icemc::Mesh> fThicknesses;
    std::map<Layer, icemc::Mesh> fSurfaceMag; ///< Here mag means Mag() of the Position, not elevation
    std::map<Layer, icemc::Mesh> fDensities;
    icemc::Mesh fSurfaceAboveGeoid;
    std::shared_ptr<TKDTreeID> fKDTree = nullptr;
    std::vector<double> fXs;
    std::vector<double> fYs;
    std::vector<double> fZs;
    
    

  }; //class Earth

  /**
   * @todo put this somewhere in a function so it gets accessed by a descriptive name
   * And maybe with a bounds check too...
   */
  constexpr double densities[3]={14000.,3400.,2900.}; // average density of each earth layer
}

#endif
