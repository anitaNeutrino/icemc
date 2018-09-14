#ifndef ICEMC_WORLD_MODEL_H
#define ICEMC_WORLD_MODEL_H

#include "Geoid.h"
#include "TVector3.h"
#include "TGraph.h"
#include "Math/Delaunay2D.h"

namespace icemc {

  class Detector;

  /**
   * @class WorldModel
   * @brief Interface to icemc for any class seeking to model the world.
   * 
   * In order to model neutrinos passing through the Earth and interacting, you 
   * need some kind of model that tells you where the Earth's mass is and where the ice is.
   * 
   * This class defines the interface icemc requires any such model to have.
   * It also defines some methods which use those data, even though this class contains none.
   */
  class WorldModel {

  public:
    virtual double SurfaceAboveGeoid(const Geoid::Position& p) const = 0;
    virtual double IceThickness(const Geoid::Position& p) const = 0;

    /** 
     * This one is actually implemented in this class.
     * If have model the surface at each point, finding the normal is easy.
     * 
     * @param p is the position at which we wish to find the surface normal
     * 
     * @return vector representing the surface normal
     */
    virtual TVector3 GetSurfaceNormal(const Geoid::Position& p) const;

    virtual inline double Surface(const Geoid::Position& p) const {
      return p.Surface() + SurfaceAboveGeoid(p);
    };

    inline double IceVolumeAtTime(double time){
      return fHorizons.Eval(time);
    }

    enum {defaultTimeStep = 60}; // seconds
    virtual double CreateHorizons(Detector* detector,  double horizonDistanceMeters, double timeStepSeconds = defaultTimeStep);
    

    virtual double IceVolumeWithinHorizon(const Geoid::Position& p,  double horizonDistance) const = 0;
  private:
    
    
    TGraph fHorizons;
  };



  /**
   * @class Mesh
   * @brief Models a surface of a layer of the Earth
   */

  class Mesh {
  public:
    Mesh();
    virtual ~Mesh();

    size_t addPoint(const Geoid::Position& p, double val);

    inline size_t N() const {
      return fPositions.size();
    }

    double eval(const Geoid::Position& p) const;

  private:
    void init() const;
    
    std::vector<Geoid::Position> fPositions;
    std::vector<double> fValues;    
    
    mutable ROOT::Math::Delaunay2D* fXUp = nullptr;
    mutable ROOT::Math::Delaunay2D* fXDown = nullptr;
    
    mutable ROOT::Math::Delaunay2D* fYUp = nullptr;
    mutable ROOT::Math::Delaunay2D* fYDown = nullptr;
    
    mutable ROOT::Math::Delaunay2D* fZUp = nullptr;
    mutable ROOT::Math::Delaunay2D* fZDown = nullptr; 
    mutable bool fDoneInit = false;

    mutable std::vector<double> xUpX, xUpY, xUpZ;
    mutable std::vector<double> yUpX, yUpY, yUpZ;
    mutable std::vector<double> zUpX, zUpY, zUpZ;
    mutable std::vector<double> xDownX, xDownY, xDownZ;
    mutable std::vector<double> yDownX, yDownY, yDownZ;
    mutable std::vector<double> zDownX, zDownY, zDownZ;

    
  };
  
}


#endif 
