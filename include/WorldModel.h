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
    // double distanceZ(const TVector3& p) const;
    // double distanceR(const TVector3& p) const;

    inline size_t N() const {
      return fXUps.size() + fXDowns.size();
    }

    double eval(const Geoid::Position& p) const;

    // inline bool above(const TVector3& p){
    //   return p.Z() >= distanceZ(p);
    // }

    // inline bool below(const TVector3& p){
    //   return p.Z() < distanceZ(p);
    // }

  private:
    void init() const;
    
    std::vector<double> fXUps;
    std::vector<double> fYUps;
    std::vector<double> fZUps;

    std::vector<double> fXDowns;
    std::vector<double> fYDowns;
    std::vector<double> fZDowns;

    std::vector<Geoid::Position> fPUps;
    std::vector<Geoid::Position> fPDowns;

    mutable ROOT::Math::Delaunay2D* fSurfUp = nullptr;
    mutable ROOT::Math::Delaunay2D* fSurfDown = nullptr;
    mutable bool fDoneInit = false;
  };
  
}


#endif 
