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

    // inline double IceVolumeAtTime(double time){
    //   return fHorizons.Eval(time);
    // }

    // enum {defaultTimeStep = 60}; // seconds
    // virtual double CreateHorizons(Detector* detector,  double horizonDistanceMeters, double timeStepSeconds = defaultTimeStep);
    

    virtual double IceVolumeWithinHorizon(const Geoid::Position& p,  double horizonDistance) const = 0;
  private:
    
    
    TGraph fHorizons;
  };



  /**
   * @class Mesh
   * @brief Models a surface of a layer of the Earth.
   * 
   * We want to know what some value is as a function of f(x,y,z)
   * It has 6 internal Delaunay triangulations, 2 per dimension: top/bottom/front/back/left/right.
   * 
   * When a value is requested with eval(const Geoid::Position& p) it finds the position on the Geoid 
   * surface, s(x,y,z), then it queries the Delaunay triangulation with it is nearest the center of.
   * This prevents edge effects in the model.
   * It may create some slightly unsmooth edge effects at the boundaries of the meshes.
   * 
   * @todo we will have multiple meshes at the same points, could make Geoid::Positions shared?
   * @todo or inherit from Delaunay2D and do all the interpolation in one place
   */

  using ROOT::Math::Delaunay2D;
  
  class Mesh {

    class Point {
    public:
      Point(const Geoid::Position& p, double v) : position(p), value(v) {}
      Geoid::Position position;
      double value;
    };
    
  public:
    Mesh();
    virtual ~Mesh();

    size_t addPoint(const Geoid::Position& p, double val);

    inline size_t N() const {
      return fPoints.size();
    }
    std::vector<Point>::const_iterator begin() const {return fPoints.begin();}
    std::vector<Point>::const_iterator end() const {return fPoints.end();}

    double eval(const Geoid::Position& p) const;
    

  private:
    void init() const;

    std::vector<Point> fPoints;

    mutable std::unique_ptr<Delaunay2D> fXUp = nullptr;
    mutable std::unique_ptr<Delaunay2D> fXDown = nullptr;

    mutable std::unique_ptr<Delaunay2D> fYUp = nullptr;
    mutable std::unique_ptr<Delaunay2D> fYDown = nullptr;

    mutable std::unique_ptr<Delaunay2D> fZUp = nullptr;
    mutable std::unique_ptr<Delaunay2D> fZDown = nullptr;
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
