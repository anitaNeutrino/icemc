#ifndef ICEMC_WORLD_MODEL_H
#define ICEMC_WORLD_MODEL_H

#include "Geoid.h"
#include "TVector3.h"
#include "TGraph.h"
#include "Math/Delaunay2D.h"
#include "TKDTree.h"


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
    virtual std::pair<Geoid::Position, double> integratePath(const Geoid::Position& start, const TVector3& direction) const = 0;
    virtual double SurfaceAboveGeoid(const Geoid::Position& p) const = 0;
    virtual double IceThickness(const Geoid::Position& p) const = 0;
    virtual double maxIceThicknessWithinDistance(const Geoid::Position& p, double distanceMeters) const = 0;    
    virtual double Density(const Geoid::Position& p) const = 0;

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
      return p.EllipsoidSurface() + SurfaceAboveGeoid(p);
    };

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
   * When a value is requested with eval(const Geoid::Position& p) it finds the position on the Geoid.
   * surface, s(x,y,z), then it queries the Delaunay triangulation it's closest to the center of.
   * This prevents edge effects in the model.
   * It may create some slightly unsmooth edge effects at the boundaries of the meshes.
   * 
   * @todo we will have multiple meshes at the same points, could make Geoid::Positions shared?
   * @todo or inherit from Delaunay2D and do all the interpolation in one place
   */

  using ROOT::Math::Delaunay2D;
  
  class Mesh : public TNamed {

  public:
    Mesh(const char* name = "Mesh", const char* title = "Mesh")
      : TNamed(name, title),
	fUpX(Hemisphere::Axis::x),
	fDownX(Hemisphere::Axis::x),	
	fUpY(Hemisphere::Axis::y),
	fDownY(Hemisphere::Axis::y),	
	fUpZ(Hemisphere::Axis::z),
	fDownZ(Hemisphere::Axis::z)
    {;}
    
    virtual ~Mesh() {;}

    size_t addPoint(const Geoid::Position& p, double val);

    inline size_t N() const {
      return fNumPoints;
    }
    void build();

    double eval(const Geoid::Position& p) const;
    double findMaxWithinDistance(const  Geoid::Position& p, double distanceMeters) const;

  private:

    struct Points {
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;
      std::vector<double> v;
      size_t add(const Geoid::Position& p, double val){
	x.push_back(p.X());
	y.push_back(p.Y());
	z.push_back(p.Z());
	v.push_back(val);
	return v.size();
      }
    };

    class Hemisphere {      
      static const int nDim = 3;
      static const UInt_t kdTreeBinSize = 1e6;
    public:
      
      enum class Axis {x, y, z};
      Hemisphere(Axis axis) : fAxis(axis) {;}
      Axis fAxis;
      std::shared_ptr<Delaunay2D> fDelaunay = nullptr;
      std::shared_ptr<TKDTreeID> fKDTree = nullptr;
      Points fPoints;

      bool build();
    };
    
    size_t fNumPoints = 0;
    Hemisphere fUpX;
    Hemisphere fDownX;
    Hemisphere fUpY;
    Hemisphere fDownY;
    Hemisphere fUpZ;
    Hemisphere fDownZ;
    
    bool fDoneBuild = false;    
  };

  
  
}


#endif 
