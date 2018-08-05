#ifndef ICEMC_WORLD_MODEL_H
#define ICEMC_WORLD_MODEL_H

#include "Geoid.h"
#include "TVector3.h"
#include "TGraph.h"

namespace icemc {

  class Detector;

  /**
   * @class WorldModel
   * @brief Interface to icemc for any class seeking to model the world.
   * 
   * In order to model neutrinos passing through the Earth and interacting, you 
   * need some kind of model that tells you where the mass is and where the ice is.
   * 
   * This class defines the interface icemc requires any such model to have.
   * It also 
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
    
    enum {defaultTimeStep = 60};
    virtual double CreateHorizons(Detector* detector,  double horizonDistanceMeters, double timeStepSeconds = defaultTimeStep);
    

  private:
    
    enum {defaultStepSize = 1000};
    virtual double IceVolumeWithinHorizon(const Geoid::Position& p,  double horizonDistance, double stepSizeMeters = defaultStepSize) const;
    
    TGraph fHorizons;
  };

}


#endif 
