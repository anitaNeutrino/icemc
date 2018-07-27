

#ifndef ICEMC_EARTH_MODEL_H
#define ICEMC_EARTH_MODEL_H

#include "Geoid.h"

namespace G = Geoid;

namespace icemc {

  /**
   * @class EarthModel
   * @brief Interface to icemc for any class seeking to model the world.
   * 
   * In order to model neutrinos passing through the Earth and interacting,
   * you need some kind of model of the Earth.
   * This class defines the interface icemc expects any world model to have
   * in order to simulate neutrinos.
   */
  class EarthModel {

  public:
    virtual double SurfaceAboveGeoid(const G::Position& p) const = 0;
    virtual double IceThickness(const G::Position& p) const = 0;

    virtual TVector3 GetSurfaceNormal(const G::Position& p) const;

    virtual inline double Surface(const G::Position& p) const {
      return p.Surface() + SurfaceAboveGeoid(p);
    };

    
  };


}
#endif 
