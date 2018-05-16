#ifndef RAY_H_
#define RAY_H_

////////////////////////////////////////////////////////////////////////////////////////////////
//class RayTracer:
////////////////////////////////////////////////////////////////////////////////////////////////

#include "TRandom3.h"
#include <iostream>
#include <cmath>

#include "vector.hh"
#include "anita.hh"
#include "position.hh"
#include "icemodel.hh"


namespace icemc {

  class Settings;
  class Anita;
  class Signal;

  /**
   * @class RayTracer
   * @brief Does optical ray tracing
   */
  class RayTracer {
    
  public:
    RayTracer();
    void Initialize();

    Vector n_exit2bn[5]; // normal vector in direction of exit point to balloon - 5 iterations, 3 directions for each 
    Vector nsurf_rfexit; // normal of the surface at the place where the rf leaves
    Vector nsurf_rfexit_db;
    Vector nrf_iceside[5];  // direction of rf [tries][3d]
    
    Vector nrf_iceside_eachboresight[5][Anita::NLAYERS_MAX][Anita::NPHI_MAX];  // direction of rf [tries][3d]
    Vector n_exit2bn_eachboresight[5][Anita::NLAYERS_MAX][Anita::NPHI_MAX];
    Position rfexit_eachboresight[5][Anita::NLAYERS_MAX][Anita::NPHI_MAX];
    
    Position rfexit_db[5];
    Position rfexit[5]; // position where the rf exits the ice- 5 iterations, 3 dimensions each
    double sum_slopeyness; // for summing the average slopeyness
    //double sum_slopeyness=0; // for summing the average slopeyness
    double slopeyx,slopeyy,slopeyz;    
    void PrintAnglesofIncidence();    
    int GetRayIceSide(const Vector &n_exit2bn, const Vector &nsurf_rfexit, double nexit, double nenter, Vector &nrf2_iceside);
    int TraceRay(const Settings *settings1,Anita *anita1,int whichiteration,double n_depth);
    int GetSurfaceNormal(const Settings *settings1, const IceModel *antarctica,Vector posnu,double &slopeyangle,int whichtry);
    int RandomizeSurface(const Settings *settings1,Position rfexit_temp,Vector posnu, const IceModel *antarctica,double &slopeyangle,int whichtry);
    void GetRFExit(const Settings *settings1,Anita *anita1,int whichray,Position posnu,Position posnu_down,Position r_bn,Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX],int whichtry, const IceModel *antarctica);
    static int WhereDoesItLeave(const Position &posnu, const Vector &ntemp, const IceModel *antarctica, Position &r_out);

    Vector xaxis;
    Vector yaxis;
    Vector zaxis;
  };

}

#endif
