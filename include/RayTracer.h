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
#include "Antarctica.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


namespace icemc {

  class Settings;
  class Anita;
  class Signal;
  class Antarctica;

  /**
   * @class RayTracer
   * @brief Does optical ray tracing
   */
  class RayTracer {

  public:
    RayTracer(const Antarctica* ice);
    virtual ~RayTracer();

    void Initialize();

    void PrintAnglesOfIncidence() const;
    
    int GetRayIceSide(const Vector &n_exit2bn, const Vector &nsurf_rfexit, double nexit, double nenter, Vector &nrf2_iceside) const;
    
    int TraceRay(const Settings *settings1, Anita *anita1, int whichiteration, double n_depth);
    
    int GetSurfaceNormal(const Settings *settings1, const Antarctica *antarctica,Vector posnu,double &slopeyangle,int whichtry);
    
    int RandomizeSurface(const Settings *settings1,Position rfexit_temp,Vector posnu, const Antarctica *antarctica,double &slopeyangle,int whichtry);

    void GetRFExit(const Settings *settings1,Anita *anita1,int whichray,Position posnu,Position posnu_down,Position r_bn,Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX],int whichtry, const Antarctica *antarctica, bool debug=false);


    
    
    // static int WhereDoesItLeaveOld(const Position &posnu, const Vector &ntemp, const Antarctica *antarctica, Position &r_out);
    static int WhereDoesItLeave(const Position &posnu, const Vector &ntemp, const Antarctica *antarctica, Position &r_out); //, bool debug = false);

    Vector findPathToDetector(const Position &posnu, const Position& detector, bool debug = false);
    static Vector refractiveBoundary(const Vector& incoming, const Vector& normal, double n_incoming, double n_outgoing, bool debug = false);
    static TCanvas* testRefractiveBoundary(double n1 = 1, double n2 = 1.31);
    
    const Vector& getSurfaceNormalAtRFExit() const {return nsurf_rfexit;}

    Vector getNormalToRFIceSide(int iter = numTries - 1) const {
      return iter >= 0 && iter < numTries ? nrf_iceside[iter] : Vector(0, 0, 0);
    }

    void initGuess(const Position& nuInteraction, const Position& detector);

    void setDebug(bool debug = true){
      fDebug = debug;
    }

  private:
    const Antarctica* fAntarctica;    

    static constexpr int numTries = 5;
    const Vector xaxis;
    const Vector yaxis;
    const Vector zaxis;

    Position fBalloonPos;
    Position fInteractionPos;
    Vector fLocalX;
    Vector fLocalY;
    Vector fLocalZ;

    mutable Vector fSurfaceNormal;
    mutable Vector fSurfacePos;
    mutable Vector fRefractedRfDir;
    mutable Vector fEndPoint;
    mutable TGraph* fMinimizerPath = nullptr;
    mutable double fBestResidual = DBL_MAX;

    ROOT::Math::Minimizer* fMinimizer = nullptr;
    ROOT::Math::Functor* fFitFunc = nullptr;

    mutable bool fDebug = false;
    bool fDoingMinimization = false;

    Position getSurfacePosition(const double* params) const;
    double evalPath(const double* params) const;
    Vector toLocal(const Vector& v, bool translate=true) const;
    void makeDebugPlots(const TString& fileName) const;

  public:
    
    Vector n_exit2bn[numTries]; // normal vector in direction of exit point to balloon - 5 iterations
    Vector nsurf_rfexit; // normal of the surface at the place where the rf leaves
    Vector nsurf_rfexit_db;
    Vector nrf_iceside[numTries];  // direction of rf [tries][3d]
    
    Vector nrf_iceside_eachboresight[numTries][Anita::NLAYERS_MAX][Anita::NPHI_MAX];
    Vector n_exit2bn_eachboresight[numTries][Anita::NLAYERS_MAX][Anita::NPHI_MAX];
    Position rfexit_eachboresight[numTries][Anita::NLAYERS_MAX][Anita::NPHI_MAX];
    
    Position rfexit_db[numTries];
    
    Position rfexit[numTries]; // position where the rf exits the ice- 5 iterations
    
    double sum_slopeyness; // for summing the average slopeyness
    
    //double sum_slopeyness=0; // for summing the average slopeyness
    double slopeyx,slopeyy,slopeyz;

    void makePlot(const char* name, const Antarctica* antarctica, const Position& nuInteraction,  const Vector& nuDirection, const Position& detector) const;    

  };
}

#endif
