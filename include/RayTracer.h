#ifndef RAY_H_
#define RAY_H_

////////////////////////////////////////////////////////////////////////////////////////////////
//class RayTracer:
////////////////////////////////////////////////////////////////////////////////////////////////

#include "TRandom3.h"
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "anita.hh"
#include "Geoid.h"
#include "Antarctica.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


namespace icemc {

  class Settings;
  class Anita;
  class Signal;
  class EarthModel;
  class Antarctica;  

  /**
   * @class RayTracer
   * @brief Does optical ray tracing
   */
  class RayTracer {

  public:
    RayTracer(const EarthModel* ice);
    virtual ~RayTracer();

    void Initialize();

    void PrintAnglesOfIncidence() const;
    
    int GetRayIceSide(const TVector3 &n_exit2bn, const TVector3 &nsurf_rfexit, double nexit, double nenter, TVector3 &nrf2_iceside) const;
    
    int TraceRay(const Settings *settings1, Anita *anita1, int whichiteration, double n_depth);
    
    int GetSurfaceNormal(const Settings *settings1, const Antarctica *antarctica,TVector3 posnu,double &slopeyangle,int whichtry);
    
    int RandomizeSurface(const Settings *settings1,Geoid::Position rfexit_temp,TVector3 posnu, const Antarctica *antarctica,double &slopeyangle,int whichtry);

    void GetRFExit(const Settings *settings1,Anita *anita1,int whichray,Geoid::Position posnu,Geoid::Position posnu_down,Geoid::Position r_bn,Geoid::Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX],int whichtry, const Antarctica *antarctica, bool debug=false);


    
    
    // static int WhereDoesItLeaveOld(const Geoid::Position &posnu, const TVector3 &ntemp, const Antarctica *antarctica, Geoid::Position &r_out);
    static int WhereDoesItLeave(const Geoid::Position &posnu, const TVector3 &ntemp, const Antarctica *antarctica, Geoid::Position &r_out); //, bool debug = false);

    TVector3 findPathToDetector(const Geoid::Position &posnu, const Geoid::Position& detector, bool debug = false);
    static TVector3 refractiveBoundary(const TVector3& incoming, const TVector3& normal, double n_incoming, double n_outgoing, bool debug = false);
    static TCanvas* testRefractiveBoundary(double n1 = 1, double n2 = 1.31);
    
    const TVector3& getSurfaceNormalAtRFExit() const {return nsurf_rfexit;}

    TVector3 getNormalToRFIceSide(int iter = numTries - 1) const {
      return iter >= 0 && iter < numTries ? nrf_iceside[iter] : TVector3(0, 0, 0);
    }

    void initGuess(const Geoid::Position& nuInteraction, const Geoid::Position& detector);

    void setDebug(bool debug = true){
      fDebug = debug;
    }

  private:
    // const Antarctica* fAntarctica;
    const EarthModel* fAntarctica;

    static constexpr int numTries = 5;
    const TVector3 xaxis;
    const TVector3 yaxis;
    const TVector3 zaxis;

    Geoid::Position fBalloonPos;
    Geoid::Position fInteractionPos;
    TVector3 fLocalX;
    TVector3 fLocalY;
    TVector3 fLocalZ;

    mutable TVector3 fSurfaceNormal;
    mutable Geoid::Position fSurfacePos;
    mutable TVector3 fRefractedRfDir;
    mutable Geoid::Position fEndPoint;
    mutable TGraph* fMinimizerPath = nullptr;
    mutable double fBestResidual = DBL_MAX;

    ROOT::Math::Minimizer* fMinimizer = nullptr;
    ROOT::Math::Functor* fFitFunc = nullptr;

    mutable bool fDebug = false;
    bool fDoingMinimization = false;

    Geoid::Position getSurfacePosition(const double* params) const;
    double evalPath(const double* params) const;
    TVector3 toLocal(const TVector3& v, bool translate=true) const;
    void makeDebugPlots(const TString& fileName) const;

  public:
    
    TVector3 n_exit2bn[numTries]; // normal vector in direction of exit point to balloon - 5 iterations
    TVector3 nsurf_rfexit; // normal of the surface at the place where the rf leaves
    TVector3 nsurf_rfexit_db;
    TVector3 nrf_iceside[numTries];  // direction of rf [tries][3d]
    
    TVector3 nrf_iceside_eachboresight[numTries][Anita::NLAYERS_MAX][Anita::NPHI_MAX];
    TVector3 n_exit2bn_eachboresight[numTries][Anita::NLAYERS_MAX][Anita::NPHI_MAX];
    Geoid::Position rfexit_eachboresight[numTries][Anita::NLAYERS_MAX][Anita::NPHI_MAX];
    
    Geoid::Position rfexit_db[numTries];
    
    Geoid::Position rfexit[numTries]; // position where the rf exits the ice- 5 iterations
    
    double sum_slopeyness; // for summing the average slopeyness
    
    //double sum_slopeyness=0; // for summing the average slopeyness
    double slopeyx,slopeyy,slopeyz;

    void makePlot(const char* name, const Antarctica* antarctica, const Geoid::Position& nuInteraction,  const TVector3& nuDirection, const Geoid::Position& detector) const;    

  };
}

#endif
