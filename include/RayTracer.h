#ifndef RAY_H_
#define RAY_H_

#include "TVector3.h"
#include "Geoid.h"
#include "Antarctica.h"
#include "LocalCoordinateSystem.h"
#include "OpticalPath.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


namespace icemc {

  class WorldModel;


  
  /**
   * @class RayTracer
   * @brief Does optical ray tracing
   */
  class RayTracer {

  public:
    RayTracer(std::shared_ptr<const WorldModel> world);
    virtual ~RayTracer();
    OpticalPath findPath(const Geoid::Position &interaction, const Geoid::Position& detector, bool debug = false);

    static TVector3 refractiveBoundary(const TVector3& incoming, const TVector3& surfaceNormal, double n_incoming, double n_outgoing, bool debug=false);
    // static TVector3 refractiveBoundaryPol(const TVector3& incoming, const TVector3& surfaceNormal, double n_incoming, double n_outgoing);

    void setDebug(bool debug = true){
      fDebug = debug;
    }

  private:
    std::shared_ptr<const WorldModel> fWorld;

    Geoid::Position fDetectorPos;
    std::shared_ptr<LocalCoordinateSystem> fLocalCoords = nullptr;
    
    Geoid::Position fInteractionPos;

    mutable TVector3 fSurfaceNormal;
    mutable Geoid::Position fSurfacePos;
    mutable TVector3 fRefractedRfDir;
    mutable Geoid::Position fEndPoint;
    mutable TGraph* fMinimizerPath = nullptr;
    mutable double fBestResidual = DBL_MAX;
    mutable OpticalPath fOpticalPath;

    ROOT::Math::Minimizer* fMinimizer = nullptr;
    ROOT::Math::Functor* fFitFunc = nullptr;

    mutable bool fDebug = false;
    bool fDoingMinimization = false;

    Geoid::Position getSurfacePosition(const double* params) const;
    double evalPath(const double* params) const;
    void makeDebugPlots(const TString& fileName) const;

  public:
    static TCanvas* testRefractiveBoundary(double n1, double n2);

  };
}

#endif
