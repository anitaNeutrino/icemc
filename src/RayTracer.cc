#include "TF1.h"
#include "TRandom3.h"

#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "Crust2.h"
#include "Antarctica.h"
#include "AskaryanRadiationModel.h"
#include "TVector3.h"
#include "Geoid.h"
#include "RayTracer.h"

#include <cmath>

#include "TCanvas.h"
#include "TPaveText.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TFile.h"
#include "TH2D.h"


icemc::RayTracer::RayTracer(std::shared_ptr<const WorldModel> world) :
  fWorld(world)
{
}

icemc::RayTracer::~RayTracer()
{
  if(fMinimizer){
    delete fMinimizer;
  }
  if(fFitFunc){
    delete fFitFunc;
  }
  if(fMinimizerPath){
    delete fMinimizerPath;
  }
}












TVector3 icemc::RayTracer::refractiveBoundary(const TVector3& incoming, const TVector3& surfaceNormal, double n_incoming, double n_outgoing, bool debug){

  // The conventions for the inputs are:
  // The incoming vector (need not be a unit) is assumed to end exactly on the surface, with orientation given by surfaceNormal (need not be unit).
  // The outoing (returned) vector will be of unit length regardless of the length of the incoming vector.
  // First: we're going to split the incoming ray into two components,
  // Parallel to the surfaceNormal and perpendicular to the surfaceNormal.

  // First, since there  are two possible directions for the surface normal,
  // let's make sure that the normal is aligned *with* the incoming ray
  const TVector3 unitNormal = surfaceNormal.Dot(incoming) > 0 ? surfaceNormal.Unit() : -surfaceNormal.Unit();

  // We only care about direction, so convert to unit vector to simplify calculations
  const TVector3 incomingUnit = incoming.Unit();

  // By definition, the dot product gives the size of the component parallel to unitNormal
  const TVector3 incomingPara = incomingUnit.Dot(unitNormal)*unitNormal;

  // and the perpendicular is whatever's not parallel
  const TVector3 incomingPerp = incomingUnit - incomingPara;

  /**
   * n_outgoing
   *                   ^
   * boundary          | surface normal
   * ----------------------------------
   * n_incoming      /|
   *                / |
   *               / *| * => theta_i
   * |incoming|=1 / \_|      (incident
   *             /    |	      angle)
   *            /     |
   *           /      |  incoming_para
   *          /       |
   *         /       _|
   *        /_______|_|
   *
   *        incoming_perp
   * 
   * Here the labels parallel/perpendicular 
   * are w.r.t. the surface NORMAL, not the boundary itself
   */

  // from the sketchy diagram above, you can hopefully see...
  const double sinThetaIncoming = incomingPerp.Mag(); // divided by incomingUnit.Mag(), which = 1;
  
  // Snell's law is that n_i sin(theta_i) = n_o sin(theta_o)
  // where _i -> incoming, _o -> outgoing
  // so sin(theta_o) = (n_i/n_o) sin(theta_i)

  const double sinThetaOutgoing = (n_incoming/n_outgoing)*sinThetaIncoming;

  // two cases here... if sinThetaOutgoing is negative or greater than 1, then it's TIR
  if(sinThetaOutgoing < 0 || sinThetaOutgoing > 1){
    // TIR, reflect the bit parallel to the surface normal
    if(debug){
      std::cout << "TIR!\n";
    }
    return (incomingPerp - incomingPara).Unit();
  }
  else{
    // refraction, bend the ray
    const double cosThetaOutgoing = TMath::Sqrt(1 - sinThetaOutgoing*sinThetaOutgoing);
    // if(debug){
    //   // std::cout << "Refract!\t" << cosThetaOutgoing << "\t" << sinThetaOutgoing << std::endl;
    //   std::cout << "Refract!\t" << TMath::ASin(sinThetaIncoming)*TMath::RadToDeg() << "\t" << TMath::ASin(sinThetaOutgoing)*TMath::RadToDeg() << std::endl;          
    // }

    // make new unit vector, which means perpendicular(parallel) components are size of sin(cos) theta_outgoing
    TVector3 refracted = cosThetaOutgoing*incomingPara.Unit() + sinThetaOutgoing*incomingPerp.Unit();
    return refracted.Unit();
  }
}



Geoid::Position icemc::RayTracer::getSurfacePosition(const double* params) const {

  // so we pick a point shifted by deltaX, and deltaY from the start position
  // Geoid::Position surfacePos = fInteractionPos + params[0]*fLocalX + params[1]*fLocalY;
  TVector3 localPos(params[0], params[1], 0);
  Geoid::Position surfacePos = fLocalCoords->localPositionToGlobal(localPos);
  double surfaceAlt = fWorld->SurfaceAboveGeoid(surfacePos);  
  surfacePos.SetAltitude(surfaceAlt); ///@todo Is this efficient yet?

  return surfacePos;
}






double icemc::RayTracer::evalPath(const double* params) const {

  // position we're aiming for from the detector...
  const Geoid::Position surfacePos = getSurfacePosition(params);

  // Gives us this initial RF direction...
  const TVector3 rfDir = (surfacePos - fDetectorPos).Unit();

  OpticalPath::Step s2; // from the surface to the balloon
  s2.start = surfacePos;  
  s2.end = fDetectorPos;
  s2.n = AskaryanRadiationModel::N_AIR;
  s2.attenuationLength = DBL_MAX; //@todo is this sensible? 

  s2.boundaryNormal = fWorld->GetSurfaceNormal(surfacePos);
  
  ///@todo get these refractive index numbers from the world model...
  const TVector3 refractedRfDir = refractiveBoundary(rfDir, s2.boundaryNormal, AskaryanRadiationModel::N_AIR, AskaryanRadiationModel::NICE, fDebug);
  const double distRemaining = (surfacePos - fInteractionPos).Mag();

  const TVector3 endPoint = surfacePos + refractedRfDir*distRemaining;

  OpticalPath::Step s1; // from the source (hopefully the end point) to the surface
  s1.start = endPoint;
  s1.end = surfacePos;
  s1.n = AskaryanRadiationModel::NICE;
  const double attenLengthIceMeters = 700; ///@todo Get this number from the world model
  s1.attenuationLength = attenLengthIceMeters;
  // order matters, think about this!
  s1.boundaryNormal = s2.boundaryNormal; ///@todo think about this...
  
  
  const TVector3 delta = (endPoint - fInteractionPos);
  double residual = delta.Mag();
  
  if(fDebug && fDoingMinimization){
    if(!fMinimizerPath){
      fMinimizerPath = new TGraph();
    }
    fMinimizerPath->SetPoint(fMinimizerPath->GetN(), params[0],  params[1]);
    std::cout << "Iter = " << fMinimizerPath->GetN()
  	      << ", dx (km) = " << 1e-3*params[0]
  	      << ", dy (km) = " << 1e-3*params[1]
  	      << ", residual (m) = " << residual
  	      // <<  ", endPoint.Mag() = " << endPoint.Mag()
  	      << "\n";
  }

  if(fDoingMinimization && residual < fBestResidual){
    fBestResidual = residual;
    fSurfaceNormal = s1.boundaryNormal;
    fSurfacePos = surfacePos;
    fEndPoint = endPoint;
    fRefractedRfDir = refractedRfDir;

    fOpticalPath.clear();

    fOpticalPath.steps.emplace_back(s1);
    fOpticalPath.steps.emplace_back(s2);
    fOpticalPath.residual = residual;    
  }
  return residual;
}



icemc::OpticalPath icemc::RayTracer::findPath(const Geoid::Position&rfStart, const Geoid::Position& rfEnd, bool debug){
  ///@todo add collision check to best fit path
  ///@todo dynamically figure out which way to fit, from lowest to highest refractive index.
  fDetectorPos = rfEnd; // RF ends up at the detector...
  fLocalCoords = std::make_shared<LocalCoordinateSystem>(fDetectorPos);
  fInteractionPos = rfStart;
  fOpticalPath.reset();
  
  
  // here we setup a new local coordinate system for the fitter to do translations in
  // The idea is to fit along the surface in terms of dx and dy

  if(!fMinimizer){
    fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
    fMinimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    // fMinimizer->SetMaxIterations(10000);  // for GSL
    fMinimizer->SetTolerance(0.01);
    fMinimizer->SetPrintLevel(fDebug ? 1 : 0);

    fFitFunc = new ROOT::Math::Functor(this, &icemc::RayTracer::evalPath, 2);
    fMinimizer->SetFunction(*fFitFunc);
  }
  
  fMinimizer->SetVariable(0, "dx_local", 0, 1);
  fMinimizer->SetVariable(1, "dy_local", 0, 1);

  fBestResidual = DBL_MAX;
  int oldErrorLevel = gErrorIgnoreLevel;
  if(fDebug==false){ // get minuit to shut the hell up
    gErrorIgnoreLevel = 1001;
  }

  fDoingMinimization = true;
  fMinimizer->Minimize();
  fDoingMinimization = false;

  gErrorIgnoreLevel = oldErrorLevel;

  if(fDebug){
    makeDebugPlots("testFitter.root");
    // exit(1);
  }

  TVector3 rfPath = (fSurfacePos - fInteractionPos).Unit();

  static int nGood = 0;
  static int nBad = 0;
  const double residualThreshold = 1;
  if(fBestResidual > residualThreshold){
    nBad++;
    if(fDebug){
      icemc::report() << severity::warning << "Outside path fitter residual tolerance!,  fBestResidual = " << fBestResidual << ", nGood = " << nGood << ", nBad = " << nBad << std::endl;
    }
    rfPath.SetXYZ(0, 0, 0);
  }
  else{
    nGood++;
    if(fDebug){
      icemc::report() << severity::info << "Successful fit! fBestResidual = " << fBestResidual << ", nGood = " << nGood << ", nBad = " << nBad << std::endl;
    }
  }

  return fOpticalPath;
}



void icemc::RayTracer::makeDebugPlots(const TString& fileName) const {
    
  TFile* fTest = new TFile(fileName, "recreate");
 
  const TVector3 toInteraction = fLocalCoords->globalPositionToLocal(fInteractionPos);
  const TVector3 toSurface = fLocalCoords->globalPositionToLocal(fSurfacePos);
  const TVector3 surfPlusRef = fLocalCoords->globalTranslationToLocal(fRefractedRfDir) + toSurface;
  const TVector3 surfPlusNorm = 1000*fLocalCoords->globalTranslationToLocal(fSurfaceNormal) + toSurface;
  std::cout << "The dot product is " << fSurfaceNormal.Unit().Dot((fSurfacePos - fDetectorPos).Unit()) << std::endl;
  // std::cout << toSurface.Unit() << "\t" << toLocal(fRefractedRfDir, false) << std::endl;    
  const TVector3 toEndPoint = fLocalCoords->globalPositionToLocal(fEndPoint);

  TGraph* grInteractionTop = new TGraph();
  grInteractionTop->SetPoint(grInteractionTop->GetN(), toInteraction.X(), toInteraction.Y());
  grInteractionTop->SetName("grInteractionTop");
  grInteractionTop->Write();
  delete grInteractionTop;

  TGraph* grInteractionSide = new TGraph();
  grInteractionSide->SetPoint(grInteractionSide->GetN(), toInteraction.X(), toInteraction.Z());
  grInteractionSide->SetName("grInteractionSide");
  grInteractionSide->Write();
  delete grInteractionSide;
    
  TGraph* grTopView = new TGraph();
  grTopView->SetPoint(0, 0, 0); // balloon at origin
  grTopView->SetPoint(1, toSurface.X(), toSurface.Y());
  grTopView->SetPoint(2, surfPlusRef.X(), surfPlusRef.Y());
  grTopView->SetPoint(3, toEndPoint.X(), toEndPoint.Y());
  grTopView->SetName("grTopView");
  grTopView->Write();
  delete grTopView;

  if(fMinimizerPath){
    fMinimizerPath->SetName("grMinuitPath");
    fMinimizerPath->Write();
  }

  TGraph* grNormalTop = new TGraph();
  grNormalTop->SetPoint(0, toSurface.X(), toSurface.Y());
  grNormalTop->SetPoint(1, surfPlusNorm.X(), surfPlusNorm.Y());
  grNormalTop->SetName("grNormalTop");
  grNormalTop->Write();
  delete grNormalTop;

  TGraph* grNormalSide = new TGraph();
  grNormalSide->SetPoint(0, toSurface.X(), toSurface.Z());
  grNormalSide->SetPoint(1, surfPlusNorm.X(), surfPlusNorm.Z());
  grNormalSide->SetName("grNormalSide");
  grNormalSide->Write();
  delete grNormalSide;

  TGraph* grSideView = new TGraph();
  grSideView->SetPoint(0, 0, 0); // balloon at origin
  grSideView->SetPoint(1, toSurface.X(), toSurface.Z());
  grSideView->SetPoint(2, surfPlusRef.X(), surfPlusRef.Z());    
  grSideView->SetPoint(3, toEndPoint.X(), toEndPoint.Z());
  grSideView->SetName("grSideView");
  grSideView->Write();
  delete grSideView;

  TGraph* grSurface = new TGraph();
  grSurface->SetName("grSurface");
  const int nPoints = 1000000;
  const int nExtra = nPoints/10;
  const double localXDistToInteraction = toInteraction.X();
  const double deltaX = localXDistToInteraction/nPoints;

  for(int d=-nExtra; d < nPoints+nExtra; d++){
    const std::array<double, 2> params {deltaX*d, 0};
    const TVector3 v = getSurfacePosition(params.data());
    const TVector3 vl = fLocalCoords->globalPositionToLocal(v);
    grSurface->SetPoint(grSurface->GetN(), vl.X(), vl.Z());
  }
  grSurface->Write();    
  delete grSurface;

  bool oldDebug = fDebug;
  fDebug = false; // this will be noisy otherwise

  // we picked our local coordinate system so that ANITA is along +x somewhere
  // so let's pick a reasonable range ...
  double x0 = fMinimizer->X()[0];
  double y0 = fMinimizer->X()[1];
    
  const double bigPlotDist = 800e3;
  const double smallPlotDist = 100;
  const int nBins = 1024;
  const std::vector<double> xMins {x0 - smallPlotDist, 0};
  const std::vector<double> xMaxs {x0 + smallPlotDist, bigPlotDist};    
  const std::vector<double> yMins {y0 - smallPlotDist, -0.5*bigPlotDist};
  const std::vector<double> yMaxs {y0 + smallPlotDist, +0.5*bigPlotDist};
    
  for(int rangeInd=0; rangeInd < xMins.size(); rangeInd++){

    TString name = TString::Format("hFitterParameterSpace_%d", rangeInd);
    TString title = "How far are we from the interaction point as a function of ray tracing at this surface position; x_{local} (meters); y_{local} (meters); Residual distance to detector (m)";
    TH2D* hFitterParameterSpace = new TH2D(name, title,
					   nBins, xMins.at(rangeInd), xMaxs.at(rangeInd),
					   nBins, yMins.at(rangeInd), yMaxs.at(rangeInd));
    for(int by=1; by < hFitterParameterSpace->GetNbinsY(); by++){
      double dy = hFitterParameterSpace->GetYaxis()->GetBinCenter(by);
      for(int bx=1; bx < hFitterParameterSpace->GetNbinsX(); bx++){
	double dx = hFitterParameterSpace->GetXaxis()->GetBinCenter(bx);
	double params[2] = {dx, dy};
	hFitterParameterSpace->SetBinContent(bx, by, evalPath(params));
      }
    }
    hFitterParameterSpace->Write();
    delete hFitterParameterSpace;
  }
    
  fTest->Write();
  fTest->Close();
  fDebug = oldDebug;
}












TCanvas* icemc::RayTracer::testRefractiveBoundary(double n1, double n2) {

  TGraph* gr = new TGraph();
  gr->SetBit(kCanDelete);
  gr->SetTitle("Angles relative to surface normal; Angle of incidence (Degrees);Angle of refraction (Degrees)");
  
  TCanvas* c = new TCanvas();
  c->Divide(2);
  c->cd(1);
  TGraph* grBoundary = new TGraph();
  grBoundary->SetTitle("Ray tracing");
  grBoundary->SetPoint(0, -1, 0);
  grBoundary->SetPoint(1, 1, 0);
  grBoundary->Draw("al");
  grBoundary->SetLineStyle(2);
  grBoundary->SetMaximum(1);
  grBoundary->SetMinimum(-1);
  grBoundary->SetBit(kCanDelete);

  TBox* bIncoming = new TBox(-1, -1, 1, 0);
  EColor b1 = n1 > n2 ? kCyan : kWhite;
  bIncoming->SetFillColorAlpha(b1, 0.5);
  bIncoming->SetFillStyle(3345);
  bIncoming->Draw("same");
  bIncoming->SetBit(kCanDelete);

  TBox* bOutgoing = new TBox(-1, 0, 1, 1);
  EColor b2 = n2 > n1 ? kCyan : kWhite;
  bOutgoing->SetFillColorAlpha(b2, 0.5);
  bOutgoing->SetFillStyle(3354);  
  bOutgoing->Draw("same");
  bOutgoing->SetBit(kCanDelete);

  TPaveText* t1 = new TPaveText(-1, -1, -0.6, -0.8);
  t1->AddText("Incoming");
  t1->AddText(TString::Format("n=%4.2lf", n1));
  t1->SetBorderSize(0);
  t1->SetFillColor(b1);
  t1->SetBit(kCanDelete);
  t1->Draw();
  
  TPaveText* t2 = new TPaveText(-1, 0.8, -0.6, 1);
  t2->AddText("Outgoing");
  t2->AddText(TString::Format("n=%4.2lf", n2));
  t2->SetBorderSize(0);
  t2->SetFillColor(b2);
  t2->SetBit(kCanDelete);  
  t2->Draw();

  const TVector3 normal(0, 0, 1);
  for(int theta_i_Deg = 1; theta_i_Deg <= 90 ; theta_i_Deg++){

    double thetaRad = TMath::DegToRad()*theta_i_Deg;
    double x = sin(thetaRad);
    double y = 0;
    double z = cos(thetaRad);
    const TVector3 incoming(x, y, z);

    TGraph* grPath = new TGraph();
    grPath->SetPoint(0, -x, -z);
    grPath->SetPoint(1, 0, 0);
    
    TVector3 outgoing = refractiveBoundary(incoming,  normal, n1, n2);    
    grPath->SetPoint(2, outgoing.X(), outgoing.Z());

    double cosThetaOut = outgoing.Dot(normal);
    gr->SetPoint(gr->GetN(), theta_i_Deg, TMath::ACos(cosThetaOut)*TMath::RadToDeg());
    grPath->SetLineColor(theta_i_Deg + 1);
    
    grPath->Draw("lsame");

    grPath->SetBit(kCanDelete);
  }

  c->cd(2);
  gr->Draw("al");
  
  return c;
}
