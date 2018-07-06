#include "TF1.h"
#include "TRandom3.h"

#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "Earth.h"
#include "Antarctica.h"
#include "AskaryanFreqsGenerator.h"
#include "vector.hh"
#include "position.hh"
#include "anita.hh"
#include "RayTracer.h"

#include <cmath>

#include "TCanvas.h"
#include "TPaveText.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TFile.h"
#include "TH2D.h"


icemc::RayTracer::RayTracer(const Antarctica* ice) :
  fAntarctica(ice),
  xaxis(1.,0.,0.),
  yaxis(0.,1.,0.),
  zaxis(0.,0.,1.)
{
  sum_slopeyness=0.;  
  slopeyx=0.;
  slopeyy=0.;
  slopeyz=0.; // these are a measure of how much the surface is sloped in the x,y and z directions
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



void icemc::RayTracer::PrintAnglesOfIncidence() const {

  std::cout << "angle of incidence (firn-air) is " << nrf_iceside[3].Angle(nsurf_rfexit)*constants::DEGRAD << "\n";
  std::cout << "angle of incidence (ice-firn) is " << nrf_iceside[4].Angle(nsurf_rfexit)*constants::DEGRAD << "\n";
}



void icemc::RayTracer::initGuess(const Position& nuInteraction, const Position& detector) {  

  // to seed the ray tracing, we come up with an obviously wrong initial guess
  // first, we assume that the direction of the RF in the ice goes vertically upwards...
  nrf_iceside[0] = nuInteraction.Unit();

  //.. so the RF exit point is directly above the interaction...
  rfexit[0] = fAntarctica->Surface(nuInteraction) * nuInteraction.Unit();

  // ... and the ray then changes direction on its exit from the ice, to go to the detector
  n_exit2bn[0] = (detector - rfexit[0]).Unit();

  // Reset iterative solver solutions ///@todo reset everything? 
  for(int i=1; i < numTries; i++){
    rfexit[i].SetXYZ(0,0,0);
    n_exit2bn[i].SetXYZ(0,0,0); 
  }

  ///@todo boresight stuff from Balloon::PickDownwardInteractionPoint
  // // do the calculation for the boresights?
  // if (settings1->BORESIGHTS) {
  //   for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
  //     for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
  // 	ray1->rfexit_eachboresight[0][ilayer][ifold] = antarctica1->Surface(interaction1->posnu) * interaction1->posnu.Unit();// this first guess rfexit is the same for all antennas too
  // 	ray1->n_exit2bn_eachboresight[0][ilayer][ifold] = (r_boresights[ilayer][ifold]- ray1->rfexit_eachboresight[0][ilayer][ifold]).Unit();
  // 	// std::cout << "ilayer, ifold, n_exit2bn are " << ilayer << "\t" << ifold << " ";
  //     }
  //   }
  // }
  
  // if (settings1->BORESIGHTS) {
  //   // this is the same for all of the antennas too
  //   for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) { // loop over layers on the payload
  //     for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
  // 	ray1->nrf_iceside_eachboresight[0][ilayer][ifold] = interaction1->posnu.Unit();
  //     } // end loop over fold
  //   } // end loop over payload layers
  // } // end if boresights
  
  
}





void icemc::RayTracer::Initialize() {
  
  for (int i = 0; i < numTries; i++) {
    n_exit2bn[i]   = Vector(0.,0.,0.); // normal vector in direction of exit point to balloon - 5 iterations, 3 directions for eac
    nrf_iceside[i] = Vector(0.,0.,0.);  // direction of rf [tries][3d]
    rfexit_db[i]   = Vector(0.,0.,0.);
    rfexit[i]      = Vector(0.,0.,0.); // position where the rf exits the ice- 5 iterations, 3 dimensions eac
    
    for (int j=0;j<Anita::NPHI_MAX;j++) {
      for (int k=0;k<Anita::NLAYERS_MAX;k++) {
	n_exit2bn_eachboresight[i][k][j]=Vector(0.,0.,0.); // normal vector in direction of exit point to each antenna boresight - 5 iterations
	nrf_iceside_eachboresight[i][k][j]=Vector(0.,0.,0.);  // direction of rf [tries][3d]
  	
	rfexit_eachboresight[i][k][j]=Vector(0.,0.,0.);
      }
    }
    
  }

  nsurf_rfexit = Vector(0.,0.,0.); // normal of the surface at the place where the rf leaves
  nsurf_rfexit_db = Vector(0.,0.,0.);
}


void icemc::RayTracer::GetRFExit(const Settings *settings1, Anita *anita1, int whichray, Position posnu, Position posnu_down, Position r_bn,
				 Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX], int whichtry, const Antarctica *antarctica, bool debug){



  if (whichray==0){
    WhereDoesItLeave(posnu, nrf_iceside[2*whichtry],antarctica, // inputs
		     rfexit[whichtry]); // output
  }  
  
  //******wufan******
  if (whichray==1){ // reflected rays
    WhereDoesItLeave(posnu_down,nrf_iceside[2*whichtry],antarctica,rfexit[whichtry]);  //from mirror point position and the direction of signals to the ice 
  //to find the exit point at the surface of the Earth.wufan 
  }

  n_exit2bn[whichtry] = (r_bn - rfexit[whichtry]).Unit();
  
  if (settings1->BORESIGHTS) { // now find rfexit and n_exit2bn for each boresight, still first iteration
    //std::cout << "first iteration.\n";
    for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
      for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
        WhereDoesItLeave(posnu,nrf_iceside_eachboresight[2*whichtry][ilayer][ifold],antarctica,
  			 rfexit_eachboresight[whichtry][ilayer][ifold]);
        n_exit2bn_eachboresight[whichtry][ilayer][ifold] = (r_boresights[ilayer][ifold] - rfexit_eachboresight[whichtry][ilayer][ifold]).Unit(); 
      }
    }
  } // end if we are doing this for each boresight
  
  // if (settings1->SLAC) {
  //   // ray comes out a little earlier because of the slope of the surface.
  //   // use law of cosines the get how much "distance" should be cut short
    
  //   double x=sin(settings1->SLACSLOPE*constants::RADDEG)*(settings1->SLACICELENGTH/2.+rfexit[0].Distance(rfexit[whichtry]))/sin(constants::PI/2.+acos(nrf_iceside[2*whichtry].Dot(nsurf_rfexit)));
    
  //   rfexit[whichtry]-=x*nrf_iceside[2*whichtry];
    
  //   if (settings1->BORESIGHTS) {
  //     for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
  //       for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {

  //         x = sin(settings1->SLACSLOPE*constants::RADDEG)*(settings1->SLACICELENGTH/2.+rfexit_eachboresight[0][ilayer][ifold].Distance(rfexit_eachboresight[whichtry][ilayer][ifold]))/sin(constants::PI/2.+acos(nrf_iceside_eachboresight[2*whichtry][ilayer][ifold].Dot(nsurf_rfexit)));

  //         rfexit_eachboresight[whichtry][ilayer][ifold]-=x*nrf_iceside_eachboresight[2*whichtry][ilayer][ifold];
          
  //         n_exit2bn_eachboresight[whichtry][ilayer][ifold] = (r_boresights[ilayer][ifold] - rfexit_eachboresight[whichtry][ilayer][ifold]).Unit();
  //       } // end loop over antennas on a layer
  //     } // end loop over layers of the payload
  //   } // end if we're keeping track of all antenna boresights
  // } // end if we're modeling the slac run
}

//###########
// icemc::RayTracer::WhereDoesItLeave() is defined in ray.hh since it is a statis member function // MS 2/1/2017


int icemc::RayTracer::RandomizeSurface(const Settings *settings1, Position rfexit_temp, Vector posnu, const Antarctica *antarctica, double &slopeyangle, int whichtry){

  double howmuch=settings1->SLOPEYSIZE;
  Position nsurf_rfexit_temp;
  
  // rotate the surface normal according to the local slope.
  if (!settings1->SLAC){
    nsurf_rfexit_temp = antarctica->GetSurfaceNormal(rfexit_temp); // find the normal to the surface taking into account the tilt from the differential heights between neighboring bins
  }
  else if (settings1->SLAC) { // if we are simulating the slac test, rotate the surface by 10 degrees away from the south pole
    Vector zaxis(0.,0.,1.);
    nsurf_rfexit_temp=(rfexit_temp.Unit()).Rotate(-settings1->SLACSLOPE*constants::RADDEG,posnu.Cross(zaxis));
  }
  
  Position nsurf_rfexit_temp_copy=nsurf_rfexit_temp;
  // tilt local surface based on slopeyness.
  if (settings1->SLOPEY) {
    
    // randomizing surface direction
    
    // 0.10=5.4, 0.2=7.4 deg mean
    
    // 02/06/07: 0.012=0.84, 0.10=6.78,3.76 deg mean
    
    double slopeyness=0;
    if (whichtry==0) { // only reset the surface slopeyness for the first try an then repeat the same for each subsequent try
      slopeyx=howmuch*gRandom->Gaus();
      slopeyy=howmuch*gRandom->Gaus();
      slopeyz=howmuch*gRandom->Gaus();
    }
    Vector ntemp2 = nsurf_rfexit_temp + slopeyx * xaxis
      + slopeyy * yaxis
      + slopeyz * zaxis;
    
      
    ntemp2 = ntemp2 / ntemp2.Mag();
    
    double rtemp= ntemp2.Dot(nsurf_rfexit_temp);
    
    if (rtemp<=1) {
      slopeyness=acos(rtemp);
      sum_slopeyness+=slopeyness;
      nsurf_rfexit_temp = ntemp2;
    }//if
    else{
      slopeyness=-999;
    }
    if (settings1->DEBUG) {
      rtemp=nsurf_rfexit_temp.Mag(); 
      if (fabs(rtemp-1.0)>=0.01) {
	std::cout << "error: nx,ny,nz_surf="<<nsurf_rfexit_temp<<std::endl;
	return 0;
      }//if
    }//if

    // to avoid FP errors, fudge the rare case of going out the poles
    if (fabs(nsurf_rfexit_temp[2])>=0.99999) {  // a hack
      nsurf_rfexit_temp = Vector(sqrt(1.-0.99999*0.99999-nsurf_rfexit_temp[1]*nsurf_rfexit_temp[1]), nsurf_rfexit_temp[1], 0.99999);
    }//if
  }
  slopeyangle=constants::DEGRAD*acos(nsurf_rfexit_temp.Dot(nsurf_rfexit_temp_copy));
  nsurf_rfexit = nsurf_rfexit_temp;
  return 1;
  
}//RandomizeSurface


// int icemc::RayTracer::GetSurfaceNormal(Antarctica *antarctica,Vector posnu,Position *rfexit) {
int icemc::RayTracer::GetSurfaceNormal(const Settings *settings1,const Antarctica *antarctica,Vector posnu,double &slopeyangle,int whichtry){

  Position rfexit_temp;

  rfexit_temp=rfexit[whichtry];

  if (!RandomizeSurface(settings1,rfexit_temp,posnu,antarctica,slopeyangle,whichtry)) {
    return 0;  //can happen only in debug mode
  }
  return 1;
}




int icemc::RayTracer::TraceRay(const Settings *settings1,Anita *anita1,int iter,double n_depth) { // iter is which iteration (1 or 2)

   
  // use snell's law to get the first guess at the 
  // direction of the rf as it leaves ice surface.
  // 0th guess was simply radially outward from interaction position
  // this now takes into account balloon position and surface normal.
    
  if (settings1->FIRN) {

    
    if (!GetRayIceSide(n_exit2bn[iter-1],nsurf_rfexit,AskaryanFreqsGenerator::N_AIR,constants::NFIRN,
		       nrf_iceside[2*iter-1])) { // nrf_iceside[1] is the rf direction in the firn

      return 0; // reject if TIR.
    } 
    // This could be throwing away events that the final guess would have kept
    //  Need to look into this and update in the future
      
    // use snell's law again to get direction of rf ray emitted from
    // interaction point.  Physically, it is a continuous bending
    // through the firn but it turns out it obeys Snell's law (Seckel)
            
    if (!GetRayIceSide(nrf_iceside[2*iter-1],nsurf_rfexit,constants::NFIRN,n_depth,
		       nrf_iceside[2*iter])) { // nrf_iceside[2] is the rf direction in the ice

      return 0;   // note to self:  need to check whether rays are totally internally reflected within ice   
    }




    // find where the refracted ray leaves.
    // rejected if the rf leaves beyond the boundary of the 
    // continent. Should never happen.
      
    if (settings1->BORESIGHTS) {
	
      for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
	for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
	    
	  if (!GetRayIceSide(n_exit2bn_eachboresight[iter-1][ilayer][ifold],nsurf_rfexit,AskaryanFreqsGenerator::N_AIR,constants::NFIRN,
			     nrf_iceside_eachboresight[2*iter-1][ilayer][ifold])) // nrf_iceside[1] is the rf direction in the firn
	    return 0; // reject if TIR.


	  //	    std::cout << "ITER " << iter-1 << " " << (n_exit2bn_eachboresight[iter-1][ilayer][ifold]) << std::endl;

	  // This could be throwing away events that the final guess would have kept
	  //  Need to look into this and update in the future
	    
	  // use snell's law again to get direction of rf ray emitted from
	  // interaction point.  Physically, it is a continuous bending
	  // through the firn but it turns out it obeys Snell's law (Seckel)
	    
	    
	  if (!GetRayIceSide(nrf_iceside_eachboresight[2*iter-1][ilayer][ifold],nsurf_rfexit,constants::NFIRN,n_depth,
			     nrf_iceside_eachboresight[2*iter][ilayer][ifold])) // nrf_iceside[2] is the rf direction in the ice
	    return 0;   // note to self:  need to check whether rays are totally internally reflected within ice   
	    
	  // find where the refracted ray leaves.
	  // rejected if the rf leaves beyond the boundary of the 
	  // continent. Should never happen.
	    
	    
	} // end phi foldings
      } // end layers
    } // end if doing this for all boresights
  } // end if we are modeling the firn
    
  else { // no firn

    if (!GetRayIceSide(n_exit2bn[iter-1],nsurf_rfexit,AskaryanFreqsGenerator::N_AIR,AskaryanFreqsGenerator::NICE,
		       nrf_iceside[2*iter-1])) // nrf_iceside[1] is the rf direction in the ice
      return 0; // reject if TIR.  
    nrf_iceside[2*iter]=nrf_iceside[2*iter-1]; // no firn so the next element is the same
      
    if (settings1->BORESIGHTS) {

      for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
	for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
	    
	  if (!GetRayIceSide(n_exit2bn_eachboresight[iter-1][ilayer][ifold],nsurf_rfexit,AskaryanFreqsGenerator::N_AIR,AskaryanFreqsGenerator::NICE,
			     nrf_iceside_eachboresight[2*iter-1][ilayer][ifold])){ // nrf_iceside[1] is the rf direction in the ice
	    return 0; // reject if TIR.
	  }	  
	  nrf_iceside_eachboresight[2*iter][ilayer][ifold]=nrf_iceside_eachboresight[2*iter-1][ilayer][ifold]; // no firn so the next element is the same

	} // number of foldings in phi
      } // number of layers
    } // end if doing this for all boresights
      
      
      
  } // end if we are not modeling the firn
    
    

  return 1;
}


int icemc::RayTracer::GetRayIceSide(const Vector &n_exit2bn,
				    const Vector &nsurf_rfexit,
				    double nexit,
				    double nenter, 
				    Vector &nrf2_iceside) const {

  // this function performs snell's law in three dimensions

  double NRATIO=nexit/nenter;

  //  std::cout << " INPUT " << n_exit2bn << " " << nsurf_rfexit << " " << nexit << " " << nenter << std::endl;
  
  double costh=(n_exit2bn.Dot(nsurf_rfexit))/(n_exit2bn.Mag() * nsurf_rfexit.Mag()); // cos(theta) of the transmission angle

  if (costh<0) {
    //    std::cout << "returning 0.  inu is " << inu << "\n";
    return 0;
  
  }
  double sinth=sqrt(1 - costh*costh);
  
  //  double sinth_i=nexit*sinth/nenter;
  
  //  double angle=asin(sinth)-asin(sinth_i);
  
  
  // std::cout << "sinth, sinth_i are " << sinth << " " << sinth_i << "\n";
  // std::cout << "th, th_i are " << asin(sinth)*DEGRAD << " " << asin(sinth_i)*DEGRAD << "\n";
  // std::cout << "angle is " << angle*DEGRAD << "\n";
  
  //  nrf2_iceside=n_exit2rx.Rotate(angle,n_exit2rx.Cross(nsurf_rfexit));


  double factor=NRATIO*costh-sqrt(1-(NRATIO*sinth*NRATIO*sinth));


  nrf2_iceside = -factor*nsurf_rfexit + NRATIO*n_exit2bn;
  nrf2_iceside = nrf2_iceside.Unit(); // normalize
  
  
    
  

  return 1;
}//GetRayIceSide











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

  const Vector normal(0, 0, 1);
  for(int theta_i_Deg = 1; theta_i_Deg <= 90 ; theta_i_Deg++){

    double thetaRad = TMath::DegToRad()*theta_i_Deg;
    double x = sin(thetaRad);
    double y = 0;
    double z = cos(thetaRad);
    const Vector incoming(x, y, z);

    TGraph* grPath = new TGraph();
    grPath->SetPoint(0, -x, -z);
    grPath->SetPoint(1, 0, 0);
    
    Vector outgoing = refractiveBoundary(incoming,  normal, n1, n2);    
    grPath->SetPoint(2, outgoing.GetX(), outgoing.GetZ());

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



icemc::Vector icemc::RayTracer::refractiveBoundary(const Vector& incoming, const Vector& surfaceNormal, double n_incoming, double n_outgoing, bool debug){

  // The conventions for the inputs are:
  // The incoming vector (need not be a unit) is assumed to end exactly on the surface, with orientation given by surfaceNormal (need not be unit).
  // The outoing (returned) vector will be of unit length regardless of the length of the incoming vector.
  // First: we're going to split the incoming ray into two components,
  // Parallel to the surfaceNormal and perpendicular to the surfaceNormal.

  // First, since there  are two possible directions for the surface normal,
  // let's make sure that the normal is aligned *with* the incoming ray
  const Vector unitNormal = surfaceNormal.Dot(incoming) > 0 ? surfaceNormal.Unit() : -surfaceNormal.Unit();

  // We only care about direction, so convert to unit vector to simplify calculations
  const Vector incomingUnit = incoming.Unit();

  // By definition, the dot product gives the size of the component parallel to unitNormal
  const Vector incomingPara = incomingUnit.Dot(unitNormal)*unitNormal;

  // and the perpendicular is whatever's not parallel
  const Vector incomingPerp = incomingUnit - incomingPara;

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
    Vector refracted = cosThetaOutgoing*incomingPara.Unit() + sinThetaOutgoing*incomingPerp.Unit();
    return refracted.Unit();
  }
}


icemc::Vector icemc::RayTracer::toLocal(const Vector& v, bool translate) const {
  // translate the earth centered vector, v, to our new local coordinate system
  int translationFactor = translate ? 1 : 0;
  const Vector relativeToNewOrigin = v - translationFactor*fBalloonPos;
  double x = relativeToNewOrigin.Dot(fLocalX);
  double y = relativeToNewOrigin.Dot(fLocalY);
  double z = relativeToNewOrigin.Dot(fLocalZ);
  return Vector(x, y, z);
}

// icemc::Position icemc::RayTracer::nudgeSurface(const double* params) const {
icemc::Position icemc::RayTracer::getSurfacePosition(const double* params) const {

  // so we pick a point shifted by deltaX, and deltaY from the start position
  // Position surfacePos = fInteractionPos + params[0]*fLocalX + params[1]*fLocalY;
  Position surfacePos = fBalloonPos + params[0]*fLocalX + params[1]*fLocalY;
  double lon = surfacePos.Lon();
  double lat = surfacePos.Lat();
  surfacePos.SetLonLatAlt(lon, lat, fAntarctica->Surface(lon, lat));

  return surfacePos;
}




double icemc::RayTracer::evalPath(const double* params) const {

  // position we're aiming for from the detector...
  const Position surfacePos = getSurfacePosition(params);

  // Gives us this initial RF direction...
  const Vector rfDir = (surfacePos - fBalloonPos).Unit();

  const Vector surfaceNormal = fAntarctica->GetSurfaceNormal(surfacePos);
  // const Vector surfaceNormal = surfacePos.Unit();

  const Vector refractedRfDir = refractiveBoundary(rfDir, surfaceNormal, AskaryanFreqsGenerator::N_AIR, AskaryanFreqsGenerator::NICE, fDebug);
  const double dist = (surfacePos - fInteractionPos).Mag();

  const Vector endPoint = surfacePos + refractedRfDir*dist;

  const Vector delta = (endPoint - fInteractionPos);
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
    fSurfaceNormal = surfaceNormal;
    fSurfacePos = surfacePos;
    fEndPoint = endPoint;
    fRefractedRfDir = refractedRfDir;
  }
  return residual;
}




icemc::Vector icemc::RayTracer::findPathToDetector(const Position& rfStart, const Position& balloon, bool debug){
  
  fInteractionPos = rfStart;
  fBalloonPos = balloon;

  // here we setup a new local coordinate system for the fitter to do translations in
  // The idea is to fit along the surface in terms of dx and dy
  
  fLocalZ = fBalloonPos.Unit(); // z is directly up from rfStart
  fLocalY = fLocalZ.Cross(fInteractionPos - fBalloonPos).Unit(); // y is perpendicual to local up AND the balloon-interaction path
  // fLocalZ = fInteractionPos.Unit(); // z is directly up from rfStart
  // fLocalY = fLocalZ.Cross(fBalloonPos - fInteractionPos).Unit(); // y is perpendicual to local up AND the balloon-interaction path
  fLocalX = fLocalY.Cross(fLocalZ).Unit(); // ... and therefore x is local horizontal in direction of balloon-interaction path

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
  
  static int nGood = 0;
  static int nBad = 0;
  if(fBestResidual > 1){
    nBad++;
    icemcLog() << icemc::warning << "Outside path fitter residual tolerance!,  fBestResidual = " << fBestResidual << ", nGood = " << nGood << ", nBad = " << nBad << std::endl;
  }
  else{
    nGood++;
    icemcLog() << icemc::info << "Successful fit! fBestResidual = " << fBestResidual << ", nGood = " << nGood << ", nBad = " << nBad << std::endl;
  }

  ///@todo add collision check to best fit path  
  
  if(fDebug){
    makeDebugPlots("testFitter.root");
    exit(1);
  }

  return fSurfacePos;
}




void icemc::RayTracer::makeDebugPlots(const TString& fileName) const {
    
  TFile* fTest= new TFile(fileName, "recreate");
 
  const Vector toInteraction = toLocal(fInteractionPos);
  const Vector toSurface = toLocal(fSurfacePos);
  const Vector surfPlusRef = toLocal(fRefractedRfDir, false) + toSurface;
  const Vector surfPlusNorm = 1000*toLocal(fSurfaceNormal, false) + toSurface;
  std::cout << "The dot product is " << fSurfaceNormal.Unit().Dot((fSurfacePos - fBalloonPos).Unit()) << std::endl;
  // std::cout << toSurface.Unit() << "\t" << toLocal(fRefractedRfDir, false) << std::endl;    
  const Vector toEndPoint = toLocal(fEndPoint);

  TGraph* grInteractionTop = new TGraph();
  grInteractionTop->SetPoint(grInteractionTop->GetN(), toInteraction.GetX(), toInteraction.GetY());
  grInteractionTop->SetName("grInteractionTop");
  grInteractionTop->Write();
  delete grInteractionTop;

  TGraph* grInteractionSide = new TGraph();
  grInteractionSide->SetPoint(grInteractionSide->GetN(), toInteraction.GetX(), toInteraction.GetZ());
  grInteractionSide->SetName("grInteractionSide");
  grInteractionSide->Write();
  delete grInteractionSide;
    
  TGraph* grTopView = new TGraph();
  grTopView->SetPoint(0, 0, 0); // balloon at origin
  grTopView->SetPoint(1, toSurface.GetX(), toSurface.GetY());
  grTopView->SetPoint(2, surfPlusRef.GetX(), surfPlusRef.GetY());
  grTopView->SetPoint(3, toEndPoint.GetX(), toEndPoint.GetY());
  grTopView->SetName("grTopView");
  grTopView->Write();
  delete grTopView;

  if(fMinimizerPath){
    fMinimizerPath->SetName("grMinuitPath");
    fMinimizerPath->Write();
  }

  TGraph* grNormalTop = new TGraph();
  grNormalTop->SetPoint(0, toSurface.GetX(), toSurface.GetY());
  grNormalTop->SetPoint(1, surfPlusNorm.GetX(), surfPlusNorm.GetY());
  grNormalTop->SetName("grNormalTop");
  grNormalTop->Write();
  delete grNormalTop;

  TGraph* grNormalSide = new TGraph();
  grNormalSide->SetPoint(0, toSurface.GetX(), toSurface.GetZ());
  grNormalSide->SetPoint(1, surfPlusNorm.GetX(), surfPlusNorm.GetZ());
  grNormalSide->SetName("grNormalSide");
  grNormalSide->Write();
  delete grNormalSide;

  TGraph* grSideView = new TGraph();
  grSideView->SetPoint(0, 0, 0); // balloon at origin
  grSideView->SetPoint(1, toSurface.GetX(), toSurface.GetZ());
  grSideView->SetPoint(2, surfPlusRef.GetX(), surfPlusRef.GetZ());    
  grSideView->SetPoint(3, toEndPoint.GetX(), toEndPoint.GetZ());
  grSideView->SetName("grSideView");
  grSideView->Write();
  delete grSideView;

  TGraph* grSurface = new TGraph();
  grSurface->SetName("grSurface");
  const int nPoints = 1000000;
  const int nExtra = nPoints/10;
  const double localXDistToInteraction = toInteraction.GetX();
  const double deltaX = localXDistToInteraction/nPoints;

  for(int d=-nExtra; d < nPoints+nExtra; d++){
    const std::array<double, 2> params {deltaX*d, 0};
    const Vector v = getSurfacePosition(params.data());
    const Vector vl = toLocal(v);
    grSurface->SetPoint(grSurface->GetN(), vl.GetX(), vl.GetZ());
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
  const vector<double> xMins {x0 - smallPlotDist, 0};
  const vector<double> xMaxs {x0 + smallPlotDist, bigPlotDist};    
  const vector<double> yMins {y0 - smallPlotDist, -0.5*bigPlotDist};
  const vector<double> yMaxs {y0 + smallPlotDist, +0.5*bigPlotDist};
    
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




// int icemc::RayTracer::WhereDoesItLeave(const Position &rfStart, const Vector &rfDirectionUnit, const Antarctica *antarctica, Position &surfaceIntersection, bool debug){

//   Position ray = rfStart;

//   double deltaRadius = antarctica->Surface(surfaceIntersection) - surfaceIntersection.Mag();
//   int iter = 0;
//   const int maxIter = 100e3;

//   // if(deltaRadius > 0){    
//   //   icemcLog() << icemc::warning << "The initial ray trace is above the surface!" << std::endl;
//   // }

//   while(deltaRadius < 0){
//     double stepSize = 1;
//     ray += stepSize * rfDirectionUnit;
//     deltaRadius = antarctica->Surface(surfaceIntersection) - surfaceIntersection.Mag();
//     iter++;
//     if(iter >= maxIter){
//       break;
//     }
//   }

//   surfaceIntersection = ray;

//   if(iter>=maxIter){
//     if(debug){
//       std::cerr << "No solution!" << "\n";
//     }
//     return 0;
//   }

//   return 1;
// }



  // if (whichray==0){
  //   WhereDoesItLeave(posnu, nrf_iceside[2*whichtry],antarctica, // inputs
  // 		     rfexit[whichtry]); // output
  // }  
  
  // //******wufan******
  // if (whichray==1){ // reflected rays
  //   WhereDoesItLeave(posnu_down,nrf_iceside[2*whichtry],antarctica,rfexit[whichtry]);  //from mirror point position and the direction of signals to the ice 
  // //to find the exit point at the surface of the Earth.wufan 
  // }

int icemc::RayTracer::WhereDoesItLeave(const Position &rfStart, const Vector &rfDirectionUnit, const Antarctica *antarctica, Position &surfaceIntersection){  
  
  double distance = 0; ///< How far along rfDirectionUnit (meters)
  double rfStartRadius = rfStart.Mag(); ///< Distance from center of earth to interaction (meters)

  // use law of cosines to get distance from interaction to exit point for the ray
  // need to solve for that distance using the quadratic formula

  // angle between posnu and rfDirectionUnit vector for law of cosines.
  double cosTheta = -1*(rfStart.Dot(rfDirectionUnit))/rfStartRadius;

  // a,b,c for quadratic formula, needed to solve for 
  double a = 1;
  double b = -2*rfStartRadius*cosTheta;
  double c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2);

  if (b*b-4*a*c<0.) {
    return 0;
  }

  // use the "+" solution because the other one is where the ray is headed downward toward the rock
  distance = 0.5*(-b + sqrt(b*b - 4*a*c));

  // now here is the exit point for the ray
  surfaceIntersection = rfStart + distance*rfDirectionUnit;

  c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
  // sometimes though the new surface is lower than the one over posnu which causes a problem.
  if (b*b-4*a*c<0.) {
    //std::cout << "inu is " << inu << "\n";  
    // try halving the distance
    distance=distance/2.;
    //std::cout << "bad.  distance 1/2 is " << distance << "\n";
    surfaceIntersection = rfStart + distance*rfDirectionUnit;
    c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
    if (b*b-4*a*c<0.) { // if we still have the problem back up more
      distance=distance/2.; // now we are at 1/4 the distance
      //std::cout << "bad.  distance 1/4 is " << distance << "\n";
      surfaceIntersection = rfStart + distance*rfDirectionUnit;
      c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
      if (b*b-4*a*c<0.) { // the problem is less then 1/4 of the way in
	distance=distance/2.; // now we are at 1/8 the distance
	//std::cout << "bad.  distance 1/8 is " << distance << "\n";
	surfaceIntersection = rfStart + distance*rfDirectionUnit;
	c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
	if (b*b-4*a*c<0.) {
	  // still have the problem so just make the distance 0
	  distance=0.; // now we are at 1/8 the distance
	  //std::cout << "bad.  distance is " << distance << "\n";
	  surfaceIntersection=antarctica->Surface(surfaceIntersection)/rfStart.Mag()*rfStart;
	}
      } // now we are at 1/8 the distance
      else {// if this surface is ok problem is between 1/4 and 1/2
	distance=distance*1.5; // now we are at 3/8 the distance
	//	std::cout << "good.  distance 3/8 is " << distance << "\n";
	surfaceIntersection = rfStart + distance*rfDirectionUnit;
	c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
	if (b*b-4.*a*c<0.) {
	  distance=distance*2./3.; // go back to 1/4
	  surfaceIntersection = rfStart + distance*rfDirectionUnit;
	  c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
	  //std::cout << "good at distance 1/4 is " << distance << "\n";
	}
      } // now we are at 3/8 the distance
		
		
    } // now we are at 1/4 the distance
    else { // if this surface at 1/2 distance is ok see if we can go a little further
      distance=distance*1.5; // now we are at 3/4 the distance
      //std::cout << "good.  distance 3/4 is " << distance << "\n";
      surfaceIntersection = rfStart + distance*rfDirectionUnit;
      c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
      if (b*b-4*a*c<0.) { // the problem is between 1/2 and 3/4 of the way in
	distance=distance*5./6.; // now we are at 5/8 the distance
	//std::cout << "bad.  distance 5/8 is " << distance << "\n";
	surfaceIntersection = rfStart + distance*rfDirectionUnit;
	c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
	if (b*b-4*a*c<0.) {
	  distance=distance*4./5.;
	  surfaceIntersection = rfStart + distance*rfDirectionUnit;
	  c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
	  //std::cout << "good at distance 3/4 is " << distance << "\n";
	}
      } // now we are at 1/8 the distance
      else {// if this surface is ok problem is between 1/4 and 1/2
	distance=distance*7./6.; // now we are at 7/8 the distance
	//std::cout << "good.  distance 7/8 is " << distance << "\n";
	surfaceIntersection = rfStart + distance*rfDirectionUnit;
	c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
	if (b*b-4*a*c<0) {
	  // now found the problem so go back to 3/4 distance
	  distance=distance*6./7.;
	  //std::cout << "good at  distance 3/4 is " << distance << "\n";
	  surfaceIntersection = rfStart + distance*rfDirectionUnit;
	  c = rfStartRadius*rfStartRadius - pow(antarctica->Surface(surfaceIntersection),2); // redo the law of cosines
	}
      } // now we are at 3/8 the distance
    } // now we are at 3/4 distance
  } // if exit point we initially found was not ok
  else {
    distance=(-1*b+sqrt(b*b-4*a*c))/2; // and quadratic formula
    surfaceIntersection = rfStart + distance*rfDirectionUnit;
  }
  return 1;
}






void icemc::RayTracer::makePlot(const char* name, const Antarctica* antarctica, const Position& nuInteraction,  const Vector& nuDirection, const Position& detector) const {

  // here we make some graphics to show what's going on...
  // everything of interest is on a plane between the detector and neutrino
  // want y axis of TGraphs to be radial outwards vector from nuInteraction
  const Vector graphYHat = nuInteraction.Unit();
  
  // want x axis of TGraphs to be perpendicular to that y-axis in plane of interaction-to-detector...
  // which is a slightly more involved calculation
  const Vector nuIntToDetUnit = (detector - nuInteraction).Unit(); // this vector lies in the plane of interest
  const Vector planeVector = graphYHat.Cross(nuIntToDetUnit).Unit(); // so this vector is perpendicular to the plane of interest...
  const Vector graphXHat = planeVector.Cross(graphYHat).Unit(); // so.. this vector is perpendicular to yhat in the plane of interest

  const double nuIntToDetDist = (detector - nuInteraction).Mag();
  
  const int nPointsSurf = 1000;
  TGraph grSurface;  
  const double altitudeOfInteraction = nuInteraction.Mag();
  const double stepSize = nuIntToDetDist/nPointsSurf;

  const double extraRange = 100;
  double maxY =  extraRange; // nuInteraction is at TGraph origin, units are meters, so start here
  double minY = -extraRange; // nuInteraction is at TGraph origin, units are meters, so start here
  
  for(int i=-100; i < nPointsSurf + 100; i++){
    Position getSurfaceHere = nuInteraction + nuIntToDetUnit*i*stepSize;
    double surface = antarctica->Surface(getSurfaceHere);
    double surfaceRelativeToNewOrigin = surface - altitudeOfInteraction;

    double x = (getSurfaceHere - nuInteraction).Dot(graphXHat);

    // std::cout << i << "\t" << x << "\t" << surfaceRelativeToNewOrigin << "\t" << surface << std::endl;
    
    grSurface.SetPoint(grSurface.GetN(), x, surfaceRelativeToNewOrigin);

    if(surfaceRelativeToNewOrigin > maxY){
      maxY = surfaceRelativeToNewOrigin + extraRange;
    }
    if(surfaceRelativeToNewOrigin < minY){
      minY = surfaceRelativeToNewOrigin - extraRange;
    }
    
  }

  
  TCanvas c1;
  grSurface.SetName("grSurface");
  grSurface.Draw("al");

  TMarker interactionPoint(0, 0, 8);
  interactionPoint.SetMarkerColor(kMagenta);
  interactionPoint.Draw("psame");  

  TMarker balloon((detector - nuInteraction).Dot(graphXHat),
		  (detector - nuInteraction).Dot(graphYHat),
		  8);
  balloon.SetMarkerColor(kCyan); 
  
  balloon.Draw("psame");

  const double arrowLength = 1000; // meters in the output plot
  TArrow a(-(arrowLength*nuDirection).Dot(graphXHat), -(arrowLength*nuDirection).Dot(graphYHat), 0, 0, 0.05, "|>");
  a.SetLineColor(interactionPoint.GetMarkerColor());
  a.SetFillColor(interactionPoint.GetMarkerColor());  
  a.Draw("same");
  
  TLegend l1(0.8, 0.8, 1, 1);
  l1.AddEntry(&grSurface, "Ice Surface", "l");
  
  std::vector<TGraph> grs;
  grs.reserve(numTries);
  for(int i=0; i < numTries; i++){

    const double bigEnoughToBeARealExitGuess = 1.1;
    if(rfexit[i].Mag() > bigEnoughToBeARealExitGuess){


      std::vector<double> x_iceside;
      std::vector<double> y_iceside;

      Vector downwards = rfexit[i];
      const double scaleFactorDownwards = 10;
      while((downwards - nuInteraction).Dot(graphYHat) > 0){

	// -ve sign very important here, otherwise 
	downwards += scaleFactorDownwards*-nrf_iceside[i];

	double x = (downwards - nuInteraction).Dot(graphXHat);
	double y = (downwards - nuInteraction).Dot(graphYHat);

	if(y > maxY){
	  maxY = y + extraRange;
	}
	if(y < minY){
	  minY = y - extraRange;
	}		

	x_iceside.push_back(x);
	y_iceside.push_back(y);	
      }

      grs.emplace_back(TGraph());
      TGraph& gr = grs.back();

      for(int j=x_iceside.size() -1; j >= 0; j--){
	gr.SetPoint(gr.GetN(), x_iceside.at(j), y_iceside.at(j));
      }
      gr.SetPoint(gr.GetN(), (rfexit[i] - nuInteraction).Dot(graphXHat), (rfexit[i] - nuInteraction).Dot(graphYHat));

      Vector upwards = rfexit[i];
      const double scaleFactorUpwards = 1000;
      while((upwards - detector).Dot(graphXHat) < 0){
	upwards += scaleFactorUpwards*n_exit2bn[i];

	double x = (upwards - nuInteraction).Dot(graphXHat);
	double y = (upwards - nuInteraction).Dot(graphYHat);

	if(y > maxY){
	  maxY = y + extraRange;
	}
	if(y < minY){
	  minY = y - extraRange;
	}	

	gr.SetPoint(gr.GetN(), x, y);
      }

      gr.SetName(TString::Format("gr_%d", i));
      gr.SetLineColor(i+2);
      gr.SetLineStyle(i+2);      
      gr.Draw("lsame");

      const Vector anitaToRfExit = rfexit[i] - detector;
      //  what's the angle between ANITA's position (vector from earth centre to ANITA), and a vector from ANITA to the RF interaction      
      const double trueTheta = detector.Angle(anitaToRfExit);
      // if the interaction was directly below ANITA, trueTheta would be 180 degrees
      // if the interaction was directly horizontal from ANITA, trueTheta would be 90 degrees
      // if the interaction was directly above ANITA, trueTheta would be 0 degrees
      // therefore to convert to payload coordinates where  theta=0 is horizontal, and +ve theta is UP, -ve theta is DOWN...
      const double trueThetaPayloadCoordinates = 0.5*TMath::Pi() - trueTheta;
      l1.AddEntry(&gr, TString::Format("rfExit[%d], #theta_{payload} = %4.2f^{#circ}", i, trueThetaPayloadCoordinates*TMath::RadToDeg()), "l");
    }
  }

  TString title = "Ray Tracing (Earth local vertical/horizontal at interaction point); Horizontal distance (m); Vertical distance (m)";
  grSurface.SetTitle(title);
  grSurface.GetXaxis()->SetNoExponent(1);
  grSurface.GetYaxis()->SetNoExponent(1);  
  grSurface.SetMaximum(maxY);
  grSurface.SetMinimum(minY);
  l1.Draw();
  c1.Modified();
  c1.Update();
  c1.SaveAs(name);

  return;
  
}
