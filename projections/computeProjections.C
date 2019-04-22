#include "TH2.h" 
#include "TVector3.h" 
#include "TMath.h"
#include "TRandom.h" 
#include "TStyle.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TView3D.h"
#include "TFile.h" 

void computeProjections(int drawobs = -1, double dcap = 1, double dobs = 1, double maxcap =30, int N = 100000000) 
{

//  gStyle->SetCanvasPreferGL(1); 

  TFile out("projections.root","RECREATE"); 
  int ncap = maxcap/dcap; 
  double min_obs = 0; 
  double max_obs = 90+maxcap; 
  double nobs = (max_obs-min_obs)/dobs + 1; 

  TH2D * h = new TH2D("projections","Projected Area of Spherical Cap onto Plane;Observation angle (deg); Cap angle (deg)", 
    nobs, min_obs-dobs/2, max_obs + dobs/2, 
    ncap, dcap/2., maxcap+dcap/2
  ); 



 for (int iobs = 1; iobs <= h->GetNbinsX(); iobs++) 
 {

    double obs = h->GetXaxis()->GetBinCenter(iobs); 
    printf("%g\n",obs); 

    //use phi = 0
    double obs_rad = obs *TMath::DegToRad(); 
    TVector3 l( sin(obs_rad), 0, cos(obs_rad)); 

    //randomly draw within a plane

    int npass = 0; 

    if (drawobs == iobs) 
    {

      TCanvas * c = new TCanvas; 
      const double maxs[] = {2,2,2}; 
      TView3D * view = (TView3D*)TView3D::CreateView(1);
      view->SetRange(-2,-2,-2,2,2,2,2); 

      double r2 =  pow(sin(maxcap*TMath::DegToRad()),2);
      TPolyMarker3D * spherical_cap = new TPolyMarker3D(2*41*41); 
      int ii = 0; 
      for (double x = -1; x<=1; x+=0.05)
      {
        for (double y = -1; y<=1; y+=0.05)
        {
          if (x*x+y*y >1) continue;

          spherical_cap->SetPoint(ii++, x,y, -(sqrt(1-x*x-y*y))); ;
          spherical_cap->SetPoint(ii++, x,y, (sqrt(1-x*x-y*y))); ;
        }
      }
      spherical_cap->Draw("psame"); 
    }

    TPolyMarker3D * pts = iobs ==drawobs ?  new TPolyMarker3D(N/10) : 0;
    TPolyMarker3D * pts2 = iobs ==drawobs ? new TPolyMarker3D(N/10) : 0;
    for (int i = 0; i < N; i++) 
    {
      double x= gRandom->Uniform(-1,1); 
      double y= gRandom->Uniform(-1,1); 

      TVector3 o(x,y,0); 

      //rotate point to observation angle 
      o.RotateUz(l);  
      o+=l; 

      // and move away so that the plane is tangent to the sphere
      

      //now find the sphere plane intersection 
      double  discr = pow(l.Dot(o),2) - ( o.Dot(o) - 1);

      if (discr <0) 
      {
        continue; //no solution
      }



      //we only take one angle to avoid double counting

      double second_term = sqrt(discr);
      double first_term = -(l.Dot(o)); 
      double d1 = first_term + second_term; 

      TVector3 x1 = o + d1*l; 
      double ang = x1.Theta()*TMath::RadToDeg(); 
      if (drawobs == iobs && ang < maxcap) 
      {
        if (i < pts2->GetN())
        {
          pts->SetPoint(i,o.X(),o.Y(),o.Z());
          pts2->SetPoint(i,x1.X(),x1.Y(),x1.Z());
        }
      }




      for (int icap = 1; icap <= h->GetNbinsY(); icap++) 
      {
        double cap = icap*dcap; 
        if (ang < cap) h->SetBinContent(iobs,icap,h->GetBinContent(iobs,icap)+1); 
      }

    }

   if (drawobs == iobs) 
   {
     pts->SetMarkerColor(2); 
     pts->Draw("psame"); 
     pts2->SetMarkerColor(3); 
     pts2->Draw("psame"); 
   }
 }
 
 h->Scale(4./N); 
  
 h->Write(); 



}
