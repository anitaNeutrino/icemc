#include "Crust2.h"
#include "TH2DAntarctica.h"
#include "TFile.h"

int main(){

  icemc::Crust2 c;

  TFile* fOut = new TFile("testCrust.root", "recreate");
  TH2DAntarctica*  h = new TH2DAntarctica();
  h->SetName("hSurfaceAboveGeoid");

  TH2DAntarctica*  h1 = new TH2DAntarctica();  
  h1->SetName("hWaterDepth");

  TH2DAntarctica*  h2 = new TH2DAntarctica();  
  h2->SetName("hIceThickness");

  TH2D* hBoundary = new TH2D("hBoundary", "hBoundary", 128, -46, -44, 128, 44, 46);
  for(int by=1; by <= hBoundary->GetNbinsY(); by++){
    const double lat = hBoundary->GetYaxis()->GetBinCenter(by);
    for(int bx=1; bx <= hBoundary->GetNbinsX(); bx++){
      const double lon = hBoundary->GetXaxis()->GetBinCenter(bx);
      Geoid::Position p;
      p.SetLonLatAlt(lon, lat, 0.);
      double s = c.SurfaceAboveGeoid(p);
      hBoundary->SetBinContent(bx, by, s);
    }
  }
  
  for(int by=1; by <= h->GetNbinsY(); by++){
    const double northing = h->GetYaxis()->GetBinCenter(by);
    for(int bx=1; bx <= h->GetNbinsX(); bx++){
      const double easting = h->GetXaxis()->GetBinCenter(bx);      
      Geoid::Position p;
      double lon, lat;
      Geoid::EastingNorthingToLonLat(easting, northing, lon, lat);
      p.SetLonLatAlt(lon, lat, 0);
      double r = c.SurfaceAboveGeoid(p);
      h->SetBinContent(bx, by, r);

      double wd = c.WaterDepth(p);
      h1->SetBinContent(bx, by, wd);
      
      double th = c.IceThickness(p);
      h2->SetBinContent(bx,  by, th);
      // std::cout << easting << "\t" << northing << "\t" << r << "\t" << p.Surface() - p.Mag() << std::endl;
    }
  }
  h->Write();  
  delete h;  
  h = nullptr;
  
  h1->Write();  
  delete h1;  
  h1 = nullptr;
  
  h2->Write();  
  delete h2;  
  h2 = nullptr;

  TH2D* h3 = new TH2D("hSurfaceAboveGeoid2", "Lon/lat", 360, -180, 180,  180, -90, 90);
  for(int by=1; by <= h3->GetNbinsY(); by++){
    const double lat = h3->GetYaxis()->GetBinCenter(by);
    for(int bx=1; bx <= h3->GetNbinsX(); bx++){
      const double lon = h3->GetXaxis()->GetBinCenter(bx);
      Geoid::Position p;
      p.SetLonLatAlt(lon, lat, 0);
      double r = c.SurfaceAboveGeoid(p);
      h3->SetBinContent(bx, by, r);
    }
  }

  fOut->Write();
  fOut->Close();
     
}
