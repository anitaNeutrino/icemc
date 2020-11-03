#include "Crust.h"
#include "TH2DAntarctica.h"
#include "TFile.h"

int main(int argc, char **argv){

  int model=2;
  if (argc==2)
    model = atoi(argv[1]);
  
  std::cout << "Using model " << model << std::endl;
  const icemc::Crust c(1,model); 
 
  std::cout << c.GetTotalIceVolume() << " m^3 of Antarctic ice" << std::endl;
  std::cout << c.GetTotalIceArea() << " m^2 of Antarctic is covered by ice" << std::endl;
  
  int nx, ny;
  nx = ny = 500;
  TFile* fOut = new TFile("testCrust.root", "recreate");
  TH2DAntarctica*  h = new TH2DAntarctica(nx, ny);
  h->SetName("hSurfaceAboveGeoid");

  TH2DAntarctica*  h1 = new TH2DAntarctica(nx, ny);  
  h1->SetName("hWaterDepth");

  TH2DAntarctica*  h2 = new TH2DAntarctica(nx, ny);  
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

  double max[3] = {0};

  double totalvolume = 0;
  
  for(int by=1; by <= h->GetNbinsY(); by++){
    const double northing = h->GetYaxis()->GetBinCenter(by);
    const double ywidth = h->GetYaxis()->GetBinWidth(by);
    for(int bx=1; bx <= h->GetNbinsX(); bx++){
      const double easting = h->GetXaxis()->GetBinCenter(bx);      
      const double xwidth = h->GetXaxis()->GetBinWidth(bx);
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

      if(th >= 0){
	totalvolume += xwidth*ywidth*th;
	//std::cout << lon << ", " << lat << ":   " << xwidth << " X " << ywidth << ", thickness=" << th << std::endl;
	if(th > max[0]){
	  max[0] = th;
	  max[1] = lon;
	  max[2] = lat;
	}
      }
    }
  }

    std::cout << "Max thickness = " << max[0] << " at " << max[1] << ", " << max[2] << std::endl;
    std::cout << "Total volume = " << totalvolume << std::endl;
  
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
