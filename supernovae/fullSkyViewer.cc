// This program gives you a view of an all-sky map
// for all different lons
// for selected astrophysical objects with
// known RA and Dec
// ---
// Uses a Mollweidge projection, see
// the Mollweide projection of Earth:
// https://www.giss.nasa.gov/tools/gprojector/help/projections/Mollweide.png

#include <iostream>
#include <libgen.h>
#include <string>
#include <map>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TMath.h" 
#include "TObjString.h" 
#include "TTimeStamp.h"
#include "TH1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "SkyMap.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCanvas.h"
#include "TStyle.h"

void fullSkyViewer()
{
  TFile *file = new TFile("./data/supernovaeTNS.root"); 
  TTree *tree = (TTree*) file->Get("tnsTree");
  UInt_t discoveryEntries = tree->GetEntries();

  // Root vars
  float ra = 0;
  float dec = 0;
  int passed = 0;

  // Branches
  tree->SetBranchAddress("ra",&ra);
  tree->SetBranchAddress("dec",&dec);

  double x[638] = {};
  double y[638] = {};

  TCanvas *c = new TCanvas("c","c",2000,1000);
  TMarker *marker = new TMarker();
  
  for(signed int lon = -180; lon < 1080; lon++)
    {
      SkyMap *skyMapOut = new SkyMap(lon);
      
      for(unsigned int i = 0; i < 638; i++)
	{
	  tree->GetEntry(i);
	  marker->SetX(ra);
	  marker->SetY(dec);

	  marker->SetMarkerColor(9);
	  marker->SetMarkerStyle(20);
	  marker->SetMarkerSize(1.5);
	  skyMapOut->addMarker(marker);
	}
      
      skyMapOut->Draw();
      c->Update();

      delete skyMapOut;
      
    }
  
  return;
}
