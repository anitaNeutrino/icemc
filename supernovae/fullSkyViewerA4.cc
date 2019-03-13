// This program gives a view of an all-sky map
// of different lons
// for the supernova observed during ANITA-4

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

void fullSkyViewerA4()
{
  TFile *file = new TFile("./data/supernovaeTNS.root"); 
  TTree *tree = (TTree*) file->Get("tnsTree");
  UInt_t discoveryEntries = tree->GetEntries();

  // Root vars
  float ra = 0;
  float dec = 0;
  int passed = 0;

  int discoveryUnixTime = 0;

  // New vars
  int a_tmin;
  int a_tmax;
  
  // A4 for now
  a_tmin = 1480707643;
  a_tmax = 1483004563;
  
  // Branches
  tree->SetBranchAddress("ra",&ra);
  tree->SetBranchAddress("dec",&dec);
  tree->SetBranchAddress("discoveryUnixTime",&discoveryUnixTime);

  double x[42] = {};
  double y[42] = {};

  TCanvas *c = new TCanvas("c","c",2000,1000);
  TMarker *marker = new TMarker();
  
  for(signed int lon = -180; lon < 1080; lon++)
    {
      SkyMap *skyMapOut = new SkyMap(lon);
      
      for(unsigned int i = 0; i < 42; i++)
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
