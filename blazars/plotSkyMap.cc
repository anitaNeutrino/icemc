/*
  Example macro to draw a skymap and associated objects.
  Use FAVA data as input, or any ra/dec data.
*/
#include <iostream>
#include <libgen.h>
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
#include "fava.h" 
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "SkyMap.h"


void plotSkyMap()
{
  // Load root file with ra/dec data
  TFile *file = new TFile("fava.root");
  TTree *tree = (TTree*) file->Get("fava");
  TCanvas *c1 = new TCanvas();

  // Use fava.h
  FAVAEntry *fava = new FAVAEntry();
  tree->SetBranchAddress("fava",&fava);
  const int maxEntries = tree->GetEntries();

  // Use TMarkers to place our objects
  TMarker *astroObject = new TMarker();
  SkyMap *test = new SkyMap();

  // Style
  astroObject->SetMarkerStyle(20);
  astroObject->SetMarkerSize(1);
  
  for (int i = 0; i < tree->GetEntries(); i++) 
    {
      tree->GetEntry(i);
      astroObject->SetX(fava->ra);
      astroObject->SetY(fava->dec);
      test->addMarker(astroObject);
    }

  cout << "Astrophysical objects plotted: " << maxEntries << endl;

  test->Draw();    
  
  return;
  
}
