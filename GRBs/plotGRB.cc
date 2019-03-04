// This program prints / draws GRBs

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

void plotGRB()
{
  TFile *file = new TFile("./data/GRBIceCube.root"); 
  TTree *tree = (TTree*) file->Get("iceCubeTree");
  UInt_t entries = tree->GetEntries();

  // Root vars
  float RA = 0;
  float dec = 0;
  std::string *name = new std::string;
  int unixTriggerTime = 0;

  // New vars
  int a_tmin;
  int a_tmax;
  int passed = 0;
  
  // A4 for now
  a_tmin = 1480707643;
  a_tmax = 1483004563;
  
  // Branches
  tree->SetBranchAddress("RA",&RA);
  tree->SetBranchAddress("dec",&dec);
  tree->SetBranchAddress("name",&name);
  tree->SetBranchAddress("unixTriggerTime",&unixTriggerTime);
  
  TMarker *astroObject = new TMarker();
  SkyMap *skyMapOut = new SkyMap(90);
  
  for(unsigned int i = 0; i < entries; i++)
    {
      tree->GetEntry(i);

      // Skip those if the discovery time falls outside of ANITA flight time
      if( (unixTriggerTime < a_tmin) || (unixTriggerTime > a_tmax) ){continue;}
      
      //cout << ra << endl;
      //cout << dec << endl;
      cout << *name << ", ";
      // Set ra, dec
      astroObject->SetX(RA);
      astroObject->SetY(dec);
      astroObject->SetMarkerStyle(20);
      astroObject->SetMarkerSize(1.5);
      astroObject->SetMarkerColor(9);

      //if(*name == "SN 2016isg"){astroObject->SetMarkerColor(2);} 
      
      skyMapOut->addMarker(astroObject);
      passed++;
	
    }
  cout << "" << endl;
  cout << passed << " GRBs with detection time coinciding with the A4 flight." << endl;
  skyMapOut->Draw("");
  
  return;
  
}
