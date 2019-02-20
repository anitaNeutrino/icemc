// This program prints / draws info confirmed supernovae

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

void plotSuperNovaeTNS()
{
  TFile *file = new TFile("./data/supernovaeTNS.root"); 
  TTree *tree = (TTree*) file->Get("tnsTree");
  UInt_t discoveryEntries = tree->GetEntries();

  // Root vars
  float ra = 0;
  float dec = 0;
  std::string *name = new std::string;
  int discoveryUnixTime = 0;

  // New vars
  int a_tmin;
  int a_tmax;
  int passed = 0;
  
  // A4 for now
  a_tmin = 1480707643;
  a_tmax = 1483004563;
  
  // Branches
  tree->SetBranchAddress("ra",&ra);
  tree->SetBranchAddress("dec",&dec);
  tree->SetBranchAddress("name",&name);
  tree->SetBranchAddress("discoveryUnixTime",&discoveryUnixTime);
  
  TMarker *astroObject = new TMarker();
  SkyMap *skyMapOut = new SkyMap(0);
  
  for(unsigned int i = 0; i < discoveryEntries; i++)
    {
      tree->GetEntry(i);

      // Skip those if the discovery time falls outside of ANITA flight time
      if( (discoveryUnixTime < a_tmin) || (discoveryUnixTime > a_tmax) ){continue;}

      //cout << ra << endl;
      //cout << dec << endl;
      cout << *name << endl;
      // Set ra, dec
      astroObject->SetX(ra);
      astroObject->SetY(dec);
      astroObject->SetMarkerStyle(20);
      astroObject->SetMarkerSize(2);
      skyMapOut->addMarker(astroObject);
      passed++;
	
    }
  cout << passed << " supernovae with discovery time coinciding with the A4 flight." << endl;
  skyMapOut->Draw("");
  
  return;
  
}
