// This program prints / draws confirmed supernovae

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
  std::string *objType = new std::string;
  int discoveryUnixTime = 0;

  // New vars
  int a_tmin;
  int a_tmax;
  int passed = 0;
  std::string desiredType;
  size_t found;

  
  // A4 for now
  a_tmin = 1480707643;
  a_tmax = 1483004563;
  
  // Branches
  tree->SetBranchAddress("ra",&ra);
  tree->SetBranchAddress("dec",&dec);
  tree->SetBranchAddress("objType",&objType);
  tree->SetBranchAddress("name",&name);
  tree->SetBranchAddress("discoveryUnixTime",&discoveryUnixTime);
  
  TMarker *astroObject = new TMarker();
  SkyMap *skyMapOut = new SkyMap(90);
  
  for(unsigned int i = 0; i < discoveryEntries; i++)
    {
      tree->GetEntry(i);

      // Skip those if the discovery time falls outside of ANITA flight time
      if( (discoveryUnixTime < a_tmin) || (discoveryUnixTime > a_tmax) ){continue;}
      if( abs(dec)>30 ){continue;}
      //desiredType = "SN II";
      //found = objType->find(desiredType);
      //if(found != 0){continue;}
      if(*name != "SN 2016jby" && *name != "SN 2016iyz"){continue;}
      
      //cout << ra << endl;
      //cout << dec << endl;
      cout << *name << ", ";
      cout << *objType << "- ";
      // Set ra, dec
      astroObject->SetX(ra);
      astroObject->SetY(dec);
      astroObject->SetMarkerStyle(20);
      astroObject->SetMarkerSize(1.5);
      astroObject->SetMarkerColor(9);

      //if(*name == "SN 2016isg"){astroObject->SetMarkerColor(2);} 
      
      skyMapOut->addMarker(astroObject);
      passed++;
	
    }
  cout << passed << " supernovae with discovery time coinciding with the A4 flight." << endl;
  skyMapOut->Draw("");
  
  return;
  
}
