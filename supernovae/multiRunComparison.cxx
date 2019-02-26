#include <iostream>
#include "TTree.h"
#include "TGraph.h"
#include "TFile.h"
#include "TAxis.h"

void multiRunComparison()
{
  std::cout << "..." << std::endl;
  const int tests = 8;
  int decCut[tests] = {2,5,10,20,30,50,70,90};
  int passed[tests] = {};
  
  for(int i = 0; i < tests; i++)
    {
      TFile *file = new TFile(TString::Format("./sourceOutput/dec%d/SimulatedAnitaTruthFile.root",decCut[i])); 
      TTree *tree = (TTree*) file->Get("truthAnitaTree");
      //tree->Scan();
      Int_t entries = tree->GetEntries();
      passed[i] = entries;
      cout << entries << " neutrinos passed for a dec cut of |" << decCut[i] << "|" << endl;
      file->Close();
    }

  TGraph *gr = new TGraph(tests,decCut,passed);
  gr->SetTitle("Neutrinos passed for objects in a declination-limited sky");
  gr->SetLineColor(4);
  gr->SetMarkerStyle(24);
  gr->SetMarkerSize(1.5);
  gr->SetMarkerColor(1);
    
  gr->GetXaxis()->SetTitle("|dec cut|");
  gr->GetYaxis()->SetTitle("# neutrinos passed");
  gr->Draw("AL*");
  
  return;
}
