#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"

void plotLightCurve()
{
  TFile *file = new TFile("3FGL.root"); 
  TTree *tree = (TTree*) file->Get("catTree");
  const Int_t maxEntries = tree->GetEntries();
  cout << "Number of total objects: " << maxEntries << endl;

  std::string * sourceName = new std::string();
  std::string * association = new std::string();
  std::vector<int> * met = 0;
  std::vector<double> * le_count = 0;

  tree->SetBranchAddress("sourceName",&sourceName);
  tree->SetBranchAddress("association",&association);
  tree->SetBranchAddress("met",&met);
  tree->SetBranchAddress("le_count",&le_count);

  Int_t chosenObject = 4; // Object to plot light curve of if you only want to plot one
  
  for (int i = chosenObject; i < chosenObject+1; i++) 
    {
      tree->GetEntry(i);
    } // just to get metSize
  
  const int metSize = met->size();
  int t[metSize], le[metSize];
  
  for (int i = chosenObject; i < chosenObject+1; i++) 
    {
      tree->GetEntry(i);
      std::cout << "Producing light curve of " << sourceName->c_str() << " --- " << association->c_str() << std::endl;
      // amount of data points is variable.
      for(int j = 0; j < met->size(); j++)
	{
	  //cout << "met = " << (*met)[j] << ", le_count = " << (*le_count)[j] << endl;
	  t[j] = (*met)[j];
	  le[j] = (*le_count)[j];
	}
    }

	TGraph *g = new TGraph(metSize,t,le);
      g->SetTitle(TString::Format("Light curve of %s - %s",sourceName->c_str(),association->c_str()));
      g->Draw();
  
      return;
  
}
