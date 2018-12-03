#include "RootOutput.h"
#include <iostream>
#include "EventGenerator.h"
#include "ConnollyEtAl2011.h"
#include "Settings.h"
#include "EnvironmentVariable.h"
#include "Tools.h"
#include "RayTracer.h"
#include "SurfaceScreen.h"
#include "Geoid.h"
#include "Antarctica.h"
#include "Report.h"




icemc::RootOutput::RootOutput(const EventGenerator* uhen, const Settings* settings, const char* outputDir, int run)
  : fOutputDir(outputDir), fRun(run), fIceFinal(NULL)
{
  std::cout << "icemc::RootOutput created in write mode." << std::endl;
  initIceFinal(uhen, settings);
}


icemc::RootOutput::RootOutput(const char* outputDir, int run)
  : fOutputDir(outputDir), fRun(run), fIceFinal(NULL)
{
  std::cout << "icemc::RootOutput created in read mode." << std::endl;
  readIceFinal();
}



void icemc::RootOutput::readIceFinal(){

  TString fileName = fOutputDir + TString::Format("/run%d/IceFinal_%d.root", fRun, fRun);
  fIceFinal = TFile::Open(fileName);

  if(!fIceFinal){
    icemc::report() << severity::error << "Can't open " << fileName << std::endl;
    return;
  }
  
  fAllTree = (TTree*) fIceFinal->Get("allTree");
  fPassTree = (TTree*) fIceFinal->Get("passTree");

  fAllTree->SetBranchAddress("eventSummary", &fEventSummary);
  fPassTree->SetBranchAddress("event", &fEvent);
}

icemc::RootOutput::~RootOutput(){

  // write, close and delete all non-NULL member files.
  const int numFiles = 1;
  TFile* fs[numFiles] = {fIceFinal};
  for(int i=0; i < numFiles;  i++){
    if(fs[i]){
      fs[i]->Write();
      fs[i]->Close();
      delete fs[i];
    }
  }
}



void icemc::RootOutput::initTree(TTree* t, const char* name, const char* title, TFile* f){
  t->SetName(name);
  t->SetTitle(title);
  t->SetDirectory(f);
}


void icemc::RootOutput::initHist(TH2* h, const char* name, const char* title,
				 int nx, double xMin, double xMax,
				 int ny, double yMin, double yMax, TFile* f){
  h->SetNameTitle(name, title);
  h->SetBins(nx, xMin, xMax, ny,  yMin, yMax);
  h->SetDirectory(f);
}


void icemc::RootOutput::initHist(TH1* h, const char* name, const char* title,
				 int nx, double xMin, double xMax, TFile* f){
  h->SetNameTitle(name, title);
  h->SetBins(nx, xMin, xMax);
  h->SetDirectory(f);
}



void icemc::RootOutput::initIceFinal(const EventGenerator* uhen2, const Settings* settings2){

  if(fIceFinal){
    icemc::report() << severity::warning << "IceFinal already initialized!"  << std::endl;
    return;
  }

  EventGenerator* uhen = const_cast<EventGenerator*>(uhen2);
  // Settings* settings = const_cast<Settings*>(settings2);

  // first the file(s)
  TString fileName = fOutputDir + TString::Format("/run%d/IceFinal_%d.root", fRun, fRun);
  fIceFinal = new TFile(fileName, "RECREATE", "ice");

  TNamed* ss = settings2->makeRootSaveableSettings();
  ss->Write();
  delete ss;
  ss = NULL;

  fAllTree = new TTree();
  initTree(fAllTree, "allTree", "allTree", fIceFinal);
  fAllTree->Branch("eventSummary", &uhen->fEventSummary);
  
  fPassTree = new TTree();
  initTree(fPassTree, "passTree", "passTree", fIceFinal);
  fPassTree->Branch("event", &uhen->fEvent);

  // these are only for read only mode,
  // disable these pointers when we're writing
  // as we are copying the EventGenerator members into the trees.
  fEventSummary = nullptr;
  fEvent = nullptr;
}








