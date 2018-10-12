#include "RootOutput.h"
#include <iostream>
#include "EventGenerator.h"
#include "ConnollyEtAl2011.h"
// #include "balloon.hh"
// #include "anita.hh"
#include "ANITA.h"
#include "Settings.h"
#include "Taumodel.hh"
#include "EnvironmentVariable.h"
#include "Tools.h"
#include "RayTracer.h"
#include "screen.hh"
#include "Geoid.h"
#include "Antarctica.h"
#include "Report.h"




icemc::RootOutput::RootOutput(const EventGenerator* uhen, const Settings* settings, const char* outputDir, int run)
  : fOutputDir(outputDir), fRun(run), fIceFinal(NULL)
{

  initIceFinal(uhen, settings);
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
  Settings* settings = const_cast<Settings*>(settings2);

  // first the file(s)
  TString fileName = fOutputDir + TString::Format("/icefinal%d.root", fRun);
  fIceFinal = new TFile(fileName, "RECREATE", "ice");

  TNamed* ss = settings2->makeRootSaveableSettings();
  ss->Write();
  delete ss;
  ss = NULL;

  // initTree(&allTree, "allTree", "allTree", fIceFinal);
  // allTree.Branch("genNu", &uhen->fGenNu);

  // initTree(&passTree, "passTree", "passTree", fIceFinal);
  // passTree.Branch("passNu", &uhen->fPassNu);

  // // histograms
  // initHist(&ref_int_coord, "ref_int_coord", "", 600, -3000, 3000, 500, -2500, 2500, fIceFinal);
  // ref_int_coord.SetMarkerSize(0.7);
  // ref_int_coord.SetMarkerColor(kBlue);

  
  // initHist(&dir_int_coord, "dir_int_coord", "", 600, -3000, 3000, 500, -2500, 2500, fIceFinal);
  // dir_int_coord.SetMarkerSize(0.7);
  // dir_int_coord.SetMarkerStyle(30);


  // initHist(&h1mybeta, "betaforall", "betaforall(deg)", 180, -15, 15, fIceFinal);
  // initHist(&h1mytheta, "mytheta", "mytheta(deg)", 180, -90, 90, fIceFinal);
  // initHist(&hundogaintoheight_e, "undogaintoheight_e", "undogaintoheight_e", 100, 0., 1., fIceFinal);
  // initHist(&hundogaintoheight_h, "undogaintoheight_h", "undogaintoheight_h", 100, 0., 1., fIceFinal);
  // initHist(&rec_diff, "rec_diff", "rec_diff", 100, -1., 1., fIceFinal);
  // initHist(&recsum_diff, "recsum_diff", "recsum_diff", 100, -1., 1., fIceFinal);
  // initHist(&rec_diff0, "rec_diff0", "rec_diff0", 100, -1., 1., fIceFinal);
  // initHist(&rec_diff1, "rec_diff1", "rec_diff1", 100, -1., 1., fIceFinal);
  // initHist(&rec_diff2, "rec_diff2", "rec_diff2", 100, -1., 1., fIceFinal);
  // initHist(&rec_diff3, "rec_diff3", "rec_diff3", 100, -1., 1., fIceFinal);
  // initHist(&prob_eachilon_bn, "prob_eachilon_bn", "prob_eachilon_bn", 180, 0., 180., fIceFinal);
  // initHist(&h6, "theta_vs_hitangle_h", "theta_vs_hitangle_h", 100, -3.14, 3.14, 100, -1.1, 1.1, fIceFinal);
  // initHist(&h10, "hitangle_e", "hitangle_e", 20, -1.6, 1.6, fIceFinal);
  // initHist(&hy, "hy", "hy", 100, 0., 1., fIceFinal);
  // initHist(&fraction_sec_muons, "fraction_sec_muons", "fraction_sec_muons", 100, 0., 1., fIceFinal);
  // initHist(&fraction_sec_taus, "fraction_sec_taus", "fraction_sec_taus", 100, 0., 1., fIceFinal);
  // initHist(&n_sec_muons, "n_sec_muons", "n_sec_muons", 100, 0., 10., fIceFinal);
  // initHist(&n_sec_taus, "n_sec_taus", "n_sec_taus", 100, 0., 10., fIceFinal);
  // initHist(&sampleweights, "sampleweights", "sampleweights", 100, -5., 0., fIceFinal);
  
  
  



}








