////////////////////////////////////////////////////////////////////////
// Simple macro to calculate trigger efficiency as a function of SNR
// To run: root -b -q plotSNR.C
// author: L. Cremonesi l.cremonesi@ucl.ac.uk
////////////////////////////////////////////////////////////////////////

void plotSNR(){

  // Change pathname to the path to your icemc runs
  string path = "/datapool/software/anitaBuildTool/build/components/icemc/outputs";

  // Change first and last run number appropriately
  int firstRun = 321;
  int lastRun  = 348;


  // Define TChains for the head and truth Trees
  TChain *tHead = new TChain("headTree");
  TChain *tTrue = new TChain("truthAnitaTree");

  // Read the trees
  for (int irun=firstRun; irun<=lastRun; irun++){
    tHead->Add(Form("%s/run%d/SimulatedAnitaHeadFile%d.root",  path.c_str(), irun, irun));
    tTrue->Add(Form("%s/run%d/SimulatedAnitaTruthFile%d.root", path.c_str(), irun, irun));

  } 

  RawAnitaHeader *header = NULL;
  TruthAnitaEvent *truth = NULL;

  tHead->SetBranchAddress("header", &header);
  tTrue->SetBranchAddress("truth", &truth);
  
  int nentries = tHead->GetEntries();

  int ntot = 40;
  double snrmin = 0.;
  double snrmax = 20.;

  // Define the numerator and denominator histograms to calculate the efficiency
  TH1D *hNum = new TH1D ("hNum", "", ntot, snrmin, snrmax);
  TH1D *hDenom = new TH1D ("hDenom", "", ntot, snrmin, snrmax);
  hNum->Sumw2();
  hDenom->Sumw2();
  
  double snr=0.;

  // Define output file
  TFile *fout = new TFile("icemcWAISeff2.root", "recreate");


  // Loop over all the entries in the TChains
  for (int ientry=0; ientry<nentries; ientry++){

    tHead->GetEntry(ientry);
    tTrue->GetEntry(ientry);

    // Get SNR value at Digitizer (change V to H in case of using HPOL)
    snr = truth->maxSNRAtDigitizerH;
    if (snr>snrmax) snr=snrmax*0.9999;
    
    hDenom->Fill(snr, 1);

    // if there is a trigger add the event to the numerator
    // NB l3TrigPattern is for VPOL
    //    l3TrigPatternH is for HPOL
    if (header->l3TrigPatternH>0)   hNum->Fill(snr, 1);
  
    // cout << header->eventNumber << " " << truth->maxSNRAtDigitizerH << " " << 1 << endl;
  }

  TEfficiency *eff = 0;

  // Define the TEfficiency from the numerator and denominator
  if (TEfficiency::CheckConsistency(*hNum, *hDenom)){
    eff = new TEfficiency(*hNum, *hDenom);
     
    TCanvas *c = new TCanvas("c");
    eff->Draw();
    
    eff->Write("icemcWAIS");

  }
}
