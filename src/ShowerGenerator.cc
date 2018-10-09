#include "Settings.h"
#include "Geoid.h"
#include "ConnollyEtAl2011.h"
#include "ShowerGenerator.h"
#include "Antarctica.h"
#include "Tools.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

#include "EnvironmentVariable.h"
#include "Report.h"

ClassImp(icemc::Shower)


icemc::ShowerGenerator::ShowerGenerator(const Settings* settings) : fSettings(settings) {
  //For Total Tau Survival probability equation
  //n.b. not in SI units.
  ////from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
  //Equation 16  &  used in Equation 30. 
  /////////////////////|Units////|Description////////////////////////
  B0=1.2*std::pow(10.,-7.); //| m^2/kg  |
  B1=0.16*std::pow(10.,-7.);//| m^2/kg  | }parameterization using a logarithmic dependence on energy for B,
  E0=std::pow(10.,10.);     //| GeV     | the tau elecromagnetic energy loss parameter.
  mT=1.777;	   //| GeV     |Mass of Tau
  cT=0.00008693; 	   //| m       |Tau Decay length (86.93 microMeters)
	                   //|         |
  Mn=1.672622E-24;   //| g       |nucleon/ proton mass in grams,also equal to 0.938 GeV. 
  A=1.;              //| none    |constant that sets the total probability to unity
  //these last two constanst from Connolly Calc 2011, used in d_dzPsurvNu().
  

  SECONDARIES=1; // include secondary interactions
  TAUDECAY=1; // include secondary interactions
  // This is just the initialization, it is set in ReadInputs

  // reading in tauola data file for tau decays
  const std::string ICEMC_SRC_DIR=icemc::EnvironmentVariable::ICEMC_SRC_DIR();
  tauolainfile.open((ICEMC_SRC_DIR+"/data/tau_decay_tauola.dat").c_str(),std::ifstream::in);
  InitTauola();
  
  TAUFRAC=.5; //fraction of tau neutrino-cc current events where the primare interaction point is the first bang   

  count_nfb=0;
  secondary_e_noncons=0;

  for(auto& a : dsdy_muon_brems) a.fill(0);
  for(auto& a : dsdy_muon_epair) a.fill(0);
  for(auto& a : dsdy_muon_pn)    a.fill(0);

  for(auto& a : y_muon_brems) a.fill(0);
  for(auto& a : y_muon_epair) a.fill(0);
  for(auto& a : y_muon_pn)    a.fill(0);

  for(auto& a : dsdy_tauon_brems)     a.fill(0);
  for(auto& a : dsdy_tauon_epair)     a.fill(0);
  for(auto& a : dsdy_tauon_pn)        a.fill(0);
  for(auto& a : dsdy_tauon_hadrdecay) a.fill(0);
  for(auto& a : dsdy_tauon_edecay)    a.fill(0);
  for(auto& a : dsdy_tauon_mudecay)   a.fill(0);

  for(auto& a : y_tauon_brems)     a.fill(0);
  for(auto& a : y_tauon_epair)     a.fill(0);
  for(auto& a : y_tauon_pn)        a.fill(0);
  for(auto& a : y_tauon_hadrdecay) a.fill(0);
  for(auto& a : y_tauon_edecay)    a.fill(0);
  for(auto& a : y_tauon_mudecay)   a.fill(0);

  int_muon_brems.fill(0);
  int_muon_epair.fill(0);
  int_muon_pn.fill(0);
  
  int_tauon_brems.fill(0);
  int_tauon_epair.fill(0);
  int_tauon_pn.fill(0);
  int_tauon_hadrdecay.fill(0);
  int_tauon_edecay.fill(0);
  int_tauon_mudecay.fill(0);

  // Read probability distributions for secondary interactions

  ReadSecondaries();

	
	
}//ShowerGenerator Constructor

// void icemc::ShowerGenerator::readData(std::string nuflavor, std::string secndryType, double (*y)[nP], double (*dsdy)[nP])
// void icemc::ShowerGenerator::readData(const std::string& nuflavor,
// 				      const std::string& secndryType,
// 				      std::array<std::array<double, nP>, nE>& y,
// 				      std::array<std::array<double, nP>, nE>& dsdy)
// 				      // double (*dsdy)[nP])  
// {

void icemc::ShowerGenerator::readData(Neutrino::Flavor flavor, Secondary secondary){
  
  const std::string ICEMC_SRC_DIR=icemc::EnvironmentVariable::ICEMC_SRC_DIR();
  const std::string ICEMC_SECONDARY_DIR=ICEMC_SRC_DIR+"/secondary/";

  std::string nuflavor;
  switch(flavor){
  case Neutrino::Flavor::mu:
    nuflavor = "muons";
    break;
  case Neutrino::Flavor::tau:
    nuflavor = "tauon";
    break;
  default:
    icemc::report() << severity::warning << "No secondary data for " << flavor << std::endl;
  }

  std::string secndryType;
  switch(secondary){
  case Secondary::brems:
    secndryType = "brems";
    break;
  case Secondary::epair:
    secndryType = "epair";
    break;
  case Secondary::pn:
    secndryType = "pn";
    break;
  case Secondary::hadrdecay:
    secndryType = "hadrdecay";
    break;
  case Secondary::edecay:
    secndryType = "edecay";
    break;
  case Secondary::mudecay:
    secndryType = "mudecay";
    break;
  }
  
  std::ifstream in;
  std::string suffix=".vec";
  if(nuflavor=="tauon"){
    suffix="_tau.vec";    
  }

  DataKey key(flavor, secondary);
  TH2F& h = fE_Y_dsdy[key];
  h.SetName(TString::Format("h_%s_%s", nuflavor.c_str(), secndryType.c_str()));
  h.SetTitle(TString::Format("%s %s;Energy (eV); y(?); dsdy(?)", nuflavor.c_str(), secndryType.c_str()));
  h.SetBins(7, 18-0.25, 21.5-0.25, 100, 0, 1);
  
  TH1F& hInt = fE_dsdy[key];
  hInt.SetBins(7, 18-0.25, 21.5-0.25);
  hInt.SetName(TString::Format("hInt_%s_%s", nuflavor.c_str(), secndryType.c_str()));
  hInt.SetTitle(TString::Format("%s %s (Integral over y);Energy (eV); total dsdy(?)", nuflavor.c_str(), secndryType.c_str()));

  TH2F& hC = fE_YCumulative_dsdy[key];
  hC.SetName(TString::Format("hCumulative_%s_%s", nuflavor.c_str(), secndryType.c_str()));
  hC.SetTitle(TString::Format("%s %s (y-CDF);Energy (eV); CDF-y (no units); dsdy(?)", nuflavor.c_str(), secndryType.c_str()));
  hC.SetBins(7, 18-0.25, 21.5-0.25, 100, 0, 1);  

  for(int index=0;index<nE;index++){
    std::stringstream senergy;
    double energy=18+0.5*index;
    int precision=(index%2==0)? 2:3;
    senergy << std::setprecision(precision) << energy;
    std::string path=ICEMC_SECONDARY_DIR+"/"+nuflavor+"/dsdy_"+secndryType+"_1e"+senergy.str()+suffix;
    TGraph gr(path.c_str());
    if(gr.GetN() != 100){
      icemc::report() << severity::warning << "Unexpected file length!"  << std::endl;
    }

    double sum=0;
    for(int i=0; i < gr.GetN(); i++){
      int bin = h.Fill(energy, gr.GetX()[i], gr.GetY()[i]);
      h.SetBinError(bin, 0);

      sum += gr.GetY()[i];      

      bin = hC.Fill(energy, gr.GetX()[i], sum);
      hC.SetBinError(bin, 0);
    }
    int bin2 = hInt.Fill(energy, sum);
    hInt.SetBinError(bin2, 0);

    // normalize cumulative distribution?
    for(int by=0;  by <= hC.GetNbinsY(); by++){
      double val = hC.GetBinContent(index+1, by);
      hC.SetBinContent(index+1, by, val/sum);
    }
    
    // in.open(path.c_str());
    // NPROB=0;
    // while(!in.eof()){
    //   in >> y[index][NPROB] >> dsdy[index][NPROB];
    //   NPROB++;
    //   if(NPROB>=nP){
    // 	// cerr << " ERROR in reading in y_muon_brem. \n";
    // 	break;
    //   }
    // }
    // in.close();
    
  }
}

void icemc::ShowerGenerator::Draw(Option_t* opt){

  std::vector<TCanvas*> cans;
  for(auto& p : fE_Y_dsdy){
    cans.emplace_back(new TCanvas());
    cans.back()->Divide(3);
    
    cans.back()->cd(1);
    p.second.Draw(opt);
    gPad->SetLogz(1);

    cans.back()->cd(2);
    fE_YCumulative_dsdy[p.first].Draw(opt);

    cans.back()->cd(3);    
    fE_dsdy[p.first].Draw();
    // gPad->SetLogz(1);
  }
}

void icemc::ShowerGenerator::ReadSecondaries() {
  // reading in data for secondary interactions
    
  icemc::report() << "Reading in data on secondary interactions." << std::endl;

  // property (y/dsdy), particle (mu, tau), secondary (brem, epair, pn... hadrdecay, edecay, mudecay)

  readData(Neutrino::Flavor::mu, Secondary::brems);
  readData(Neutrino::Flavor::mu, Secondary::epair);
  readData(Neutrino::Flavor::mu, Secondary::pn);

  readData(Neutrino::Flavor::tau, Secondary::brems);
  readData(Neutrino::Flavor::tau, Secondary::epair);
  readData(Neutrino::Flavor::tau, Secondary::pn);
  readData(Neutrino::Flavor::tau, Secondary::hadrdecay);
  readData(Neutrino::Flavor::tau, Secondary::edecay);
  readData(Neutrino::Flavor::tau, Secondary::mudecay);

  // readData("muons", "brems", fY[{Neutrino::Flavor::mu, Secondary::brems}], fDsdy[{Neutrino::Flavor::mu, Secondary::brems}]);
  // readData("muons", "epair", fY[{Neutrino::Flavor::mu, Secondary::epair}], fDsdy[{Neutrino::Flavor::mu, Secondary::epair}]);
  // readData("muons", "pn",    fY[{Neutrino::Flavor::mu, Secondary::pn}],    fDsdy[{Neutrino::Flavor::mu, Secondary::pn}]);
  		    
  // readData("tauon", "brems", fY[{Neutrino::Flavor::tau, Secondary::brems}], fDsdy[{Neutrino::Flavor::tau, Secondary::brems}]);
  // readData("tauon", "epair", fY[{Neutrino::Flavor::tau, Secondary::epair}], fDsdy[{Neutrino::Flavor::tau, Secondary::epair}]);
  // readData("tauon", "pn",    fY[{Neutrino::Flavor::tau, Secondary::pn}],    fDsdy[{Neutrino::Flavor::tau, Secondary::pn}]);
  // readData("tauon", "hadrdecay", fY[{Neutrino::Flavor::tau, Secondary::hadrdecay}],    fDsdy[{Neutrino::Flavor::tau, Secondary::hadrdecay}]);
  // readData("tauon", "edecay", fY[{Neutrino::Flavor::tau, Secondary::edecay}],    fDsdy[{Neutrino::Flavor::tau, Secondary::edecay}]);
  // readData("tauon", "mudecay",fY[{Neutrino::Flavor::tau, Secondary::mudecay}],    fDsdy[{Neutrino::Flavor::tau, Secondary::mudecay}]);

  // readData("muons", "brems", y_muon_brems, dsdy_muon_brems);
  // readData("muons", "epair", y_muon_epair, dsdy_muon_epair);
  // readData("muons", "pn",    y_muon_pn,    dsdy_muon_pn);
  		    
  // readData("tauon", "brems",     y_tauon_brems,     dsdy_tauon_brems);
  // readData("tauon", "epair",     y_tauon_epair,     dsdy_tauon_epair);
  // readData("tauon", "pn",        y_tauon_pn,        dsdy_tauon_pn);
  // readData("tauon", "hadrdecay", y_tauon_hadrdecay, dsdy_tauon_hadrdecay);
  // readData("tauon", "edecay",    y_tauon_edecay,    dsdy_tauon_edecay);
  // readData("tauon", "mudecay",   y_tauon_mudecay,   dsdy_tauon_mudecay);
  

  //cout << "NPROB=" << NPROB << ",  nP=" << nP << endl;

  // auto integrate = [&] (const std::array<std::array<double, nP>, nE>& toIntegrate, std::array<double, nE>& result)->void
  //   {     
  //    for(int j=0; j < nE;  j++){
  //       result.at(j) = std::accumulate(toIntegrate.begin(), toIntegrate.end(), 0.0);
  //    }
  //    return;
  //   };

  // integrate(dsdy_muon_brems,      int_muon_brems);
  // integrate(dsdy_muon_epair,      int_muon_epair);
  // integrate(dsdy_muon_pn,         int_muon_pn);
  // integrate(dsdy_tauon_brems,     int_tauon_brems);
  // integrate(dsdy_tauon_epair,     int_tauon_epair);
  // integrate(dsdy_tauon_pn,        int_tauon_pn);
  // integrate(dsdy_tauon_brems,     int_tauon_brems);
  // integrate(dsdy_tauon_epair,     int_tauon_epair);
  // integrate(dsdy_tauon_pn,        int_tauon_pn);
  // integrate(dsdy_tauon_hadrdecay, int_tauon_hadrdecay);
  // integrate(dsdy_tauon_edecay,    int_tauon_edecay);
  // integrate(dsdy_tauon_mudecay,   int_tauon_mudecay);
  
  // for(int j=0;j<nE;j++) {
  //   // integrating prob. distributions.
  //   // int_muon_brems[j]=std::accumulate(dsdy_muon_brems[j],dsdy_muon_brems[j]+nP,0.);//very important to keep the initial value the same type as the elements type
  //   // int_muon_epair[j]=std::accumulate(dsdy_muon_epair[j],dsdy_muon_epair[j]+nP,0.);
  //   // int_muon_pn[j]=std::accumulate(dsdy_muon_pn[j],dsdy_muon_pn[j]+nP,0.);
  //   // int_tauon_brems[j]=std::accumulate(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+nP,0.);
  //   // int_tauon_epair[j]=std::accumulate(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+nP,0.);
  //   // int_tauon_pn[j]=std::accumulate(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+nP,0.);
  //   // int_tauon_hadrdecay[j]=std::accumulate(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+nP,0.);
  //   // int_tauon_edecay[j]=std::accumulate(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+nP,0.);
  //   // int_tauon_mudecay[j]=std::accumulate(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+nP,0.);
    

  //   // maximum value of prob. dist.
  //   max_muon_brems = *std::max_element(dsdy_muon_brems[j],dsdy_muon_brems[j]+nP);
  //   //cout << "max_muon_brems=" << max_muon_brems << endl;//fenfang
  //   max_muon_epair = *std::max_element(dsdy_muon_epair[j],dsdy_muon_epair[j]+nP);
  //   max_muon_pn = *std::max_element(dsdy_muon_pn[j],dsdy_muon_pn[j]+nP);   
  //   max_tauon_brems = *std::max_element(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+nP);
  //   max_tauon_epair = *std::max_element(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+nP);
  //   max_tauon_pn = *std::max_element(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+nP);
  //   max_tauon_hadrdecay = *std::max_element(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+nP);
  //   max_tauon_edecay = *std::max_element(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+nP);
  //   max_tauon_mudecay = *std::max_element(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+nP);
     
  //   // minimum value of prob. dist.
  //   min_muon_brems=Tools::dMinNotZero(dsdy_muon_brems[j],nP);
  //   min_muon_epair=Tools::dMinNotZero(dsdy_muon_epair[j],nP);
  //   min_muon_pn=Tools::dMinNotZero(dsdy_muon_pn[j],nP);   
  //   min_tauon_brems=Tools::dMinNotZero(dsdy_tauon_brems[j],nP);
  //   min_tauon_epair=Tools::dMinNotZero(dsdy_tauon_epair[j],nP);
  //   min_tauon_pn=Tools::dMinNotZero(dsdy_tauon_pn[j],nP);
  //   min_tauon_hadrdecay=Tools::dMinNotZero(dsdy_tauon_hadrdecay[j],nP);
  //   min_tauon_edecay=Tools::dMinNotZero(dsdy_tauon_edecay[j],nP);
  //   min_tauon_mudecay=Tools::dMinNotZero(dsdy_tauon_mudecay[j],nP);
     
  //   if (min_muon_brems<=0){
  //     std::cout << "Minimum probability is <=0!\n";
  //   }    
  //   std::partial_sum(dsdy_muon_brems[j],dsdy_muon_brems[j]+nP,y_cumulative_muon_brems[j]);
  //   std::partial_sum(dsdy_muon_epair[j],dsdy_muon_epair[j]+nP,y_cumulative_muon_epair[j]);
  //   std::partial_sum(dsdy_muon_pn[j],dsdy_muon_pn[j]+nP,y_cumulative_muon_pn[j]);
  //   std::partial_sum(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+nP,y_cumulative_tauon_brems[j]);
  //   std::partial_sum(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+nP,y_cumulative_tauon_epair[j]);
  //   std::partial_sum(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+nP,y_cumulative_tauon_pn[j]);
  //   std::partial_sum(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+nP,y_cumulative_tauon_hadrdecay[j]);
  //   std::partial_sum(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+nP,y_cumulative_tauon_mudecay[j]);
  //   std::partial_sum(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+nP,y_cumulative_tauon_edecay[j]);
     
  //   for (int i=0;i<nP;i++) {
  //     y_cumulative_muon_brems[j][i]      /= y_cumulative_muon_brems[j][nP-1];
  //     y_cumulative_muon_epair[j][i]      /= y_cumulative_muon_epair[j][nP-1];
  //     y_cumulative_muon_pn[j][i]         /= y_cumulative_muon_pn[j][nP-1];
  //     y_cumulative_tauon_brems[j][i]     /= y_cumulative_tauon_brems[j][nP-1];
  //     y_cumulative_tauon_epair[j][i]     /= y_cumulative_tauon_epair[j][nP-1];
  //     y_cumulative_tauon_pn[j][i]        /= y_cumulative_tauon_pn[j][nP-1];
  //     y_cumulative_tauon_hadrdecay[j][i] /= y_cumulative_tauon_hadrdecay[j][nP-1];
  //     y_cumulative_tauon_mudecay[j][i]   /= y_cumulative_tauon_mudecay[j][nP-1];
  //     y_cumulative_tauon_edecay[j][i]    /= y_cumulative_tauon_edecay[j][nP-1];
  //   } //for

  // }
  std::cout<<"Finished reading secondary interaction data.\n"; 
} //end method ReadSecondaries


void icemc::ShowerGenerator::doShower(Neutrino::Flavor nuflavor, double plepton, double &em_secondaries_max, double &had_secondaries_max,int &n_interactions, TH1F *hy) {

  em_secondaries_max=0.;
  had_secondaries_max=0.;
  
  double energy = plepton;

  // ok it it looks like we're doing a lookup table with i.
  // this could probably be improved.
  int i=(int)((log10(plepton)-18.)*2.);

  if (i>6){
    i=6;
  }
  else if (i<0){
    i=0;
  }

  // int n_brems,n_epair,n_pn; // number of interactions of each type.
  // int index_y; // index along the horizontal axis of ped's plots
  // double rnd1=1000.;
  // double rnd2=1000.;  // random numbers for throwing at dart board
  double y = 0; // inelasticity
 
  // std::string whichtype; // which type of interaction corresponds to that index
  Secondary whichtype;
  
  if (nuflavor==Neutrino::Flavor::mu){
    /**
     * Pick poisson distribution for each, uses i.
     */
    std::vector<Secondary>secondaries;
    
    for(auto sec : {Secondary::brems, Secondary::epair, Secondary::pn}){
      DataKey key(nuflavor, sec);
      TH1F& h_E_dsdy = fE_dsdy[key];
      
      int nPoisson = pickPoisson(h_E_dsdy.Interpolate(energy));
      for(int k=0; k < nPoisson; k++){
      	secondaries.emplace_back(sec);
      }
    }
    // then we can sort them to get our particles in a random order
    // std::sort(secondaries.begin(),  secondaries.end());
    
    std::random_shuffle(secondaries.begin(), secondaries.end(),
			[&](int z){ return int(pickUniform(z))%z; } );

    // secondaryCount[Secondary::brems];// = 
    // secondaryCount[Secondary::epair] = pickPoisson(int_muon_epair[i]);
    // secondaryCount[Secondary::pn] = pickPoisson(int_muon_pn[i]);

    // int n_brems = pickPoisson(int_muon_brems[i]); // pick number of brem interactions
    // int n_epair = pickPoisson(int_muon_epair[i]); // # of pair production
    // int n_pn    = pickPoisson(int_muon_pn[i]); // # photonuclear interactions
    // int n_total = n_brems + n_epair + n_pn;

    
    for(auto secondary : secondaries){
      const double* d = nullptr;

      switch(secondary){
      case Secondary::brems: d = y_cumulative_muon_brems[i].data();  break;
      case Secondary::epair: d = y_cumulative_muon_epair[i].data();  break;
      case Secondary::pn:    d = y_cumulative_muon_pn[i].data();     break;
      }
      // double rnd1 = pickUniform(secondaries.size());

      DataKey key(nuflavor, secondary);
      fE_YCumulative_dsdy[key];
      // this is stupid...
      // Picky(d, NPROB, rnd1, y);

      if (y*plepton>std::max(em_secondaries_max,had_secondaries_max)) {  // if this is the largest interaction for this event so far
	if (secondary == Secondary::brems || secondary == Secondary::epair) {  // save it
	  em_secondaries_max=y*plepton;
	}
	if (secondary == Secondary::pn) {
	  had_secondaries_max=y*plepton;
	}
      }      
    }
    
    const int n_total = secondaries.size();
    const int n_brems = std::count(secondaries.begin(), secondaries.end(), Secondary::brems);
    const int n_epair = std::count(secondaries.begin(), secondaries.end(), Secondary::epair);
    // we generated n_total particles, loop through them
    for (int j=0; j<n_total; j++) {

      double rnd1 = pickUniform(n_total);
      
      if (rnd1<=(double)n_brems/(double)(n_total)){
	whichtype= Secondary::brems;
      }
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_total)){
	whichtype= Secondary::epair;
      }
      else{
	whichtype= Secondary::pn;
      }
      // rnd2=1000.;  // random numbers for throwing at dart board
      // index_y=0;

      if (whichtype== Secondary::brems) {	
	rnd1 = pickUniform();
	Picky(&y_cumulative_muon_brems[i][0],NPROB,rnd1,y);
      }
      else if (whichtype== Secondary::epair) {
	rnd1 = pickUniform();
	Picky(&y_cumulative_muon_epair[i][0],NPROB,rnd1,y);
      }
      else if (whichtype== Secondary::pn) {
	rnd1 = pickUniform();
	Picky(&y_cumulative_muon_pn[i][0],NPROB,rnd1,y);
      }
     
      if (y*plepton>std::max(em_secondaries_max,had_secondaries_max)) {  // if this is the largest interaction for this event so far
	if (whichtype== Secondary::brems || whichtype== Secondary::epair) {  // save it
	  em_secondaries_max=y*plepton;
	}
	if (whichtype== Secondary::pn) {
	  had_secondaries_max=y*plepton;
	}
      }
    } // loop over secondary interactions
  } // end if it was a muon neutrino

  if (nuflavor == Neutrino::Flavor::tau) {
    int n_brems = pickPoisson(int_tauon_brems[i]);
    int n_epair = pickPoisson(int_tauon_epair[i]);
    int n_pn =    pickPoisson(int_tauon_pn[i]);

    n_interactions += (n_brems+n_epair+n_pn); // increment number of secondary interactions.

    for (int j=0;j<n_brems+n_epair+n_pn;j++) { // loop over secondary interactions. 
      
      double rnd1 = pickUniform();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn)){
	whichtype= Secondary::brems;
      }
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn)){
	whichtype= Secondary::epair;
      }
      else{
	whichtype= Secondary::pn;
      }
      
      rnd1=1000.;
      // rnd2=1000.;  // random numbers for throwing at dart board
      // index_y=0;

      if (whichtype== Secondary::brems) {  // bremstrahlung interaction
	rnd1 = pickUniform();
	Picky(&y_cumulative_tauon_brems[i][0],NPROB,rnd1,y);
      }
      if (whichtype== Secondary::epair) { // pair production
	rnd1 = pickUniform();
	Picky(&y_cumulative_tauon_epair[i][0],NPROB,rnd1,y);
      }
      if (whichtype== Secondary::pn) {
	rnd1 = pickUniform();
	Picky(&y_cumulative_tauon_pn[i][0],NPROB,rnd1,y);
      }

      if (fSettings->HIST==1 && !fSettings->ONLYFINAL && hy->GetEntries()<fSettings->HIST_MAX_ENTRIES){
	hy->Fill(y);
      }

      if (y*plepton>std::max(em_secondaries_max,had_secondaries_max)) { // if this is the biggest secondary signal yet,
	if (whichtype== Secondary::brems || whichtype== Secondary::epair){ // save it.
	  em_secondaries_max=y*plepton;
	}
	if (whichtype== Secondary::pn){
	  had_secondaries_max=y*plepton;
	}
      }
    }
   

    if (TAUDECAY) {
      n_interactions++; // increment number of interactions, for plotting.

      double rnd1 = pickUniform();
      if (rnd1<0.65011){  // pick which type of decay it is.
	whichtype = Secondary::hadrdecay;
      }
      if (rnd1>=0.65011 && rnd1<0.8219){
	whichtype = Secondary::mudecay;
      }
      if (rnd1>=0.8219){
	whichtype = Secondary::edecay;
      }
      
      rnd1=1000.;
      // rnd2=1000.;  // random numbers for throwing at dart board
      // index_y=0;     
      
      if (whichtype== Secondary::hadrdecay) { // hadronic decay
	rnd1 = pickUniform();
	Picky(&y_cumulative_tauon_hadrdecay[i][0], NPROB, rnd1, y);	
      }
      else if (whichtype== Secondary::edecay) { // e decay	
	rnd1 = pickUniform();
	Picky(&y_cumulative_tauon_edecay[i][0], NPROB, rnd1, y);
      }
      else if (whichtype== Secondary::mudecay) { // mu decay
	rnd1 = pickUniform();
	Picky(&y_cumulative_tauon_mudecay[i][0], NPROB, rnd1, y);
      }
      
     
      if (y*plepton>std::max(em_secondaries_max, had_secondaries_max)) {  // if this is the biggest interaction yet,    
	if (whichtype== Secondary::edecay){ // save it.
	  em_secondaries_max=y*plepton;
	}
	if (whichtype== Secondary::hadrdecay){
	  had_secondaries_max=y*plepton;
	}
      } //if
    } //if (TAUDECAY)
  } //if (nutau)

} //GetShowerGenerator



icemc::Shower icemc::ShowerGenerator::generate(const Neutrino& nu){
  Shower s;
  
  // doShower(nu.flavor,nu.energy,em_secondaries_max,had_secondaries_max,s.nInteractions,hy);   
  return s;
}


icemc::Shower icemc::ShowerGenerator::GetEMFrac(Neutrino::Flavor nuflavor,
						Neutrino::Interaction::Current current,
						// const string& nuflavor,
						// const string& current,
						const std::string& taudecay,
						double y,
						TH1F *hy,
						double pnu,
						int inu,
						// double& emfrac,
						// double& hadfrac,
						// int& n_interactions,
						int taumodes1) {


  Shower s;
  s.pnu = pnu;

  if (current==Neutrino::Interaction::Current::Charged){
    plepton=(1.-y)*pnu;
  }
  else{
    plepton=0.;
  }


  const double negligible = 1e-10;
  if (nuflavor==Neutrino::Flavor::e && current==Neutrino::Interaction::Current::Charged) {
    s.emFrac=1.-y;
    s.hadFrac=y;
  }
  else if(nuflavor==Neutrino::Flavor::mu && current==Neutrino::Interaction::Current::Charged) {
    s.emFrac=negligible;
    s.hadFrac=y;
  }
  else if(nuflavor==Neutrino::Flavor::tau && current==Neutrino::Interaction::Current::Charged) {
    // behaves like a muon
    if(taumodes1 ==1){//taumodes==1; tau created somewhere in rock and decays at posnu.
      this->pickEMFracDB(s.emFrac,s.hadFrac);
    }
    else if (taumodes1 == 0){
      s.emFrac=negligible;
      s.hadFrac=y;
    }
  }
  else if (current==Neutrino::Interaction::Current::Neutral) {
    s.emFrac=negligible;
    s.hadFrac=y;
  }


  em_secondaries_max = s.emFrac; // initialize search for maximum signal among primary, secondary interactions.
  had_secondaries_max = s.hadFrac;

  if (SECONDARIES==1 && current==Neutrino::Interaction::Current::Charged) {

    while (1) {

      // find how much em and hadronic energies comes from secondary interactions.
      // keep picking until you get a bunch of secondary interactions that conserve energy
      doShower(nuflavor,plepton,em_secondaries_max,had_secondaries_max,s.nInteractions,hy); 

      if (em_secondaries_max+had_secondaries_max<=plepton*(1.+1.E-5)) // if conserves energy, break.
	break;
      else {
	secondary_e_noncons++; //Record how many times we come up with something that doesn't conserve energy
	em_secondaries_max=s.emFrac;
	had_secondaries_max=s.hadFrac;
      } //else
    } //while(1)

    if ((em_secondaries_max+had_secondaries_max)>(s.sumFrac())*pnu) { // if maximum signal from secondaries is larger than
                                                                         // signal from primary interaction
      s.emFrac=em_secondaries_max/pnu; // then use that one.
      s.hadFrac=had_secondaries_max/pnu;
      if (s.emFrac <= negligible){
	s.emFrac=negligible;
      }
      if (s.hadFrac<= negligible){
	s.hadFrac=negligible;
      }
    } //if
  } //if (charged current, secondaries on)

  if (nuflavor==Neutrino::Flavor::mu && current==Neutrino::Interaction::Current::Charged && s.nInteractions==0){
    icemc::report() << severity::warning << "Look at this one.  inu is " << inu << "\n";
  }  

  if ((y<0 || y>1) && y != -999.) {
    icemc::report() << severity::error <<  "Illegal value of y, y =" << y << "\n";
  }
  
  if (s.sumFrac()>1.00001) {
    icemc::report() << severity::error << "emFrac,hadfrac=" << s.emFrac << "," << s.hadFrac << ": sum = " << s.sumFrac() << "\n";
    icemc::report() << "nuflavor,taudecay=" << nuflavor << " " << taudecay << "\n";
  }
  
  return s;

} //GetEMFrac


//----------------------------------------------------------
//InitTauola()
//Initializes the tau decay information

void icemc::ShowerGenerator::InitTauola() {
  for(int k=0;k<5;k++)
    tauolainfile >> tauola[0][k];
  for(int i=1;i<N_TAUOLA;i++)
    for(int j=0;j<6;j++)
      tauolainfile >> tauola[i][j];

  return;
}//InitTauola


void icemc::ShowerGenerator::GetTauDecay(Neutrino::Flavor nuflavor, Neutrino::Interaction::Current current, std::string& taudecay, double& emfrac_db, double& hadfrac_db) {
 
  if (!(nuflavor==Neutrino::Flavor::tau || current==Neutrino::Interaction::Current::Charged || interestedintaus)){
    return;
  }
  
  // if nu_tau choose tau decay type
  
  double rnd = pickUniform();
  int decay = static_cast<int>(rnd*(N_TAUOLA-2)+1);
  
  hadfrac_db=tauola[decay][3];
  emfrac_db=tauola[decay][4];
  
  if(tauola[decay][1]!=0)
    taudecay="m";
  else
    taudecay="e";
 

  if(taudecay=="m")
    secondbang=false; //for all muon decays, the interaction point chosen is the neutrino interaction since we don't detect the decay if
  //the tau decays into a muon.
  else {
    double rnd=pickUniform();
    if(rnd>TAUFRAC) {
      secondbang=true;
      count_nfb++;
    } else
      secondbang=false;
  }
  

} //GetTauDecay

//-----------------------------------------------------
//GetEMFracDB()
//Gets the emfrac_db and hadfrac_db for a doublebang

void icemc::ShowerGenerator::pickEMFracDB(double& emfrac_db, double& hadfrac_db) {


  double rnd = pickUniform();
  int decay = static_cast<int>(rnd*(N_TAUOLA-2)+1);

  hadfrac_db=tauola[decay][3];
  emfrac_db=tauola[decay][4];

  return;
}//GetEMFracDB

//------------------------------------------------------
//GetDBViewAngle()
//Gets the viewangle of the second bang

double icemc::ShowerGenerator::GetDBViewAngle(const TVector3 &refr, const TVector3 &nnu) {

  /// @todo FIX CHANGECOORD
  // return ((nnu.ChangeCoord(refr)).Angle(z_axis));
  return nnu.Angle(TVector3(0, 0, 1));

}//GetDBViewAngle

//------------------------------------------------------
//GetFirstBang()
//Gets the position of the first bang when the interaction point is the tau decay point

//  void icemc::ShowerGenerator::GetFirstBang(const Geoid::Position &r_in, const TVector3 &nnu, Geoid::Position &posnu, double len_int_kgm2, double chord, double &nuentrancelength) {
  
//   double weightbang;
//   double junk1;
//   double junk2;
//   double junk3;
//   int junk4,junk5,junk6;
//   double myair=0;

//   TVector3 r_out = r_in + chord*nnu;

//   antarctica->Getchord(len_int_kgm2,r_in,r_out,
// 		  junk1,weightbang,junk2,myair,junk3,junk4,junk5,junk6);
//   double r1,r2;
//   if(weightbang>.999)
//     r2=pickUniform()*chord;
//   else {
//     do {
//       r1=pickUniform();
//       r2=pickUniform()*chord;
//       r_out = r_in + r2*nnu;
//       antarctica->Getchord(len_int_kgm2,r_in,r_out,
// 		      junk1,weightbang,junk2,myair,junk3,junk4,junk5,junk6);
//     }
//     while(r1>1-weightbang);
//   }
//   posnu = r_in + r2*nnu;
//   nuentrancelength=r2;

//   return;
// }//GetFirstBang

//---------------------------------------------------------
//NFBWeight()
//Gets the weight of the tau decay for second bang events
double icemc::ShowerGenerator::NFBWeight(double ptau, double taulength) {
  
  double gamma=ptau/constants::MTAU;
  double D=constants::TAUDECAY_TIME*constants::CLIGHT*gamma;

  return exp(-taulength/D);

}

void icemc::ShowerGenerator::Picky(const double *y_cumulative,int NPROB,double rnd,double& y) const {
  for (int i=0;i<NPROB;i++) {
    if (y_cumulative[i] <= rnd && y_cumulative[i+1] > rnd) {
      y=(double)i/(double)NPROB;
      continue; // once you found the right bin, stop looping.
    } //if
  } //for
} //Picky

