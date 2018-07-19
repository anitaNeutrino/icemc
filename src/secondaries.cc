#include "TVector3.h"
#include "TRandom3.h"
#include "Settings.h"
#include "TVector3.h"
#include "GeoidModel.h"
#include "Primaries.h"
#include "secondaries.hh"
#include "Antarctica.h"
#include "Tools.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

#include "TH1F.h"
#include "Constants.h"
#include "TTreeIndex.h"
#include "TChain.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h" 
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"

#include "EnvironmentVariable.h"
#include "IcemcLog.h"

ClassImp(icemc::ShowerProperties)

const std::string ICEMC_SRC_DIR=icemc::EnvironmentVariable::ICEMC_SRC_DIR();
const std::string ICEMC_SECONDARY_DIR=ICEMC_SRC_DIR+"/secondary/";


icemc::Secondaries::Secondaries() {
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
  
  flavors[0]="nue";
  flavors[1]="numu";
  flavors[2]="nutau"; // the gps path of the anita-lite flight

  SECONDARIES=1; // include secondary interactions
  TAUDECAY=1; // include secondary interactions
  // This is just the initialization, it is set in ReadInputs


  // reading in tauola data file for tau decays
  tauolainfile.open((ICEMC_SRC_DIR+"/data/tau_decay_tauola.dat").c_str(),std::ifstream::in);
  InitTauola();
  
  TAUFRAC=.5; //fraction of tau neutrino-cc current events where the primare interaction point is the first bang   

  


  count_nfb=0;
  secondary_e_noncons=0;

  for (int i=0;i<7;i++) {
    Tools::Zero(dsdy_muon_brems[i],NPROB_MAX);
    Tools::Zero(dsdy_muon_epair[i],NPROB_MAX);
    Tools::Zero(dsdy_muon_pn[i],NPROB_MAX);
    
    Tools::Zero(y_muon_brems[i],NPROB_MAX);
    Tools::Zero(y_muon_epair[i],NPROB_MAX);
    Tools::Zero(y_muon_pn[i],NPROB_MAX);
    
    Tools::Zero(dsdy_tauon_brems[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_epair[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_pn[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_hadrdecay[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_edecay[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_mudecay[i],NPROB_MAX);
  


    Tools::Zero(y_tauon_brems[i],NPROB_MAX);
    Tools::Zero(y_tauon_epair[i],NPROB_MAX);
    Tools::Zero(y_tauon_pn[i],NPROB_MAX);
    Tools::Zero(y_tauon_hadrdecay[i],NPROB_MAX);
    Tools::Zero(y_tauon_edecay[i],NPROB_MAX);
    Tools::Zero(y_tauon_mudecay[i],NPROB_MAX);
  } //for (Tools::Zeroing)

  Tools::Zero(int_muon_brems,7);
  Tools::Zero(int_muon_epair,7);
  Tools::Zero(int_muon_pn,7);
  
  Tools::Zero(int_tauon_brems,7);
  Tools::Zero(int_tauon_epair,7);
  Tools::Zero(int_tauon_pn,7);
  Tools::Zero(int_tauon_hadrdecay,7);
  Tools::Zero(int_tauon_edecay,7);
  Tools::Zero(int_tauon_mudecay,7);

  // Read probability distributions for secondary interactions

  ReadSecondaries();

	
	
}//Secondaries Constructor

void icemc::Secondaries::readData(std::string nuflavor,std::string secndryType, double (*y)[NPROB_MAX], double (*dsdy)[NPROB_MAX])
{
  
  std::stringstream senergy;
  
  std::ifstream ifile;
  std::string suffix=".vec";
  if(nuflavor=="tauon")
    suffix="_tau.vec";
  
  for(int index=0;index<7;index++)
    {senergy.str("");
      double energy=18+0.5*index;
      int precision=(index%2==0)?2:3;
      senergy << std::setprecision(precision) << energy;
      std::string path=ICEMC_SECONDARY_DIR+"/"+nuflavor+"/dsdy_"+secndryType+"_1e"+senergy.str()+suffix;
      //cout << "openning file " << path.c_str() << endl;
      ifile.open(path.c_str());
      NPROB=0;
      while(!ifile.eof())
	{
	  ifile >> y[index][NPROB] >> dsdy[index][NPROB];
	  NPROB++;
	  if(NPROB>=NPROB_MAX)
	    {
	      // cerr << " ERROR in reading in y_muon_brem. \n";
	      break;
	    }
	 
	}
      ifile.close();
    }
  
}

void icemc::Secondaries::ReadSecondaries() {
  // reading in data for secondary interactions
  
  std::cout<<"Reading in data on secondary interactions.\n";

  readData("muons","brems",y_muon_brems,dsdy_muon_brems);
  readData("muons","epair",y_muon_epair,dsdy_muon_epair);
  readData("muons","pn",y_muon_pn,dsdy_muon_pn);
  readData("tauon","brems",y_tauon_brems,dsdy_tauon_brems);
  readData("tauon","epair",y_tauon_epair,dsdy_tauon_epair);
  readData("tauon","pn",y_tauon_pn,dsdy_tauon_pn);
  readData("tauon","hadrdecay",y_tauon_hadrdecay,dsdy_tauon_hadrdecay);
  readData("tauon","edecay",y_tauon_edecay,dsdy_tauon_edecay);
  readData("tauon","mudecay",y_tauon_mudecay,dsdy_tauon_mudecay);
  //cout << "NPROB=" << NPROB << ",  NPROB_MAX=" << NPROB_MAX << endl;
  for(int j=0;j<7;j++) {
    // integrating prob. distributions.
    int_muon_brems[j]=std::accumulate(dsdy_muon_brems[j],dsdy_muon_brems[j]+NPROB_MAX,0.);//very important to keep the initial value the same type as the elements type
    int_muon_epair[j]=std::accumulate(dsdy_muon_epair[j],dsdy_muon_epair[j]+NPROB_MAX,0.);
    int_muon_pn[j]=std::accumulate(dsdy_muon_pn[j],dsdy_muon_pn[j]+NPROB_MAX,0.);
    int_tauon_brems[j]=std::accumulate(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+NPROB_MAX,0.);
    int_tauon_epair[j]=std::accumulate(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+NPROB_MAX,0.);
    int_tauon_pn[j]=std::accumulate(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+NPROB_MAX,0.);
    int_tauon_hadrdecay[j]=std::accumulate(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+NPROB_MAX,0.);
    int_tauon_edecay[j]=std::accumulate(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+NPROB_MAX,0.);
    int_tauon_mudecay[j]=std::accumulate(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+NPROB_MAX,0.);
    
    // maximum value of prob. dist.
    max_muon_brems=*std::max_element(dsdy_muon_brems[j],dsdy_muon_brems[j]+NPROB_MAX);
    //cout << "max_muon_brems=" << max_muon_brems << endl;//fenfang
    max_muon_epair=*std::max_element(dsdy_muon_epair[j],dsdy_muon_epair[j]+NPROB_MAX);
    max_muon_pn=*std::max_element(dsdy_muon_pn[j],dsdy_muon_pn[j]+NPROB_MAX);   
    max_tauon_brems=*std::max_element(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+NPROB_MAX);
    max_tauon_epair=*std::max_element(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+NPROB_MAX);
    max_tauon_pn=*std::max_element(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+NPROB_MAX);
    max_tauon_hadrdecay=*std::max_element(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+NPROB_MAX);
    max_tauon_edecay=*std::max_element(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+NPROB_MAX);
    max_tauon_mudecay=*std::max_element(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+NPROB_MAX);
     
    // minimum value of prob. dist.
    min_muon_brems=Tools::dMinNotZero(dsdy_muon_brems[j],NPROB_MAX);
    min_muon_epair=Tools::dMinNotZero(dsdy_muon_epair[j],NPROB_MAX);
    min_muon_pn=Tools::dMinNotZero(dsdy_muon_pn[j],NPROB_MAX);   
    min_tauon_brems=Tools::dMinNotZero(dsdy_tauon_brems[j],NPROB_MAX);
    min_tauon_epair=Tools::dMinNotZero(dsdy_tauon_epair[j],NPROB_MAX);
    min_tauon_pn=Tools::dMinNotZero(dsdy_tauon_pn[j],NPROB_MAX);
    min_tauon_hadrdecay=Tools::dMinNotZero(dsdy_tauon_hadrdecay[j],NPROB_MAX);
    min_tauon_edecay=Tools::dMinNotZero(dsdy_tauon_edecay[j],NPROB_MAX);
    min_tauon_mudecay=Tools::dMinNotZero(dsdy_tauon_mudecay[j],NPROB_MAX);
     
    if (min_muon_brems<=0)
      std::cout << "Minimum probability is <=0!\n";
    
    std::partial_sum(dsdy_muon_brems[j],dsdy_muon_brems[j]+NPROB_MAX,y_cumulative_muon_brems[j]);
    std::partial_sum(dsdy_muon_epair[j],dsdy_muon_epair[j]+NPROB_MAX,y_cumulative_muon_epair[j]);
    std::partial_sum(dsdy_muon_pn[j],dsdy_muon_pn[j]+NPROB_MAX,y_cumulative_muon_pn[j]);
    std::partial_sum(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+NPROB_MAX,y_cumulative_tauon_brems[j]);
    std::partial_sum(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+NPROB_MAX,y_cumulative_tauon_epair[j]);
    std::partial_sum(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+NPROB_MAX,y_cumulative_tauon_pn[j]);
    std::partial_sum(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+NPROB_MAX,y_cumulative_tauon_hadrdecay[j]);
    std::partial_sum(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+NPROB_MAX,y_cumulative_tauon_mudecay[j]);
    std::partial_sum(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+NPROB_MAX,y_cumulative_tauon_edecay[j]);
     
    for (int i=0;i<NPROB_MAX;i++) {
      y_cumulative_muon_brems[j][i]      /= y_cumulative_muon_brems[j][NPROB_MAX-1];
      y_cumulative_muon_epair[j][i]      /= y_cumulative_muon_epair[j][NPROB_MAX-1];
      y_cumulative_muon_pn[j][i]         /= y_cumulative_muon_pn[j][NPROB_MAX-1];
      y_cumulative_tauon_brems[j][i]     /= y_cumulative_tauon_brems[j][NPROB_MAX-1];
      y_cumulative_tauon_epair[j][i]     /= y_cumulative_tauon_epair[j][NPROB_MAX-1];
      y_cumulative_tauon_pn[j][i]        /= y_cumulative_tauon_pn[j][NPROB_MAX-1];
      y_cumulative_tauon_hadrdecay[j][i] /= y_cumulative_tauon_hadrdecay[j][NPROB_MAX-1];
      y_cumulative_tauon_mudecay[j][i]   /= y_cumulative_tauon_mudecay[j][NPROB_MAX-1];
      y_cumulative_tauon_edecay[j][i]    /= y_cumulative_tauon_edecay[j][NPROB_MAX-1];
    } //for

  }
  std::cout<<"Finished reading secondary interaction data.\n"; 
  
 
} //end method ReadSecondaries


void icemc::Secondaries::GetSecondaries(const Settings *settings1, NuFlavor nuflavor, double plepton, double &em_secondaries_max, double &had_secondaries_max,int &n_interactions, TH1F *hy) {

  em_secondaries_max=0.;
  had_secondaries_max=0.;

  int i=(int)((log10(plepton)-18.)*2.);
  if (i>6){
    i=6;
  }
  else if (i<0){
    i=0;
  }

  int n_brems,n_epair,n_pn; // number of interactions of each type.
  // int index_y; // index along the horizontal axis of ped's plots
  double rnd1=1000.;
  // double rnd2=1000.;  // random numbers for throwing at dart board
  double y = 0; // inelasticity
 
  std::string whichtype; // which type of interaction corresponds to that index
  


  if (nuflavor==NuFlavor::mu) {
    n_brems=gRandom->Poisson(int_muon_brems[i]); // pick number of brem interactions
    n_epair=gRandom->Poisson(int_muon_epair[i]); // # of pair production
    n_pn=gRandom->Poisson(int_muon_pn[i]); // # photonuclear interactions   
    
    n_interactions+=(n_brems+n_epair+n_pn);	


    for (int j=0;j<n_brems+n_epair+n_pn;j++) {
      rnd1=gRandom->Rndm();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn))
	whichtype="brems";
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn))
	whichtype="epair";
      else
	whichtype="pn";

      rnd1=1000.;
      // rnd2=1000.;  // random numbers for throwing at dart board
      // index_y=0;

      if (whichtype=="brems") {	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_brems[i],NPROB,rnd1,y);
      }
      else if (whichtype=="epair") {	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_epair[i],NPROB,rnd1,y);	
      }
      else if (whichtype=="pn") {
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_pn[i],NPROB,rnd1,y);
      }
     
      if (y*plepton>std::max(em_secondaries_max,had_secondaries_max)) {  // if this is the largest interaction for this event so far
	if (whichtype=="brems" || whichtype=="epair") {  // save it
	  em_secondaries_max=y*plepton;

	}
	if (whichtype=="pn") 
	  had_secondaries_max=y*plepton;
	 
		
      }
    } // loop over secondary interactions
  } // end if it was a muon neutrino
  if (nuflavor==NuFlavor::tau) {
    n_brems=gRandom->Poisson(int_tauon_brems[i]);
    n_epair=gRandom->Poisson(int_tauon_epair[i]);
    n_pn=gRandom->Poisson(int_tauon_pn[i]);

    n_interactions+=(n_brems+n_epair+n_pn); // increment number of secondary interactions.

    for (int j=0;j<n_brems+n_epair+n_pn;j++) { // loop over secondary interactions. 
      
      rnd1=gRandom->Rndm();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn))
	whichtype="brems";
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn))
	whichtype="epair";
      else
	whichtype="pn";
  
      rnd1=1000.;
      // rnd2=1000.;  // random numbers for throwing at dart board
      // index_y=0;

      if (whichtype=="brems") {  // bremstrahlung interaction
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_brems[i],NPROB,rnd1,y);
      }
      if (whichtype=="epair") { // pair production
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_epair[i],NPROB,rnd1,y);
      }
      if (whichtype=="pn") {
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_pn[i],NPROB,rnd1,y);
      }

      if (settings1->HIST==1 && !settings1->ONLYFINAL && hy->GetEntries()<settings1->HIST_MAX_ENTRIES)
	hy->Fill(y);
      if (y*plepton>std::max(em_secondaries_max,had_secondaries_max)) { // if this is the biggest secondary signal yet,
	if (whichtype=="brems" || whichtype=="epair") // save it.
	  em_secondaries_max=y*plepton;
	if (whichtype=="pn")
	  had_secondaries_max=y*plepton;
      }
    }
   

    if (TAUDECAY) {
      n_interactions++; // increment number of interactions, for plotting.

      rnd1=gRandom->Rndm();
      if (rnd1<0.65011)  // pick which type of decay it is.
	whichtype="hadrdecay";
      if (rnd1>=0.65011 && rnd1<0.8219)
	whichtype="mudecay";
      if (rnd1>=0.8219)
	whichtype="edecay";
           
      rnd1=1000.;
      // rnd2=1000.;  // random numbers for throwing at dart board
      // index_y=0;     
      
      if (whichtype=="hadrdecay") { // hadronic decay
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_hadrdecay[i],NPROB,rnd1,y);	
      }
      else if (whichtype=="edecay") { // e decay	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_edecay[i],NPROB,rnd1,y);
      }
      else if (whichtype=="mudecay") { // mu decay
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_mudecay[i],NPROB,rnd1,y);
      }
      
     
      if (y*plepton>std::max(em_secondaries_max, had_secondaries_max)) {  // if this is the biggest interaction yet,    
	if (whichtype=="edecay") // save it.
	  em_secondaries_max=y*plepton;
	if (whichtype=="hadrdecay")
	  had_secondaries_max=y*plepton;
      } //if     
    } //if (TAUDECAY)
  } //if (nutau)

} //GetSecondaries



icemc::ShowerProperties icemc::Secondaries::GetEMFrac(const Settings *settings1,
						      NuFlavor nuflavor,
						      CurrentType current,
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


  ShowerProperties sp;
  sp.pnu = pnu;

  if (current==CurrentType::Charged){
    plepton=(1.-y)*pnu;
  }
  else{
    plepton=0.;
  }


  const double negligible = 1e-10;
  if (nuflavor==NuFlavor::e && current==CurrentType::Charged) {
    sp.emFrac=1.-y;
    sp.hadFrac=y;
  }
  else if(nuflavor==NuFlavor::mu && current==CurrentType::Charged) {
    sp.emFrac=negligible;
    sp.hadFrac=y;
  }
  else if(nuflavor==NuFlavor::tau && current==CurrentType::Charged) {
    // behaves like a muon
    if(taumodes1 ==1){//taumodes==1; tau created somewhere in rock and decays at posnu.
      this->GetEMFracDB(sp.emFrac,sp.hadFrac);
    }
    else if (taumodes1 == 0){
      sp.emFrac=negligible;
      sp.hadFrac=y;
    }
  }
  else if (current==CurrentType::Neutral) {
    sp.emFrac=negligible;
    sp.hadFrac=y;
  }


  em_secondaries_max = sp.emFrac; // initialize search for maximum signal among primary, secondary interactions.
  had_secondaries_max = sp.hadFrac;

  if (SECONDARIES==1 && current==CurrentType::Charged) {

    while (1) {

      // find how much em and hadronic energies comes from secondary interactions.
      // keep picking until you get a bunch of secondary interactions that conserve energy
      GetSecondaries(settings1,nuflavor,plepton,em_secondaries_max,had_secondaries_max,sp.nInteractions,hy); 

      if (em_secondaries_max+had_secondaries_max<=plepton*(1.+1.E-5)) // if conserves energy, break.
	break;
      else {
	secondary_e_noncons++; //Record how many times we come up with something that doesn't conserve energy
	em_secondaries_max=sp.emFrac;
	had_secondaries_max=sp.hadFrac;
      } //else
    } //while(1)

    if ((em_secondaries_max+had_secondaries_max)>(sp.sumFrac())*pnu) { // if maximum signal from secondaries is larger than
                                                                         // signal from primary interaction
      sp.emFrac=em_secondaries_max/pnu; // then use that one.
      sp.hadFrac=had_secondaries_max/pnu;
      if (sp.emFrac <= negligible){
	sp.emFrac=negligible;
      }
      if (sp.hadFrac<= negligible){
	sp.hadFrac=negligible;
      }
    } //if
  } //if (charged current, secondaries on)

  if (nuflavor==NuFlavor::mu && current==CurrentType::Charged && sp.nInteractions==0){
    icemcLog() << icemc::warning << "Look at this one.  inu is " << inu << "\n";
  }  

  if ((y<0 || y>1) && y != -999.) {
    icemcLog() << icemc::error <<  "Illegal value of y, y =" << y << "\n";
  }
  
  if (sp.sumFrac()>1.00001) {
    icemcLog() << icemc::error << "emFrac,hadfrac=" << sp.emFrac << "," << sp.hadFrac << ": sum = " << sp.sumFrac() << "\n";
    icemcLog() << "nuflavor,taudecay=" << nuflavor << " " << taudecay << "\n";
  }
  
  return sp;

} //GetEMFrac


//----------------------------------------------------------
//InitTauola()
//Initializes the tau decay information

void icemc::Secondaries::InitTauola() {
  for(int k=0;k<5;k++)
    tauolainfile >> tauola[0][k];
  for(int i=1;i<N_TAUOLA;i++)
    for(int j=0;j<6;j++)
      tauolainfile >> tauola[i][j];

  return;
}//InitTauola


void icemc::Secondaries::GetTauDecay(NuFlavor nuflavor, CurrentType current, std::string& taudecay, double& emfrac_db, double& hadfrac_db) {
 
  if (!(nuflavor==NuFlavor::tau || current==CurrentType::Charged || interestedintaus)){
    return;
  }
  
  // if nu_tau choose tau decay type
  
  double rnd = gRandom->Rndm();
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
    double rnd=gRandom->Rndm();
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

void icemc::Secondaries::GetEMFracDB(double& emfrac_db, double& hadfrac_db) const {


  double rnd = gRandom->Rndm();
  int decay = static_cast<int>(rnd*(N_TAUOLA-2)+1);

  hadfrac_db=tauola[decay][3];
  emfrac_db=tauola[decay][4];

  return;
}//GetEMFracDB

//------------------------------------------------------
//GetDBViewAngle()
//Gets the viewangle of the second bang

double icemc::Secondaries::GetDBViewAngle(const TVector3 &refr, const TVector3 &nnu) {

  /// @todo FIX CHANGECOORD
  // return ((nnu.ChangeCoord(refr)).Angle(z_axis));
  return nnu.Angle(TVector3(0, 0, 1));

}//GetDBViewAngle

//------------------------------------------------------
//GetFirstBang()
//Gets the position of the first bang when the interaction point is the tau decay point

//  void icemc::Secondaries::GetFirstBang(const GeoidModel::Position &r_in, const TVector3 &nnu, GeoidModel::Position &posnu, double len_int_kgm2, double chord, double &nuentrancelength) {
  
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
//     r2=gRandom->Rndm()*chord;
//   else {
//     do {
//       r1=gRandom->Rndm();
//       r2=gRandom->Rndm()*chord;
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
double icemc::Secondaries::NFBWeight(double ptau, double taulength) {
  
  double gamma=ptau/constants::MTAU;
  double D=constants::TAUDECAY_TIME*constants::CLIGHT*gamma;

  return exp(-taulength/D);

}
void icemc::Secondaries::Picky(double *y_cumulative,int NPROB,double rnd,double& y) {
  for (int i=0;i<NPROB;i++) {
    if (y_cumulative[i]<=rnd && y_cumulative[i+1]>rnd) {
      y=(double)i/(double)NPROB;
      continue; // once you found the right bin, stop looping.
    } //if
  } //for
} //Picky

