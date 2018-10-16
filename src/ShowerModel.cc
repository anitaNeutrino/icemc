#include "Settings.h"
#include "Geoid.h"
#include "ConnollyEtAl2011.h"
#include "ShowerModel.h"
#include "Antarctica.h"
#include "Tools.h"
#include "EnvironmentVariable.h"
#include "Report.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

ClassImp(icemc::Shower)

icemc::ShowerModel::~ShowerModel(){
}


icemc::ShowerModel::ShowerModel(const Settings* settings) : fSettings(settings) {
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
  

  SECONDARIES=fSettings->SECONDARIES;
  TAUDECAY=fSettings->TAUDECAY;

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

	
	
}//ShowerModel Constructor

// void icemc::ShowerModel::readData(std::string nuflavor, std::string secndryType, double (*y)[nP], double (*dsdy)[nP])
// void icemc::ShowerModel::readData(const std::string& nuflavor,
// 				      const std::string& secndryType,
// 				      std::array<std::array<double, nP>, nE>& y,
// 				      std::array<std::array<double, nP>, nE>& dsdy)
// 				      // double (*dsdy)[nP])  
// {

void icemc::ShowerModel::readData(Neutrino::Flavor flavor, Secondary secondary){
  
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

void icemc::ShowerModel::Draw(Option_t* opt){

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

void icemc::ShowerModel::ReadSecondaries() {
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

  // }
  std::cout<<"Finished reading secondary interaction data.\n"; 
} //end method ReadSecondaries


void icemc::ShowerModel::doShower(Neutrino::Flavor nuflavor, Energy plepton, Energy &em_secondaries_max, Energy &had_secondaries_max,int &n_interactions) {
  // @todo currently disabled as it's not totally necessary and I don't understand what it used to do because it was written so fucking badly
  return;

  em_secondaries_max.setZero();
  had_secondaries_max.setZero();
  
  Energy energy = plepton;
  if (nuflavor==Neutrino::Flavor::mu || nuflavor==Neutrino::Flavor::tau){
    /**
     * Pick poisson distribution for each, uses i.
     */
    std::vector<Secondary>secondaries;
    
    for(auto sec : {Secondary::brems, Secondary::epair, Secondary::pn}){
      DataKey key(nuflavor, sec);
      TH1F& h_E_dsdy = fE_dsdy[key];
      
      int nPoisson = pickPoisson(h_E_dsdy.Interpolate(energy.in(Energy::Unit::eV)));
      for(int k=0; k < nPoisson; k++){
      	secondaries.emplace_back(sec);
      }
    }
    std::random_shuffle(secondaries.begin(), secondaries.end(),
			[&](int z){ return int(pickUniform(z))%z; } );

    for(auto secondary : secondaries){
      // double rnd1 = pickUniform(secondaries.size());

      DataKey key(nuflavor, secondary);
      fE_YCumulative_dsdy[key];
    }
  }    
} //GetShowerModel



icemc::Shower icemc::ShowerModel::generate(const Neutrino& nu) {
  Shower s = GetEMFrac(nu.flavor, nu.interaction.current, nu.interaction.y,  nu.energy);
  s.axis = nu.path.direction.Unit();
  // doShower(nu.flavor,nu.energy,em_secondaries_max,had_secondaries_max,s.nInteractions);   
  return s;
}


icemc::Shower icemc::ShowerModel::GetEMFrac(Neutrino::Flavor nuflavor,
					    Neutrino::Interaction::Current current,
					    double y,
					    Energy pnu) {
  // int taumodes1) {


  Shower s;
  s.pnu = pnu;

  // if (current==Neutrino::Interaction::Current::Charged){
  //   plepton=(1.-y)*pnu;
  // }
  // else{
  //   plepton.setZero();
  // }


  const double negligible = 1e-10;
  if (nuflavor==Neutrino::Flavor::e && current==Neutrino::Interaction::Current::Charged) {
    s.emFrac = 1.-y;
    s.hadFrac = y;
  }
  else if(nuflavor==Neutrino::Flavor::mu && current==Neutrino::Interaction::Current::Charged) {
    s.emFrac = negligible;
    s.hadFrac = y;
  }
  else if(nuflavor==Neutrino::Flavor::tau && current==Neutrino::Interaction::Current::Charged) {
    // behaves like a muon
    // if(taumodes1 ==1){//taumodes==1; tau created somewhere in rock and decays at posnu.
    //   this->pickEMFracDB(s.emFrac,s.hadFrac);
    // }
    // else if (taumodes1 == 0){
    s.emFrac=negligible;
    s.hadFrac=y;
    // }
  }
  else if (current==Neutrino::Interaction::Current::Neutral) {
    s.emFrac = negligible;
    s.hadFrac = y;
  }

  // em_secondaries_max.setZero(); // initialize search for maximum signal among primary, secondary interactions.
  // had_secondaries_max.setZero();

  ///@todo restore a better version of this
  // if (SECONDARIES==1 && current==Neutrino::Interaction::Current::Charged) {

  //   while (1) {

  //     // find how much em and hadronic energies comes from secondary interactions.
  //     // keep picking until you get a bunch of secondary interactions that conserve energy
  //     doShower(nuflavor,plepton,em_secondaries_max,had_secondaries_max,s.nInteractions); 

  //     if (em_secondaries_max+had_secondaries_max<=plepton*(1.+1.E-5)) // if conserves energy, break.
  // 	break;
  //     else {
  // 	secondary_e_noncons++; //Record how many times we come up with something that doesn't conserve energy
  // 	em_secondaries_max.setZero();
  // 	had_secondaries_max.setZero();
  //     } //else
  //   } //while(1)

  //   if ((em_secondaries_max+had_secondaries_max)>(s.sumFrac())*pnu) { // if maximum signal from secondaries is larger than
  //                                                                        // signal from primary interaction
  //     s.emFrac=em_secondaries_max/pnu; // then use that one.
  //     s.hadFrac=had_secondaries_max/pnu;
  //     if (s.emFrac <= negligible){
  // 	s.emFrac=negligible;
  //     }
  //     if (s.hadFrac<= negligible){
  // 	s.hadFrac=negligible;
  //     }
  //   } //if
  // } //if (charged current, secondaries on)

  // if (nuflavor==Neutrino::Flavor::mu && current==Neutrino::Interaction::Current::Charged && s.nInteractions==0){
  //   icemc::report() << severity::warning << "Look at this one.  inu is " << inu << "\n";
  // }  

  if (s.sumFrac()>1.00001) {
    icemc::report() << severity::error << "emFrac,hadfrac=" << s.emFrac << "," << s.hadFrac << ": sum = " << s.sumFrac() << "\n";
  }
  
  return s;
} //GetEMFrac


//----------------------------------------------------------
//InitTauola()
//Initializes the tau decay information

void icemc::ShowerModel::InitTauola() {
  for(int k=0;k<5;k++)
    tauolainfile >> tauola[0][k];
  for(int i=1;i<N_TAUOLA;i++)
    for(int j=0;j<6;j++)
      tauolainfile >> tauola[i][j];

  return;
}//InitTauola


void icemc::ShowerModel::GetTauDecay(Neutrino::Flavor nuflavor, Neutrino::Interaction::Current current, std::string& taudecay, double& emfrac_db, double& hadfrac_db) {
 
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

void icemc::ShowerModel::pickEMFracDB(double& emfrac_db, double& hadfrac_db) {


  double rnd = pickUniform();
  int decay = static_cast<int>(rnd*(N_TAUOLA-2)+1);

  hadfrac_db=tauola[decay][3];
  emfrac_db=tauola[decay][4];

  return;
}//GetEMFracDB

//------------------------------------------------------
//GetDBViewAngle()
//Gets the viewangle of the second bang

double icemc::ShowerModel::GetDBViewAngle(const TVector3 &refr, const TVector3 &nnu) {

  /// @todo FIX CHANGECOORD
  // return ((nnu.ChangeCoord(refr)).Angle(z_axis));
  return nnu.Angle(TVector3(0, 0, 1));

}//GetDBViewAngle

//------------------------------------------------------
//GetFirstBang()
//Gets the position of the first bang when the interaction point is the tau decay point

//  void icemc::ShowerModel::GetFirstBang(const Geoid::Position &r_in, const TVector3 &nnu, Geoid::Position &posnu, double len_int_kgm2, double chord, double &nuentrancelength) {
  
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
double icemc::ShowerModel::NFBWeight(double ptau, double taulength) {
  
  double gamma=ptau/constants::MTAU;
  double D=constants::TAUDECAY_TIME*constants::CLIGHT*gamma;

  return exp(-taulength/D);

}

void icemc::ShowerModel::Picky(const double *y_cumulative,int NPROB,double rnd,double& y) const {
  for (int i=0;i<NPROB;i++) {
    if (y_cumulative[i] <= rnd && y_cumulative[i+1] > rnd) {
      y=(double)i/(double)NPROB;
      continue; // once you found the right bin, stop looping.
    } //if
  } //for
} //Picky

