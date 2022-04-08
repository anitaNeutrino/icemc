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


    //@todo tried switching y (on Y-axis) and probability (content of bin) so that we could interpolate between X=Energy and Y=Probability to get a y-value, but didn't get it to work
    // double total=0;
    // for(int i=0; i < gr.GetN(); i++){
    //   int bin = h.Fill(energy, gr.GetX()[i], gr.GetY()[i]);
    //   h.SetBinError(bin, 0);

    //   total += gr.GetY()[i];      
    // }
    // int bin2 = hInt.Fill(energy, total);
    // hInt.SetBinError(bin2, 0);

    // double sum=0;
        // for(int i=0; i < gr.GetN(); i++){
    //   sum += gr.GetY()[i];      

    //   int bin = hC.GetBin(energy, sum/total);
    //   hC.Fill(bin, gr.GetX()[i]);
    //   hC.SetBinError(bin, 0);
    // }
    
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

    cans.back()->SaveAs(".pdf");
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


// Finds a shower that conserved energy
icemc::Shower icemc::ShowerModel::doShower(const Shower& s, const Neutrino::Flavor nuflavor, const double y) {
  //return;

  Shower shower = s;
  
  double em_secondaries_max = 0;
  double had_secondaries_max = 0;

  Energy plepton = (1-y)*shower.pnu;

  // Create a secondary shower for this interaction
  do {
    em_secondaries_max = s.emFrac;
    had_secondaries_max = s.hadFrac;
  
    /**
     * Pick poisson distribution for each, uses i.
     */
    std::vector<Secondary>secondaries;
    
    for(auto sec : {Secondary::brems, Secondary::epair, Secondary::pn}){
      DataKey key(nuflavor, sec);
      TH1F& h_E_dsdy = fE_dsdy[key];

      int nPoisson = pickPoisson(h_E_dsdy.Interpolate(log10(plepton.in(Energy::Unit::eV))));
      for(int k=0; k < nPoisson; k++){
	secondaries.emplace_back(sec);
      }
    }
    
    std::random_shuffle(secondaries.begin(), secondaries.end(),
			[&](int z){ return int(pickUniform(z))%z; } );

    if (nuflavor==Neutrino::Flavor::tau && TAUDECAY){
      // If we're interested in tau decay, add one to the end of the secondary interactions
      for (int n=0; n<10; n++)
	secondaries.emplace_back(pickTauDecayType());
    }

    shower.nInteractions += secondaries.size();
    // Go over each of these secondary particles and check if it creates a bigger interaction than anything before
    for(auto secondary : secondaries){
      
      DataKey key(nuflavor, secondary);
      TH2F& h_E_YCumulative_dsdy = fE_YCumulative_dsdy[key];
      
      double rn = pickUniform();
      //double newY = h_E_YCumulative_dsdy.Interpolate(log10(plepton.in(Energy::Unit::eV)), rn);
      int index = 2*log10(plepton.in(Energy::Unit::EeV))+1;
      double newY = 0;
      //@todo better way to seek through this TH2?
      if(h_E_YCumulative_dsdy.GetBinContent(index, 1)<=rn){ // If first bin is already larger than random probability, y is 0
	for(int by=0;  by <= h_E_YCumulative_dsdy.GetNbinsY(); by++){
	  if (h_E_YCumulative_dsdy.GetBinContent(index, by+1)<=rn && h_E_YCumulative_dsdy.GetBinContent(index, by+2)>rn){
	    newY = (double)by/(double)h_E_YCumulative_dsdy.GetNbinsY();
	    break;
	  }
	}
      }
      
      // If this secondary particle has the largest interaction so far
      if (newY*plepton.in(Energy::Unit::eV) > std::max(em_secondaries_max, had_secondaries_max)){
	if (secondary==Secondary::brems || secondary==Secondary::epair){
	  em_secondaries_max = newY*plepton.in(Energy::Unit::eV);
	}
	if (secondary==Secondary::pn){
	  had_secondaries_max = newY*plepton.in(Energy::Unit::eV);
	}
      }
    }
  } while(em_secondaries_max+had_secondaries_max > plepton.in(Energy::Unit::eV)*(1.+1.E-5)); // Until we get one that conserves energy


  
  // if maximum signal from secondaries is larger than signal from primary interaction
  if ((em_secondaries_max+had_secondaries_max)>(shower.sumFrac())*shower.pnu.in(Energy::Unit::eV)) {
    shower.secondary = true;

    const double negligible = 1e-10; ///@todo add to constants?
    shower.emFrac=em_secondaries_max/shower.pnu.in(Energy::Unit::eV);
    if (shower.emFrac < negligible){
      shower.emFrac=negligible;
    }
    
    shower.hadFrac=had_secondaries_max/shower.pnu.in(Energy::Unit::eV);
    if (shower.hadFrac < negligible){
      shower.hadFrac=negligible;
    }
  }

  return shower;
}



icemc::Shower icemc::ShowerModel::generate(const Neutrino& nu, const Interaction& i) {
  Shower s = GetEMFrac(nu.flavor, i.current, i.y,  nu.energy);
  if (s.secondary)
    std::cout << "Generated shower using secondary signal\n";
  s.axis = nu.path.direction.Unit();  
  return s;
}


icemc::Shower icemc::ShowerModel::GetEMFrac(Neutrino::Flavor nuflavor,
					    Interaction::Current current,
					    double y,
					    Energy pnu) {
  Shower s;
  s.pnu = pnu;

  const double negligible = 1e-10; ///@todo make a constant?
  if (nuflavor==Neutrino::Flavor::e && current==Interaction::Current::Charged) {
    s.emFrac = 1 - y; // if it turns into an electron, it will get stopped quickly and all the energy ends up in the shower
    s.hadFrac = y;
  }
  else if(nuflavor==Neutrino::Flavor::mu && current==Interaction::Current::Charged) {
    s.emFrac = negligible; // if it turns into a muon, it lives for a long time, and goes much further longer than the shower size?
    s.hadFrac = y;
  }
  else if(nuflavor==Neutrino::Flavor::tau && current==Interaction::Current::Charged) {
    ///@todo model taus better
    // behaves like a muon

    // if(taumodes1 ==1){//taumodes==1; tau created somewhere in rock and decays at posnu.
    //   this->pickEMFracDB(s.emFrac,s.hadFrac);
    // }
    // else if (taumodes1 == 0){
    s.emFrac=negligible; //?
    s.hadFrac=y;
    // }
  }
  else if (current==Interaction::Current::Neutral) {
    s.emFrac = negligible;
    s.hadFrac = y;
  }

  //@Askaryan parameterization
  //s.hadFrac = y;
  //s.emFrac = negligible;
  
  //@todo Restructure this -- needs to actually change values stored in shower
  if (SECONDARIES==1 && current==Interaction::Current::Charged && (nuflavor==Neutrino::Flavor::mu || nuflavor==Neutrino::Flavor::tau)) {
    
    // find how much em and hadronic energies comes from secondary interactions.
    s = doShower(s, nuflavor, y); 
    
  } //if (charged current, secondaries on)

  if (nuflavor==Neutrino::Flavor::mu && current==Interaction::Current::Charged && s.nInteractions==0){
    icemc::report() << severity::warning << "Look at this one, mu neutrino with CC interaction and zero interactions\n";
  }  
  
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


void icemc::ShowerModel::GetTauDecay(Neutrino::Flavor nuflavor, Interaction::Current current, std::string& taudecay, double& emfrac_db, double& hadfrac_db) {
 
  if (!(nuflavor==Neutrino::Flavor::tau || current==Interaction::Current::Charged || interestedintaus)){
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

// void icemc::ShowerModel::Picky(const double *y_cumulative,int NPROB,double rnd,double& y) const {
//   for (int i=0;i<NPROB;i++) {
//     if (y_cumulative[i] <= rnd && y_cumulative[i+1] > rnd) {
//       y=(double)i/(double)NPROB;
//       continue; // once you found the right bin, stop looping.
//     } //if
//   } //for
// } //Picky



icemc::ShowerModel::Secondary icemc::ShowerModel::pickTauDecayType(){
  //@todo verify these, maybe there are better values now
  double rn = pickUniform();
  if (rn<0.65011){
    return Secondary::hadrdecay;
  }
  else if (rn>0.8219){
    return Secondary::edecay;
  }
  else{
    return Secondary::mudecay;
  }
}
