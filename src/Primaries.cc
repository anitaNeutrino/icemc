#include "TRandom3.h"
#include "Constants.h"
#include "TVector3.h"
#include "Geoid.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include "Settings.h"
#include "Crust2.h"
#include "Antarctica.h"
#include "Primaries.h"
#include "Report.h"
#include "RayTracer.h"
#include <cmath>

#include "TH2D.h"
#include "TCanvas.h"


#include "Inelasticity.h"




icemc::Primaries::Primaries(const Settings* settings) : fSettings(settings){//constructor

  // This is for parametrizations in Connolly et al. 2011  
  //in the form of [i][j] where i is neutrino type(nu_nubar) and j is current type, "nc" vs "cc".
  //[nu_nubar][currentint]
  //[0=nu, 1=nubar][0=neutral current, 1=charged current]
  //[0][0]->[nu][neutral current]
  //[0][1]->[nu][charged current]
  //[1][0]->[nubar][neutral current]
  //[1][1]->[nubar][charged current]
  
  //[nu][neutral current]
  auto nu_nc = std::make_pair(Neutrino::L::Matter, Neutrino::Current::Neutral);
  c0[nu_nc] = -1.826;
  c1[nu_nc] = -17.31;
  c2[nu_nc] = -6.448; 
  c3[nu_nc] =  1.431;
  c4[nu_nc] = -18.61;
  
  //[nu][charged current]
  auto nu_cc = std::make_pair(Neutrino::L::Matter, Neutrino::Current::Charged);
  c0[nu_cc] = -1.826;
  c1[nu_cc] = -17.31;
  c2[nu_cc] = -6.406; 
  c3[nu_cc] =  1.431;
  c4[nu_cc] = -17.91;
  
  //[nubar][neutral current]
  auto nubar_nc = std::make_pair(Neutrino::L::AntiMatter, Neutrino::Current::Neutral);
  c0[nubar_nc] = -1.033;
  c1[nubar_nc] = -15.95;
  c2[nubar_nc] = -7.296; 
  c3[nubar_nc] =  1.569;
  c4[nubar_nc] = -18.30;
  
  //[nubar][charged current]
  auto nubar_cc = std::make_pair(Neutrino::L::AntiMatter, Neutrino::Current::Charged);
  c0[nubar_cc] = -1.033;
  c1[nubar_cc] = -15.95;
  c2[nubar_cc] = -7.247;
  c3[nubar_cc] =  1.569;
  c4[nubar_cc] = -17.72;
  
  for(auto l : {Neutrino::L::Matter, Neutrino::L::AntiMatter}){ 
    for(auto cc : {Neutrino::Current::Neutral, Neutrino::Current::Charged}){
      std::stringstream name;
      name << "fSigma_" << l << "_" << cc;

      auto p = std::make_pair(l,cc);
      auto it_b = fSigma.emplace(p, TF1(name.str().c_str(),"pow(10, [1]+[2]*log(x-[0])+[3]*pow(log(x-[0]),2)+[4]/log(x-[0]))", 4., 21.));//check bounds. they're in log10 GeV)
      auto it = it_b.first;
      it->second.SetParameters(c0[p], c1[p], c2[p], c3[p], c4[p]);
    }
  }
  // m_csigma=new TCanvas("m_csigma","m_csigma title",1000, 700);
  // m_hsigma=new TH2D("hsigma","title hsigma", 600, 7., 12., 600, -40., -30.);
  
  // m_hsigma->SetTitle("log10 (pnu) vs.log10 Cross Section Sigma");
  // m_hsigma->GetXaxis()->SetTitle("Log10(Ev/ GeV)");
  // m_hsigma->GetYaxis()->SetTitle("log10(Cross Section/ m^2)");
  
  // m_hsigma->Draw("scat");
  // m_hsigma->SetMarkerStyle(7);
  // m_hsigma->SetMarkerSize(3);
  
  // again y distributions from Connolly et al. 2011
  // m_myY = std::unique_ptr<Y>(new icemc::Y());
  
  //From Table V. Connolly Calc 2011.
  //A_low[4];//same for any [i]nu_nubar and [j]currentint.
  A_low[0]=0.0;
  A_low[1]=0.0941;
  A_low[2]=4.72;
  A_low[3]=0.456;
  
  //high y///////////////////
  //[0=nu, 1=nubar][0=neutral current, 1=charged current]
  //[nu_bar][currentint];
  
  
  A0_high[nu_nc] = -0.005;
  A1_high[nu_nc] = 0.23;
  A2_high[nu_nc] = 3.0;
  A3_high[nu_nc] = 1.7;
  
  A0_high[nu_cc] = -0.008;
  A1_high[nu_cc] = 0.26;
  A2_high[nu_cc] = 3.0;
  A3_high[nu_cc] = 1.7;
  
  A0_high[nubar_nc] = -0.005;
  A1_high[nubar_nc] = 0.23;
  A2_high[nubar_nc] = 3.0;
  A3_high[nubar_nc] = 1.7;
  
  A0_high[nubar_cc] = -0.0026;
  A1_high[nubar_cc] = 0.085;
  A2_high[nubar_cc] = 4.1;
  A3_high[nubar_cc] = 1.7;
  
  b0=2.55;
  b1=-0.0949; //C2_low=b0+b1*epsilon;
  
  // 0=Reno
  // 1=Connolly et al. 2011
  mine[0] = 1.2E15;
  mine[1] = 1.E4;// minimum energy for cross section parametrizations
  maxe[0] = 1.E21;
  maxe[1] = 1.E21; // use the same upper limit for reno as for connolly et al.
}

// double icemc::Primaries::Getyweight(double pnu,double y,int nu_nubar, Neutrino::Current currentint) {
double icemc::Primaries::Getyweight(double pnu,double y, Neutrino::L leptonNumber, Neutrino::Current currentint) {  
  return m_myY.Getyweight(pnu,y,leptonNumber,currentint);
}


double icemc::Primaries::pickY(double pnu,Neutrino::L leptonNumber,Neutrino::Current currentint) {
  return m_myY.pickY(fSettings,pnu,leptonNumber,currentint);
}
// double icemc::Primaries::pickY(const Settings *settings1,double pnu,int nu_nubar,Neutrino::Current currentint) {
//   return m_myY->pickY(settings1,pnu,nu_nubar,currentint);
// }


icemc::Primaries::~Primaries(){//default deconstructor
}//deconstructor



TCanvas* icemc::Primaries::plotSigma(Neutrino::L l, Neutrino::Current c, int nSamples) {
  TH2D* h = new TH2D("hsigma","title hsigma", 600, 7., 12., 600, -40., -30.);

  for(int i=0; i < nSamples; i++){

    double pnu = pickUniform(mine[fSettings->SIGMAPARAM], maxe[fSettings->SIGMAPARAM]);
    double pnuGeV=pnu/1.E9;//Convert eV to GeV.
    double epsilon=log10(pnuGeV);

    double sigma, len;
    GetSigma(pnu, sigma, len, l, c);

    h->Fill(epsilon, log10(sigma));
  }
  h->SetBit(kCanDelete);

  auto p = std::make_pair(l, c);

  auto can = new TCanvas();
  h->Draw("colz");
  
  fSigma[p].Draw("lsame");
  
  return can;
}




// int icemc::Primaries::GetSigma(double pnu, double& sigma,double &len_int_kgm2, const Settings *settings1, int nu_nubar, Neutrino::Current current){
int icemc::Primaries::GetSigma(double pnu, double& sigma,double &len_int_kgm2, Neutrino::L leptonNumber, Neutrino::Current current){
  
  // int currentint = static_cast<int>(current);
  // int nu_nubar = leptonNumber == Neutrino::L::Matter ? 0 : 1;
  // calculate cross section
  if (pnu<mine[fSettings->SIGMAPARAM] || pnu>maxe[fSettings->SIGMAPARAM]) {
    icemc::report() <<  severity::error << "Need a parameterization for this energy region.\n";
    return 0;
  }
  else {
    if(fSettings->SIGMAPARAM==0){ // Reno
      // fit to cross sections calculated by M.H. Reno using the same method as Gandhi et al, but with the CTEQ6-DIS parton distribution functions instead of the CTEQ4-DIS distribution functions
      sigma=(2.501E-39)*pow(pnu/1.E9,0.3076)*fSettings->SIGMA_FACTOR; // 10^18 eV - 10^21 eV(use this one for ANITA)
      //sigma=(1.2873E-39)*pow(pnu/1.E9,0.33646)*SIGMA_FACTOR; // 10^17 eV - 10^20 eV (use this one for SalSA)
    }//old code
    else if (fSettings->SIGMAPARAM==1) {//Connolly et al.
      double pnuGeV=pnu/1.E9;//Convert eV to GeV.
      double epsilon=log10(pnuGeV);
      sigma=fSettings->SIGMA_FACTOR*(fSigma[{leptonNumber, current}].Eval(epsilon))/1.E4;//convert cm to meters. multiply by (1m^2/10^4 cm^2).
      
      // if(m_hsigma->GetEntries()<2000){
      //   m_hsigma->Fill(epsilon, log10(sigma));
      // }
    }
  }
  // interaction length in kg/m^2
  
  len_int_kgm2=constants::M_NUCL/sigma; // kg/m^2
  return 1;
} //GetSigma




//! pick a neutrino type, flavor ratio 1:1:1
icemc::Neutrino::Flavor icemc::Primaries::pickFlavor() {
  double r = pickUniform(0, 3);
  if (r <= 1){  
    return Neutrino::Flavor::e;
  }
  else if(r <= 2){
    return Neutrino::Flavor::mu;
  }
  else if(r <= 3) { 
    return Neutrino::Flavor::tau;
  }
  else{
    auto t = Neutrino::Flavor::tau;
    report() << severity::error << "Random number too large to pick neutrino flavor, returning " << t << std::endl;
    return t;
  }
}


















icemc::Interaction::Interaction(Primaries *primary1, const Settings *settings1) {

  noway=0;
  wheredoesitleave_err=0;
  neverseesice=0;

  wheredoesitenterice_err=0;
  toohigh=0;
  toolow=0;

  iceinteraction=0;
  dtryingdirection=0.;
  dnutries=0.;

  weight_nu=0;
  weight_nu_prob=0;

  // setNuFlavor(primary1,settings1);
  current = pickCurrent(); //setCurrent();
}


void icemc::Interaction::PickAnyDirection() {
  double rndlist[2];
  gRandom->RndmArray(2,rndlist);
  
  costheta_nutraject=2*rndlist[0]-1;

  // pick a neutrino azimuthal angle
  phi_nutraject=2*constants::PI*rndlist[1];
  
  // check that these give the right result
  double thetanu=acos(costheta_nutraject);
  
  double sinthetanu=sin(thetanu);
  
  // find direction vector of neutrino
  // **** are cosine and sine flipped?
  nnu.SetX(sinthetanu*cos(phi_nutraject));
  nnu.SetY(sinthetanu*sin(phi_nutraject));
  nnu.SetZ(costheta_nutraject);
}




int icemc::Interaction::PickDownwardInteractionPoint(const Geoid::Position&r_bn, const Settings *settings1, const Antarctica *antarctica1) {

  if(settings1->UNBIASED_SELECTION==1){
    if(antarctica1->PickUnbiased(this)){ // pick neutrino direction and interaction point
      dtryingdirection=1.;
      iceinteraction=1;
    }
    else{
      iceinteraction=0;
    }
  }
  else{
    iceinteraction=1;

    // If we require neutrinos from a particular position
    // we generate that cartesian position here

    static Geoid::Position specific_position;

    if (settings1->SPECIFIC_NU_POSITION) 
      {
	specific_position.SetLonLatAlt(settings1->SPECIFIC_NU_POSITION_LONGITUDE,
				       settings1->SPECIFIC_NU_POSITION_LATITUDE,
				       settings1->SPECIFIC_NU_POSITION_ALTITUDE);
      }

    do{
      ///@todo ibnposition
      posnu = antarctica1->pickInteractionPosition(r_bn);
    } while(settings1->SPECIFIC_NU_POSITION &&  (posnu - specific_position).Mag() > settings1->SPECIFIC_NU_POSITION_DISTANCE);    
  }

  // // first pass at rf exit point,  straight above the interaction point!
  // // ray1->rfexit[0] = antarctica1->Surface(interaction1->posnu) * interaction1->posnu.Unit(); 
  // ray1->initGuess(posnu, r_bn);

  // double r_down = 2*(antarctica1->Surface(posnu)-antarctica1->IceThickness(posnu))-posnu.Mag();
  // posnu_down = r_down * posnu.Unit();
  // //position of the mirror point of interaction
  
  // //interaction1->posnu is downward interaction1->posnu.
  // // distance=interaction1->posnu.Distance(r_bn);

  // // depth of interaction
  // // gets distance between interaction and exit point, this time it's same as depth
  // // because our first guess at exit point is radially outward from interaction.
  // // negative means below surface
  // altitude_int=-1*ray1->rfexit[0].Distance(posnu);
  // altitude_int_mirror=-1*ray1->rfexit[0].Distance(posnu_down);//get depth of mirror point

  // r_fromballoon[0]=r_bn.Distance(posnu);

  // //distance from the mirror point to the balloon because it is equal to the path that signals pass
  // r_fromballoon[1]=r_bn.Distance(posnu_down);

  // // if the angle from the initial exit point guess to the balloon is more than 90 degree
  // // then it's over the horizon, so return 1
  // if (ray1->n_exit2bn[0].Angle(posnu) > constants::PI/2){
  //   return 1;
  // }
  
  return 0;
}//PickDownwardInteractionPoint



/**
 * Choose CC or NC: get from ratios in Ghandi etal paper,
 * updated for the CTEQ6-DIS parton distribution functions (M.H. Reno, personal communication).
 * Need to add capability of using ratios from Connolly et al.
 */

icemc::Neutrino::Current icemc::Interaction::pickCurrent() {
  double rnd = pickUniform();
  if (rnd<=0.6865254){ // 10^18 eV - 10^21 eV (use this one for ANITA)
//if (rnd<=0.6893498) // 10^17 eV - 10^20 eV (use this one for SalSA)
    return Neutrino::Current::Charged;//"cc";
  }
  else{
    return Neutrino::Current::Neutral;//"nc";  
  }
} //GetCurrent



int icemc::Interaction::getPdgCode() const {
  int pdgcode = -1;
  if (nuflavor==Neutrino::Flavor::e){
    pdgcode = 12;
  }
  else if (nuflavor==Neutrino::Flavor::mu){
    pdgcode = 14;
  }
  else if (nuflavor==Neutrino::Flavor::tau){
    pdgcode = 16;
  }
  return pdgcode;
}


