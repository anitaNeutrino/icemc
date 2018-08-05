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
#include "IcemcLog.h"
#include "RayTracer.h"
#include <cmath>

#include "TH2D.h"
#include "TCanvas.h"


#include "Inelasticity.h"





icemc::Primaries::Primaries(){//constructor

  // This is for parametrizations in Connolly et al. 2011  
  //in the form of [i][j] where i is neutrino type(nu_nubar) and j is current type, "nc" vs "cc".
  //[nu_nubar][currentint]
  //[0=nu, 1=nubar][0=neutral current, 1=charged current]
  //[0][0]->[nu][neutral current]
  //[0][1]->[nu][charged current]
  //[1][0]->[nubar][neutral current]
  //[1][1]->[nubar][charged current]
  
  //[nu][neutral current]
  c0[0][0]=-1.826;
  c1[0][0]=-17.31;
  c2[0][0]=-6.448; 
  c3[0][0]=1.431;
  c4[0][0]=-18.61;
  
  //[nu][charged current]
  c0[0][1]=-1.826;
  c1[0][1]=-17.31;
  c2[0][1]=-6.406; 
  c3[0][1]=1.431;
  c4[0][1]=-17.91;
  
  //[nubar][neutral current]	
  c0[1][0]=-1.033;
  c1[1][0]=-15.95;
  c2[1][0]= -7.296; 
  c3[1][0]=1.569;
  c4[1][0]=-18.30;
  
  //[nubar][charged current]
  c0[1][1]=-1.033;
  c1[1][1]=-15.95;
  c2[1][1]=-7.247; 
  c3[1][1]=1.569;
  c4[1][1]=-17.72;
  
  char ch[50];
  std::string stmp;
  std::string sbase="fsigma";
  for(int i=0; i<=1;i++){ // nu, nubar
    for(int j=0; j<=1; j++){ // nc, cc
      sprintf(ch,"%d%d",i,j);
      stmp=ch;	
      m_fsigma[i][j]=new TF1((sbase+stmp).c_str(),"pow(10, [1]+[2]*log(x-[0])+[3]*pow(log(x-[0]),2)+[4]/log(x-[0]))", 4., 21.);//check bounds. they're in log10 GeV.
      //x=log10(pnu/GeV).
      m_fsigma[i][j]->SetParameters(c0[i][j], c1[i][j], c2[i][j], c3[i][j], c4[i][j]);
      //"fsigma00"->[nu][neutral current]
      //"fsigma01"->[nu][charged current]
      //"fsigma10"->[nubar][neutral current]
      //"fsigma11"->[nubar][charged current]
    }		
  }
  m_csigma=new TCanvas("m_csigma","m_csigma title",1000, 700);
  m_hsigma=new TH2D("hsigma","title hsigma", 600, 7., 12., 600, -40., -30.);
  
  m_hsigma->SetTitle("log10 (pnu) vs.log10 Cross Section Sigma");
  m_hsigma->GetXaxis()->SetTitle("Log10(Ev/ GeV)");
  m_hsigma->GetYaxis()->SetTitle("log10(Cross Section/ m^2)");
  
  m_hsigma->Draw("scat");
  m_hsigma->SetMarkerStyle(7);
  m_hsigma->SetMarkerSize(3);
  
  // again y distributions from Connolly et al. 2011
  m_myY = new icemc::Y();
  
  //From Table V. Connolly Calc 2011.
  //A_low[4];//same for any [i]nu_nubar and [j]currentint.
  A_low[0]=0.0;
  A_low[1]=0.0941;
  A_low[2]=4.72;
  A_low[3]=0.456;
  //high y///////////////////
  //[0=nu, 1=nubar][0=neutral current, 1=charged current]
  //[nu_bar][currentint];
  A0_high[0][0]=-0.005;
  A1_high[0][0]=0.23;
  A2_high[0][0]=3.0;
  A3_high[0][0]=1.7;
  
  A0_high[0][1]=-0.008;
  A1_high[0][1]=0.26;
  A2_high[0][1]=3.0;
  A3_high[0][1]=1.7;
  
  A0_high[1][0]=-0.005;
  A1_high[1][0]=0.23;
  A2_high[1][0]=3.0;
  A3_high[1][0]=1.7;
  
  A0_high[1][1]=-0.0026;
  A1_high[1][1]=0.085;
  A2_high[1][1]=4.1;
  A3_high[1][1]=1.7;
  
  b0=2.55;
  b1=-0.0949; //C2_low=b0+b1*epsilon;
  
  run_old_code=0;//for GetSigma() & Gety() runs of the old code if 1, else runs current code.

  // 0=Reno
  // 1=Connolly et al. 2011
  mine[0]=1.2E15;
  mine[1]=1.E4;// minimum energy for cross section parametrizations
  maxe[0]=1.E21;
  maxe[1]=1.E21; // use the same upper limit for reno as for connolly et al.
}


double icemc::Primaries::Getyweight(double pnu,double y,int nu_nubar, Neutrino::CurrentType currentint) {
  return m_myY->Getyweight(pnu,y,nu_nubar,currentint);
}


double icemc::Primaries::pickY(const Settings *settings1,double pnu,int nu_nubar,Neutrino::CurrentType currentint) {
  return m_myY->pickY(settings1,pnu,nu_nubar,currentint);
}


icemc::Primaries::~Primaries(){//default deconstructor
  m_hsigma->Draw("same");
  m_csigma->Print("sigmaCrossSection.pdf");
  delete m_hsigma;
  delete m_myY;
  for(int i=0; i<=1;i++){ // nu, nubar
    for(int j=0; j<=1; j++){ // nc, cc
      delete m_fsigma[i][j];
    }
  }
}//deconstructor


int icemc::Primaries::GetSigma(double pnu,double& sigma,double &len_int_kgm2,const Settings *settings1,int nu_nubar,Neutrino::CurrentType current){
  int currentint = static_cast<int>(current);
  // calculate cross section
  if (pnu<mine[settings1->SIGMAPARAM] || pnu>maxe[settings1->SIGMAPARAM]) {
    icemcLog() <<  icemc::error << "Need a parameterization for this energy region.\n";
    return 0;
  }
  else {
   
    //nu=0, nubar=1
    if(nu_nubar!=0 && nu_nubar!=1){   
      std::cout<<"nu_nubar is not defined correctly!\n";
      return 0;
    }
    if (current!=Neutrino::CurrentType::Charged && current!=Neutrino::CurrentType::Neutral){//default "cc"
      std::cout<<"Current is not cc or nc!\n";
      return 0;
    }
    
    if(settings1->SIGMAPARAM==0){ // Reno
      // fit to cross sections calculated by M.H. Reno using the same method as Gandhi et al, but with the CTEQ6-DIS parton distribution functions instead of the CTEQ4-DIS distribution functions
      sigma=(2.501E-39)*pow(pnu/1.E9,0.3076)*settings1->SIGMA_FACTOR; // 10^18 eV - 10^21 eV(use this one for ANITA)
      //sigma=(1.2873E-39)*pow(pnu/1.E9,0.33646)*SIGMA_FACTOR; // 10^17 eV - 10^20 eV (use this one for SalSA)
    }//old code
    else if (settings1->SIGMAPARAM==1) {//Connolly et al.
      double pnuGeV=pnu/1.E9;//Convert eV to GeV.
      double epsilon=log10(pnuGeV);
      sigma=settings1->SIGMA_FACTOR*(m_fsigma[nu_nubar][currentint]->Eval(epsilon))/1.E4;//convert cm to meters. multiply by (1m^2/10^4 cm^2).
      
      if(m_hsigma->GetEntries()<2000){
        m_hsigma->Fill(epsilon, log10(sigma));
      }
    }//else current code
  }//if
  // interaction length in kg/m^2
  
  len_int_kgm2=constants::M_NUCL/sigma; // kg/m^2
  return 1;
} //GetSigma




//! pick a neutrino type, flavor ratio 1:1:1
icemc::Neutrino::Flavor icemc::Primaries::GetNuFlavor() const {
  Neutrino::Flavor nuflavor = Neutrino::Flavor::e;

  double rnd=gRandom->Rndm();

  if (rnd<=(1./3.)) {  
    nuflavor=Neutrino::Flavor::e;
  } //if
  else if(rnd<=(2./3.)) { 
    nuflavor=Neutrino::Flavor::mu;
  } //else if
  else if(rnd<=(1.)) { 
    nuflavor=Neutrino::Flavor::tau;
  } //else if
  else{
    std::cout << "unable to pick nu flavor\n";
  }
  return nuflavor;
} //GetNuFlavor


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

  setNuFlavor(primary1,settings1);
  setCurrent();
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


void  icemc::Interaction::setNuFlavor(const Primaries *primary1, const Settings *settings1) {
  // pick the neutrino flavor,  type of tau decay when relevant,
  //  lpm energy.
  nuflavor=primary1->GetNuFlavor();
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
      posnu = antarctica1->PickInteractionLocation(r_bn);
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

icemc::Neutrino::CurrentType icemc::Interaction::GetCurrent() {
  Neutrino::CurrentType current;
  double rnd=gRandom->Rndm();
  if (rnd<=0.6865254){ // 10^18 eV - 10^21 eV (use this one for ANITA)
//if (rnd<=0.6893498) // 10^17 eV - 10^20 eV (use this one for SalSA)
    current = Neutrino::CurrentType::Charged;//"cc";
  }
  else{
    current = Neutrino::CurrentType::Neutral;//"nc";  
  }
  return current;
} //GetCurrent

void icemc::Interaction::setCurrent() {
  // pick whether it is neutral current
  // or charged current
  current=this->GetCurrent();
}//setCurrent


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


