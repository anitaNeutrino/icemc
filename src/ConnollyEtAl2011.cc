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
#include "ConnollyEtAl2011.h"
#include "Report.h"
#include "RayTracer.h"
#include <cmath>

#include "TH2D.h"
#include "TCanvas.h"


#include "Inelasticity.h"


icemc::ConnollyEtAl2011::ConnollyEtAl2011(const Settings* settings)
  : fSettings(settings), fY(settings)
{//constructor

  // This is for parametrizations in Connolly et al. 2011  
  //in the form of [i][j] where i is neutrino type(nu_nubar) and j is current type, "nc" vs "cc".
  //[nu_nubar][currentint]
  //[0=nu, 1=nubar][0=neutral current, 1=charged current]
  //[0][0]->[nu][neutral current]
  //[0][1]->[nu][charged current]
  //[1][0]->[nubar][neutral current]
  //[1][1]->[nubar][charged current]
  
  //[nu][neutral current]
  auto nu_nc = std::make_pair(Neutrino::L::Matter, Neutrino::Interaction::Current::Neutral);
  c0[nu_nc] = -1.826;
  c1[nu_nc] = -17.31;
  c2[nu_nc] = -6.448; 
  c3[nu_nc] =  1.431;
  c4[nu_nc] = -18.61;
  
  //[nu][charged current]
  auto nu_cc = std::make_pair(Neutrino::L::Matter, Neutrino::Interaction::Current::Charged);
  c0[nu_cc] = -1.826;
  c1[nu_cc] = -17.31;
  c2[nu_cc] = -6.406; 
  c3[nu_cc] =  1.431;
  c4[nu_cc] = -17.91;
  
  //[nubar][neutral current]
  auto nubar_nc = std::make_pair(Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Neutral);
  c0[nubar_nc] = -1.033;
  c1[nubar_nc] = -15.95;
  c2[nubar_nc] = -7.296; 
  c3[nubar_nc] =  1.569;
  c4[nubar_nc] = -18.30;
  
  //[nubar][charged current]
  auto nubar_cc = std::make_pair(Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Charged);
  c0[nubar_cc] = -1.033;
  c1[nubar_cc] = -15.95;
  c2[nubar_cc] = -7.247;
  c3[nubar_cc] =  1.569;
  c4[nubar_cc] = -17.72;
  
  for(auto l : {Neutrino::L::Matter, Neutrino::L::AntiMatter}){ 
    for(auto cc : {Neutrino::Interaction::Current::Neutral, Neutrino::Interaction::Current::Charged}){
      std::stringstream name;
      name << "fSigma_" << l << "_" << cc;

      auto p = std::make_pair(l,cc);
      auto it_b = fSigma.emplace(p, TF1(name.str().c_str(),"pow(10, [1]+[2]*log(x-[0])+[3]*pow(log(x-[0]),2)+[4]/log(x-[0]))", 4., 21.));//check bounds. they're in log10 GeV)
      auto it = it_b.first;
      it->second.SetParameters(c0[p], c1[p], c2[p], c3[p], c4[p]);
    }
  }

  // again y distributions from Connolly et al. 2011
  // fY = std::unique_ptr<Y>(new icemc::Y());
  
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
  // mine[0] = 1.2E15;
  // mine[1] = 1.E4;
  // maxe[0] = 1.E21;
  // maxe[1] = 1.E21;

  fMinEnergy = Energy(1.E4, Energy::Unit::eV); // minimum energy for cross section parametrizations
  fMaxEnergy = Energy(1.E21, Energy::Unit::eV); // use the same upper limit for reno as for connolly et al.
  
}

// double icemc::ConnollyEtAl2011::Getyweight(double pnu,double y,int nu_nubar, Neutrino::Interaction::Current currentint) {
double icemc::ConnollyEtAl2011::Getyweight(Energy pnu,double y, Neutrino::L leptonNumber, Neutrino::Interaction::Current currentint) {  
  return fY.Getyweight(pnu,y,leptonNumber,currentint);
}


double icemc::ConnollyEtAl2011::pickY(double pnu,Neutrino::L leptonNumber,Neutrino::Interaction::Current currentint) {
  return fY.pickY(pnu,leptonNumber,currentint);
}
// double icemc::ConnollyEtAl2011::pickY(const Settings *settings1,double pnu,int nu_nubar,Neutrino::Interaction::Current currentint) {
//   return fY->pickY(settings1,pnu,nu_nubar,currentint);
// }







// int icemc::ConnollyEtAl2011::GetSigma(double pnu, double& sigma,double &len_int_kgm2, const Settings *settings1, int nu_nubar, Neutrino::Interaction::Current current){
double icemc::ConnollyEtAl2011::getSigma(Energy energy, Neutrino::L leptonNumber, Neutrino::Interaction::Current current) const {
  
  // int currentint = static_cast<int>(current);
  // int nu_nubar = leptonNumber == Neutrino::L::Matter ? 0 : 1;
  // calculate cross section
  if(!validEnergy(energy)){
    icemc::report() << severity::error << "Need a parameterization for this energy region, energy = "
		    << energy << " but min = " << fMinEnergy << ", max = " << fMaxEnergy << std::endl;
    return -1;
  }
  double sigma = 0;
  if(fSettings->SIGMAPARAM==0){ // Reno
      // fit to cross sections calculated by M.H. Reno using the same method as Gandhi et al, but with the CTEQ6-DIS parton distribution functions instead of the CTEQ4-DIS distribution functions
    sigma=(2.501E-39)*pow(energy.in(Energy::Unit::GeV),0.3076)*fSettings->SIGMA_FACTOR; // 10^18 eV - 10^21 eV(use this one for ANITA)
    //sigma=(1.2873E-39)*pow(pnu/1.E9,0.33646)*SIGMA_FACTOR; // 10^17 eV - 10^20 eV (use this one for SalSA)
  }//old code
  else if (fSettings->SIGMAPARAM==1) {//Connolly et al.
    double pnuGeV=energy.in(Energy::Unit::GeV);//Convert eV to GeV.
    double epsilon=log10(pnuGeV);

    auto f_it = fSigma.find({leptonNumber, current});
    const TF1& f = f_it->second;
    sigma=fSettings->SIGMA_FACTOR*(f.Eval(epsilon))/1.E4;//convert cm to meters. multiply by (1m^2/10^4 cm^2).
  }
  // interaction length in kg/m^2
  
  // len_int_kgm2=constants::M_NUCL/sigma; // kg/m^2
  return sigma;
} //GetSigma
