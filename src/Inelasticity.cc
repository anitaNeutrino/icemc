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

///////////////// Y //////////////
icemc::Y::Y() { // Constructor

  /**
   * The Y class contains all of the parameterizations for generating
   * inelasticity distributions
   * We are following Connolly et al. (2011) but any model can be added.
   */
  
  ffrac = std::unique_ptr<TF1>(new TF1("ffrac","[0]*sin([1]*(x-[2]))",7.,12.)); // This is the fraction of the distribution in the low y region given by Equation 18. 
  
  ffrac->FixParameter(0,0.128); // These parameters are the same for all interaction types
  ffrac->FixParameter(1,-0.197);
  ffrac->FixParameter(2,21.8);

  std::string sbase="C1_high";
  std::string which;
  for(auto l : {Neutrino::L::Matter, Neutrino::L::AntiMatter}){
    for(auto c : {Neutrino::Current::Charged, Neutrino::Current::Neutral}){
      std::stringstream name;
      name  << "fC1_high_" << l << "_" <<  c;
      auto lc = std::make_pair(l, c);
      fC1_high.emplace(lc, TF1(name.str().c_str(),"[0]-[1]*(-exp(-(x-[2])/[3]))",7.,12.)); // parameterization of parameter C1 in the high y region according to Equation 16
    }
  }

  // parameter A_0 in Table V for the high y region
  fC1_high[{Neutrino::L::AntiMatter, Neutrino::Current::Charged}].FixParameter(0,-0.0026);//nubar, CC
  fC1_high[{Neutrino::L::Matter,     Neutrino::Current::Charged}].FixParameter(0,-0.008); //nu,    CC
  fC1_high[{Neutrino::L::AntiMatter, Neutrino::Current::Neutral}].FixParameter(0,-0.005); //nubar, NC
  fC1_high[{Neutrino::L::Matter,     Neutrino::Current::Neutral}].FixParameter(0,-0.005); //nu,    NC

  // parameter A_1 in Table V for the high y region
  fC1_high[{Neutrino::L::AntiMatter, Neutrino::Current::Charged}].FixParameter(1,0.085); // nubar, CC
  fC1_high[{Neutrino::L::Matter,     Neutrino::Current::Charged}].FixParameter(1,0.26); // nu, CC
  fC1_high[{Neutrino::L::AntiMatter, Neutrino::Current::Neutral}].FixParameter(1,0.23); // nubar, NC
  fC1_high[{Neutrino::L::Matter,     Neutrino::Current::Neutral}].FixParameter(1,0.23); // nu, NC

  // parameter A_2 in Table V for the high y region
  fC1_high[{Neutrino::L::AntiMatter, Neutrino::Current::Charged}].FixParameter(2,4.1); // nubar, CC
  fC1_high[{Neutrino::L::Matter,     Neutrino::Current::Charged}].FixParameter(2,3.0); // nu, CC   
  fC1_high[{Neutrino::L::AntiMatter, Neutrino::Current::Neutral}].FixParameter(2,3.0); // nubar, NC
  fC1_high[{Neutrino::L::Matter,     Neutrino::Current::Neutral}].FixParameter(2,3.0); // nu, NC

  // parameter A_3 in Table V for the high y region.  This parameter is the same for all four interaction types
  for(auto& f : fC1_high){
    f.second.FixParameter(3, 1.7);
  }
  
  // for (int i=0;i<2;i++) { // nu, nubar
  //   for (int j=0;j<2;j++) { // CC, NC
  //     fC1_high[i][j]->FixParameter(3,1.7);
  //   }
  // }

  fC1_low = std::unique_ptr<TF1>(new TF1("C1_low","[0]-[1]*(-exp(-(x-[2])/[3]))",7.,12.)); // parameterization of parameter C1 in the low y region according to Equation 16.
  // This parameterization is the same for all interaction types.
  
  fC1_low->FixParameter(0,0.);
  fC1_low->FixParameter(1,0.0941);
  fC1_low->FixParameter(2,4.72);
  fC1_low->FixParameter(3,0.456);
  
  fC2 = std::unique_ptr<TF1>(new TF1("C2","[0]+[1]*x",7.,12.)); // parameterization of parameter C2 in the low y region according to Equation 17.
  // This parameterization is the same for all interaction types.
  fC2->FixParameter(0,2.55);
  fC2->FixParameter(1,-9.49E-2);

  // For picking inelasticity in low y region according to Equation 14.
  fy0_low = std::unique_ptr<TF3>(new TF3("fy0_low","x+(z*([1]-x)^(-1./y+1)+(1-z)*([0]-x)^(-1./y+1))^(y/(y-1))")); // x=C_1, y=C_2, z=R
  fy0_low->SetParameter(0,ymin_low);  // y_min
  fy0_low->SetParameter(1,ymax_low); // y_max

  // For picking inelasticity in high y region according to Equation 15.
  fy0_high = std::unique_ptr<TF2>(new TF2("fy0_high","([1]-x)^y/([0]-x)^(y-1.)+x")); // x=C_1, y=R
  fy0_high->SetParameter(0,ymin_high); // y_min
  fy0_high->SetParameter(1,ymax_high); // y_max
}//Y Constructor


//! Pick an inelasticity y according to the model chosen
// double icemc::Y::pickY(const Settings *settings1,double pnu,int nu_nubar,Neutrino::Current currentint) {
double icemc::Y::pickY(const Settings *settings1,double pnu,Neutrino::L leptonNumber,Neutrino::Current currentint) {  
  if(settings1->YPARAM==0){
    return pickYGandhietal();
  }//old Gety
  else { //use prescription in Connolly et al.2011
    leptonNumber=Neutrino::L::Matter; ///@todo ? 
    double elast_y=pickYConnollyetal2011(leptonNumber,currentint,pnu);
    return elast_y;   
  }//current Gety
} //Gety


//! THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM Ghandhi,Reno,Quigg,Sarcevic  hep-ph/9512364 (the curves are not in their later article).  There is also a slow energy dependence.
double icemc::Y::pickYGandhietal() {
  double rnd;
  double x = 0;
  // generate according to Ghandi fig. 6 
  	// adjust exponent until looks like the curve
  	//  and has right mean.
  	//  (Note this is not the fcn, but the inverse of the integral...)
  rnd = pickUniform(1); //gRandom->Rndm(1); // (0,1)
  //  cout << "R1, R2, rnd are " << R1 << " " << R2 << " " << rnd << "\n";
  x=pow(-log(R1+rnd*R2),2.5);
  return x;   
}


// double icemc::Y::pickYConnollyetal2011(int NU, Neutrino::Current CURRENT,double pnu) {
double icemc::Y::pickYConnollyetal2011(Neutrino::L leptonNumber, Neutrino::Current CURRENT,double pnu) {

  // Select a y according to recipe in Connolly et al. (2011)
  //pnu is in eV.
  double epsilon = log10(pnu/1.E9);
  // pick a y region 
  //double R1=Rand3Y.Rndm(); // choose our first random number
  double r1 = pickUniform(); //gRandom->Rndm();
 
  int iyregion=0; // 0 for high y region 
  if (r1<ffrac->Eval(epsilon)){ // Is it going to be in low y region?
    iyregion=1; // 1 for low y region
  } 
  double C1_this;
  if (iyregion==0){ // high y region
    C1_this = fC1_high[{leptonNumber, CURRENT}].Eval(epsilon); // C1 for this event
  }
  else{ // low y region
    C1_this = fC1_low->Eval(epsilon); // C1 for this event
  }
  
  double C2_this =  fC2->Eval(epsilon); // C2 for this event
  
  // pick another random number
  double r2 = pickUniform(); //gRandom->Rndm();//  double r2=Rand3Y.Rndm();
  double y0 = 0.;
 
  if (iyregion==0){  // high y region
    y0 = fy0_high->Eval(C1_this,r2); // pick y0 according to Equation 15
  }
  else if (iyregion==1){  // low y region
    y0 = fy0_low->Eval(C1_this,C2_this,r2); // pick y0 according to Equation 14
  }
  return y0;
}//pickY


// double icemc::Y::Getyweight(double pnu, double y, int nu_nubar, Neutrino::Current currentint){
double icemc::Y::Getyweight(double pnu, double y, Neutrino::L leptonNumber, Neutrino::Current current){

  // int nu_nubar = leptonNumber == Neutrino::L::Matter ? 0 : 1;
  //from Connolly Calc 2011, Equations 9, 10, 11, 16, and 17.
  // double dy=0.;//default
  //Ev, cc or nc, nu or nubar.
  
  double C0_highbar, C0_lowbar,C0_high, C0_low;//these C0's are normalization factors.
  double dNdy=0.;//default
  double U, W, B, T;//are added in to help with readability of equations.
  double C1_low, C2, C1_high;
  double weighty;
  double epsilon=log10(pnu/1.E9);
  
  C2 = fC2->Eval(epsilon);//Eq(17)
  C1_low = fC1_low->Eval(epsilon);//Eq(16) (Low region)
  C1_high = fC1_high[{leptonNumber, current}].Eval(epsilon);//Eq(16)(High region) 
  
   
  if(leptonNumber==Neutrino::L::Matter){//nu_nubar==0) {
    U = 1-1/C2;
    W = fabs( (ymax_high-C1_high)/(ymin_high-C1_high));
    B = (pow(ymax_low-C1_low, 1/C2)/(ymax_low-C1_high));
    T = B*((pow(ymax_low-C1_low, U)-pow(ymin_low-C1_low, U) )/U)+log(W);
    C0_high = 1/T;	
    C0_low = C0_high*(pow(ymax_low-C1_low, 1/C2))/(ymax_low-C1_high);
    
    if(y<ymax_low){//Eq(9)
      // dy=0.00002;
      dNdy=C0_low/pow(y-C1_low, 1/C2);//Eq(10)
    }
    else if(y>=ymax_low && y<1.){//Eq(9)
      // dy=0.001;
      dNdy=C0_high/(y-C1_high);//Eq(10)
    }
    else{
      dNdy=0.;
      icemc::report() << severity::warning <<"y value is outside of the domain of y.\n";
    }
  }
  // else if(nu_nubar==1){
  else if(leptonNumber==Neutrino::L::AntiMatter){
    U=1-1/C2;
    W=fabs( (ymax_high-C1_high)/(ymin_high-C1_high));
    B=(pow(ymax_low-C1_low, 1/C2)/(ymax_low-C1_high));
    T=B*((pow(ymax_low-C1_low, U)-pow(ymin_low-C1_low, U) )/U)+log(W);
    C0_highbar=1/T;	
    C0_lowbar=C0_highbar*(pow(ymax_low-C1_low, 1/C2))/(ymax_low-C1_high);
   
    if(y<ymax_low){
      // dy=0.00002;
      dNdy=C0_lowbar/pow(y-C1_low, 1/C2);
    }
    else if(y>=ymax_low && y<1.){
      // dy=0.001;
      dNdy=C0_highbar/(y-C1_high);
    }
    else{
      dNdy=0;
      icemc::report() << severity::warning << "y value is outside of the domain of y.\n";
    }
  }
  else{
    icemc::report() << severity::warning << "Nu_nubar is not defined!\n";
  }		
  weighty=dNdy;
  return weighty;
}//Getyweight
