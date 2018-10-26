#include <sstream>
#include "ConnollyEtAl2011.h"
#include "Report.h"
#include "Settings.h"

icemc::ConnollyEtAl2011::CrossSectionModel::CrossSectionModel(const icemc::Settings* settings)
  : icemc::CrossSectionModel(settings)
{//constructor

  for(auto l : {Neutrino::L::Matter, Neutrino::L::AntiMatter}){
    for(auto c : {Neutrino::Interaction::Current::Neutral, Neutrino::Interaction::Current::Charged}){
      std::stringstream name;
      name << "fSigma_" << l << "_" << c;

      auto p = std::make_pair(l,c);
      auto mapIter_bool = fSigma.emplace(p, TF1(name.str().c_str(),"pow(10, [1]+[2]*log(x-[0])+[3]*pow(log(x-[0]),2)+[4]/log(x-[0]))", 4., 21.)); //check bounds. they're in log10 
      auto mapIter = mapIter_bool.first;
      mapIter->second.SetParameters(C0.at(p),
				    C1.at(p),
				    C2.at(p),
				    C3.at(p),
				    C4.at(p));
    }
  }
  
  // 0=Reno
  // 1=Connolly et al. 2011
  // mine[0] = 1.2E15;
  // mine[1] = 1.E4;
  // maxe[0] = 1.E21;
  // maxe[1] = 1.E21;

  fMinEnergy = 1.E4*Energy::Unit::eV; // minimum energy for cross section parametrisations
  fMaxEnergy = 1.E21*Energy::Unit::eV; // use the same upper limit for reno as for Connolly et al.  
}




double icemc::ConnollyEtAl2011::CrossSectionModel::getSigma(icemc::Energy energy, icemc::Neutrino::L leptonNumber, icemc::Neutrino::Interaction::Current current) const {
  
  // int currentint = static_cast<int>(current);
  // int nu_nubar = leptonNumber == Neutrino::L::Matter ? 0 : 1;
  // calculate cross section
  if(!validEnergy(energy)){
    icemc::report() << severity::error << "Need a parameterization for this energy region, energy = "
		    << energy << " but min = " << fMinEnergy << ", max = " << fMaxEnergy << std::endl;
    return -1;
  }
  double epsilon=log10(energy.in(Energy::Unit::GeV));
  const TF1& sigmaVsEpsilon = fSigma.at({leptonNumber, current});
  double sigma = sigmaVsEpsilon.Eval(epsilon);
  sigma /= 1.E4;//convert cm to meters. multiply by (1m^2/10^4 cm^2).
  sigma *= fSettings ? fSettings->SIGMA_FACTOR : 1; // any constant factor from settings
  return sigma;
}
































icemc::ConnollyEtAl2011::YGenerator::YGenerator(const Settings* s) :
  fSettings(s),
  fFracLowHigh("fFracLowHigh", equation18.c_str(), epsLow, epsHigh), // This is the fraction of the distribution in the low y region given by Equation 18.   
  fC1_low("fC1_low", equation16.c_str(), epsLow, epsHigh),
  fC2("fC2", equation17.c_str(), epsLow, epsHigh), 
  fy0_low("fy0_low", equation14.c_str()),
  fy0_high("fy0_high", equation15.c_str())
{

  fFracLowHigh.FixParameter(0, F0);
  fFracLowHigh.FixParameter(1, F1);
  fFracLowHigh.FixParameter(2, F2);

  const std::map<std::pair<icemc::Neutrino::L, icemc::Neutrino::Interaction::Current>, EColor> fCols =
    {{{icemc::Neutrino::L::AntiMatter, icemc::Neutrino::Interaction::Current::Charged}, kRed},
     {{icemc::Neutrino::L::Matter,     icemc::Neutrino::Interaction::Current::Charged}, kBlue },
     {{icemc::Neutrino::L::AntiMatter, icemc::Neutrino::Interaction::Current::Neutral}, kGreen },
     {{icemc::Neutrino::L::Matter,     icemc::Neutrino::Interaction::Current::Neutral}, kMagenta }};
  
  for(auto l : {Neutrino::L::Matter, Neutrino::L::AntiMatter}){
    for(auto c : {Neutrino::Interaction::Current::Charged,
		  Neutrino::Interaction::Current::Neutral}){
      std::stringstream ss;
      ss  << "fC1_high_" << l << "_" <<  c;
      std::string name = ss.str();
      std::replace(name.begin(), name.end(), ':', '_');
      auto p = std::make_pair(l, c);
      auto mapIter_bool = fC1_high.emplace(p, TF1(name.c_str(), equation16.c_str(), 7.,12.)); // parameterization of parameter C1 in the high y region according to Equation 16
      auto mapIter = mapIter_bool.first;
      TF1& f = mapIter->second;
      f.FixParameter(0, A0_high.at(p));
      f.FixParameter(1, A1_high.at(p));
      f.FixParameter(2, A2_high.at(p));
      f.FixParameter(3, A3_high.at(p));

      ss << "; log_{10}(E) (GeV); C1 (unitless)";
      std::string title = ss.str();
      std::replace(title.begin(), title.end(), '_', ' ');
      f.SetTitle(title.c_str());
      f.SetLineColor(fCols.at(p));
      f.SetLineWidth(3);
      f.SetFillColor(0);

      // std::cout << l << "\t" << c << std::endl;
      // for(int p=0; p < 4; p++){
      // 	std::cout << f.GetParameter(p) << std::endl;
      // }
      // std::cout << std::endl;
    }
  }
  
  fC1_low.FixParameter(0,A0_low);
  fC1_low.FixParameter(1,A1_low);
  fC1_low.FixParameter(2,A2_low);
  fC1_low.FixParameter(3,A3_low);

  fC2.FixParameter(0,B0);
  fC2.FixParameter(1,B1);

  // For picking inelasticity in low y region according to Equation 14.
  fy0_low.FixParameter(0,ymin_low);  // y_min
  fy0_low.FixParameter(1,ymax_low);  // y_max

  // For picking inelasticity in high y region according to Equation 15.
  fy0_high.SetParameter(0,ymin_high); // y_min
  fy0_high.SetParameter(1,ymax_high); // y_max
}




// void icemc::ConnollyEtAl2011::YGenerator::Draw(Option_t* opt){
TH1D* icemc::ConnollyEtAl2011::YGenerator::plot(Energy energy, Neutrino::L l, Neutrino::Interaction::Current c, int n) {

  std::stringstream ss;
  ss << "h_"  << energy << "_" << l << "_" << c;
  std::string s = ss.str();
  std::replace( s.begin(), s.end(), ' ', '_');  
  std::replace( s.begin(), s.end(), ':', '_');
  TH1D* h = new TH1D(s.c_str(), s.c_str(), n/100, 0, 1);
  for(int i=0; i < n; i++){
    double y = pickY(energy, l, c);
    h->Fill(y);
  }
  return h;
}



//! Pick an inelasticity y according to the model chosen
// double icemc::YGenerator::pickY(const Settings *settings1,double pnu,int nu_nubar,Neutrino::Interaction::Current currentint) {
double icemc::ConnollyEtAl2011::YGenerator::pickY(Energy energy,
						  Neutrino::L leptonNumber,
						  Neutrino::Interaction::Current current) {
  // Select a y according to recipe in Connolly et al. (2011)
  double epsilon = log10(energy.in(Energy::Unit::GeV));

  // pick a y region 
  auto region = pickUniform() < fFracLowHigh.Eval(epsilon) ? YRegion::Low : YRegion::High;

  double C1_this = 0;

  switch(region){
  case YRegion::High:
    C1_this = fC1_high[{leptonNumber, current}].Eval(epsilon);
    break;
  case YRegion::Low:
    C1_this = fC1_low.Eval(epsilon);
    break;
  }
  
  double C2_this =  fC2.Eval(epsilon); // C2 for this event
  double R = pickUniform();
  double y = 0.;
  switch(region){
  case YRegion::High:
    {
      y = fy0_high.Eval(C1_this,R); // pick y0 according to Equation 15
      if(TMath::IsNaN(y)){
	icemc::report() << "Warning got NaN for inelasticity with " << energy << ", " << leptonNumber << "\t" << current << std::endl;
      }
    }
    break;  
  case YRegion::Low:  // low y region
    y = fy0_low.Eval(C1_this,C2_this,R); // pick y0 according to Equation 14
    break;
  }
  return y;

    
  // leptonNumber=Neutrino::L::Matter; ///@todo ? 
}




// double icemc::YGenerator::Getyweight(double pnu, double y, int nu_nubar, Neutrino::Interaction::Current currentint){
double icemc::ConnollyEtAl2011::YGenerator::Getyweight(Energy pnu, double y, Neutrino::L leptonNumber, Neutrino::Interaction::Current current){

  // int nu_nubar = leptonNumber == Neutrino::L::Matter ? 0 : 1;
  //from Connolly Calc 2011, Equations 9, 10, 11, 16, and 17.
  // double dy=0.;//default
  //Ev, cc or nc, nu or nubar.
  
  double C0_highbar, C0_lowbar,C0_high, C0_low;//these C0's are normalization factors.
  double dNdy=0.;//default
  double U, W, B, T;//are added in to help with readability of equations.
  double C1_low, C2, C1_high;
  double weighty;
  double epsilon=log10(pnu.in(Energy::Unit::GeV));
  
  C2 = fC2.Eval(epsilon);//Eq(17)
  C1_low = fC1_low.Eval(epsilon);//Eq(16) (Low region)
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
