#include "Inelasticity.h"
#include "ConnollyEtAl2011.h"
#include "Settings.h"

#include <fstream>
#include "TRandom3.h"
int main(){
  using namespace icemc;

  Settings s;
  std::ofstream os;
  s.ReadInputs("inputs.anita3.conf", os);

  ConnollyEtAl2011::YGenerator yGen(&s);

  TRandom3 randy;
  double epsMin = 7;
  double epsMax = 12;
  const int n = 1000;
  const double deltaEps = (epsMax - epsMin)/(n-1);
  for(int i=0; i < n; i++){
    double epsilon = epsMin + deltaEps*i; //randy.Uniform(7, 12);
    double energy_eV = pow(10.0, 9+epsilon);
    auto energy = energy_eV*Energy::Unit::eV;

    std::vector<double> ys; ys.reserve(4);
    for(auto l : {Neutrino::L::Matter, Neutrino::L::AntiMatter}){
      for(auto c : {Neutrino::Interaction::Current::Neutral, Neutrino::Interaction::Current::Charged}){
	double y = yGen.pickY(energy, l, c);
	ys.push_back(y);
      }
    }
    std::cout << energy;
    for(auto y : ys){std::cout << "\t" << y;}
    std::cout << std::endl;
  }
  
}
