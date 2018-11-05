#include "ShowerModel.h"
#include "CommandLineOptions.h"
#include "Settings.h"
#include "Neutrino.h"

using namespace icemc;

int main(int argc, char* argv[]){

  Settings settings;
  CommandLineOptions clOpts(argc, argv, settings);

  ShowerModel sg(&settings);

  for(auto flavor : {Neutrino::Flavor::e, Neutrino::Flavor::mu, Neutrino::Flavor::tau}) {
    for(auto current : {Interaction::Current::Charged, Interaction::Current::Neutral}){
      Energy energy = 1*Energy::Unit::EeV;
      while(energy <= 1*Energy::Unit::ZeV){
	double y = 1;
	Shower s = sg.GetEMFrac(flavor, current, y, energy);
	std::cout << s.emFrac << "\t" << s.hadFrac << std::endl;

	Energy dE = pow(10, floor(TMath::Log10(energy.in(Energy::Unit::eV))))*Energy::Unit::eV;
	energy += dE;
      }
    }
  }
  return 0;
}
