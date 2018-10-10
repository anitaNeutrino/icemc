#include "OldShowerGenerator.h"
#include "ShowerGenerator.h"
#include "Settings.h"
#include "Neutrino.h"

using namespace icemc;

int main(int argc, char* argv[]){

  Settings settings;
  CommandLineOptions clOpts(argc, argv, settings);

  ShowerGenerator sg(&settings);
  OldShowerGenerator osg(&settings);

  for(auto flavor : {Neutrino::Flavor::e, Neutrino::Flavor::mu, Neutrino::Flavor::tau}) {
    for(auto current : {Neutrino::Interaction::Current::Charged, Neutrino::Interaction::Current::Neutral}){
      Energy energy = 1*Energy::Unit::EeV;
      while(energy <= 1*Energy::Unit::ZeV){
	std::string taudecay = "";
	double y = 1;
	TH1F *hy = nullptr;
	int inu = 0;
	int taumodes1 = 0;
	OldShower os = osg.GetEMFrac(flavor, current, taudecay, y, hy, energy.in(Energy::Unit::eV),  inu, taumodes1);
	std::cout << flavor  <<  "\t" << current  << "\t" << energy << "\t\t\t\t" << os.emFrac << "\t" << os.hadFrac << "\t\t\t\t";

	Shower s = sg.GetEMFrac(flavor, current, y, energy);
	std::cout << s.emFrac << "\t" << s.hadFrac << std::endl;

	Energy dE = pow(10, floor(TMath::Log10(energy.in(Energy::Unit::eV))))*Energy::Unit::eV;
	energy += dE;
      }
    }
  }
  return 0;
}
