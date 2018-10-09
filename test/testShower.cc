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

  double energy = 1e18;
  while(energy <= 1e21){

    Neutrino::Flavor flavor = Neutrino::Flavor::e;
    Neutrino::Interaction::Current current = Neutrino::Interaction::Current::Charged;
    std::string taudecay = "";
    double y = 1;
    TH1F *hy = nullptr;
    double pnu = energy;
    int inu = 0;
    int taumodes1  = 0;
    OldShower os = osg.GetEMFrac(flavor, current, taudecay, y, hy, pnu,  inu, taumodes1);
    std::cout << pnu <<  "\t\t\t\t" << os.emFrac << "\t" << os.hadFrac << "\t\t\t\t";

    Shower s = sg.GetEMFrac(flavor, current, taudecay, y, hy, pnu,  inu, taumodes1);
    std::cout << s.emFrac << "\t" << s.hadFrac << std::endl;

    energy += pow(10, floor(TMath::Log10(energy)));
  }
  return 0;
}
