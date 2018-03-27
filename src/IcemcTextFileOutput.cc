#include "IcemcTextFileOutput.h"

#include <sstream>
#include <string>

icemc::TextFileOutput::TextFileOutput(const char* outputDir, int run)
  : fOutputDir(outputDir), fRun(run)
{
  std::stringstream ss;
  ss << run;
  std::string run_num;
  ss >> run_num;
  
  std::string fileName = fOutputDir+"/nu_position"+run_num+".txt";
  nu_out.open(fileName,  std::ios::app); //Positions,  direction of momentum,  and neutrino type for Ryan.

  fileName = fOutputDir+"/veff"+run_num+".txt";
  veff_out.open(fileName,  std::ios::app);//to output only energy and effective volume to veff.txt
  
  fileName = fOutputDir+"/distance"+run_num+".txt";
  distanceout.open(fileName, std::ios::app);

  fileName = fOutputDir+"/debug"+run_num+".txt";
  outfile.open(fileName, std::ios::out);

  fileName = fOutputDir+"/forbrian"+run_num+".txt";
  forbrian.open(fileName, std::ios::out);

  fileName = fOutputDir+"/al_voltages_direct"+run_num+".dat";
  al_voltages_direct.open(fileName, std::ios::out); //added djg ------provide anita-lite voltages and noise from MC for anita-lite analysis

  fileName = fOutputDir+"/events"+run_num+".txt";
  eventsthatpassfile.open(fileName, std::ios::out);

  fileName = fOutputDir+"/numbers"+run_num+".txt";
  fnumbers.open(fileName, std::ios::out); // debugging

  fileName = fOutputDir+"/output"+run_num+".txt";
  foutput.open(fileName,  std::ios::app);

  fileName = fOutputDir+"/slacviewangles"+run_num+".dat";
  fslac_viewangles.open(fileName, std::ios::out); // this outputs numbers that we need for analyzing slac data

  fileName = fOutputDir+"/slac_hitangles"+run_num+".dat";  
  fslac_hitangles.open(fileName,  std::ios::out);
  
}

