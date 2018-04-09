#include "IcemcLog.h"

icemc::Logger& icemc::Log(const char* outputdir, int run){
  static icemc::Logger log(outputdir, run);
  return log;
}


icemc::Logger::Logger(const char* outputDir, int run)
  : fOutputDir(outputDir), fRun(run), fMustReset(false), fUseStdErr(false), fUseColorCodes(true)
{


}


icemc::Logger::~Logger(){
  std::cout << getColorReset() << std::flush;
  std::cerr << getColorReset() << std::flush;
}



void icemc::Logger::openLogFiles(const char* outputDir, int run){

  // Update member vars if non-default arguments passed
  if(outputDir){
    fOutputDir = outputDir;
  }
  if(run > -1){
    fRun = run;
  }

  // convert the run to a string
  std::stringstream ss;
  ss << fRun;
  std::string run_num;
  ss >> run_num;


  // open the files
  
  std::string fileName = fOutputDir + "/nu_position" + run_num + ".txt";
  nu_out.open(fileName,  std::ios::app); //Positions,  direction of momentum,  and neutrino type for Ryan.

  fileName = fOutputDir + "/veff" + run_num + ".txt";
  veff_out.open(fileName,  std::ios::app);//to output only energy and effective volume to veff.txt
  
  fileName = fOutputDir + "/distance" + run_num + ".txt";
  distanceout.open(fileName, std::ios::app);

  fileName = fOutputDir + "/debug" + run_num + ".txt";
  outfile.open(fileName, std::ios::out);

  fileName = fOutputDir + "/forbrian" + run_num + ".txt";
  forbrian.open(fileName, std::ios::out);

  fileName = fOutputDir + "/al_voltages_direct" + run_num + ".dat";
  al_voltages_direct.open(fileName, std::ios::out); //added djg ------provide anita-lite voltages and noise from MC for anita-lite analysis

  fileName = fOutputDir + "/events" + run_num + ".txt";
  eventsthatpassfile.open(fileName, std::ios::out);

  fileName = fOutputDir + "/numbers" + run_num + ".txt";
  fnumbers.open(fileName, std::ios::out); // debugging

  fileName = fOutputDir + "/output" + run_num + ".txt";
  foutput.open(fileName,  std::ios::app);

  fileName = fOutputDir + "/slacviewangles" + run_num + ".dat";
  fslac_viewangles.open(fileName, std::ios::out); // this outputs numbers that we need for analyzing slac data

  fileName = fOutputDir + "/slac_hitangles" + run_num + ".dat";
  fslac_hitangles.open(fileName,  std::ios::out);
  
}

icemc::Logger& icemc::Logger::message(icemc::severity s){

  const char* red     = fUseColorCodes ? "\x1b[31m" : "";
  const char* blue    = fUseColorCodes ? "\x1b[34m" : "";
  const char* magenta = fUseColorCodes ? "\x1b[35m" : "";

  fMustReset = true;
  switch(s){
  case info:
    fUseStdErr = false;
    getStream() << blue << "Info! ";
    foutput << "Info! ";
    break;
  case warning:
    fUseStdErr = true;
    getStream() << magenta << "Warning! ";
    foutput << "Warning! ";
    break;
  case error:
    fUseStdErr = true;
    getStream() << red << "Error! ";
    foutput << "Error! ";
    break;
  }
  return *this;
}
