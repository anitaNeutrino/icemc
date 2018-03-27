#ifndef ICEMC_TEXT_FILE_OUTPUT_H
#define ICEMC_TEXT_FILE_OUTPUT_H

#include <fstream>

namespace icemc{

  class TextFileOutput{
  public:

    TextFileOutput(const char* outputDir, int run);
    std::string fOutputDir;
    
    std::ofstream nu_out;
    std::ofstream veff_out;
    std::ofstream distanceout;
    std::ofstream outfile;
    std::ofstream forbrian;
    std::ofstream al_voltages_direct;
    std::ofstream eventsthatpassfile;
    std::ofstream fnumbers;
    std::ofstream foutput;
    std::ofstream fslac_viewangles;
    std::ofstream fslac_hitangles;

  private:
    std::string fOutputdir;
    int fRun;
  };

  


}



#endif
