#include "CommandLineOpts.h"
#include "Settings.h"

#include <iostream>
#include <unistd.h> // for getopt
#include <sstream>



// Prettify warnings 
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"




icemc::CommandLineOpts::CommandLineOpts(int argc, char* argv[], Settings& settings) :
  outputdir("."), input(""), run_no(0)
{

  setExecName(argv[0]);
  
  startNu = 0;
  nnu_tmp = 0;
  exp_tmp = 0;  
  trig_thresh = 0.;
  
  char clswitch; // command line switch
  std::string mkdirCommand("mkdir -p ");
  bool got_t = false;
  bool got_n = false;
  bool got_r = false;
  bool got_i = false;
  bool got_o = false;
  bool got_e = false;
  bool got_x = false;
  bool got_h = false;
  if (argc>1) {
    while ((clswitch = getopt(argc, argv, "ht:i:o:r:n:e:x:")) != EOF && !got_h) { // abort early if -h passed
      switch(clswitch) {
      case 'h':
	if(!got_h){ // print help message
	  got_h = true;
	}
	break;
      case 'n':
	if(!got_n){
	  got_n = true;	
	  nnu_tmp = atoi(optarg);
	  std::cout << ANSI_COLOR_GREEN << "Changed number of simulated neutrinos to " << nnu_tmp << ANSI_COLOR_RESET << std::endl;
	}
        break;
      case 't':
	if(!got_t){
	  got_t = true;	
	  trig_thresh = atof(optarg);
	  std::cout << ANSI_COLOR_GREEN << "Set trigger threshold to " << trig_thresh << ANSI_COLOR_RESET << std::endl;
	}
        break;
      case 'i':
	if(!got_i){
	  got_i = true;	
	  input = optarg;
	  std::cout << ANSI_COLOR_GREEN << "Changed input file to: " << input << ANSI_COLOR_RESET << std::endl;
	}
        break;
      case 'o':
	if(!got_o){
	  got_o = true;
	  outputdir = optarg;
	  std::cout << ANSI_COLOR_GREEN << "Changed output directory to: " << outputdir << ANSI_COLOR_RESET << std::endl;
	  mkdirCommand += outputdir;
	  system(mkdirCommand.c_str());
	}
        break;
      case 'e':
	if(!got_e){
	  got_e = true;
	  exp_tmp = atof(optarg);
	  std::cout << ANSI_COLOR_GREEN << "Changed neutrino energy exponent to " << exp_tmp << ANSI_COLOR_RESET << std::endl;
	}
	break;
      case 'x':
	if(!got_x){
	  got_x = true;
	  startNu = atoi(optarg);
	  std::cout << ANSI_COLOR_GREEN << "Running icemc for just 1 event with eventNumber : " << startNu << ANSI_COLOR_RESET << std::endl;
	}
	break;
      case 'r':
	if(!got_r){
	  got_r = true;	
	  run_num = optarg;
	  std::stringstream convert(run_num);
	  convert>>run_no;
	}
        break;
      }
    }
  }

  // post parsing logic...
  if(got_x){ // change the total so we only do 1 neutrino
    nnu_tmp = startNu+1;
  }

  /// @todo validate this logic for whether the command line options were sufficient
  are_good = true;
  if(got_h){
    // we just want to print the help message
    are_good = false;
  }
  if(!got_i){ // always need input file
    are_good = false;
  }
  else if(!got_r || !got_i){
    are_good = false;
  }
  
  if(!are_good){ // print a pretty error message
    if(!got_h){
      fprintf(stderr, "%sError in %s, insufficient command line arguments passed%s\n",    ANSI_COLOR_RED, __PRETTY_FUNCTION__,  ANSI_COLOR_RESET);
    }
    fprintf(stderr, "%s%s basic options:%s\n",                                            ANSI_COLOR_BLUE, executable.c_str(), ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-h%s print help message, do nothing else                  %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);
    fprintf(stderr, "%s%s basic options:%s\n",                                            ANSI_COLOR_BLUE, executable.c_str(), ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-n%s number of neutrinos to simulate                      %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);    
    fprintf(stderr, "  %s-r%s run number                                           %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-i%s config input file                                    %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-o%s output directory                                     %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-e%s neutrino exponent (e.g. 20 for 1e20 eV neutrinos)    %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);
    fprintf(stderr, "%s%s advanced options:%s\n",                                         ANSI_COLOR_MAGENTA, executable.c_str(), ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-x%s simulate just one event with this eventNumber        %s\n", ANSI_COLOR_RED, ANSI_COLOR_MAGENTA, ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-t%s trigger threshold                                    %s\n", ANSI_COLOR_RED, ANSI_COLOR_MAGENTA, ANSI_COLOR_RESET);
  }
  else{

    // update the settings object with the extracted info    
    std::ofstream foutput("temp.txt"); // need to replace this with the proper output at some point...
    settings.ReadInputs(input.c_str(),  foutput); //, NNU, RANDOMISEPOL);

    if (exp_tmp!=0){
      // notify?
      settings.EXPONENT = exp_tmp;
    }
    if(got_o){
      // notify?
      settings.outputDir = outputdir;
    }


    settings.SEED += run_no; // uniquify per-run seed by adding run number
  }

}


void icemc::CommandLineOpts::setExecName(char* argv0){

  std::string fullPath(argv0);
  executable = fullPath;

  std::size_t found = fullPath.rfind("/");
  if (found!=std::string::npos){
    executable = fullPath.substr(found+1);
    // std::cout << found << "\t" << fullPath << "\t" << executable << std::endl;
  }
}
