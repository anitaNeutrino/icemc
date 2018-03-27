#include "CommandLineOpts.h"
#include "Settings.h"

#include <iostream>
#include <unistd.h> // for getopt
#include <sstream>




// Prettify because, why not?
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"




icemc::CommandLineOpts::CommandLineOpts(int argc, char* argv[]){

  setExecName(argv[0]);

  nnu_tmp = 0;
  exp_tmp = 0;
  trig_thresh = 0.;
  startNu = 0;
  outputdir = "."; // use pwd
  
  char clswitch; // command line switch
  std::string mkdirCommand("mkdir -p ");
  bool got_t = false;
  bool got_n = false;
  bool got_r = false;
  bool got_i = false;
  bool got_o = false;
  bool got_e = false;
  bool got_x = false;
  
  if (argc>1) {
    while ((clswitch = getopt(argc, argv, "t:i:o:r:n:e:x:")) != EOF) {
      switch(clswitch) {
      case 'n':
	nnu_tmp = atoi(optarg);
	got_n = true;
	std::cout << "Changed number of simulated neutrinos to " << nnu_tmp << std::endl;
        break;
      case 't':
	trig_thresh = atof(optarg);
	got_t = true;	
        break;
      case 'i':
        input = optarg;
	got_i = true;
        std::cout << "Changed input file to: " << input << std::endl;
        break;
      case 'o':
        outputdir = optarg;
        std::cout << "Changed output directory to: " << outputdir << std::endl;
	mkdirCommand += outputdir;
        system(mkdirCommand.c_str());
	got_o = true;
        break;
      case 'e':
	exp_tmp = atof(optarg);
	got_e = true;
	std::cout << "Changed neutrino energy exponent to " << exp_tmp << std::endl;
	break;
      case 'x':
	startNu = atoi(optarg);
	got_x = true;
	std::cout << "Running icemc for just 1 event with eventNumber : " << startNu << std::endl;
	break;
      case 'r':
        run_num=optarg;
	std::stringstream convert(run_num);
        convert>>run_no;
        break;
      }
    }
  }



  // Logic for whether the command line options were sufficient
  are_good = true;
  if(!got_i){ // always need input file
    are_good = false;
  }
  else if(got_x){ // 
    run_no = 0;
    nnu_tmp = startNu+1;
    are_good = true;
  }
  else if(!got_r || !got_i){
    are_good = false;    
  }

  if(!are_good){
    fprintf(stderr, "%sError in %s, insufficient command line arguments passed%s\n",  ANSI_COLOR_RED, __PRETTY_FUNCTION__,  ANSI_COLOR_RESET);
    fprintf(stderr, "%s%s basic options:%s\n",                                        ANSI_COLOR_BLUE, executable.c_str(), ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-n%s number of neutrinos to simulate                      %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-r%s run number                                           %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-i%s config input file                                    %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-o%s output directory                                     %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-e%s neutrino exponent (e.g. 20 for 1e20 eV neutrinos)    %s\n", ANSI_COLOR_RED, ANSI_COLOR_BLUE, ANSI_COLOR_RESET);
    fprintf(stderr, "%s%s advanced options:%s\n",                                     ANSI_COLOR_MAGENTA, executable.c_str(), ANSI_COLOR_RESET);
    fprintf(stderr, "  %s-x%s simulate just one event with this eventNumber        %s\n", ANSI_COLOR_RED, ANSI_COLOR_MAGENTA, ANSI_COLOR_RESET);
    // Yeah, at the time of separating this out into it's own class, I don't know what this switch does...
    fprintf(stderr, "  %s-t%s turn on some kind of trigger threshold?              %s\n", ANSI_COLOR_RED, ANSI_COLOR_MAGENTA, ANSI_COLOR_RESET);
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
