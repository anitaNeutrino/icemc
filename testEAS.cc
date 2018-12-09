/* -*- C++ -*-.*********************************************************************************************
   Author: Linda Cremonesi
   Email: l.cremonesi@ucl.ac.uk

   Description:
   Template to generate triggers coming from EAS
***********************************************************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <ctype.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <time.h>
#include "TTreeIndex.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
#include "Math/InterpolationTypes.h"
#include "Math/Interpolator.h"
#include "TGaxis.h"
#include "signal.h"

#include "EnvironmentVariable.h"

//#include "rx.hpp"
#include "Constants.h"
#include "Settings.h"
#include "position.hh"

#include "earthmodel.hh"
#include "Tools.h"
#include "vector.hh"
#include "roughness.hh"
#include "anita.hh"
#include "balloon.hh"
#include "icemodel.hh"
// #include "trigger.hh"
#include "Spectra.h"
#include "signal.hh"
#include "secondaries.hh"
#include "ray.hh"
#include "counting.hh"
#include "Primaries.h"
#include "Taumodel.hh"
#include "screen.hh"
#include "GlobalTrigger.h"
#include "ChanTrigger.h"

#include <string>
#include <sstream>

#if __cplusplus > 199711L
#include <type_traits>
#endif

#include <typeinfo>

#ifdef ANITA_UTIL_EXISTS
#include "UsefulAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "RawAnitaHeader.h"
#include "Adu5Pat.h"
#include "FFTtools.h"
UsefulAnitaEvent*     realEvPtr    = NULL;
RawAnitaHeader*       rawHeaderPtr = NULL;
Adu5Pat*         Adu5PatPtr   = NULL;
#ifdef ANITA3_EVENTREADER
#include "TruthAnitaEvent.h"
TruthAnitaEvent*      truthEvPtr   = NULL;
#endif
#endif
#include "TApplication.h"
#include <fftw3.h>
#include "hot-loop.h"
#include "cr-ft.h"
#include "hot-module-test.h"

namespace bvv {
  
  //! Fill amplitude and phase arrays given FT of a ZHS waveform.
  int interpolate_zhs_ft
  (
   double (&amp)[Anita::NFREQ]   /*!< [out] ft amplitude. */,
   double (&phase)[2 * Anita::NFREQ] /*!< [out] ft phase. */,
   Anita *anita /*!< [in] Anita class with the highest frequency bin. */,
   cr_ft_state *cr_ft_result  /*!< [in] a struct returned from from the zhs ft module. */
   )
  {
    // cr_ft_result->FftPhi, cr_ft_result->vis_nbins / 2 + 1
    for (int i = 0; i < Anita::NFREQ; i++) {
      double target_freq = anita->freq[i];
        // double target_freq = (anita->FREQ_HIGH  - anita->FREQ_LOW) / (double) Anita::NFREQ * i + anita->FREQ_LOW;
      double target_freq_src_units = target_freq / cr_ft_result->dfreq;
      int int_target_freq_src_units = target_freq_src_units + 0.5;
      // Hack alert!:
      // amp[i] = cr_ft_result->FftRho[int_target_freq_src_units];
      amp[i] = 1;
    }
    for (int i = 0; i < Anita::NFREQ * 2; i++) {
      double target_freq = anita->freq_forfft[i * 2];
      // double target_freq = anita->FREQ_HIGH / (2 * Anita::NFREQ) * i;
      double target_freq_src_units = target_freq / cr_ft_result->dfreq;
      int int_target_freq_src_units = target_freq_src_units + 0.5;
      phase[i] = cr_ft_result->FftPhi[int_target_freq_src_units];
    }
      
    return 0;
  }
  }


extern const double pi = atan(1)*4;

// cpp'ed version of https://stackoverflow.com/questions/9210528/split-string-with-delimiters-in-c
vector<string> str_split(char* a_str, const char a_delim)
{
  char delim[2];
  delim[0] = a_delim;
  delim[1] = 0;

  vector<string> result;

  char* token = strtok(a_str, delim);

  while (token) {
    result.push_back(token);
    token = strtok(0, delim);
  }

  return result;
}


// https://stackoverflow.com/questions/3501338/c-read-file-line-by-line:
// Notice that index "_with_line_ind" starts from 1:
#define WITH_LINES(_with_lines_name, _with_lines_ind, _with_lines_tokens, _with_lines_codeblock...) \
  do \
    {									\
      vector<string> _with_lines_tokens;				\
      FILE * _with_lines_fp;						\
      char * _with_lines_line = NULL;					\
      size_t _with_lines_len = 0;					\
      ssize_t _with_lines_read;						\
      int _with_lines_ind = 0;						\
      _with_lines_fp = fopen(_with_lines_name, "r");						\
      if (_with_lines_fp == NULL)					\
	exit(EXIT_FAILURE);						\
      while ((_with_lines_read = getline(&_with_lines_line, &_with_lines_len, _with_lines_fp)) != -1) { \
        _with_lines_ind++;						\
	_with_lines_tokens = str_split(_with_lines_line,' ');		\
	_with_lines_codeblock						\
	  }								\
      fclose(_with_lines_fp);						\
      if (_with_lines_line)						\
	free(_with_lines_line);						\
    }									\
  while(0)

#define SET3_IND(_SET3_IND_TARGET_ROOT_NAME, _SET3_IND_NAME_SPECIFIER1, _SET3_IND_NAME_SPECIFIER2, _SET3_IND_NAME_SPECIFIER3, _SET3_IND_TARGET_CATEGORY, _SET3_IND_SRC_IND) \
do									\
{									\
 _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER1##_SET3_IND_TARGET_CATEGORY = _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER1[_SET3_IND_SRC_IND]; \
 _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER2##_SET3_IND_TARGET_CATEGORY = _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER2[_SET3_IND_SRC_IND]; \
 _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER3##_SET3_IND_TARGET_CATEGORY = _SET3_IND_TARGET_ROOT_NAME##_SET3_IND_NAME_SPECIFIER3[_SET3_IND_SRC_IND]; \
 } \
while(0)

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

Taumodel* TauPtr = NULL;

const string ICEMC_SRC_DIR = EnvironmentVariable::ICEMC_SRC_DIR();

ClassImp(RX);

using namespace std;

class EarthModel;
class Position;


bool FWHM(long int n, double *x, double *y, double &xmin, double &xmax, int &ind_maxval, double &xmaxval, double &ymaxval, double threshold_rel = 0.5){
  // xmin, xmax: interval where y > ymaxval * threshold_rel;
  long int ind_max = TMath::LocMax(n, y);
  ind_maxval = ind_max;
  double max_val = y[ind_max]; 
  ymaxval = max_val;
  xmaxval = x[ind_max];
  double threshold_abs = max_val * threshold_rel;
  bool found_xmin = false;
  bool found_xmax = false;
  for (int i = 0; i <= ind_max; i++){
    double val = y[i];
    if (val > threshold_abs){
      xmin = (x[i] - x[i - 1]) / (y[i] - y[i - 1]) * (threshold_abs - y[i - 1]) + x[i - 1];
      printf("ind_a: %d\n", i);
      found_xmin = true;
      break;
    }
  }
  if (!found_xmin) return false;

  for (int i = n-1; i >= ind_max; i--){
    double val = y[i];
    if (val > threshold_abs){
      xmax = (x[i] - x[i + 1]) / (y[i] - y[i + 1]) * (threshold_abs - y[i + 1]) + x[i + 1];
      printf("ind_b: %d\n", i);
      found_xmax = true;
      break;
    }
  }
  return (found_xmax && found_xmin);
}

#ifdef ANITA_UTIL_EXISTS
int GetIceMCAntfromUsefulEventAnt(Settings *settings1,  int UsefulEventAnt);
#ifdef R_EARTH
#undef R_EARTH
#endif
#endif

bool ABORT_EARLY = false;    // This flag is set to true when interrupt_signal_handler() is called


void interrupt_signal_handler(int sig);

int ZhsTimeN = -999;
double ZhsTimeStart = -999;
double ZhsTimeDelta = -999;
vector<double> ZhsTimeArr;
vector<double> ZhsTimeE(25000);
vector<double> ZhsAlpha(25000);
Settings* global_settings1;
Anita *global_anita1;
Balloon *global_bn1;
double NrFT[ANITA_TIME_SAMPLES];

int main(int argc,  char **argv) {

  extern void InitGui();
  // VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
  // TROOT HSroot("HSroot", "Some very smart title", initfuncs);
  // TApplication theApp("App", 0, 0);
  TApplication *theApp __attribute__ ((unused)) = new TApplication("tapp", NULL, NULL);
  // theApp = theApp;
  gROOT->SetStyle("Plain");

  string stemp;

  Settings* settings1 = new Settings();
  global_settings1 = settings1;

  string input="inputs.txt";
  string sim_inp="/nfs/data_disks/herc0a/users/bugaev/ANITA/SIMS/Event_4212/timefresnel-root.dat";
  string run_num;//current run number as string
  int run_no = 0;//current run number as integer
  TString outputdir;

  
  if( (argc%3!=1)&&(argc%2!=1)) {
    cout << "Syntax for options: -i inputfile -o outputdir -r run_number\n";
    return 0;
  }
  int nnu_tmp=0;
  double exp_tmp=0;
  double trig_thresh=0.;
  char clswitch; // command line switch
  if (argc>1) {
    while ((clswitch = getopt(argc, argv, "t:i:o:r:n:e:s:")) != EOF) {
      switch(clswitch) {
      case 'n':
	nnu_tmp=atoi(optarg);
	cout << "Changed number of simulated neutrinos to " << nnu_tmp << endl;
        break;
      case 't':
	trig_thresh=atof(optarg);
        break;
      case 'i':
        input=optarg;
        cout << "Changed input file to: " << input << endl;
        break;
      case 'o':
        outputdir=optarg;
        cout << "Changed output directory to: " << outputdir.Data() << endl;
        stemp="mkdir " + string(outputdir.Data());
        system(stemp.c_str());
        break;
      case 'e':
	exp_tmp=atof(optarg);
	cout << "Changed neutrino energy exponent to " << exp_tmp << endl;
	break;
      case 'r':
        {	
          run_num=optarg;
          stringstream convert(run_num);
          convert>>run_no;
          break;
        }
      case 's':
        sim_inp=optarg;
        cout << "Changed input simulation directory to: " << sim_inp << endl;
        break;
      } // end switch
    } // end while
  } // end if arg>1

  
  vector<double> Ex, Ey, Ez; 
  WITH_LINES(
             sim_inp.c_str(),
             ind,
             tokens,
             if (ind > 19) {
               Ex.push_back(atof(tokens[11].c_str()));
               Ey.push_back(atof(tokens[12].c_str()));
               Ez.push_back(atof(tokens[13].c_str()));
               ZhsTimeArr.push_back(atof(tokens[5].c_str()));
             }
             );

  ZhsTimeN = Ex.size();
  double E2max = -1;
  double Ex_max = 1e-99, Ey_max = 1e-99, Ez_max = 1e-99;
  int Emax_ind = -1;
  for (int i = 0; i < ZhsTimeN; i++){
    double E2 = Ex[i] * Ex[i] + Ey[i] * Ey[i] + Ez[i] * Ez[i];
    if (E2 > E2max) {
      E2max = E2;
      Emax_ind = i;
    }
  };
  assert(Emax_ind > -1);

  double Emax = sqrt(E2max);

  SET3_IND(E, x, y, z, _max, Emax_ind);


  // Do not buffer output:
  // setbuf(stdout, NULL);
  ZhsTimeE.resize(ZhsTimeN);
  ZhsAlpha.resize(ZhsTimeN);
  for (int i = 0; i < ZhsTimeN; i++){
    double E2 = Ex[i] * Ex[i] + Ey[i] * Ey[i] + Ez[i] * Ez[i];
    double E = sqrt(E2);
    double EEmax = (Ex[i] * Ex_max + Ey[i] * Ey_max + Ez[i] * Ez_max);
    double Eproj = EEmax / Emax;
    ZhsTimeE[i] = Eproj;
    double cosEEmax; // = abs(E) > 1e-99 ? Eproj / E : 999.0;
    double EEmaxAlpha = 190.0; // = cosEEmax != 999.0 ? acos(cosEEmax) * 180 / pi: 999;

    if (abs(E) > 0) {
      cosEEmax = Eproj / E;
      if (cosEEmax >= 1.0) {
        EEmaxAlpha = 0;
      }
      else if (cosEEmax <= -1.0) {
        EEmaxAlpha = 180;
      }
      else {
        EEmaxAlpha = std::acos(cosEEmax) * 180 / pi;
      }
    }
      
    ZhsAlpha[i] = EEmaxAlpha;
    // printf("i: %d/%d, Time: %11.8e, ZhsTimeE[i]: %11.4e, Eproj/E: %11.8e, EEmaxAlpha: %11.4e\n", i, ZhsTimeN - 1, ZhsTimeArr.at(i), ZhsTimeE.at(i), Eproj/E, EEmaxAlpha);
  }

  // Projections on the Emax:

  printf("The maximal field Emax: %11.4e\n", Emax);
  printf("The maximal index of Emax: %d\n", Emax_ind);
  printf("The maximal components: (%11.4e, %11.4e, %11.4e)\n", Ex_max, Ey_max, Ez_max);
  ZhsTimeStart = ZhsTimeArr[0]; 
  ZhsTimeDelta = ZhsTimeArr[1] - ZhsTimeArr[0];

  gStyle->SetOptTitle(0);

  WITH_LINES(
             "./GDB/cut_lo_hi/ft/icemc/__plot1d.dat",
             // "./GDB/ift/inter/__plot1d.dat",
             ind,
             tokens,
             NrFT[ind - 1] = atof(tokens[0].c_str());
             );
  for (int ind = 0; ind < ANITA_TIME_SAMPLES; ind++) {
    NrFT[ind] = NrFT[ind] * 5.5e+8; // The coefficient is an empirical value derived from eyeballing using gnuplot.
    //    // cout << ind << " NrFT: " << NrFT[ind] << endl;
  }
  // cout << "Before the call to hot_loop" << endl;

  struct cr_ft_state *cr_ft_result = (struct cr_ft_state *) hot_loop("/nfs/data_disks/herc0a/users/bugaev/ANITA/anitaBuildTool/components/icemc/cr-ft.so", false /* bInteractive */);
  // struct hot_test_state *hot_test_result = (struct hot_test_state *) hot_loop("/nfs/data_disks/herc0a/users/bugaev/ANITA/anitaBuildTool/components/icemc/hot-module-test.so", true /* bInteractive */);
  // exit(0);
  
  double vmmhz[Anita::NFREQ];                        //  V/m/MHz at balloon (after all steps)
  // given the angle you are off the Cerenkov cone,  the fraction of the observed e field that comes from the em shower
   
  Anita *anita1=new Anita();// right now this constructor gets banding info
  global_anita1 = anita1;

  cout << "testEAS reached the end of the development block" << endl;
  //  exit(0);

  settings1->SEED=settings1->SEED +run_no;
  cout <<"seed is " << settings1->SEED << endl;

  TRandom *rsave = gRandom;
  TRandom3 *Rand3 = new TRandom3(settings1->SEED);//for generating random numbers
  gRandom=Rand3;


  Balloon *bn1=new Balloon(); // instance of the balloon
  global_bn1 = bn1;
  Secondaries *sec1=new Secondaries();
  // Primaries *primary1=new Primaries();
  Signal *sig1=new Signal();
  Ray *ray1=new Ray(); // create new instance of the ray class
  Counting *count1=new Counting();
  GlobalTrigger *globaltrig1;
  // Taumodel *taus1 = new Taumodel();
  // input parameters

  stemp=string(outputdir.Data())+"/output"+run_num+".txt";
  ofstream foutput(stemp.c_str(),  ios::app);

  int NNU;
  double RANDOMISEPOL=0.;
  int inu=0;
  int neutrinos_passing_all_cuts=0;
  double weight=0; 
  int discones_passing=0;  // number of discones that pass
 
  settings1->ReadInputs(input.c_str(),  foutput, NNU, RANDOMISEPOL);
  settings1->ApplyInputs(anita1,  sec1,  sig1,  bn1,  ray1);
  
  settings1->SEED=settings1->SEED + run_no;
  gRandom->SetSeed(settings1->SEED);

  bn1->InitializeBalloon();
  anita1->Initialize(settings1, foutput, inu, outputdir);
  bvv::interpolate_zhs_ft(vmmhz, anita1->v_phases, anita1, cr_ft_result);
  // for (int i = 0; i < Anita::NFREQ; i++) {
  //   cout << "i: " << i << " vmmhz: " << vmmhz[i] << " phase: " << phase[i] << endl;
  // }

  if (trig_thresh!=0)
    anita1->powerthreshold[4]=trig_thresh;
  if (nnu_tmp!=0)
    NNU=nnu_tmp;
  if (exp_tmp!=0)
    settings1->EXPONENT=exp_tmp;

  // create new instance of the screen class
  // only using the specular case as we are ignoring the ice
  Screen *panel1 = new Screen(0);              
  
  time_t raw_start_time = time(NULL);
  struct tm * start_time = localtime(&raw_start_time);

  cout << "Date and time at start of run are: " << asctime (start_time) << "\n";


  Tools::Zero(anita1->NRX_PHI, Anita::NLAYERS_MAX);
  for (int i=0;i<Anita::NLAYERS_MAX;i++) {
    Tools::Zero(anita1->PHI_EACHLAYER[i], Anita::NPHI_MAX);
  }
  Tools::Zero(anita1->PHI_OFFSET, Anita::NLAYERS_MAX);
  Tools::Zero(anita1->THETA_ZENITH, Anita::NLAYERS_MAX);
  Tools::Zero(anita1->LAYER_VPOSITION, Anita::NLAYERS_MAX);
  Tools::Zero(anita1->RRX, Anita::NLAYERS_MAX);

  double volts_rx_rfcm_lab_e_all[48][512];
  double volts_rx_rfcm_lab_h_all[48][512];

 
  Vector n_eplane = const_z;
  Vector n_hplane = -const_y;
  Vector n_normal = const_x;

  // Antenna #30's horizontal polarization: 
  Vector n_pol = Vector(+0.0274856, +0.998238,  +0.0525939); // direction of polarization
  // Direction opposite to antenna #30's normal:
  Vector direction2bn = -1 * Vector(-0.918422,  +0.0459897, -0.392921); // direction from EAS to balloon
  Vector n_pol_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX]; // direction of polarization of signal seen at each antenna
  Vector direction2bn_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX]; // direction from EAS to balloon

  // variable declarations for functions GetEcompHcompEvector and GetEcompHcompkvector - oindree
  double e_component=0; // E comp along polarization
  double h_component=0; // H comp along polarization
  double n_component=0; // normal comp along polarization

  double e_component_kvector=0; // component of e-field along the rx e-plane
  double h_component_kvector=0; // component of the e-field along the rx h-plane
  double n_component_kvector=0; // component of the e-field along the normal


  double hitangle_e, hitangle_h;       // angle the ray hits the antenna wrt e-plane, h-plane
  double hitangle_e_all[Anita::NANTENNAS_MAX];         // hit angles rel. to e plane stored for each antenna
  double hitangle_h_all[Anita::NANTENNAS_MAX];         // hit angles rel. to h plane stored for each antenna
  
  double sourceLon;
  double sourceAlt;
  double sourceLat;


  int l3trig[Anita::NPOL];  // 16 bit number which says which phi sectors pass L3 V-POL
  // For each trigger layer,  which "clumps" pass L2.  16 bit,  16 bit and 8 bit for layers 1 & 2 and nadirs
  int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];
  //For each trigger layer,  which antennas pass L1.  16 bit,  16 bit and 8 bit and layers 1,  2 and nadirs
  int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];

  // these are declared here so that they can be stuck into trees
  int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //counting how many pass trigger requirement

  int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX]; //counting how many pass trigger requirement

  int nchannels_triggered /* __attribute__ ((unused)) */ = 0; // total number of channels triggered
  int nchannels_perrx_triggered[48] __attribute__ ((unused)); // total number of channels triggered


  Tools::Zero(count1->npass, 2); // sums events that pass,  without weights

#ifdef ANITA_UTIL_EXISTS

  string outputAnitaFile =string(outputdir.Data())+"/SimulatedAnitaEventFile"+run_num+".root";
  TFile *anitafileEvent = new TFile(outputAnitaFile.c_str(), "RECREATE");

  TTree *eventTree = new TTree("eventTree", "eventTree");
  eventTree->Branch("event",             &realEvPtr           );
  eventTree->Branch("run",               &run_no,   "run/I"   );
  eventTree->Branch("weight",            &weight,   "weight/D");

  outputAnitaFile =string(outputdir.Data())+"/SimulatedAnitaHeadFile"+run_num+".root";
  TFile *anitafileHead = new TFile(outputAnitaFile.c_str(), "RECREATE");

  TTree *headTree = new TTree("headTree", "headTree");
  headTree->Branch("header",  &rawHeaderPtr           );
  headTree->Branch("weight",  &weight,      "weight/D");

  outputAnitaFile =string(outputdir.Data())+"/SimulatedAnitaGpsFile"+run_num+".root";
  TFile *anitafileGps = new TFile(outputAnitaFile.c_str(), "RECREATE");

  UInt_t eventNumber;
  TTree *adu5PatTree = new TTree("adu5PatTree", "adu5PatTree");
  adu5PatTree->Branch("pat",          &Adu5PatPtr                   );
  adu5PatTree->Branch("eventNumber",  &eventNumber,  "eventNumber/I");
  adu5PatTree->Branch("weight",       &weight,       "weight/D"     );

  AnitaGeomTool *AnitaGeom1 = AnitaGeomTool::Instance();

#ifdef ANITA3_EVENTREADER
 
  // Set AnitaVersion so that the right payload geometry is used
  AnitaVersion::set(settings1->ANITAVERSION);

  outputAnitaFile =string(outputdir.Data())+"/SimulatedAnitaTruthFile"+run_num+".root";
  TFile *anitafileTruth = new TFile(outputAnitaFile.c_str(), "RECREATE");

  TString icemcgitversion = TString::Format("%s", EnvironmentVariable::ICEMC_VERSION(outputdir));  
  printf("ICEMC GIT Repository Version: %s\n", icemcgitversion.Data());
  unsigned int timenow = time(NULL);

  TTree *configAnitaTree = new TTree("configIcemcTree", "Config file and settings information");
  configAnitaTree->Branch("gitversion",   &icemcgitversion  );
  configAnitaTree->Branch("startTime",    &timenow          );
  // configAnitaTree->Branch("settings",  &settings1                    );
  configAnitaTree->Fill();

  
  TTree *truthAnitaTree = new TTree("truthAnitaTree", "Truth Anita Tree");
  truthAnitaTree->Branch("truth",     &truthEvPtr                   );
#endif

#endif
  //end ROOT variable definitions
  ///////////////////////////////////////////////////////////////////////


  IceModel *antarctica = new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10, settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0, settings1->WEIGHTABSORPTION);

  // fills arrays according to antenna specs
  anita1->GetBeamWidths(settings1); // this is used if GAINS set to 0
  // Antenna measured gain vs. frequency
  anita1->ReadGains(); // this is used if GAINS set to 1
  anita1->Set_gain_angle(settings1, sig1->NMEDIUM_RECEIVER);

  
  // sets position of balloon and related quantities
  // these are all passed as pointers
  // theta,  phi,  altitude of balloon
  // position of balloon,  altitude and position of surface of earth (relative to the center of the earth) under balloon
  bn1->SetDefaultBalloonPosition(antarctica);

  anita1->SetNoise(settings1, bn1, antarctica);
  //pathtree->Fill(); //Added by Stephen for verification of path

  // find the maximum distance the interaction could be from the balloon and still be within the horizon.
  antarctica->GetMAXHORIZON(bn1);

  // calculate the volume of ice seen by the balloon for all balloon positions
  antarctica->CreateHorizons(settings1, bn1, bn1->theta_bn, bn1->phi_bn, bn1->altitude_bn, foutput);
  cout << "Done with CreateHorizons.\n";

  // builds payload based on read inputs
  anita1->GetPayload(settings1,  bn1);

  if (settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME==5) {
    Vector plusz(0., 0., 1.);
    bn1->PickBalloonPosition(plusz, antarctica, settings1, anita1);
    anita1->calculate_all_offsets();
    double angle_theta=16.;
    double angle_phi=0.;
    Vector x = Vector(cos(angle_theta * RADDEG) * cos((angle_phi+11.25) * RADDEG),
                      cos(angle_theta * RADDEG) * sin((angle_phi+11.25) * RADDEG),
                      sin(angle_theta * RADDEG));
    anita1->GetArrivalTimes(x,bn1,settings1);
    cout << "end of getarrivaltimes\n";
  }

  time_t raw_loop_start_time = time(NULL);
  cout<<"Starting loop over events.  Time required for setup is "<<(int)((raw_loop_start_time - raw_start_time)/60)<<":"<< ((raw_loop_start_time - raw_start_time)%60)<<endl;


  //  TRandom r(settings1->SEED); // use seed set as input

  signal(SIGINT,  interrupt_signal_handler);     // This function call allows icemc to gracefully abort and write files as usual rather than stopping abruptly.

  int passes_thisevent __attribute__ ((unused)) =0;
  int passestrigger __attribute__ ((unused)) = 0;
  int count_total=0;
  int count_rx=0;
  int antNum=0;
  
  // begin looping over NNU neutrinos doing the things
  for (inu = 0; inu < NNU; inu++) {

    if (NNU >= 100) {
      if (inu % (NNU / 100) == 0)
        cout << inu << " neutrinos. " << (double(inu)/double(NNU)) * 100 << "% complete.\n";
    }
    else
      cout << inu << " neutrinos.  " << (double(inu) / double(NNU)) * 100 << "% complete.\n";
    

    eventNumber=(UInt_t)(run_no)*NNU+inu;
    
    // Set seed of all random number generators to be dependent on eventNumber
    gRandom->SetSeed(eventNumber+6e7);
    TRandom3 r(eventNumber+7e8);
    anita1->fRand->SetSeed(eventNumber+8e9);
    
    anita1->passglobtrig[0]=0;
    anita1->passglobtrig[1]=0;
    passes_thisevent=0;
    passestrigger=0;
    count_total++;
    // initializing the voltage seen by each polarization of each antenna
    bn1->dtryingposition=0;

    // for (int i=0; i<Anita::NFREQ;i++) {
    //   vmmhz[i] = 0.; // the full signal with all factors accounted for (1/r,  atten. etc.)
    //   // vmmhz_em[i]=0.; // for keeping track of just the em component of the shower
    // } //Zero the vmmhz array - helpful for banana plots,  shouldn't affect anything else - Stephen
    
    // Picks the balloon position and at the same time sets the masks and thresholds
    bn1->PickBalloonPosition(antarctica,  settings1,  inu,  anita1,  r.Rndm());
    
    // BR: Here calculate the direction to the balloon
    // And the polarization
    // direction2bn = something something something
    // n_pol        = something something something
    // if you decide to evaluate a different direction per antenna, then you should define also
    // direction2bn_eachboresight = something something something
    // n_pol_eachboresight        = something something something

    
    
    // make a global trigger object (but don't touch the electric fences)
    globaltrig1 = new GlobalTrigger(settings1, anita1);
    
    Tools::Zero(anita1->arrival_times, Anita::NLAYERS_MAX*Anita::NPHI_MAX);
    
    if(settings1->BORESIGHTS)
      anita1->GetArrivalTimesBoresights(direction2bn_eachboresight);
    else
      anita1->GetArrivalTimes(direction2bn,bn1,settings1);
    
    anita1->rx_minarrivaltime=Tools::WhichIsMin(anita1->arrival_times, settings1->NANTENNAS);
    
    
    globaltrig1->volts_rx_rfcm_trigger.assign(16,  vector <vector <double> >(3,  vector <double>(0)));
    anita1->rms_rfcm_e_single_event = 0;


    // std::map<std::string, void *> *hot_antenna_result = (std::map<std::string, void *> *) hot_loop("/nfs/data_disks/herc0a/users/bugaev/ANITA/anitaBuildTool/components/icemc/hot-antenna.so", false /* bInteractive */);

    count_rx=0;
    for (int ilayer=0; ilayer < settings1->NLAYERS; ilayer++) { // loop over layers on the payload
      for (int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) { // ifold loops over phi
	
	ChanTrigger *chantrig1 = new ChanTrigger();
	chantrig1->InitializeEachBand(anita1);

	//reset screen parameters (even for no roughness) for the new event
	panel1->ResetParameters();

	panel1->SetNvalidPoints(1);

	// BR: Add here the electric field at the payload before it gets through the antennas
	// For each antenna you need to define the electric field Vmmhz by doing
        for (int k=0;k<Anita::NFREQ;k++) {
          panel1->AddVmmhz_freq(vmmhz[k]);
        }
        panel1->AddDelay( 0. );
	// Use this to add the direction
	panel1->AddVec2bln(direction2bn);
	// Use this to add direction of polarization
	panel1->AddPol(n_pol);
        panel1->AddWeight( 1. );
        panel1->SetWeightNorm( 1. );

	// set the position of the source
	sourceLon=sourceLat=sourceAlt=0;
	
	bn1->GetAntennaOrientation(settings1,  anita1,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);

	// for this (hitangle_h_all[count_rx]=hitangle_h;) and histogram fill, use specular case
	//although the GetEcomp..() functions are called in ConvertInputWFtoAntennaWF() to calculate the actual waveforms
	if (!settings1->BORESIGHTS) {
	  bn1->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal,  direction2bn, e_component_kvector,  h_component_kvector,  n_component_kvector);
	  bn1->GetEcompHcompEvector(settings1,  n_eplane,  n_hplane,  n_pol,  e_component,  h_component,  n_component);
	}
	else{ // i.e. if BORESIGHTS is true
	  bn1->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal,  ray1->n_exit2bn_eachboresight[2][ilayer][ifold],  e_component_kvector,  h_component_kvector,  n_component_kvector);
	  bn1->GetEcompHcompEvector(settings1,  n_eplane,  n_hplane,  n_pol_eachboresight[ilayer][ifold], e_component,  h_component,  n_component);
	}
	bn1->GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector, hitangle_e, hitangle_h);
	// store hitangles for plotting
	hitangle_h_all[count_rx]=hitangle_h;
	hitangle_e_all[count_rx]=hitangle_e;
	// for debugging

	
	antNum = anita1->GetRxTriggerNumbering(ilayer, ifold);

	chantrig1->ApplyAntennaGain(settings1, anita1, bn1, panel1, antNum, n_eplane, n_hplane, n_normal);
	
	chantrig1->TriggerPath(settings1, anita1, antNum);
	
	chantrig1->DigitizerPath(settings1, anita1, antNum);
	
	chantrig1->TimeShiftAndSignalFluct(settings1, anita1, ilayer, ifold, volts_rx_rfcm_lab_e_all,  volts_rx_rfcm_lab_h_all);
	

	double thresholds[2][5];
	//+++++//+++++//+++++//+++++//+++++//+++++//+++++
	chantrig1->WhichBandsPass(settings1, anita1, globaltrig1, bn1, ilayer, ifold,  0, 0, 0, thresholds);
	
		
	delete chantrig1;

	count_rx++;
      } //loop through the phi-fold antennas
    }  //loop through the layers of antennas
    
    
    anita1->rms_rfcm_e_single_event = sqrt(anita1->rms_rfcm_e_single_event / (anita1->HALFNFOUR * settings1->NANTENNAS));

    for (int irx=0;irx<settings1->NANTENNAS;irx++) {
      nchannels_perrx_triggered[irx]=globaltrig1->nchannels_perrx_triggered[irx];
    }

    nchannels_triggered=Tools::iSum(globaltrig1->nchannels_perrx_triggered, settings1->NANTENNAS); // find total number of antennas that were triggered.
    cout << "nchannels_triggered: " << nchannels_triggered << endl;

    //////////////////////////////////////
    //       EVALUATE GLOBAL TRIGGER    //
    //          FOR VPOL AND HPOL       //
    //////////////////////////////////////


    int thispasses[Anita::NPOL]={0,0};

    globaltrig1->PassesTrigger(settings1, anita1, discones_passing, 2, l3trig, l2trig, l1trig, settings1->antennaclump, loctrig, loctrig_nadironly, inu,
			       thispasses);

    for (int i=0;i<2;i++) {
      for (int j=0;j<16;j++) {
	for (int k=0;k<anita1->HALFNFOUR;k++) {
	  count1->nl1triggers[i][0]+=anita1->l1trig_anita3and4_inanita[i][j][k];
	}
      }
    }

    ///////////////////////////////////////
    //       Require that it passes      //
    //            global trigger         //
    ///////////////////////////////////////
    // for Anita-lite,  Anita Hill, just L1 requirement on 2 antennas. This option is currently disabled
    // Save events that generate an RF trigger or that are part of the min bias sample
    // Minimum bias sample: save all events that we could see at the payload
    // Independentely from the fact that they generated an RF trigger

    if ( (thispasses[0]==1 && anita1->pol_allowed[0]==1)
	 || (thispasses[1]==1 && anita1->pol_allowed[1]==1)
	 || (settings1->MINBIAS==1)) {

      //	cout << inu << endl;

      anita1->passglobtrig[0]=thispasses[0];
      anita1->passglobtrig[1]=thispasses[1];


      // keep track of events passing trigger
      count1->npassestrigger[0]++;
      // tags this event as passing
      passestrigger=1;


#ifdef ANITA_UTIL_EXISTS
      realEvPtr   = new UsefulAnitaEvent();
      rawHeaderPtr  = new RawAnitaHeader();
      Adu5PatPtr  = new Adu5Pat();

      Adu5PatPtr->latitude= bn1->latitude;
      Adu5PatPtr->longitude=bn1->longitude;
      Adu5PatPtr->altitude=bn1->altitude;
      Adu5PatPtr->realTime=bn1->realTime_flightdata;
      Adu5PatPtr->heading = bn1->heading;
      Adu5PatPtr->pitch = bn1->pitch;
      Adu5PatPtr->roll = bn1->roll;
      Adu5PatPtr->run = run_no;

      int fNumPoints = 260;
      for (int i = 0; i < 96; i++){
	for (int j = 0; j < 260; j++){
	  realEvPtr->fVolts[i][j] = 0.;
	  realEvPtr->fTimes[i][j] = 0.;
	}
      }

      for (int iant = 0; iant < settings1->NANTENNAS; iant++){
	//int IceMCAnt = GetIceMCAntfromUsefulEventAnt(anita1,  AnitaGeom1,  iant);
	int IceMCAnt = GetIceMCAntfromUsefulEventAnt(settings1,  iant);
	int UsefulChanIndexH = AnitaGeom1->getChanIndexFromAntPol(iant,  AnitaPol::kHorizontal);
	int UsefulChanIndexV = AnitaGeom1->getChanIndexFromAntPol(iant,  AnitaPol::kVertical);
	realEvPtr->fNumPoints[UsefulChanIndexV] = fNumPoints;
	realEvPtr->fNumPoints[UsefulChanIndexH] = fNumPoints;
	realEvPtr->chanId[UsefulChanIndexV] = UsefulChanIndexV;
	realEvPtr->chanId[UsefulChanIndexH] = UsefulChanIndexH;

        cout << "Antenna " << iant << " ";
	for (int j = 0; j < fNumPoints; j++) {
	  // convert seconds to nanoseconds
	  realEvPtr->fTimes[UsefulChanIndexV][j] = j * anita1->TIMESTEP * 1.0E9;
	  realEvPtr->fTimes[UsefulChanIndexH][j] = j * anita1->TIMESTEP * 1.0E9;
	  // convert volts to millivolts
	  realEvPtr->fVolts[UsefulChanIndexH][j] =  volts_rx_rfcm_lab_h_all[IceMCAnt][j+128]*1000;
	  realEvPtr->fCapacitorNum[UsefulChanIndexH][j] = 0;
          cout << "H:" << realEvPtr->fVolts[UsefulChanIndexH][j] << " ";
	  realEvPtr->fVolts[UsefulChanIndexV][j] =  volts_rx_rfcm_lab_e_all[IceMCAnt][j+128]*1000;
	  realEvPtr->fCapacitorNum[UsefulChanIndexV][j] = 0;
          cout << "V:" << realEvPtr->fVolts[UsefulChanIndexV][j] << " ";
	}//end int j
        cout << endl;
      }// end int iant

      realEvPtr->eventNumber = eventNumber;

      rawHeaderPtr->eventNumber = eventNumber;
      rawHeaderPtr->surfSlipFlag = 0;
      rawHeaderPtr->errorFlag = 0;

      if (settings1->MINBIAS==1)
	rawHeaderPtr->trigType = 8; // soft-trigger
      else
	rawHeaderPtr->trigType = 1; // RF trigger

      rawHeaderPtr->run = run_no;
      // put the vpol only as a placeholder - these are only used in Anita-2 anyway
      rawHeaderPtr->upperL1TrigPattern = l1trig[0][0];
      rawHeaderPtr->lowerL1TrigPattern = l1trig[0][1];
      rawHeaderPtr->nadirL1TrigPattern = l1trig[0][2];

      rawHeaderPtr->upperL2TrigPattern = l2trig[0][0];
      rawHeaderPtr->lowerL2TrigPattern = l2trig[0][1];
      rawHeaderPtr->nadirL2TrigPattern = l2trig[0][2];

      if (settings1->WHICH<9){
	rawHeaderPtr->phiTrigMask  = (short) anita1->phiTrigMask;
	rawHeaderPtr->l3TrigPattern = (short) l3trig[0];
      }

      rawHeaderPtr->calibStatus = 31;
      rawHeaderPtr->realTime = bn1->realTime_flightdata;
      Adu5PatPtr->latitude= bn1->latitude;
      Adu5PatPtr->longitude=bn1->longitude;
      Adu5PatPtr->altitude=bn1->altitude;
      Adu5PatPtr->realTime=bn1->realTime_flightdata;
      Adu5PatPtr->heading = bn1->heading;
      Adu5PatPtr->pitch = bn1->pitch;
      Adu5PatPtr->roll = bn1->roll;
      Adu5PatPtr->run = run_no;

#ifdef ANITA3_EVENTREADER
      if (settings1->WHICH==9 || settings1->WHICH==10) {
	rawHeaderPtr->setTrigPattern((short) l3trig[0], AnitaPol::kVertical);
	rawHeaderPtr->setTrigPattern((short) l3trig[1], AnitaPol::kHorizontal);
	rawHeaderPtr->setMask( (short) anita1->l1TrigMask,  (short) anita1->phiTrigMask,  AnitaPol::kVertical);
	rawHeaderPtr->setMask( (short) anita1->l1TrigMaskH, (short) anita1->phiTrigMaskH, AnitaPol::kHorizontal);
      }

      truthEvPtr        = new TruthAnitaEvent();
      truthEvPtr->eventNumber      = eventNumber;
      truthEvPtr->realTime         = bn1->realTime_flightdata;
      truthEvPtr->run              = run_no;
      truthEvPtr->nuMom            = 0; //pnu;
      truthEvPtr->nu_pdg           = 0; //pdgcode;
      truthEvPtr->e_component      = e_component;
      truthEvPtr->h_component      = h_component;
      truthEvPtr->n_component      = n_component;
      truthEvPtr->e_component_k    = e_component_kvector;
      truthEvPtr->h_component_k    = h_component_kvector;
      truthEvPtr->n_component_k    = n_component_kvector;
      truthEvPtr->sourceLon        = sourceLon;
      truthEvPtr->sourceLat        = sourceLat;
      truthEvPtr->sourceAlt        = sourceAlt;
      truthEvPtr->weight           = weight;
      for (int i=0;i<3;i++){
	truthEvPtr->balloonPos[i]  = bn1->r_bn[i];
	truthEvPtr->balloonDir[i]  = bn1->n_bn[i];
	truthEvPtr->nuPos[i]       = 0; //interaction1->posnu[i];
	truthEvPtr->nuDir[i]       = 0; //interaction1->nnu[i];
      }
      for (int i=0;i<5;i++){
	for (int j=0;j<3;j++){
	  truthEvPtr->rfExitNor[i][j] = 0;//ray1->n_exit2bn[i][j];
	  truthEvPtr->rfExitPos[i][j] = 0;//ray1->rfexit[i][j];
	}
      }
      for (int i=0;i<48;i++){
	truthEvPtr->hitangle_e[i]  = hitangle_e_all[i];
	truthEvPtr->hitangle_h[i]  = hitangle_h_all[i];
      }
      if(settings1->ROUGHNESS){
	for (int i=0;i<Anita::NFREQ;i++)
	  truthEvPtr->vmmhz[i]       = panel1->GetVmmhz_freq(i);
      }
      truthAnitaTree->Fill();
      delete truthEvPtr;
#endif

      headTree->Fill();
      eventTree->Fill();
      adu5PatTree->Fill();

      delete realEvPtr;
      delete rawHeaderPtr;
      delete Adu5PatPtr;
#endif

      neutrinos_passing_all_cuts++;


      passes_thisevent=1; // flag this event as passing
      cout << "PRODUCED TRIGGER" << endl;
    } // end if passing global trigger conditions
    else {
      passes_thisevent=0; // flag this event as not passing
      cout << "Event produced no triggers." << endl;
    }// end else event does not pass trigger

    ///////////////////////////////////////
    //
    // WE GET HERE REGARDLESS OF WHETHER THE TRIGGER PASSES
    //
    /////////////

    delete globaltrig1;

    // keeping track of intermediate counters,  incrementing by weight1.
    // weight1 was not yet determined when integer counters were incremented.

  
    //looping over two types of rays - upgoing and downgoing.
    if (ABORT_EARLY){
      std::cout << "\n***********************************************************";
      std::cout << "\n* SIGINT received,  aborting loop over events early.";
      std::cout << "\n* Stopped after event " << inu << " instead of " << NNU;
      std::cout << "\n* Any output which relied on NNU should be corrected for.";
      std::cout << "\n***********************************************************\n";
      foutput << "\n***********************************************************";
      foutput << "\n* SIGINT received,  aborting loop over events early.";
      foutput << "\n* Stopped after event " << inu << " instead of " << NNU;
      foutput << "\n* Any output which relied on NNU should be corrected for.";
      foutput << "\n***********************************************************\n";
      break;
    }
  }//end NNU neutrino loop

  // Finished with individual neutrinos now ...
  //roughout.close();

  gRandom=rsave;
  delete Rand3;


#ifdef ANITA_UTIL_EXISTS

  anitafileEvent->cd();
  eventTree->Write("eventTree");
  anitafileEvent->Close();
  delete anitafileEvent;

  anitafileHead->cd();
  headTree->Write("headTree");
  anitafileHead->Close();
  delete anitafileHead;

  anitafileGps->cd();
  adu5PatTree->Write("adu5PatTree");
  anitafileGps->Close();
  delete anitafileGps;

#ifdef ANITA3_EVENTREADER
  anitafileTruth->cd();
  truthAnitaTree->Write("truthAnitaTree");
  anitafileTruth->Close();
  delete anitafileTruth;
#endif

#endif



  delete anita1;
  return 0;

} //END MAIN PROGRAM



void interrupt_signal_handler(int sig){
  signal(sig,  SIG_IGN);
  ABORT_EARLY = true;
  return;
}
//end interrupt_signal_handler()

#ifdef ANITA_UTIL_EXISTS
//int GetIceMCAntfromUsefulEventAnt(Anita *anita1,  AnitaGeomTool *AnitaGeom1,  int UsefulEventAnt){
int GetIceMCAntfromUsefulEventAnt(Settings *settings1,  int UsefulEventAnt){
  //int layer_temp = IceMCLayerPosition[UsefulEventIndex][0];
  //int position_temp = IceMCLayerPosition[UsefulEventIndex][1];
  //int IceMCIndex = anita1->GetRx(layer_temp,  position_temp);
  int IceMCAnt = UsefulEventAnt;
  if ((settings1->WHICH==9 || settings1->WHICH==10) && UsefulEventAnt<16) {
    IceMCAnt = (UsefulEventAnt%2==0)*UsefulEventAnt/2 + (UsefulEventAnt%2==1)*(UsefulEventAnt/2+8);
  }

  return IceMCAnt;
}
//end GetIceMCAntfromUsefulEventAnt()


#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
