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
#include "signal.h"

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
#include "SimulatedSignal.h"
#include "EnvironmentVariable.h"

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

Taumodel* TauPtr = NULL;

const string ICEMC_SRC_DIR = EnvironmentVariable::ICEMC_SRC_DIR();

ClassImp(RX);

using namespace std;

class EarthModel;
class Position;

/************MOVED FROM shared.hh and shared.cc*****************/
// These need to be moved elsewhere.
// They were all global variables inside of shared.hh.
// Most of these should not be global variables.
const int NVIEWANGLE=100; // number of viewing angles to look at the signal,  on Seckel's request
// int irays; // counts rays (for roughness) // LC: commented out because not used

int inu; // counts neutrinos as they are generated
UInt_t eventNumber;
double eventsfound_beforetrigger=0.;
double eventsfound_crust=0; //number of events that only traverse the crust
double eventsfound_weightgt01=0; // summing weights > 0.1
double eventsfound_belowhorizon=0; // how many are below horizon
double eventsfound=0;   // how many events found
double eventsfound_prob=0;   // how many events found,  probabilities
double sum[3]; // sum of weight for events found for 3 flavors
// These numbers are from Feldman and Cousins wonderful paper:  physics/9711021
double poissonerror_minus[21] = {0.-0.00, 1.-0.37, 2.-0.74, 3.-1.10, 4.-2.34, 5.-2.75, 6.-3.82, 7.-4.25, 8.-5.30, 9.-6.33, 10.-6.78, 11.-7.81, 12.-8.83, 13.-9.28, 14.-10.30, 15.-11.32, 16.-12.33, 17.-12.79, 18.-13.81, 19.-14.82, 20.-15.83};
double poissonerror_plus[21] = {1.29-0., 2.75-1., 4.25-2., 5.30-3., 6.78-4., 7.81-5., 9.28-6., 10.30-7., 11.32-8., 12.79-9., 13.81-10., 14.82-11., 16.29-12., 17.30-13., 18.32-14., 19.32-15., 20.80-16., 21.81-17., 22.82-18., 23.82-19., 25.30-20.};
const int NBINS=10; // keep track of the number of events found,  binned
// by weights
double MIN_LOGWEIGHT=-3;
double MAX_LOGWEIGHT=-1;
int index_weights=0; // which bin the weight falls into
double eventsfound_binned[NBINS];
double eventsfound_binned_e[NBINS];
double eventsfound_binned_mu[NBINS];
double eventsfound_binned_tau[NBINS];
double km3sr = 0;                      // total km3sr
double km3sr_e=0;                       // to calculate km3sr for electrons
double km3sr_mu=0;                       // to calculate km3sr for muons
double km3sr_tau=0;                       // to calculate km3sr for taus
double error_plus=0;          // keeping track of asymmetric error bars
double error_e_plus=0;
double error_mu_plus=0;
double error_tau_plus=0;
double error_minus=0;
double error_e_minus=0;
double error_mu_minus=0;
double error_tau_minus=0;
int ierr=0; // integer for returning error from functions.
double gain_dipole=2.15;  // antenna gain (nominal value for dipole is 2.15)
double changle_deg=0; // same,  in degrees.
// inputs
int NNU;        // number of neutrinos
int whichray=0; // indexes the rays that we look at (actually just used for ice but we use it in GetNuFlavor so keep it here)
double RANDOMISEPOL=0.;
/************MOVED FROM shared.hh and shared.cc*****************/


double volume_thishorizon; // for plotting volume within the horizon of the balloon
int realtime_this;  // for plotting real unix time
double longitude_this; // for plotting longitude
double latitude_this; // for plotting latitude
double altitude_this; // for plotting altitude
double heading_this=0.;// for plotting heading
double gps_offset=0;

// inputs

double pnu=pow(10., 20);   //!< energy of neutrinos

double MEANX=0;
double MEANY=0.;

double SIGNALRADIUS=2.; // in degrees


double bwslice_vnoise_thislayer[4];// for filling tree6b,  noise for each bandwidth on each layer
int passes_thisevent=0; // this event passes
int unmasked_thisevent=0; // this event is unmasked

int discones_passing;  // number of discones that pass
int NDISCONES=8;
double heff_discone=0; // effective height of a discone antenna
double polarfactor_discone=0.;// factor to account for antenna beam pattern.
double thislambda=0;// for finding wavelength at each frequency
double volts_discone=0.;// for finding voltage at each discone
double vnoise_discone=0.; // noise on each discone


double BW_DISCONES=300.E6-120.E6; // bandwidth of the discones

// ray tracing
double fresnel1=0;  //!< net fresnel factor on field at ice-firn interface
double fresnel1_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX];  // for slac simulation
double fresnel2=0;  //!< net fresnel factor on field at firn-air
double mag1=0;  // magnification factor on field at ice-firn interface
double mag1_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX];// for slac simulation
double mag2=0;  // magnification factor on field in firn-air
double t_coeff_pokey, t_coeff_slappy;
double rflength=0;  // distance from interaction to ice-air exit point

double e_comp_max1=0;
double h_comp_max1=0;
double e_comp_max2=0;
double h_comp_max2=0;
double e_comp_max3=0;
double h_comp_max3=0;

double diffexit=0;  // checking exit point between_MAX interations
double diff_3tries=0;
double diffnorm=0;  // checking angle of surf normal between iterations
double diffrefr=0;  // checking angle of refr between iterations
double costheta_inc=0;  // cos angle of incidence wrt surface normal
double costheta_exit=0; // theta of exit point wrt earth (costheta=1 at south pole)
double theta_rf_atbn; // polar angle of the signal as seen by perfect eyes at the balloon.
double theta_rf_atbn_measured; //polar angle of the signal as measured at the balloon (just the above variable smeared by 0.5 degrees)

double costhetanu=-1000; // costheta of neutrino direction wrt earth (costheta=1 at south pole)

// neutrino path
double theta_in=0; // theta where neutrino enters earth (radians, south pole=0)
double lat_in=0; // latitude where neutrino enters earth (degrees, south pole=-90)
double nearthlayers=0; // how many layers (core, mantle, crust) does nnu traverse
double weight_prob=0.;//event weight,  including probability it interacts somewhere in ice along its path
double weight1=0; // event weight,  just based on absorption in earth,  see note
double weight=0.; // total event weight (product of the previous 2)
double logweight=0.;// log of the previous number
double len_int=0;// interaction length in m
double pieceofkm2sr=0; // Use this for making plots comparing different cross sections.  The integral of a plot from a run will be the total Area*sr of the detector.  That way it is proportional to the cross section and the integral is something meaningful to people.

double CUTONWEIGHTS=1.E-10; // cut out events with small enough weight that they don't matter,  to save time

// counting variables
int count_inthisloop1=0; // for debugging
int count_inthisloop2=0;
int count_inthisloop3=0;
double averaging_thetas1=0; // for debugging
double averaging_thetas2=0;
double averaging_thetas3=0;
int count_total=0; // number of neutrinos looped over.  should always equal inu

// i.e.,  is theoretically visible by balloon

int count_asktrigger=0;
int count_asktrigger_nfb=0;
int count_pass=0;  // how many total trigger channels pass (4 bandwidth slices*2 pol * nrx)

double count_passestrigger_w=0; // same as above,  but sum weights
int passestrigger=0; // 1=this event passes trigger, 0=does not
// so that we know to increment count_passestrigger_w at the end of the event
int allcuts[2]={0, 0}; // index is which ray (upward or downward)
// 1=this ray for this event passes all cuts,  0=does not
double allcuts_weighted[2]={0, 0}; // same as above but weighted
double allcuts_weighted_polarization[3]={0, 0, 0}; // same as above but divided into [0] vpol, [1] hpol, [2] both vpol and hpol


//signal has a chance to pass after accounting for 1/r
int count_chanceofsurviving=0; // based on neutrino direction,  has a chance of making it through earth.

int count_chanceinhell0=0; // based on 1/r,  and best case attenuation
// (after 2nd guess at rf exit point)
// signal has a chance of passing

double count_chanceinhell2_w=0; //same as above,  but sum weights
int chanceinhell2=0; // 1=this event has chance to pass, 0=does not-
// so that we know to increment count_chanceinhell2_w by weight1 at the end.

int count_chordgoodlength=0; // Incremented if neutrino path through earth is more than 1m
int count_d2goodlength=0; // neutrino path through ice is more than 1m
int count_rx=0; // counting antennas we loop through them

double sum_frac[3]; // fraction of passing events that are e, mu, tau adding weights
double sum_frac_db[3]; // same for double bangs

const int NBINS_DISTANCE=28; // keep track of number that pass as a function of distance to make Peter's plot
double eventsfound_binned_distance[NBINS_DISTANCE] = {0.};  // binning cumulative events found vs. distance
int index_distance=0; // index for filling array above
double km3sr_distance[NBINS_DISTANCE] = {0.}; // result of conversion of events to sensitivity
double error_distance_plus[NBINS_DISTANCE] = {0.}; // errors on above
double error_distance_minus[NBINS_DISTANCE] = {0.};
int eventsfound_binned_distance_forerror[NBINS_DISTANCE][NBINS] = {{0}}; // for propagation of errors

//taus
double km3sr_db = 0;
double km3sr_nfb=0;
double ptau=0;
int count_passestrigger_nfb=0;
double percent_increase_db=0;
double percent_increase_nfb=0;
double percent_increase_total=0;
double error_nfb=0;
double error_km3sr_nfb=0;
double error_percent_increase_nfb=0;

Vector n_exit2bn_db[5];
Vector nrf_iceside_db[5];  // direction of rf [tries][3d]
double n_exit_phi;  //phi angle of the ray from the surface to the balloon
int count_dbexitsice=0;

// int count_interacts_in_earth=0;
double eventsfound_nfb_binned[NBINS]; // counting events without first bang

// rf parameters
double heff_max=0.62639; // maximum value of the effective height based on antenna specs

// event geometry
double scalefactor_distance=0; // 1/r scalefactor
double scalefactor_attenuation=0; //scalefactor due to attenuation in ice
double MAX_ATTENLENGTH=1671;
double maxtaper=0; // this is just for plotting - maximum you are ever off cerenkov cone while
//an event is detectable
double dviewangle_deg=0; //!< deviation from the cherenkov angle

double forseckel[NVIEWANGLE][Anita::NFREQ];// Per Seckel's request,  get strength of signal across frequencies for different viewing angles.
double viewangles[NVIEWANGLE];
double GetAirDistance(double altitude_bn,  double beta); // given beta=angle wrt horizontal that the ray hits the balloon,  calculate distance that the ray traveled in air,  including curvature of earth

//Input files
int err=0; // errors from GetDirection function

double djunk; // junk variable

//For verification plots - added by Stephen
int max_antenna0 = -1;  //antenna with the peak voltage,  top layer
int max_antenna1 = -1;  //antenna with the peak voltage,  middle layer
int max_antenna2 = -1;  //antenna with the peak voltage,  bottom layer
double max_antenna_volts0 = 0; //Voltage on the antenna with maximum signal,  top layer
double max_antenna_volts0_em = 0; //Component of voltage from em shower on the antenna with maximum signal,  top layer

double max_antenna_volts1 = 0; //Voltage on the antenna with maximum signal,  middle layer
double max_antenna_volts2 = 0; //Voltage on the antenna with maximum signal,  bottom layer

double rx0_signal_eachband[2][5];
double rx0_threshold_eachband[2][5];
double rx0_noise_eachband[2][5];
int rx0_passes_eachband[2][5];

double voltagearray[Anita::NLAYERS_MAX*Anita::NPHI_MAX]; //Records max voltages on each antenna for one neutrino

Vector ant_max_normal0; //Vector normal to the face of the antenna with the maximum signal for a single neutrino,  top layer
Vector ant_max_normal1; //Vector normal to the face of the antenna with the maximum signal for a single neutrino,  middle layer
Vector ant_max_normal2; //Vector normal to the face of the antenna with the maximum signal for a single neutrino,  bottom layer
double vmmhz1m_visible = 0; //Actual V/m/Mhz at 1m
int freq_bins = Anita::NFREQ; //Because the compiler objected to using the const directly
double total_kgm2 = 0; // output of Getchord
double nnu_array[3];
double r_in_array[3];
double nsurf_rfexit_array[3];
double nsurf_rfexit_db_array[3];
double r_bn_array[3];
double n_bn_array[3];
double posnu_array[3];
double nrf_iceside_array[5][3];
double nrf_iceside_db_array[5][3];
double ant_max_normal0_array[3];
double ant_max_normal1_array[3];
double ant_max_normal2_array[3];
double n_pol_array[3];
double n_exit2bn_array[5][3];
double r_enterice_array[3];
double n_exit2bn_db_array[5][3];
double rfexit_array[5][3];
double rfexit_db_array[5][3];

int times_crust_entered_det=0;  //Counter for total times each Earth layer is entered for detected neutrinos only
int times_mantle_entered_det=0;
int times_core_entered_det=0;
int crust_entered=0;
int mantle_entered=0;
int core_entered=0;
int neutrinos_passing_all_cuts=0;
double sum_weights=0;
//End verification plot block



double justNoise_trig[2][48][512];
double justSignal_trig[2][48][512];
double justNoise_dig[2][48][512];
double justSignal_dig[2][48][512];

// functions

// set up array of viewing angles for making plots for seckel
void SetupViewangles(Signal *sig1);

void GetAir(double *col1); // get air column as a function of theta- only important for black hole studies
double GetThisAirColumn(Settings*,  Position r_in,  Vector nnu, Position posnu,  double *col1,  double& cosalpha, double& mytheta,  double& cosbeta0, double& mybeta);

double ScaleVmMHz(double vmmhz1m_max,  const Position &posnu,  const Position &r_bn);


double IsItDoubleBang(double exitlength,  double plepton);

int WhereIsSecondBang(const Position& posnu,  const Vector& nnu,  double nuexitlength,  double pnu,  IceModel *antarctica1,  const Position& r_bn, Position &posnu2,  Position &rfexit_db,  Vector &n_exit2bn_db);
void GetCurrent(string& current);


double GetAverageVoltageFromAntennasHit(Settings *settings1,  int *nchannels_perrx_triggered,  double *voltagearray,  double& volts_rx_sum);


Vector GetPolarization(const Vector &nnu,  const Vector &nrf2_iceside);

void Attenuate(IceModel *antartica1, Settings *settings1,  double& vmmhz_max,  double rflength,  const Position &posnu);

void Attenuate_down(IceModel *antarctica1,  Settings *settings1,  double& vmmhz_max,  const Position &rfexit2,  const Position &posnu,  const Position &posnu_down);

void IsAbsorbed(double chord_kgm2,  double len_int_kgm2,  double& weight);


int GetRayIceSide(const Vector &n_exit2rx,  const Vector &nsurf_rfexit,  double nexit,  double nenter,  Vector &nrf2_iceside);

double GetViewAngle(const Vector &nrf2_iceside,  const Vector &nnu);
int TIR(const Vector &n_surf,  const Vector &nrf2_iceside,  double N_IN,  double N_OUT);

void IntegrateBands(Anita *anita1,  int k,  Screen *panel1,  double *freq,  double scalefactor,  double *sumsignal);

void Integrate(Anita *anita1,  int j,  int k,  double *vmmhz,  double *freq,  double scalefactor,  double sumsignal);

void interrupt_signal_handler(int);  // This catches the Control-C interrupt,  SIGINT

bool ABORT_EARLY = false;    // This flag is set to true when interrupt_signal_handler() is called

void WriteNeutrinoInfo(Position&,  Vector&,  Position&,  double,  string,  string,  double,  ofstream &nu_out);

void CloseTFile(TFile *hfile);

int Getmine(double*,  double*,  double*,  double*);

void Getearth(double*,  double*,  double*,  double*);


#ifdef ANITA_UTIL_EXISTS
int GetIceMCAntfromUsefulEventAnt(Settings *settings1,  int UsefulEventAnt);
#ifdef R_EARTH
#undef R_EARTH
#endif
#endif


double thresholdsAnt[48][2][5];
double thresholdsAntPass[48][2][5];


//do a threshold scan
double threshold_start=-1.;
double threshold_end=-6.;
const int NTHRESHOLDS=20;
double threshold_step=(threshold_end-threshold_start)/(double)NTHRESHOLDS;

double npass_v_thresh[NTHRESHOLDS]={0.};
double denom_v_thresh[NTHRESHOLDS]={0.};
double npass_h_thresh[NTHRESHOLDS]={0.};
double denom_h_thresh[NTHRESHOLDS]={0.};
double thresholds[NTHRESHOLDS];

int main(int argc,  char **argv) {
  //--------------------------------------------------------------
  //  MC Anita
  //
  // 12/01/03
  //
  //--------------------------------------------------------------


  // // for comparing with peter
  // double sumsignal[5]={0.};

  string stemp;

  Settings* settings1 = new Settings();

  string input= ICEMC_SRC_DIR + "/inputs.conf";
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
    while ((clswitch = getopt(argc, argv, "t:i:o:r:n:e:")) != EOF) {
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
        cout << "Changed output directory to: " << string(outputdir.Data()) << endl;
        stemp="mkdir -p " + string(outputdir.Data());
        system(stemp.c_str());
        break;
      case 'r':
        run_num=optarg;
        stringstream convert(run_num);
        convert>>run_no;
        break;
      } // end switch
    } // end while
  } // end if arg>1


  settings1->SEED=settings1->SEED +run_no;
  cout <<"seed is " << settings1->SEED << endl;

  TRandom *rsave = gRandom;
  TRandom3 *Rand3 = new TRandom3(settings1->SEED);//for generating random numbers
  gRandom=Rand3;

  stemp=string(outputdir.Data())+"/nu_position"+run_num+".txt";
  ofstream nu_out(stemp.c_str(),  ios::app); //Positions,  direction of momentum,  and neutrino type for Ryan.

  stemp=string(outputdir.Data())+"/veff"+run_num+".txt";
  ofstream veff_out(stemp.c_str(),  ios::app);//to output only energy and effective volume to veff.txt

  stemp=string(outputdir.Data())+"/distance"+run_num+".txt";
  ofstream distanceout(stemp.c_str());

  stemp=string(outputdir.Data())+"/debug"+run_num+".txt";
  fstream outfile(stemp.c_str(), ios::out);

  stemp=string(outputdir.Data())+"/forbrian"+run_num+".txt";
  fstream forbrian(stemp.c_str(), ios::out);

  stemp=string(outputdir.Data())+"/al_voltages_direct"+run_num+".dat";
  fstream al_voltages_direct(stemp.c_str(), ios::out); //added djg ------provide anita-lite voltages and noise from MC for anita-lite analysis

  stemp=string(outputdir.Data())+"/events"+run_num+".txt";
  ofstream eventsthatpassfile(stemp.c_str());

  stemp=string(outputdir.Data())+"/numbers"+run_num+".txt";
  ofstream fnumbers(stemp.c_str()); // debugging

  stemp=string(outputdir.Data())+"/output"+run_num+".txt";
  ofstream foutput(stemp.c_str(),  ios::app);

  stemp=string(outputdir.Data())+"/slacviewangles"+run_num+".dat";
  ofstream fslac_viewangles(stemp.c_str()); // this outputs numbers that we need for analyzing slac data

  stemp=string(outputdir.Data())+"/slac_hitangles"+run_num+".dat";
  ofstream fslac_hitangles(stemp.c_str()); // this outputs numbers that we need for analyzing slac data

  Balloon *bn1=new Balloon(); // instance of the balloon
  Anita *anita1=new Anita();// right now this constructor gets banding info
  Secondaries *sec1=new Secondaries();
  Signal *sig1=new Signal();
  Ray *ray1=new Ray(); // create new instance of the ray class
  Counting *count1=new Counting();

  // input parameters
  settings1->ReadInputs(input.c_str(),  foutput, NNU, RANDOMISEPOL);
  settings1->ApplyInputs(anita1,  sec1,  sig1,  bn1,  ray1);
  sig1->Initialize();

  settings1->SEED=settings1->SEED + run_no;
  gRandom->SetSeed(settings1->SEED);

  bn1->InitializeBalloon();
  anita1->Initialize(settings1, foutput, inu, outputdir);

  if (nnu_tmp!=0)
    NNU=nnu_tmp;
  if (exp_tmp!=0)
    settings1->EXPONENT=exp_tmp;

  settings1->SIGNAL_FLUCT=1;
  settings1->ZEROSIGNAL=1;
  
  time_t raw_start_time = time(NULL);
  struct tm * start_time = localtime(&raw_start_time);

  cout << "Date and time at start of run are: " << asctime (start_time) << "\n";

  // zeroing global variables.
  Tools::Zero(sum_frac, 3);
  Tools::Zero(sum_frac_db, 3);
  Tools::Zero(anita1->NRX_PHI, Anita::NLAYERS_MAX);
  for (int i=0;i<Anita::NLAYERS_MAX;i++) {
    Tools::Zero(anita1->PHI_EACHLAYER[i], Anita::NPHI_MAX);
  }
  Tools::Zero(anita1->PHI_OFFSET, Anita::NLAYERS_MAX);
  Tools::Zero(anita1->THETA_ZENITH, Anita::NLAYERS_MAX);
  Tools::Zero(anita1->LAYER_VPOSITION, Anita::NLAYERS_MAX);
  Tools::Zero(anita1->RRX, Anita::NLAYERS_MAX);


  // frequency binning
  double vmmhz[Anita::NFREQ];                        //  V/m/MHz at balloon (after all steps)

  // given the angle you are off the Cerenkov cone,  the fraction of the observed e field that comes from the em shower
  double vmmhz_em[Anita::NFREQ];
  double vmmhz_min_thatpasses=1000;
 

  string taudecay;                   // tau decay type: e, m, h

  double volts_rx_rfcm_lab_e_all[48][512];
  double volts_rx_rfcm_lab_h_all[48][512];

  // variable declarations for functions GetEcompHcompEvector and GetEcompHcompkvector - oindree
  double e_component=0; // E comp along polarization
  double h_component=0; // H comp along polarization
  double n_component=0; // normal comp along polarization

  double e_component_kvector=0; // component of e-field along the rx e-plane
  double h_component_kvector=0; // component of the e-field along the rx h-plane
  double n_component_kvector=0; // component of the e-field along the normal


  // Vector n_eplane = const_z;
  // Vector n_hplane = -const_y;
  // Vector n_normal = const_x;

  Vector ant_normal; //Vector normal to the face of the antenna

  // double hitangle_e, hitangle_h;       // angle the ray hits the antenna wrt e-plane, h-plane
  double hitangle_e_all[Anita::NANTENNAS_MAX];         // hit angles rel. to e plane stored for each antenna
  double hitangle_h_all[Anita::NANTENNAS_MAX];         // hit angles rel. to h plane stored for each antenna

  double eventsfound_db=0; // same,  for double bang
  double eventsfound_nfb=0; // for taus


  double sourceLon=-999;
  double sourceAlt=-999;
  double sourceLat=-999;

  Vector n_nutraject_ontheground; //direction of the neutrino from the person standing on the ground just below the balloon.
  Vector n_pol; // direction of polarization
  // Vector n_pol_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX]; // direction of polarization of signal seen at each antenna
  Vector n_pol_db; // same,  double bangs

  int l3trig[Anita::NPOL];  // 16 bit number which says which phi sectors pass L3 V-POL
  // For each trigger layer,  which "clumps" pass L2.  16 bit,  16 bit and 8 bit for layers 1 & 2 and nadirs
  int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];
  //For each trigger layer,  which antennas pass L1.  16 bit,  16 bit and 8 bit and layers 1,  2 and nadirs
  int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];



  eventsfound=0.; // sums weights for events that pass

  Tools::Zero(count1->npass, 2); // sums events that pass,  without weights
  Tools::Zero(sum, 3);
  eventsfound_db=0;
  eventsfound_nfb=0;

  Tools::Zero(eventsfound_binned, NBINS);
  Tools::Zero(eventsfound_binned_e, NBINS);
  Tools::Zero(eventsfound_binned_mu, NBINS);
  Tools::Zero(eventsfound_binned_tau, NBINS);
  Tools::Zero(eventsfound_nfb_binned, NBINS);

  
  Tools::Zero(anita1->arrival_times[0], Anita::NLAYERS_MAX*Anita::NPHI_MAX);
  Tools::Zero(anita1->arrival_times[1], Anita::NLAYERS_MAX*Anita::NPHI_MAX);

  //we pick both the interaction point and its corresponding mirror point


  int pdgcode=-999;

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

  TTree *adu5PatTree = new TTree("adu5PatTree", "adu5PatTree");
  adu5PatTree->Branch("pat",          &Adu5PatPtr                   );
  adu5PatTree->Branch("eventNumber",  &eventNumber,  "eventNumber/I");
  adu5PatTree->Branch("weight",       &weight,       "weight/D"     );

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

  AnitaGeomTool *AnitaGeom1 = AnitaGeomTool::Instance();
  
#endif

  
  IceModel *antarctica = new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10, settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0, settings1->WEIGHTABSORPTION);
  cout << "area of the earth's surface covered by antarctic ice is " << antarctica->ice_area << "\n";

  
  // sets position of balloon and related quantities
  // these are all passed as pointers
  // theta,  phi,  altitude of balloon
  // position of balloon,  altitude and position of surface of earth (relative to the center of the earth) under balloon
  bn1->SetDefaultBalloonPosition(antarctica);

  
  // builds payload based on read inputs
  anita1->GetPayload(settings1,  bn1);

  // if (settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME==5) {
  Vector plusz(0., 0., 1.);
  bn1->PickBalloonPosition(plusz, antarctica, settings1, anita1);
  
  time_t raw_loop_start_time = time(NULL);
  cout<<"Starting loop over events.  Time required for setup is "<<(int)((raw_loop_start_time - raw_start_time)/60)<<":"<< ((raw_loop_start_time - raw_start_time)%60)<<endl;

 
  // for averaging balloon altitude and distance from center of earth
  // for comparing with Peter
  double average_altitude=0.;
  double average_rbn=0.;

  //  TRandom r(settings1->SEED); // use seed set as input

  signal(SIGINT,  interrupt_signal_handler);     // This function call allows icemc to gracefully abort and write files as usual rather than stopping abruptly.


  int antNum;

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
    unmasked_thisevent=1;
    vmmhz_min_thatpasses=1000; // initializing.  want to find the minumum voltage that passes a

   
    count_pass=0;
    passestrigger=0;
    chanceinhell2=0;
    // sec1->secondbang=false;
    count_total++;
    // initializing the voltage seen by each polarization of each antenna
    bn1->dtryingposition=-999;
    for (int i=0; i<Anita::NFREQ;i++) {
      vmmhz[i] = 0.; // the full signal with all factors accounted for (1/r,  atten. etc.)
      vmmhz_em[i]=0.; // for keeping track of just the em component of the shower
    } //Zero the vmmhz array - helpful for banana plots,  shouldn't affect anything else - Stephen

    
    // Picks the balloon position and at the same time sets the masks and thresholds
    bn1->PickBalloonPosition(antarctica,  settings1,  inu,  anita1,  r.Rndm());
      
    // find average balloon altitude and distance from center of earth for
    // making comparisons with Peter
    average_altitude+=bn1->altitude_bn/(double)NNU;
    average_rbn+=bn1->r_bn.Mag()/(double)NNU;

    realtime_this=bn1->realTime_flightdata;
    longitude_this=bn1->longitude;
    latitude_this=bn1->latitude;
    altitude_this=bn1->altitude;
    heading_this=bn1->heading;


    for (int ilayer=0; ilayer < settings1->NLAYERS; ilayer++) { // loop over layers on the payload
      for (int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) { // ifold loops over phi
          
	ChanTrigger *chantrig1 = new ChanTrigger();
	chantrig1->InitializeEachBand(anita1);

	antNum = anita1->GetRxTriggerNumbering(ilayer, ifold);
	  

	
#ifdef ANITA_UTIL_EXISTS
	if (settings1->SIGNAL_FLUCT && (settings1->NOISEFROMFLIGHTDIGITIZER || settings1->NOISEFROMFLIGHTTRIGGER) )
	  chantrig1->getNoiseFromFlight(settings1, anita1, antNum);
#endif
	//	  chantrig1->ApplyAntennaGain(settings1, anita1, bn1, panel1, antNum, n_eplane, n_hplane, n_normal);
	
	// chantrig1->injectImpulseAfterAntenna(anita1, antNum);

	
	memset(chantrig1->volts_rx_forfft, 0, sizeof(chantrig1->volts_rx_forfft));

	chantrig1->TriggerPath(settings1, anita1, antNum, bn1);
	
	chantrig1->DigitizerPath(settings1, anita1, antNum, bn1);
	
	chantrig1->TimeShiftAndSignalFluct(settings1, anita1, ilayer, ifold, volts_rx_rfcm_lab_e_all,  volts_rx_rfcm_lab_h_all);

	chantrig1->saveTriggerWaveforms(anita1, justSignal_trig[0][antNum], justSignal_trig[1][antNum], justNoise_trig[0][antNum], justNoise_trig[1][antNum]);
	chantrig1->saveDigitizerWaveforms(anita1, justSignal_dig[0][antNum], justSignal_dig[1][antNum], justNoise_dig[0][antNum], justNoise_dig[1][antNum]);
 
	// chantrig1->WhichBandsPass(settings1, anita1, globaltrig1, bn1, ilayer, ifold,  viewangle-sig1->changle, emfrac, hadfrac, thresholdsAnt[antNum]);

	//	cout << inu << " " << ilayer << " " << ifold << endl;
	delete chantrig1;
	       
      } //loop through the phi-fold antennas
    }  //loop through the layers of antennas


    if ( true ) {

	      
#ifdef ANITA_UTIL_EXISTS
      realEvPtr   = new UsefulAnitaEvent();
      rawHeaderPtr  = new RawAnitaHeader();
      Adu5PatPtr  = new Adu5Pat();

      Adu5PatPtr->latitude  = bn1->latitude;
      Adu5PatPtr->longitude = bn1->longitude;
      Adu5PatPtr->altitude  = bn1->altitude;
      Adu5PatPtr->realTime  = bn1->realTime_flightdata;
      Adu5PatPtr->heading   = bn1->heading;
      Adu5PatPtr->pitch     = bn1->pitch;
      Adu5PatPtr->roll      = bn1->roll;
      Adu5PatPtr->run       = run_no;


      memset(realEvPtr->fNumPoints, 0, sizeof(realEvPtr->fNumPoints) );
      memset(realEvPtr->fVolts,     0, sizeof(realEvPtr->fVolts)     );
      memset(realEvPtr->fTimes,     0, sizeof(realEvPtr->fTimes)     );
      
      int fNumPoints = 260;

      for (int iant = 0; iant < settings1->NANTENNAS; iant++){
	//int IceMCAnt = GetIceMCAntfromUsefulEventAnt(anita1,  AnitaGeom1,  iant);
	int IceMCAnt = GetIceMCAntfromUsefulEventAnt(settings1,  iant);
	int UsefulChanIndexH = AnitaGeom1->getChanIndexFromAntPol(iant,  AnitaPol::kHorizontal);
	int UsefulChanIndexV = AnitaGeom1->getChanIndexFromAntPol(iant,  AnitaPol::kVertical);
	realEvPtr->fNumPoints[UsefulChanIndexV] = fNumPoints;
	realEvPtr->fNumPoints[UsefulChanIndexH] = fNumPoints;
	realEvPtr->chanId[UsefulChanIndexV] = UsefulChanIndexV;
	realEvPtr->chanId[UsefulChanIndexH] = UsefulChanIndexH;

	for (int j = 0; j < fNumPoints; j++) {
	  // convert seconds to nanoseconds
	  realEvPtr->fTimes[UsefulChanIndexV][j] = j * anita1->TIMESTEP * 1.0E9;
	  realEvPtr->fTimes[UsefulChanIndexH][j] = j * anita1->TIMESTEP * 1.0E9;
	  // convert volts to millivolts
	  realEvPtr->fVolts[UsefulChanIndexH][j] =  volts_rx_rfcm_lab_h_all[IceMCAnt][j+64]*1000;
	  realEvPtr->fVolts[UsefulChanIndexV][j] =  volts_rx_rfcm_lab_e_all[IceMCAnt][j+64]*1000;
	  realEvPtr->fCapacitorNum[UsefulChanIndexH][j] = 0;
	  realEvPtr->fCapacitorNum[UsefulChanIndexV][j] = 0;
	}//end int j
      }// end int iant

      realEvPtr->eventNumber = eventNumber;

      rawHeaderPtr->eventNumber = eventNumber;
      rawHeaderPtr->surfSlipFlag = 0;
      rawHeaderPtr->errorFlag = 0;

      rawHeaderPtr->trigType = 8; // soft-trigger

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
      rawHeaderPtr->triggerTime = bn1->realTime_flightdata;
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
      truthEvPtr->nuMom            = pnu;
      truthEvPtr->nu_pdg           = pdgcode;
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
	truthEvPtr->balloonPos[i]  = -999;
	truthEvPtr->balloonDir[i]  = -999;
	truthEvPtr->nuPos[i]       = -999;
	truthEvPtr->nuDir[i]       = -999;
      }
      // for (int i=0;i<5;i++){
      // 	for (int j=0;j<3;j++){
      // 	  truthEvPtr->rfExitNor[i][j] = ray1->n_exit2bn[i][j];
      // 	  truthEvPtr->rfExitPos[i][j] = ray1->rfexit[i][j];
      // 	}
      // }
      for (int i=0;i<48;i++){
	truthEvPtr->hitangle_e[i]  = hitangle_e_all[i];
	truthEvPtr->hitangle_h[i]  = hitangle_h_all[i];
      }
      if(settings1->ROUGHNESS){
	for (int i=0;i<Anita::NFREQ;i++)
	  truthEvPtr->vmmhz[i]       = -999;
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

      sum_weights+=weight;
      neutrinos_passing_all_cuts++;

    } // end if passing global trigger conditions
    else {
      passes_thisevent=0; // flag this event as not passing
    }// end else event does not pass trigger

    // delete globaltrig1;

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



  time_t raw_end_time = time(NULL);
  struct tm * end_time = localtime(&raw_end_time);
  cout << "Date and time at end of run are: " << asctime (end_time) << "\n";
  cout<<"Total time elapsed is "<<(int)((raw_end_time - raw_start_time)/60)<<":"<< ((raw_end_time - raw_start_time)%60)<<endl;

  foutput << "\nTotal time elapsed in run is " <<(int)((raw_end_time - raw_start_time)/60)<<":"<< ((raw_end_time - raw_start_time)%60)<<endl;

  //  delete anita1;
  return 0;

} //END MAIN PROGRAM









//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//
//  Auxiliary functions
//
//
//



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
