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
      case 'e':
	exp_tmp=atof(optarg);
	cout << "Changed neutrino energy exponent to " << exp_tmp << endl;
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
  GlobalTrigger *globaltrig1;

  // input parameters
  settings1->ReadInputs(input.c_str(),  foutput, NNU, RANDOMISEPOL);
  settings1->ApplyInputs(anita1,  sec1,  sig1,  bn1,  ray1);

  settings1->SEED=settings1->SEED + run_no;
  gRandom->SetSeed(settings1->SEED);

  bn1->InitializeBalloon();
  anita1->Initialize(settings1, foutput, inu, outputdir);

  if (nnu_tmp!=0)
    NNU=nnu_tmp;
  if (exp_tmp!=0)
    settings1->EXPONENT=exp_tmp;

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


  //added djg ////////////////////////////////////////////////////////
  al_voltages_direct<<"antenna #"<<"   "<<"volts chan 1"<<"   "<<"volts chan 2"<<"    "<<"volts chan 3"<<"    "<<"volts chan 4"<<"    "<<"noise chan 1"<<"    "<<"noise chan 2"<<"    "<<"noise chan 3"<<"   "<<"noise chan 4"<<"  "<<"weight"<<endl;
  ////////////////////////////////////////////////////////////////////

  // ray tracing
  double viewangle=0;
  double offaxis=0; // viewangle-changle,  for plotting

  int beyondhorizon = 0;  //Switch tells if neutrino interacts beyond the balloon's horizon. (0: inside horizon,  1: outside horizon)

  // frequency binning
  double vmmhz1m_max=0; // maximum V/m/MHz at 1m from Jaime (highest frequency)
  double vmmhz_lowfreq=0.; // V/m/MHz after 1/r,  attenuation at the lowest freq.
  double vmmhz[Anita::NFREQ];                        //  V/m/MHz at balloon (after all steps)

  // given the angle you are off the Cerenkov cone,  the fraction of the observed e field that comes from the em shower
  double vmmhz_em[Anita::NFREQ];
  double vmmhz_min_thatpasses=1000;
  double vmmhz_min=0;   // minimum of the above array
  double vmmhz_max=0;                        // maximum of the above array
  double deltheta_em_max, deltheta_had_max;     // maximum value of above array angular distribution
 
  // shower properties
  double emfrac, hadfrac, sumfrac;               // em and had fractions

  string taudecay;                   // tau decay type: e, m, h

  double elast_y=0;                   // inelasticity
  double volts_rx_max=0; // max voltage seen on an antenna - just for debugging purposes
  double volts_rx_ave=0; // ave voltage seen on an antenna,  among hit antennas
  double volts_rx_sum=0; // ave voltage seen on an antenna,  among hit antennas

  double volts_rx_max_highband; // max voltage seen on an antenna - just for debugging purposes
  double volts_rx_max_lowband; // max voltage seen on an antenna - just for debugging purposes
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

  // positions in earth frame
  double horizcoord=0; // x component of interaction position
  double vertcoord=0; // y component of interaction position

  double dist_int_bn_2d=0; // 2-d (surface) distance between interaction and balloon

  double cosalpha; // angle between nu momentum and surface normal at earth entrance point
  double mytheta; //!< alpha minus 90 degrees
  double cosbeta0; // angle between nu momentum and surface normal at interaction point.
  double mybeta; //!< beta minus 90 degrees
  double nuexitice=0;
  double theta_pol_measured; // theta of the polarization as measured at the payload (for now we don't correct for the 10 degree cant)

  double ptaui=0;
  double ptauf =0;


  double sourceLon;
  double sourceAlt;
  double sourceLat;
  double sourceMag;

  Vector n_nutraject_ontheground; //direction of the neutrino from the person standing on the ground just below the balloon.
  Vector n_pol; // direction of polarization
  // Vector n_pol_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX]; // direction of polarization of signal seen at each antenna
  Vector n_pol_db; // same,  double bangs

  int l3trig[Anita::NPOL];  // 16 bit number which says which phi sectors pass L3 V-POL
  // For each trigger layer,  which "clumps" pass L2.  16 bit,  16 bit and 8 bit for layers 1 & 2 and nadirs
  int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];
  //For each trigger layer,  which antennas pass L1.  16 bit,  16 bit and 8 bit and layers 1,  2 and nadirs
  int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];

  // these are declared here so that they can be stuck into trees
  int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //counting how many pass trigger requirement

  int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX]; //counting how many pass trigger requirement

  int nchannels_triggered = 0; // total number of channels triggered
  int nchannels_perrx_triggered[48]; // total number of channels triggered

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

  //we pick both the interaction point and its corresponding mirror point



  //to determine where the cut should be.
  stemp=string(outputdir.Data())+"/icefinal"+run_num+".root";
  TFile *hfile = new TFile(stemp.c_str(), "RECREATE", "ice");
 

  TTree *finaltree = new TTree("passing_events", "passing_events"); // finaltree filled for all events that pass
  finaltree->Branch("inu", &inu, "inu/I");
  finaltree->Branch("vmmhz_min", &vmmhz_min, "vmmhz_min/D");
  finaltree->Branch("vmmhz_max", &vmmhz_max, "vmmhz_max/D");
  finaltree->Branch("thresholdsAnt", &thresholdsAnt, "thresholdsAnt[48][2][5]/D");
  finaltree->Branch("thresholdsAntPass", &thresholdsAntPass, "thresholdsAntPass[48][2][5]/D");
  finaltree->Branch("deadTime", &anita1->deadTime, "deadTime/D");
  finaltree->Branch("horizcoord", &horizcoord, "horizcoord/D");
  finaltree->Branch("vertcoord", &vertcoord, "vertcoord/D");
  finaltree->Branch("horizcoord_bn", &bn1->horizcoord_bn, "horizcoord_bn/D");
  finaltree->Branch("vertcoord_bn", &bn1->vertcoord_bn, "vertcoord_bn/D");
  finaltree->Branch("r_bn", &r_bn_array, "r_bn_array[3]/D");
  finaltree->Branch("n_bn", &n_bn_array, "n_bn_array[3]/D");
  finaltree->Branch("longitude_bn", &longitude_this, "longitude_bn/D");
  finaltree->Branch("heading_bn", &heading_this, "heading_bn/D");
  finaltree->Branch("gps_offset", &gps_offset, "gps_offset/D");
  // this one is just weight due to earth absorption
  finaltree->Branch("weight1", &weight1, "weight1/D");
  // this is the total weight - the one you want to use!
  finaltree->Branch("weight", &weight, "weight/D");
  finaltree->Branch("logweight", &logweight, "logweight/D");
  finaltree->Branch("posnu", &posnu_array, "posnu_array[3]/D");
  finaltree->Branch("nnu", &nnu_array, "nnu_array[3]/D");
  finaltree->Branch("n_exit2bn", &n_exit2bn_array, "n_exit2bn_array[5][3]/D");
  finaltree->Branch("n_exit_phi", &n_exit_phi, "n_exit_phi/D");
  finaltree->Branch("rfexit", &rfexit_array, "rfexit_array[5][3]/D");

  finaltree->Branch("pnu", &pnu, "pnu/D");
  finaltree->Branch("elast_y", &elast_y, "elast_y/D");
  finaltree->Branch("emfrac", &emfrac, "emfrac/D");
  finaltree->Branch("hadfrac", &hadfrac, "hadfrac/D");
  finaltree->Branch("sumfrac", &sumfrac, "sumfrac/D");
  finaltree->Branch("nuexitice",  &nuexitice,  "nuexitice/D");
  finaltree->Branch("l3trig", &l3trig, "l3trig[2]/I");
  finaltree->Branch("l2trig", &l2trig, "l2trig[2][3]/I");
  finaltree->Branch("l1trig", &l1trig, "l1trig[2][3]/I");
  finaltree->Branch("phiTrigMask", &anita1->phiTrigMask, "phiTrigMask/s");
  finaltree->Branch("phiTrigMaskH", &anita1->phiTrigMaskH, "phiTrigMaskH/s");
  finaltree->Branch("l1TrigMask", &anita1->l1TrigMask, "l1TrigMask/s");
  finaltree->Branch("l1TrigMaskH", &anita1->l1TrigMaskH, "l1TrigMaskH/s");
  finaltree->Branch("max_antenna0", &max_antenna0, "max_antenna0/I");
  finaltree->Branch("max_antenna1", &max_antenna1, "max_antenna1/I");
  finaltree->Branch("max_antenna2", &max_antenna2, "max_antenna2/I");

  finaltree->Branch("viewangle", &viewangle, "viewangle/D");
  finaltree->Branch("offaxis", &offaxis, "offaxis/D");
  finaltree->Branch("rx0_signal_eachband", &rx0_signal_eachband, "rx0_signal_eachband[2][5]/D");
  finaltree->Branch("rx0_threshold_eachband", &rx0_threshold_eachband, "rx0_threshold_eachband[2][5]/D");
  finaltree->Branch("rx0_noise_eachband", &rx0_noise_eachband, "rx0_noise_eachband[2][5]/D");
  finaltree->Branch("rx0_passes_eachband", &rx0_passes_eachband, "rx0_passes_eachband[2][5]/I");
  finaltree->Branch("e_component", &e_component, "e_component/D");
  finaltree->Branch("h_component", &h_component, "h_component/D");
  finaltree->Branch("dist_int_bn_2d", &dist_int_bn_2d, "dist_int_bn_2d/D");

  finaltree->Branch("cosalpha", &cosalpha, "cosalpha/D");
  finaltree->Branch("mytheta", &mytheta, "mytheta/D");
  finaltree->Branch("cosbeta0", &cosbeta0, "cosbeta0/D");
  finaltree->Branch("mybeta", &mybeta, "mybeta/D");

  //Begin block added by Stephen for verification plots
  finaltree->Branch("fresnel1", &fresnel1, "fresnel1/D");
  finaltree->Branch("fresnel2", &fresnel2, "fresnel2/D");
  finaltree->Branch("mag1", &mag1, "mag1/D");
  finaltree->Branch("mag2", &mag2, "mag2/D");
  finaltree->Branch("t_coeff_pokey", &t_coeff_pokey, "t_coeff_pokey/D");
  finaltree->Branch("t_coeff_slappy", &t_coeff_slappy, "t_coeff_slappy/D");

  finaltree->Branch("hitangle_e_all", &hitangle_e_all, "hitangle_e_all[48]/D");
  finaltree->Branch("hitangle_h_all", &hitangle_h_all, "hitangle_h_all[48]/D");

  finaltree->Branch("e_comp_max1", &e_comp_max1, "e_comp_max1/D");
  finaltree->Branch("h_comp_max1", &h_comp_max1, "h_comp_max1/D");
  finaltree->Branch("e_comp_max2", &e_comp_max2, "e_comp_max2/D");
  finaltree->Branch("h_comp_max2", &h_comp_max2, "h_comp_max2/D");
  finaltree->Branch("e_comp_max3", &e_comp_max3, "e_comp_max3/D");
  finaltree->Branch("h_comp_max3", &h_comp_max3, "h_comp_max3/D");
  finaltree->Branch("max_antenna_volts0", &max_antenna_volts0, "max_antenna_volts0/D");
  finaltree->Branch("max_antenna_volts0_em", &max_antenna_volts0_em, "max_antenna_volts0_em/D");
  finaltree->Branch("max_antenna_volts1", &max_antenna_volts1, "max_antenna_volts1/D");
  finaltree->Branch("max_antenna_volts2", &max_antenna_volts2, "max_antenna_volts2/D");
  finaltree->Branch("triggers", &nchannels_perrx_triggered, "nchannels_perrx_triggered[48]/I");
  finaltree->Branch("nchannels_triggered", &nchannels_triggered, "nchannels_triggered/I");
  finaltree->Branch("volts_rx_max", &volts_rx_max, "volts_rx_max/D");
  finaltree->Branch("volts_rx_ave", &volts_rx_ave, "volts_rx_ave/D");
  finaltree->Branch("volts_rx_sum", &volts_rx_sum, "volts_rx_sum/D");

  finaltree->Branch("volts_rx_max_highband", &volts_rx_max_highband, "volts_rx_max_highband/D");
  finaltree->Branch("volts_rx_max_lowband", &volts_rx_max_lowband, "volts_rx_max_lowband/D");
  finaltree->Branch("theta_pol_measured", &theta_pol_measured, "theta_pol_measured/D");
  finaltree->Branch("theta_rf_atbn", &theta_rf_atbn, "theta_rf_atbn/D");
  finaltree->Branch("theta_rf_atbn_measured", &theta_rf_atbn_measured, "theta_rf_atbn_measured/D");
  finaltree->Branch("voltage", &voltagearray, "voltagearray[48]/D");
  finaltree->Branch("vmmhz1m_max", &vmmhz1m_max, "vmmhz1m_max/D");
  finaltree->Branch("vmmhz_lowfreq", &vmmhz_lowfreq, "vmmhz_lowfreq/D");

  finaltree->Branch("deltheta_em_max", &deltheta_em_max, "deltheta_em_max/D");
  finaltree->Branch("deltheta_had_max", &deltheta_had_max, "deltheta_had_max/D");
  finaltree->Branch("r_enterice", &r_enterice_array, "r_enterice_array[3]/D");
  finaltree->Branch("n_exit2bn_db", &n_exit2bn_db_array, "n_exit2bn_db_array[5][3]/D");

  finaltree->Branch("rfexit_db", &rfexit_db_array, "rfexit_db_array[5][3]/D");
  finaltree->Branch("r_in", &r_in_array, "r_in_array[3]/D");
  finaltree->Branch("nsurf_rfexit", &nsurf_rfexit_array, "nsurf_rfexit_array[3]/D");
  finaltree->Branch("nsurf_rfexit_db", &nsurf_rfexit_db_array, "nsurf_rfexit_db_array[3]/D");

  finaltree->Branch("ant_normal0", &ant_max_normal0_array, "ant_max_normal0_array[3]/D");
  finaltree->Branch("ant_normal1", &ant_max_normal1_array, "ant_max_normal1_array[3]/D");
  finaltree->Branch("ant_normal2", &ant_max_normal2_array, "ant_max_normal2_array[3]/D");
  finaltree->Branch("vmmhz1m_visible", &vmmhz1m_visible, "vmmhz1m_visible/D");
  finaltree->Branch("freq_bins", &freq_bins, "freq_bins/I");
  finaltree->Branch("vmmhz", &vmmhz, "vmmhz[freq_bins]/D");

  finaltree->Branch("n_pol", &n_pol_array, "n_pol_array[3]/D");
  finaltree->Branch("vmmhz_min_thatpasses", &vmmhz_min_thatpasses, "vmmhz_min_thatpasses/D");

  finaltree->Branch("anita1->PHI_OFFSET", &anita1->PHI_OFFSET, "anita1->PHI_OFFSET/D");
  finaltree->Branch("igps", &bn1->igps, "igyps/I");
  finaltree->Branch("volts_rx_rfcm_lab_e_all", &volts_rx_rfcm_lab_e_all, "volts_rx_rfcm_lab_e_all[48][512]/D");
  finaltree->Branch("volts_rx_rfcm_lab_h_all", &volts_rx_rfcm_lab_h_all, "volts_rx_rfcm_lab_h_all[48][512]/D");
  finaltree->Branch("ptaui", &ptaui, "ptaui/D");
  finaltree->Branch("ptauf", &ptauf, "ptauf/D");
  finaltree->Branch("sourceLon", &sourceLon, "sourceLon/D");
  finaltree->Branch("sourceLat", &sourceLat, "sourceLat/D");
  finaltree->Branch("sourceAlt", &sourceAlt, "sourceAlt/D");
  finaltree->Branch("sourceMag", &sourceMag, "sourceMag/D");

  double rms_rfcm_e;
  double rms_rfcm_h;
  double rms_lab_e;
  double rms_lab_h;

  double avgfreq_rfcm[Anita::NFREQ];
  double avgfreq_rfcm_lab[Anita::NFREQ];
  double freq[Anita::NFREQ];

  int pdgcode=-999;
  UInt_t eventNumber;

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
  
  TTree *triggerSettingsTree = new TTree("triggerSettingsTree", "Trigger settings");
  triggerSettingsTree->Branch("dioderms", anita1->bwslice_dioderms_fullband_allchan, "dioderms[2][48][7]/D");
  triggerSettingsTree->Branch("diodemean", anita1->bwslice_diodemean_fullband_allchan, "diodemean[2][48][7/D");
  triggerSettingsTree->Fill();

  TTree *truthAnitaTree = new TTree("truthAnitaTree", "Truth Anita Tree");
  truthAnitaTree->Branch("truth",     &truthEvPtr                   );
#endif

  AnitaGeomTool *AnitaGeom1 = AnitaGeomTool::Instance();
  
#endif

  
  IceModel *antarctica = new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10, settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0, settings1->WEIGHTABSORPTION);
  cout << "area of the earth's surface covered by antarctic ice is " << antarctica->ice_area << "\n";

  // fills arrays according to antenna specs
  anita1->GetBeamWidths(settings1); // this is used if GAINS set to 0
  // Antenna measured gain vs. frequency
  anita1->ReadGains(); // this is used if GAINS set to 1
  anita1->Set_gain_angle(settings1, sig1->NMEDIUM_RECEIVER);
  if(settings1->WHICH == 2 || settings1->WHICH == 6) anita1->SetDiffraction(); // for the upper ring

  
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

 
  // for averaging balloon altitude and distance from center of earth
  // for comparing with Peter
  double average_altitude=0.;
  double average_rbn=0.;

  //  TRandom r(settings1->SEED); // use seed set as input

  signal(SIGINT,  interrupt_signal_handler);     // This function call allows icemc to gracefully abort and write files as usual rather than stopping abruptly.

  // Setting gps offset for plotting direction wrt north
  if (settings1->WHICH==7){
    gps_offset=atan2(-0.7042,0.71)*DEGRAD;
  } else if(settings1->WHICH==8){
    gps_offset=atan2(-0.7085,0.7056)*DEGRAD;
  } else if (settings1->WHICH==9 || settings1->WHICH==10){
    gps_offset=45;
  } else gps_offset=0;

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
    anita1->inu = inu;
    
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
    sec1->secondbang=false;
    count_total++;
    // initializing the voltage seen by each polarization of each antenna
    bn1->dtryingposition=0;
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

   
    // pick random point in ice.
    // also get initial guess shower exit position
    // unit vector from shower exit to balloon.
    // get both positions of interaction point and its mirror point.
    //-------------------------------------------------------
    beyondhorizon = 0;

    // make a global trigger object (but don't touch the electric fences)
    globaltrig1 = new GlobalTrigger(settings1, anita1);

    Tools::Zero(anita1->arrival_times[0], Anita::NLAYERS_MAX*Anita::NPHI_MAX);
    Tools::Zero(anita1->arrival_times[1], Anita::NLAYERS_MAX*Anita::NPHI_MAX);

    anita1->calculateDelaysForEfficiencyScan();
    
    globaltrig1->volts_rx_rfcm_trigger.assign(16,  vector <vector <double> >(3,  vector <double>(0)));
    anita1->rms_rfcm_e_single_event = 0;

    for (int ilayer=0; ilayer < settings1->NLAYERS; ilayer++) { // loop over layers on the payload
      for (int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) { // ifold loops over phi

	antNum = anita1->GetRxTriggerNumbering(ilayer, ifold);
	
	ChanTrigger *chantrig1 = new ChanTrigger();
	chantrig1->InitializeEachBand(anita1);
	  
	//	  chantrig1->ApplyAntennaGain(settings1, anita1, bn1, panel1, antNum, n_eplane, n_hplane, n_normal);

	
#ifdef ANITA_UTIL_EXISTS
	if (settings1->SIGNAL_FLUCT && (settings1->NOISEFROMFLIGHTDIGITIZER || settings1->NOISEFROMFLIGHTTRIGGER) )
	  chantrig1->getNoiseFromFlight(anita1, antNum);
  
	if(!settings1->APPLYIMPULSERESPONSETRIGGER) chantrig1->injectImpulseAfterAntenna(anita1, antNum);
#endif
	
	chantrig1->TriggerPath(settings1, anita1, antNum, bn1);
	  
	chantrig1->DigitizerPath(settings1, anita1, antNum, bn1);
	
	chantrig1->WhichBandsPass(settings1, anita1, globaltrig1, bn1, ilayer, ifold,  viewangle-sig1->changle, emfrac, hadfrac, thresholdsAnt[antNum]);

	chantrig1->TimeShiftAndSignalFluct(settings1, anita1, ilayer, ifold, volts_rx_rfcm_lab_e_all,  volts_rx_rfcm_lab_h_all);

	chantrig1->saveTriggerWaveforms(anita1, justSignal_trig[0][antNum], justSignal_trig[1][antNum], justNoise_trig[0][antNum], justNoise_trig[1][antNum]);
	chantrig1->saveDigitizerWaveforms(anita1, justSignal_dig[0][antNum], justSignal_dig[1][antNum], justNoise_dig[0][antNum], justNoise_dig[1][antNum]);
	
	//	cout << inu << " " << ilayer << " " << ifold << endl;
	delete chantrig1;
	       
      } //loop through the phi-fold antennas
    }  //loop through the layers of antennas

    anita1->rms_rfcm_e_single_event = sqrt(anita1->rms_rfcm_e_single_event / (anita1->HALFNFOUR * settings1->NANTENNAS));


    for (int irx=0;irx<settings1->NANTENNAS;irx++) {
      nchannels_perrx_triggered[irx]=globaltrig1->nchannels_perrx_triggered[irx];
    }

    nchannels_triggered=Tools::iSum(globaltrig1->nchannels_perrx_triggered, settings1->NANTENNAS); // find total number of antennas that were triggered.
    volts_rx_ave=GetAverageVoltageFromAntennasHit(settings1, globaltrig1->nchannels_perrx_triggered, voltagearray, volts_rx_sum);

   
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
	  count1->nl1triggers[i][whichray]+=anita1->l1trig_anita3and4_inanita[i][j][k];
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
	 || (settings1->TRIGTYPE==0 && count_pass>=settings1->NFOLD)
	 || (settings1->MINBIAS==1)) {
      if (bn1->WHICHPATH==4)
	cout << "This event passes.\n";

      //	cout << inu << endl;

      anita1->passglobtrig[0]=thispasses[0];
      anita1->passglobtrig[1]=thispasses[1];

      //calculate the phi angle wrt +x axis of the ray from exit to balloon
      n_exit_phi = Tools::AbbyPhiCalc(ray1->n_exit2bn[2][0], ray1->n_exit2bn[2][1]);

      // tags this event as passing
      passestrigger=1;

      crust_entered=0; //These are switches that let us tell how far a given neutrino penetrated.  Clear them before entering Getchord.
      mantle_entered=0;
      core_entered=0;


      finaltree->Fill();
	      
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
      for (int i=0;i<5;i++){
	for (int j=0;j<3;j++){
	  truthEvPtr->rfExitNor[i][j] = ray1->n_exit2bn[i][j];
	  truthEvPtr->rfExitPos[i][j] = ray1->rfexit[i][j];
	}
      }
      for (int i=0;i<48;i++){
	truthEvPtr->hitangle_e[i]  = hitangle_e_all[i];
	truthEvPtr->hitangle_h[i]  = hitangle_h_all[i];
      }
      if(settings1->ROUGHNESS){
	for (int i=0;i<Anita::NFREQ;i++)
	  truthEvPtr->vmmhz[i]       = -999;
      }

      memset(truthEvPtr->SNRAtTrigger,       0, sizeof(truthEvPtr->SNRAtTrigger)       );
      memset(truthEvPtr->fSignalAtTrigger,   0, sizeof(truthEvPtr->fSignalAtTrigger)   );
      memset(truthEvPtr->fNoiseAtTrigger,    0, sizeof(truthEvPtr->fNoiseAtTrigger)    );
      memset(truthEvPtr->SNRAtDigitizer,     0, sizeof(truthEvPtr->SNRAtDigitizer)     );
      memset(truthEvPtr->thresholds,         0, sizeof(truthEvPtr->thresholds)         );
      memset(truthEvPtr->fDiodeOutput,       0, sizeof(truthEvPtr->fDiodeOutput)       );
	    
      truthEvPtr->maxSNRAtTriggerV=0;
      truthEvPtr->maxSNRAtTriggerH=0;
      truthEvPtr->maxSNRAtDigitizerV=0;
      truthEvPtr->maxSNRAtDigitizerH=0;

      for (int iant = 0; iant < settings1->NANTENNAS; iant++){
	int UsefulChanIndexH = AnitaGeom1->getChanIndexFromAntPol(iant,  AnitaPol::kHorizontal);
	int UsefulChanIndexV = AnitaGeom1->getChanIndexFromAntPol(iant,  AnitaPol::kVertical);

	truthEvPtr->SNRAtTrigger[UsefulChanIndexV] = Tools::calculateSNR(justSignal_trig[0][iant], justNoise_trig[0][iant]);
	truthEvPtr->SNRAtTrigger[UsefulChanIndexH] = Tools::calculateSNR(justSignal_trig[1][iant], justNoise_trig[1][iant]);
	      
	if (truthEvPtr->SNRAtTrigger[UsefulChanIndexV]>truthEvPtr->maxSNRAtTriggerV) truthEvPtr->maxSNRAtTriggerV=truthEvPtr->SNRAtTrigger[UsefulChanIndexV];
	if (truthEvPtr->SNRAtTrigger[UsefulChanIndexH]>truthEvPtr->maxSNRAtTriggerH) truthEvPtr->maxSNRAtTriggerH=truthEvPtr->SNRAtTrigger[UsefulChanIndexH];

	truthEvPtr->SNRAtDigitizer[UsefulChanIndexV] = Tools::calculateSNR(justSignal_dig[0][iant], justNoise_dig[0][iant]);
	truthEvPtr->SNRAtDigitizer[UsefulChanIndexH] = Tools::calculateSNR(justSignal_dig[1][iant], justNoise_dig[1][iant]);
	      
	if (truthEvPtr->SNRAtDigitizer[UsefulChanIndexV]>truthEvPtr->maxSNRAtDigitizerV) truthEvPtr->maxSNRAtDigitizerV=truthEvPtr->SNRAtDigitizer[UsefulChanIndexV];
	if (truthEvPtr->SNRAtDigitizer[UsefulChanIndexH]>truthEvPtr->maxSNRAtDigitizerH) truthEvPtr->maxSNRAtDigitizerH=truthEvPtr->SNRAtDigitizer[UsefulChanIndexH];

	      
	truthEvPtr->thresholds[UsefulChanIndexV] = thresholdsAnt[iant][0][4];
	truthEvPtr->thresholds[UsefulChanIndexH] = thresholdsAnt[iant][1][4];
	int irx = iant;
	if (iant<16){
	  if (iant%2) irx = iant/2;
	  else        irx = iant/2 + 1;
	}
	      
	for (int j = 0; j < fNumPoints; j++) {
	  truthEvPtr->fTimes[UsefulChanIndexV][j]             = j * anita1->TIMESTEP * 1.0E9;
	  truthEvPtr->fTimes[UsefulChanIndexH][j]             = j * anita1->TIMESTEP * 1.0E9;
		
	  truthEvPtr->fSignalAtTrigger[UsefulChanIndexV][j]   = justSignal_trig[0][iant][j]*1000;
	  truthEvPtr->fSignalAtTrigger[UsefulChanIndexH][j]   = justSignal_trig[1][iant][j]*1000;
	  truthEvPtr->fNoiseAtTrigger[UsefulChanIndexV][j]    = justNoise_trig[0][iant][j]*1000;
	  truthEvPtr->fNoiseAtTrigger[UsefulChanIndexH][j]    = justNoise_trig[1][iant][j]*1000;
	  truthEvPtr->fSignalAtDigitizer[UsefulChanIndexV][j] = justSignal_dig[0][iant][j]*1000;
	  truthEvPtr->fSignalAtDigitizer[UsefulChanIndexH][j] = justSignal_dig[1][iant][j]*1000;
	  truthEvPtr->fNoiseAtDigitizer[UsefulChanIndexV][j]  = justNoise_dig[0][iant][j]*1000;
	  truthEvPtr->fNoiseAtDigitizer[UsefulChanIndexH][j]  = justNoise_dig[1][iant][j]*1000;
		
	  truthEvPtr->fDiodeOutput[UsefulChanIndexV][j]       = anita1->timedomain_output_allantennas[0][irx][j];
	  truthEvPtr->fDiodeOutput[UsefulChanIndexH][j]       = anita1->timedomain_output_allantennas[1][irx][j];
	}//end int j
	      
      }// end int iant

      
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
      anita1->tdata->Fill();
    } // end if passing global trigger conditions
    else {
      passes_thisevent=0; // flag this event as not passing
    }// end else event does not pass trigger

    delete globaltrig1;

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
  configAnitaTree->Write("configAnitaTree");
  truthAnitaTree->Write("truthAnitaTree");
  triggerSettingsTree->Write("triggerSettingsTree");
  anitafileTruth->Close();
  delete anitafileTruth;
#endif

#endif


  anita1->rms_rfcm[0] = sqrt(anita1->rms_rfcm[0] / (double)anita1->count_getnoisewaveforms)*1000.;
  anita1->rms_rfcm[1] = sqrt(anita1->rms_rfcm[1] / (double)anita1->count_getnoisewaveforms)*1000.;
  anita1->rms_lab[0] = sqrt(anita1->rms_lab[0] / (double)anita1->count_getnoisewaveforms)*1000.;
  anita1->rms_lab[1] = sqrt(anita1->rms_lab[1] / (double)anita1->count_getnoisewaveforms)*1000.;

  cout << "RMS noise in rfcm e-pol is " << anita1->rms_rfcm[0] << " mV.\n";
  cout << "RMS noise in rfcm h-pol is " << anita1->rms_rfcm[1] << " mV.\n";
  cout << "RMS noise in lab e-pol is " << anita1->rms_lab[0] << "mV.\n";
  cout << "RMS noise in lab h-pol is " << anita1->rms_lab[1] << "mV.\n";
  for (int i=0;i<Anita::NFREQ;i++) {
    anita1->avgfreq_rfcm[i]/=(double)anita1->count_getnoisewaveforms;
    anita1->avgfreq_rfcm_lab[i]/=(double)anita1->count_getnoisewaveforms;
  }

  rms_rfcm_e=anita1->rms_rfcm[0];
  rms_rfcm_h=anita1->rms_rfcm[1];
  rms_lab_e=anita1->rms_lab[0];
  rms_lab_h=anita1->rms_lab[1];
  for (int i=0;i<Anita::NFREQ;i++) {
    avgfreq_rfcm[i]=anita1->avgfreq_rfcm[i];
    avgfreq_rfcm_lab[i]=anita1->avgfreq_rfcm_lab[i];
    freq[i]=anita1->freq[i];
  }


  // for each neutrino flavor,  fraction each contributes to sensitivity.
  sum_frac[0]=sum[0]/eventsfound;
  sum_frac[1]=sum[1]/eventsfound;
  sum_frac[2]=sum[2]/eventsfound;

  // for taus.
  sum_frac_db[0]=sum[0]/(eventsfound+eventsfound_db+eventsfound_nfb);
  sum_frac_db[1]=sum[1]/(eventsfound+eventsfound_db+eventsfound_nfb);
  sum_frac_db[2]=(sum[2]+eventsfound_db+eventsfound_nfb)/(eventsfound+eventsfound_db+eventsfound_nfb);
  //if (tree17->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && HIST==1)
  //tree17->Fill();


  cout << "closing file.\n";
  CloseTFile(hfile);

  time_t raw_end_time = time(NULL);
  struct tm * end_time = localtime(&raw_end_time);
  cout << "Date and time at end of run are: " << asctime (end_time) << "\n";
  cout<<"Total time elapsed is "<<(int)((raw_end_time - raw_start_time)/60)<<":"<< ((raw_end_time - raw_start_time)%60)<<endl;

  foutput << "\nTotal time elapsed in run is " <<(int)((raw_end_time - raw_start_time)/60)<<":"<< ((raw_end_time - raw_start_time)%60)<<endl;
  anita1->fdata->Write();
  anita1->fdata->Close();

  delete anita1;
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

void IntegrateBands(Anita *anita1, int k, Screen *panel1, double *freq, double scalefactor, double *sumsignal) {
  for (int j=0;j<5;j++) {
    // if this frequency is in this bandwidth slice
    for (int jpt=0; jpt<panel1->GetNvalidPoints(); jpt++){
      if (anita1->bwslice_min[j]<=freq[k] && anita1->bwslice_max[j]>freq[k])
        sumsignal[j]+=panel1->GetVmmhz_freq(jpt*Anita::NFREQ + k)*(freq[k+1]-freq[k])*scalefactor;
    }
  }
}
//end IntegrateBands()


void Integrate(Anita *anita1, int j, int k, double *vmmhz, double *freq, double scalefactor, double sumsignal) {
  // if this frequency is in this bandwidth slice
  if (anita1->bwslice_min[j]<=freq[k] && anita1->bwslice_max[j]>freq[k])
    sumsignal+=vmmhz[k]*(freq[k+1]-freq[k])*scalefactor;
}
//end Integrate()


void WriteNeutrinoInfo(Position &posnu,  Vector &nnu,  Position &r_bn,  double altitude_int,  string nuflavor,  string current,  double elast_y,  ofstream &nu_out) {
  nu_out << "\n" << inu << "\t" << posnu[0] << " " << posnu[1] << " " << posnu[2] << "\t" << altitude_int << "\t" << nnu[0] << " " << nnu[1] << " " << nnu[2] << "\t" << r_bn[0] << " " << r_bn[1] << " " << r_bn[2] << "\t" << nuflavor << "\t" << current << "\t" << elast_y << "\n\n";
}
//end WriteNeutrinoInfo()



double GetAirDistance(double altitude_bn, double beta) { // given beta=angle wrt horizontal that the ray hits the balloon,  calculate distance that the ray traveled in air,  including curvature of earth
  return EarthModel::R_EARTH*acos((altitude_bn+EarthModel::R_EARTH)/EarthModel::R_EARTH*(1-sin(beta)*sin(beta))+1/EarthModel::R_EARTH*sin(beta)*sqrt((altitude_bn+EarthModel::R_EARTH)*(altitude_bn+EarthModel::R_EARTH)*sin(beta)*sin(beta)-2*EarthModel::R_EARTH*altitude_bn-altitude_bn*altitude_bn));
}
//end GetAirDistance()


double GetAverageVoltageFromAntennasHit(Settings *settings1, int *nchannels_perrx_triggered, double *voltagearray, double& volts_rx_sum) {
  double sum=0;
  int count_hitantennas=0;
  for (int i=0;i<settings1->NANTENNAS;i++) {
    if (nchannels_perrx_triggered[i]>=3) {
      sum+=voltagearray[i];
      count_hitantennas++;
    } //if
  } //for
  volts_rx_sum = sum;
  sum = sum/(double)count_hitantennas;
  return sum;
}
//end GetAverageVoltageFromAntennasHit()


Vector GetPolarization(const Vector &nnu, const Vector &nrf2_iceside) {
  // Want to find a unit vector in the same plane as nnu and n_refr,
  // but perpendicular to n_refr,  pointing away from nnu.

  // cross nnu with n_refr to get the direction of the B field.
  Vector n_bfield = nnu.Cross(nrf2_iceside);
  // cross b-field with nrf2_iceside to get the polarization vector.
  Vector n_pol = n_bfield.Cross(nrf2_iceside);
  n_pol = n_pol.Unit();
  // check and make sure E-field is pointing in the right direction.
  if (nnu*nrf2_iceside>0 && n_pol*nnu>0){
    cout << "error in GetPolarization.  Event is " << inu << "\n";
  }
  return n_pol;
}
//end GetPolarization()


void CloseTFile(TFile *hfile) {
  hfile->cd();
  hfile->Write();
  hfile->Close();
}
//end CloseTFile()


double IsItDoubleBang(double exitlength, double plepton) {
  double gamma=plepton/MTAU;
  return 1-exp(-1*exitlength/(TAUDECAY_TIME*CLIGHT*gamma));
}
//end IsItDoubleBang()


int WhereIsSecondBang(const Position &posnu, const Vector &nnu, double nuexitlength, double pnu, IceModel *antarctica1, const Position &r_bn, Position &posnu2, Position &rfexit_db, Vector &n_exit2bn_db) {
  double rnd1=0;
  double rnd2=2;
  double gamma=pnu/MTAU;

  if (exp(-1*nuexitlength/(TAUDECAY_TIME*CLIGHT*gamma))>0.999){
    rnd1=gRandom->Rndm()*nuexitlength;
  }
  else {
    while (rnd2>1-exp(-1*rnd1/(TAUDECAY_TIME*CLIGHT*gamma))) {
      rnd1=gRandom->Rndm()*nuexitlength;
      rnd2=gRandom->Rndm();
    } //while
  } //else
  posnu2 = posnu + rnd1*nnu;
  rfexit_db = antarctica1->Surface(posnu2)*posnu2.Unit();

  // unit vector pointing to antenna from exit point.
  n_exit2bn_db = (r_bn - rfexit_db) / r_bn.Distance(rfexit_db);

  double cosangle=(n_exit2bn_db * posnu2) / posnu2.Mag();
  if (cosangle<0){
    return 0;
  }
  return 1;
}
//end WhereIsSecondBang()


//the following is  a new function only for reflected case.
void Attenuate_down(IceModel *antarctica1, Settings *settings1, double& vmmhz_max, const Position &rfexit2, const Position &posnu, const Position &posnu_down) {
  double ATTENLENGTH=700;
  if(!settings1->VARIABLE_ATTEN){
    ATTENLENGTH=antarctica1->EffectiveAttenuationLength(settings1, posnu, 1);
  }

  int position_in_iceshelves=antarctica1->IceOnWater(posnu);
  int position_in_rossexcept=antarctica1->RossExcept(posnu);
  // int position_in_ross = antarctica->RossIceShelf(posnu);
  // int position_in_ronne = antarctica->RonneIceShelf(posnu);
  double dtemp=posnu_down.Distance(rfexit2)/ATTENLENGTH;

  if (dtemp<20) {
    // if(position_in_ross || position_in_ronne) {
    if(position_in_iceshelves && (!position_in_rossexcept)){
      // scalefactor_attenuation=0.310227766*exp(-dtemp);
      // vmmhz_max=vmmhz_max*exp(-dtemp)*0.310227766;//10% of power reflected
      scalefactor_attenuation=0.71*exp(-dtemp);
      vmmhz_max=vmmhz_max*0.71*exp(-dtemp);//50% of power reflected. -3dB
    } //end if
    else if(position_in_rossexcept){
      scalefactor_attenuation=0.1*exp(-dtemp);
      vmmhz_max=0.1*vmmhz_max*exp(-dtemp);//1% of power reflected. -20dB
    }//end else if
    else {
      scalefactor_attenuation=sqrt(0.001)*exp(-dtemp);
      vmmhz_max=sqrt(0.001)*vmmhz_max*exp(-dtemp);//0.1% of power reflected.-30dB
    } //else
  } //if
  else {
    scalefactor_attenuation=0;
    vmmhz_max=0;
  } //else
}
//end Attenuate_down()


void Attenuate(IceModel *antarctica1, Settings *settings1, double& vmmhz_max,  double rflength, const Position &posnu) {
  double ATTENLENGTH=700;  // constant attenuation length for now.
  if (!settings1->VARIABLE_ATTEN){
    ATTENLENGTH = antarctica1->EffectiveAttenuationLength(settings1, posnu, 0);
  }

  double dtemp=(rflength/ATTENLENGTH);
  if(!settings1->ROUGHNESS){
    if (dtemp<20) {
      scalefactor_attenuation=exp(-dtemp);
      vmmhz_max=vmmhz_max*exp(-dtemp);
    } //if
    else {
      scalefactor_attenuation=0;
      vmmhz_max=0;
    } //else
  }
  else{       // use a larger allowable value in case of roughness
    if (dtemp<10000) {
      scalefactor_attenuation=exp(-dtemp);
      vmmhz_max=vmmhz_max*exp(-dtemp);
    } //if
    else {
      scalefactor_attenuation=0;
      vmmhz_max=0;
    } //else
  }
}
//end Attenuate()


void IsAbsorbed(double chord_kgm2, double len_int_kgm2, double &weight1) {
  // see if neutrino is absorbed
  //  weighting works,  but not to much purpose since nu's always
  //   interact at these energies.
  double rtemp;

  rtemp=chord_kgm2/len_int_kgm2;
  if (rtemp<=20)
    weight1=exp(-rtemp);
  else
    weight1=0;
}
//end IsAbsorbed()


void GetSmearedIncidentAngle(Vector &specular, Vector &nrf_iceside, Vector &n_exit2bn, double SMEARINCIDENTANGLE){
  //  void GetSmearedIncidentAngle(Vector &specular, Vector &nsurf_rfexit, Vector &nrf_iceside, Vector &n_exit2bn, double SMEARINCIDENTANGLE, double theta_inc_smeared) {
  // Smear the incident angle for roughness studies
  specular+=nrf_iceside; // specular is the ray that we got from Snell's law
  Vector parallel_to_surface; // find vector parallel to surface to rotate the vector around
  parallel_to_surface+=n_exit2bn; // want to cross specular with n_exit2bn
  parallel_to_surface.Cross(specular);
  nrf_iceside.Rotate(SMEARINCIDENTANGLE*(2*gRandom->Rndm()-1.), parallel_to_surface); // smear the incident ray
  //   theta_inc_smeared=acos(nrf_iceside.Dot(nsurf_rfexit));
}
//end GetSmearedIncidentAngle()


int GetRayIceSide(const Vector &n_exit2rx,  const Vector &nsurf_rfexit, double nexit,  double nenter,  Vector &nrf2_iceside) {
  // this function performs snell's law in three dimensions
  double costh=0;
  double NRATIO=nexit/nenter;
  costh=(n_exit2rx*nsurf_rfexit)/(n_exit2rx.Mag() * nsurf_rfexit.Mag()); // cos(theta) of the transmission angle

  if (costh<0) {
    //cout << "returning 0.  inu is " << inu << "\n";
    return 0;
  }
  double sinth=sqrt(1 - costh*costh);
  double factor=NRATIO*costh-sqrt(1-(NRATIO*sinth*NRATIO*sinth));
  nrf2_iceside = -factor*nsurf_rfexit + NRATIO*n_exit2rx;
  nrf2_iceside = nrf2_iceside.Unit(); // normalize
  return 1;
}
//end GetRayIceSide()



double ScaleVmMHz(double vmmhz1m_max, const Position &posnu1, const Position &r_bn) {
  double dtemp= r_bn.Distance(posnu1);
  vmmhz1m_max= vmmhz1m_max/dtemp;
  scalefactor_distance=1/dtemp;
  //cout << "dtemp is " << dtemp << "\n";
  return vmmhz1m_max;
}
//end ScaleVmMHz()


void SetupViewangles(Signal *sig1) {
  double viewangle_max=90.*RADDEG;
  double viewangle_min=30.*RADDEG;
  for (int i=0;i<NVIEWANGLE-2;i++) {
    viewangles[i]=viewangle_max-(viewangle_max-viewangle_min)/(double)(NVIEWANGLE-2)*(double)i;
  }
  viewangles[NVIEWANGLE-2]=acos(1/sig1->N_DEPTH);
  viewangles[NVIEWANGLE-1]=90.*RADDEG;
}
//end SetupViewAngles()


double GetThisAirColumn(Settings* settings1,  Position r_in, Vector nnu, Position posnu,  double *col1,  double& cosalpha, double& mytheta, double& cosbeta0, double& mybeta) {
  double myair=0; // this is the output
  // it is the column of air in kg/m^2
  cosalpha=(r_in * nnu) / r_in.Mag(); // cosangle that the neutrino enters the earth wrt surface normal at its entrry point
  mytheta=(double)(acos(cosalpha)*DEGRAD)-90.; // turn this into an angle

  //------------------added on Dec 8------------------------
  if (settings1->ATMOSPHERE) {
    int index11=int(mytheta*10.); // which index this theta corresponds to
    int index12=index11+1;
    // find column of air at this theta
    myair=(col1[index11]+(col1[index12]-col1[index11])*(mytheta*10.-double(index11)))*10.;//unit is kg/m^2
  }
  else
    myair=0.;//don't include effect of atmosphere

  //cout<<"mytheta="<<mytheta<<"; myair="<<myair<<endl;
  //------------------added on Dec 8------------------------
  cosbeta0= (posnu * nnu) / posnu.Mag(); // cos angle of neutrino wrt person standing over the interaction point
  mybeta=(double)(acos(cosbeta0)*DEGRAD)-90.; // turn that into a theta
  return myair;
}
//end GetThisAirColumn()


void GetAir(double *col1) {
  double nothing;
  ifstream air1(ICEMC_SRC_DIR+"/data/atmosphere.dat"); // length of chord in air vs. theta (deg)
  //where theta is respect to "up"
  // binned in 0.1 degrees
  for(int iii=0;iii<900;iii++) {
    air1>>nothing>>col1[iii];
  } // read in chord lengths
}
//end GetAir()


int TIR(const Vector &n_surf, const Vector &nrf2_iceside,  double N_IN, double N_OUT) {
  double test=sin(nrf2_iceside.Angle(n_surf))*N_IN/N_OUT;
  if(test>=1)
    return 1;
  else
    return 0;

  return 0;
}
//end TIR()


double GetViewAngle(const Vector &nrf2_iceside, const Vector &nnu) {
  // get viewing angle of shower
  double dtemp=nrf2_iceside*nnu;
  if (dtemp>=1 && dtemp<1.02)
    dtemp=0.999999;
  if (dtemp<=-1 && dtemp>-1.02)
    dtemp=-0.9999999;

  return acos(dtemp);
}
//end ())etViewAngle


TStyle* RootStyle() {
  TStyle *RootStyle = new TStyle("Root-Style", "The Perfect Style for Plots ;-)");

#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                  // save the global style reference
  gStyle = RootStyle;
#endif
  // otherwise you need to call TROOT::SetStyle("Root-Style")

  // Paper size
  RootStyle->SetPaperSize(TStyle::kUSLetter);

  // Canvas
  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderSize(10);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH      (600);
  RootStyle->SetCanvasDefW      (600);
  RootStyle->SetCanvasDefX      (10);
  RootStyle->SetCanvasDefY      (10);

  // Pads
  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderSize  (10);
  RootStyle->SetPadBorderMode  (0);
  //  RootStyle->SetPadBottomMargin(0.13);
  RootStyle->SetPadBottomMargin(0.16);
  RootStyle->SetPadTopMargin   (0.08);
  RootStyle->SetPadLeftMargin  (0.18);
  RootStyle->SetPadRightMargin (0.05);
  RootStyle->SetPadGridX       (0);
  RootStyle->SetPadGridY       (0);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);

  // Frames
  RootStyle->SetFrameFillStyle ( 0);
  RootStyle->SetFrameFillColor ( 0);
  RootStyle->SetFrameLineColor ( 1);
  RootStyle->SetFrameLineStyle ( 0);
  RootStyle->SetFrameLineWidth ( 2);
  RootStyle->SetFrameBorderSize(10);
  RootStyle->SetFrameBorderMode( 0);


  // Histograms
  RootStyle->SetHistFillColor(0);
  RootStyle->SetHistFillStyle(1);
  RootStyle->SetHistLineColor(1);
  RootStyle->SetHistLineStyle(0);
  RootStyle->SetHistLineWidth(2);

  // Functions
  RootStyle->SetFuncColor(1);
  RootStyle->SetFuncStyle(0);
  RootStyle->SetFuncWidth(1);

  //Legends
  RootStyle->SetStatBorderSize(2);
  RootStyle->SetStatFont      (42);
  //  RootStyle->SetOptStat       (111111);
  RootStyle->SetOptStat       (0);
  RootStyle->SetStatColor     (0);
  RootStyle->SetStatX         (0.93);
  RootStyle->SetStatY         (0.90);
  RootStyle->SetStatFontSize  (0.07);
  //  RootStyle->SetStatW         (0.2);
  //  RootStyle->SetStatH         (0.15);

  // Labels,  Ticks,  and Titles
  RootStyle->SetTickLength ( 0.015, "X");
  RootStyle->SetTitleSize  ( 0.10, "X");
  RootStyle->SetTitleOffset( 1.20, "X");
  RootStyle->SetTitleBorderSize(0);
  //  RootStyle->SetTitleFontSize((double)3.);
  RootStyle->SetLabelOffset( 0.015, "X");
  RootStyle->SetLabelSize  ( 0.050, "X");
  RootStyle->SetLabelFont  ( 42   , "X");
  RootStyle->SetTickLength ( 0.015, "Y");
  RootStyle->SetTitleSize  ( 0.10, "Y");
  RootStyle->SetTitleOffset( 0.600, "Y");
  RootStyle->SetLabelOffset( 0.015, "Y");
  RootStyle->SetLabelSize  ( 0.050, "Y");
  RootStyle->SetLabelFont  ( 42   , "Y");
  RootStyle->SetTitleFont  (42, "XY");
  RootStyle->SetTitleColor  (1);

  // Options
  RootStyle->SetOptFit     (0);
  RootStyle->SetMarkerStyle(20);
  RootStyle->SetMarkerSize(0.4);

  //  cout << ">> Style initialized with the Root Style!" << endl;
  //  cout << ">> " << modified << endl << endl;
  return RootStyle;
}
//end RootStyle()




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
