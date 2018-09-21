
#include "EventGenerator.h"

// system includes
#include "signal.h"

// ROOT includes
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TVector3.h"

// libAntarcticaRoot includes
#include "Geoid.h"

// icemc includes
#include "RootOutput.h"
#include "Constants.h"
#include "Settings.h"
#include "Crust2.h"
#include "Tools.h"
#include "Roughness.h"
#include "ANITA.h"
#include "Antarctica.h"
#include "Spectra.h"
#include "AskaryanFreqsGenerator.h"
#include "AskaryanFreqs.h"
#include "secondaries.hh"
#include "RayTracer.h"
#include "Primaries.h"
#include "Taumodel.hh"
#include "screen.hh"
#include "GlobalTrigger.h"
#include "ChanTrigger.h"
#include "SimulatedSignal.h"
#include "GeneratedNeutrino.h"
#include "EnvironmentVariable.h"
#include "Report.h"
#include "FTPair.h"
#include "Neutrino.h"
#include "NeutrinoFactory.h"
#include "DiffuseFlux.h"
#include "RNG.h"

#include <string>
#include <sstream>

#if __cplusplus > 199711L
#define isnan std::isnan 
#include <type_traits>
#endif

#include <typeinfo>
#include <fenv.h> 



bool ABORT_EARLY = false;    // This flag is set to true when interrupt_signal_handler() is called

void icemc::EventGenerator::interrupt_signal_handler(int sig){
  signal(sig,  SIG_IGN);
  ABORT_EARLY = true;
  static int numSignals = 0;
  numSignals++;
  const int maxSignals = 3;
  if(numSignals >= maxSignals){
    std::cerr << "Got " <<  numSignals << " interrupts, quitting more aggressively!" << std::endl;
    raise(SIGINT);
  }
  return;
}




/** 
 * Constructor 
 * 
 * @todo properly zero member variables
 */
icemc::EventGenerator::EventGenerator() : fNeutrinoPath(NULL), interaction1(NULL), fDetector(NULL), fTauPtr(NULL), fGenNu(NULL),  fPassNu(NULL)
{
  pnu = pow(10., 20);   //!< energy of neutrinos
  // inu = 0;

  // This function call allows icemc to gracefully abort and write files as usual rather than stopping abruptly.  
  signal(SIGINT, icemc::EventGenerator::interrupt_signal_handler);
}


/** 
 * Destructor
 * 
 */
icemc::EventGenerator::~EventGenerator()
{
  signal(SIGINT, SIG_DFL); /// unset signal handler
  
  if(fNeutrinoPath){
    delete fNeutrinoPath;
  }
  if(interaction1){
    delete interaction1;
  }
  if(fDetector){
    delete fDetector;
  }
  if(fTauPtr){
    delete fTauPtr;
  }
  if(fGenNu){
    delete fGenNu;
  }
  if(fPassNu){
    delete fPassNu;
  }  
}




void icemc::EventGenerator::Summarize(const Settings *settings1,  Anita* anita1,  Source::Spectra *nuSpectra, const AskaryanFreqsGenerator *askFreqGen, Primaries *primary1, double pnu, double eventsfound, double eventsfound_db, double eventsfound_nfb, double sigma, double* sum, double volume, double ice_area, double& km3sr, double& km3sr_e, double& km3sr_mu, double& km3sr_tau, TString outputdir) {

  double rate_v_thresh[NTHRESHOLDS];
  double errorup_v_thresh[NTHRESHOLDS];
  double errordown_v_thresh[NTHRESHOLDS];
  double rate_h_thresh[NTHRESHOLDS];
  double errorup_h_thresh[NTHRESHOLDS];
  double errordown_h_thresh[NTHRESHOLDS];
  double zeroes[NTHRESHOLDS];
  Tools::Zero(zeroes, NTHRESHOLDS);

  // plot result of threshold scan
  for (int i=0;i<NTHRESHOLDS;i++) {
    rate_v_thresh[i]=npass_v_thresh[i]/denom_v_thresh[i];
    // std::cout << "i,  npass_v_thresh are " << i << "\t" << npass_v_thresh[i] << " " << denom_v_thresh[i] << " " << rate_v_thresh[i] << "\n";
    if (npass_v_thresh[i]<=20) {
      // errorup_v_thresh[i]=poissonerror_plus[(int)npass_v_thresh[i]]/denom_v_thresh[i];
      // errordown_v_thresh[i]=poissonerror_minus[(int)npass_v_thresh[i]]/denom_v_thresh[i];
      const int nv = (int)npass_v_thresh[i];
      double err_plus=0, err_minus = 0;
      icemc::constants::getPoissonError(nv, err_plus, err_minus);
      errorup_v_thresh[i]=err_plus/denom_v_thresh[i];
      errordown_v_thresh[i]=err_minus/denom_v_thresh[i];
      // std::cout << errorup_v_thresh[i] << " " << errordown_v_thresh[i] << std::endl;
    }//end if
    if (settings1->WHICH==Payload::Anita3 || settings1->WHICH==Payload::Anita4){ // Anita-3
      rate_h_thresh[i]=npass_h_thresh[i]/denom_h_thresh[i];
      if (npass_h_thresh[i]<=20) {
	const int nh = (int)npass_h_thresh[i];
	double err_plus = 0, err_minus = 0;
	icemc::constants::getPoissonError(nh, err_plus, err_minus);
        errorup_h_thresh[i]=err_plus/denom_h_thresh[i];
        errordown_h_thresh[i]=err_minus/denom_h_thresh[i];
        // errorup_h_thresh[i]=poissonerror_plus[(int)npass_h_thresh[i]]/denom_h_thresh[i];
        // errordown_h_thresh[i]=poissonerror_minus[(int)npass_h_thresh[i]]/denom_h_thresh[i];
      }//end if
    }//end if WHICH==9

  }//end for NTHRESHOLDS

  std::string stemp=std::string(outputdir.Data())+"/thresholds.root";
  TFile *fthresholds=new TFile(stemp.c_str(), "RECREATE");
  TCanvas *cthresh=new TCanvas("cthresh", "cthresh", 880, 800);
  cthresh->SetLogy();

  TGraph *gnpass=new TGraph(NTHRESHOLDS, thresholds, npass_v_thresh);
  gnpass->SetName("npass");
  TGraph *gdenom=new TGraph(NTHRESHOLDS, thresholds, denom_v_thresh);
  gdenom->SetName("denom");

  TGraphAsymmErrors* g = new TGraphAsymmErrors(NTHRESHOLDS, thresholds, rate_v_thresh, zeroes, zeroes, errorup_v_thresh, errordown_v_thresh);
  g->SetName("rate");

  g->SetLineWidth(2);
  g->SetMarkerStyle(21);
  g->Draw("ape");

  stemp = std::string(outputdir.Data())+"/thresholds.eps";
  cthresh->Print((TString)stemp);
  g->Write();
  gdenom->Write();
  gnpass->Write();


  if (settings1->WHICH==Payload::Anita3 || settings1->WHICH==Payload::Anita4){ // Anita-3 or Anita-4
    TGraph *gnpassH=new TGraph(NTHRESHOLDS, thresholds, npass_h_thresh);
    gnpassH->SetName("npassH");
    TGraph *gdenomH=new TGraph(NTHRESHOLDS, thresholds, denom_h_thresh);
    gdenomH->SetName("denomH");
    TGraphAsymmErrors *gH=new TGraphAsymmErrors(NTHRESHOLDS, thresholds, rate_h_thresh, zeroes, zeroes, errorup_h_thresh, errordown_h_thresh);
    gH->SetName("rate");
    gH->SetLineWidth(2);
    gH->SetMarkerStyle(21);
    gH->Draw("ape");
    stemp = std::string(outputdir.Data())+"/thresholds_HPOL.eps";
    cthresh->Print((TString)stemp);

    gH->Write();
    gdenomH->Write();
    gnpassH->Write();
  }//end if WHICH==9 or WHICH==10

  fthresholds->Write();
  fthresholds->Close();

  //  double ses;                          // single-event sensitivity
  double km2sr;                        // aperture km**2-sr

  icemc::report() << "\n";

  // write out summary
  icemc::report() << "Generated " << NNU << " neutrinos with energy " << pnu << "\n";
  icemc::report() << "Number of (weighted) neutrinos that pass is: " << eventsfound << "\n";
  icemc::report() << "Number of (weighted) neutrinos that pass,  multiplied by prob. of interacting in the ice,  is: " << eventsfound_prob << "\n";
  icemc::report() << "Number of weighted direct,  reflected that pass is: " << allcuts_weighted[0] << "\t" << allcuts_weighted[1] << "\n";
  icemc::report() << "Number of (weighted) neutrinos that pass (with weight>0.001) is: " << eventsfound_weightgt01 << "\n";
  icemc::report() << "Number of (weighted) neutrinos that only traverse the crust is " << eventsfound_crust << " -> " << eventsfound_crust/eventsfound*100 << "%\n\n";
  icemc::report() << "Volume of ice is " << volume << "\n";
  icemc::report() << "Value of 4*pi*pi*r_earth*r_earth in km^2 " << 4*constants::PI*constants::PI*(Geoid::R_EARTH*Geoid::R_EARTH/1.E6) << "\n";


  // single event sensitivity for 150 MHz array per year
  //  everything is normalized to the active volume
  //  assume dF/dE propto E**-2
  //  assume 2pi sr or 4pi sr depending on whether below horizon is used
  //  assume 3.16e7 seconds/year
  //------------------------------------------------

  // ses=1/(A_eff * sr * t)
  //ses=1.0/(sigma*volume*RHOSALT*(1./M_NUCL)*3.16E7);  // 1 per year


  double nevents = 0;
  //double nweights=0;
  double nevents_db=0;
  double nevents_nfb=0;


  nevents+=eventsfound;
  nevents_db+=eventsfound_db;
  nevents_nfb+=eventsfound_nfb;
  error_plus=0;
  error_e_plus=0;
  error_mu_plus=0;
  error_tau_plus=0;
  error_minus=0;
  error_e_minus=0;
  error_mu_minus=0;
  error_tau_minus=0;

  error_nfb=0;
  error_km3sr_nfb=0;
  error_percent_increase_nfb=0;


  // loop over bins in weights
  for (int i=0;i<NBINS;i++) {
    if (eventsfound_binned[i]<=20) {

      // What used to be here...
      // error_plus+=pow(poissonerror_plus[(int)eventsfound_binned[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      // error_e_plus+=pow(poissonerror_plus[(int)eventsfound_binned_e[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      // error_mu_plus+=pow(poissonerror_plus[(int)eventsfound_binned_mu[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      // error_tau_plus+=pow(poissonerror_plus[(int)eventsfound_binned_tau[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      // error_minus+=pow(poissonerror_minus[(int)eventsfound_binned[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      // error_e_minus+=pow(poissonerror_minus[(int)eventsfound_binned_e[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      // error_mu_minus+=pow(poissonerror_minus[(int)eventsfound_binned_mu[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      // error_tau_minus+=pow(poissonerror_minus[(int)eventsfound_binned_tau[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      
      const int n = (int)eventsfound_binned[i];
      const int n_e = (int)eventsfound_binned_e[i];
      const int n_mu = (int)eventsfound_binned_mu[i];
      const int n_tau = (int)eventsfound_binned_tau[i];
      double err_plus = 0, err_minus = 0, err_e_plus = 0, err_e_minus = 0, err_mu_plus = 0, err_mu_minus = 0, err_tau_plus = 0, err_tau_minus = 0;
      icemc::constants::getPoissonError(n,     err_plus,     err_minus);
      icemc::constants::getPoissonError(n_e,   err_e_plus,   err_e_minus);
      icemc::constants::getPoissonError(n_mu,  err_mu_plus,  err_mu_minus);
      icemc::constants::getPoissonError(n_tau, err_tau_plus, err_tau_minus);
      error_plus+=pow(err_plus*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_e_plus+=pow(err_e_plus*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_mu_plus+=pow(err_mu_plus*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_tau_plus+=pow(err_tau_plus*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_minus+=pow(err_minus*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_e_minus+=pow(err_e_minus*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_mu_minus+=pow(err_mu_minus*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_tau_minus+=pow(err_tau_minus*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
    }
    else {
      error_plus     += eventsfound_binned[i]*pow(pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_e_plus   += eventsfound_binned_e[i]*pow(pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_mu_plus  += eventsfound_binned_mu[i]*pow(pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_tau_plus += eventsfound_binned_tau[i]*pow(pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_minus     = error_plus;
      error_e_minus   = error_e_plus;
      error_mu_minus  = error_mu_plus;
      error_tau_minus = error_tau_plus;
    }

    error_nfb += eventsfound_nfb_binned[i]*pow(pow(10., ((double)i+0.5)/(double)NBINS*4.+-5.), 2);
    for (int j=0;j<NBINS_DISTANCE;j++) {
      if (eventsfound_binned_distance_forerror[j][i]<=20) {

	// What was here before...
        // error_distance_plus[j]+=pow(poissonerror_plus[(int)eventsfound_binned_distance_forerror[j][i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
        // error_distance_minus[j]+=pow(poissonerror_minus[(int)eventsfound_binned_distance_forerror[j][i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
		
	double poissonErrPlus = 0, poissonErrMinus = 0;
	const int n = (int)eventsfound_binned_distance_forerror[j][i];
	icemc::constants::getPoissonError(n, poissonErrPlus, poissonErrMinus);
        error_distance_plus[j]+=pow(poissonErrPlus*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
        error_distance_minus[j]+=pow(poissonErrMinus*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);

      }
    }//end for NBINS_DISTANCE
  }//end for NBINS


  error_plus=sqrt(error_plus);
  error_e_plus=sqrt(error_e_plus);
  error_mu_plus=sqrt(error_mu_plus);
  error_tau_plus=sqrt(error_tau_plus);

  error_minus=sqrt(error_minus);
  error_e_minus=sqrt(error_e_minus);
  error_mu_minus=sqrt(error_mu_minus);
  error_tau_minus=sqrt(error_tau_minus);


  error_nfb=sqrt(error_nfb);
  for (int j=0;j<NBINS_DISTANCE;j++) {
    error_distance_plus[j]=sqrt(error_distance_plus[j]);
    error_distance_minus[j]=sqrt(error_distance_minus[j]);
  }


  // account for efficiency
  if (NNU != 0 && nevents!=0) {
    km3sr=volume*pow(1.E-3, 3)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr*nevents/(double)NNU/settings1->SIGMA_FACTOR;

    icemc::report() << nevents << " events passed out of " << NNU << "\n";

    error_plus=volume*pow(1.E-3, 3)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr*error_plus/(double)NNU/settings1->SIGMA_FACTOR;
    error_minus=volume*pow(1.E-3, 3)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr*error_minus/(double)NNU/settings1->SIGMA_FACTOR;

    percent_increase_db=km3sr_db/km3sr*100;
    percent_increase_nfb=km3sr_nfb/km3sr*100;

    error_percent_increase_nfb=sqrt(pow(error_km3sr_nfb/km3sr_nfb, 2)+pow(error_plus/km3sr, 2))*percent_increase_nfb;

    percent_increase_total=percent_increase_db+percent_increase_nfb;

    // km3sr_e = (pow(1.e-3, 3))*volume*constants::sr*(sum[0]/(double)count1->nnu_e)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20/settings1->SIGMA_FACTOR;

    icemc::report() << sum[0]/(double)nevents*100. << "% are electron neutrinos\n";

    // error_e_plus=volume*pow(1.E-3, 3)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr*error_e_plus/(double)count1->nnu_e/settings1->SIGMA_FACTOR;
    // error_e_minus=volume*pow(1.E-3, 3)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr*error_e_minus/(double)count1->nnu_e/settings1->SIGMA_FACTOR;

    // km3sr_mu = (pow(1.e-3, 3))*volume*constants::sr*(sum[1]/(double)count1->nnu_mu)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20/settings1->SIGMA_FACTOR;

    icemc::report() << sum[1]/(double)nevents*100. << "% are muon neutrinos\n";

    // error_mu_plus=volume*pow(1.E-3, 3)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr*error_mu_plus/(double)count1->nnu_mu/settings1->SIGMA_FACTOR;
    // error_mu_minus=volume*pow(1.E-3, 3)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr*error_mu_minus/(double)count1->nnu_mu/settings1->SIGMA_FACTOR;

    // km3sr_tau = (pow(1.e-3, 3))*volume*constants::sr*(sum[2]/(double)count1->nnu_tau)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20/settings1->SIGMA_FACTOR;

    icemc::report() << sum[2]/(double)nevents*100. << "% are tau neutrinos\n";

    // error_tau_plus=volume*pow(1.E-3, 3)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr*error_tau_plus/(double)count1->nnu_tau/settings1->SIGMA_FACTOR;
    // error_tau_minus=volume*pow(1.E-3, 3)*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr*error_tau_minus/(double)count1->nnu_tau/settings1->SIGMA_FACTOR;

    double sum_km3sr=0;
    double sum_km3sr_error_plus=0;
    double sum_km3sr_error_minus=0;
    for (int i=0;i<NBINS_DISTANCE;i++) {
      sum_km3sr+= (pow(1.e-3, 3))*volume*constants::sr*eventsfound_binned_distance[i]*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20/(double)NNU/settings1->SIGMA_FACTOR;
      km3sr_distance[i]=sum_km3sr;
      sum_km3sr_error_plus+=pow(pow(1.E-3, 3)*volume*error_distance_plus[i]*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr/(double)NNU, 2)/settings1->SIGMA_FACTOR;
      error_distance_plus[i]=sqrt(sum_km3sr_error_plus);
      sum_km3sr_error_minus+=pow(pow(1.E-3, 3)*volume*error_distance_minus[i]*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*constants::sr/(double)NNU, 2)/settings1->SIGMA_FACTOR;
      error_distance_minus[i]=sqrt(sum_km3sr_error_minus);
      icemc::report().distanceout << 700000./(double)NBINS_DISTANCE*(double)i << "\t" << km3sr_distance[i] << "\t" << error_distance_plus[i] << "\t" << error_distance_minus[i] << "\n";
    } //for

    icemc::report() << "Total volume * solid angle is \t\t\t\t" << km3sr << " + " << error_plus << " - " << error_minus << " km^3 str\n";
    icemc::report() << "Total volume * solid angle for electron neutrinos is \t" << km3sr_e << " + " << error_e_plus << " - " << error_e_minus << " km^3 str\n";
    icemc::report() << "Total volume * solid angle for muon neutrinos is \t" << km3sr_mu << " + " << error_mu_plus << " - " << error_mu_minus << " km^3 str\n";
    icemc::report() << "Total volume * solid angle for tau neutrinos is \t" << km3sr_tau << " + " << error_tau_plus << " - " << error_tau_minus << " km^3 str\n";

    if ( nuSpectra->IsMonoenergetic() )  {
      std::cout << "Cross section is " << sigma << "m^2\n";
      double len_int=1.0/(sigma*askFreqGen->RHOH20*(1./constants::M_NUCL)*1000);
      icemc::report() << "Interaction length is " << len_int << "\n";
      km2sr=km3sr/len_int;
      icemc::report() << "Total area x steradians using km3sr/len_int is \t\t\t\t" << km2sr << " km^2 str\n\n";
      km2sr=ice_area/(1.E6)*constants::PI*eventsfound_prob/(double)NNU;
      icemc::report() << "Total area x steradians using 4*PI*R_EARTH^2*eff. is \t" << km2sr << " km^2 str\n\n";
      icemc::report() << "These are not the same because we are not throwing all directions on all points of the surface.  Believe the first one as an approximation,  we are working on this for high cross sections.\n";
      // ses=(pnu/1.E9)/(km2sr*3.16E7);
    }//end if IsMonoenergetic

  }//end if NNU!=0 and nevents!=0


  icemc::report() << "\t\t\t\t\t\t\tprobability of passing \t# events passing\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "No way this neutrino will see any ice \t\t\t" << (double)count1->noway[0]/(double)count_total << "\t" <<
  //   (double)count1->noway[1]/(double)count_total << "\t\t" <<
  //   count1->noway[0] << "\t" << count1->noway[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Wheredoesitleave in PickUnbiased gives an error \t" << (double)count1->wheredoesitleave_err[0]/(double)count1->noway[0] << "\t" <<
  //   (double)count1->wheredoesitleave_err[1]/(double)count1->noway[1] << "\t\t" <<
  //   count1->wheredoesitleave_err[0] << "\t" << count1->wheredoesitleave_err[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "This neutrino direction never sees ice \t\t\t" << (double)count1->neverseesice[0]/(double)count1->wheredoesitleave_err[0] << "\t" <<
  //   (double)count1->neverseesice[1]/(double)count1->wheredoesitleave_err[1] << "\t\t" <<
  //   count1->neverseesice[0] << "\t" << count1->neverseesice[1] << "\n";


  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "WhereDoesItEnterIce in PickUnbiased gives an error \t\t\t" << (double)count1->wheredoesitenterice_err[0]/(double)count1->neverseesice[0] << "\t" <<
  //   (double)count1->wheredoesitenterice_err[1]/(double)count1->neverseesice[1] << "\t\t" <<
  //   count1->wheredoesitenterice_err[0] << "\t" << count1->wheredoesitenterice_err[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Interaction point too high \t\t\t\t" << (double)count1->toohigh[0]/(double)count1->wheredoesitenterice_err[0] << "\t" <<
  //   (double)count1->toohigh[1]/(double)count1->wheredoesitenterice_err[1] << "\t\t" <<
  //   count1->toohigh[0] << "\t" << count1->toohigh[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Interaction point too low \t\t\t\t" << (double)count1->toolow[0]/(double)count1->toohigh[0] << "\t" <<
  //   (double)count1->toolow[1]/(double)count1->toohigh[1] << "\t\t" <<
  //   count1->toolow[0] << "\t" << count1->toolow[1] << "\n";




  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "There is an interaction in ice \t\t\t\t" << (double)count1->iceinteraction[0]/(double)count1->toolow[0] << "\t" <<
  //   (double)count1->iceinteraction[1]/(double)count1->toolow[1] << "\t\t" <<
  //   count1->iceinteraction[0] << "\t" << count1->iceinteraction[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "In horizon \t\t\t\t\t\t" << (double)count1->inhorizon[0]/(double)count1->iceinteraction[0] << "\t" <<
  //   (double)count1->inhorizon[1]/(double)count1->iceinteraction[1] << "\t\t" <<
  //   count1->inhorizon[0] << "\t" << count1->inhorizon[1] << "\n";



  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "From surface to balloon,  ray not intersected by earth \t" << (double)count1->nraypointsup1[0]/(double)count1->inhorizon[0] << "\t" <<
  //   (double)count1->nraypointsup1[1]/(double)count1->inhorizon[1] << "\t\t" <<
  //   count1->nraypointsup1[0] << "\t" << count1->nraypointsup1[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "After 1/r scaling and best case attenuation, \n\tMaximum signal is detectable\t\t\t" << (double)count1->nnottoosmall[0]/(double)count1->nraypointsup1[0] << "\t" <<
  //   (double)count1->nnottoosmall[1]/(double)count1->nraypointsup1[1] << "\t\t"
  // 	  << count1->nnottoosmall[0] << "\t" << count1->nnottoosmall[1] << "\n";


  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Viewing angle lt 90 degrees\t\t\t" << (double)count1->nviewangle_lt_90[0]/(double)count1->nnottoosmall[0] << "\t" <<
  //   (double)count1->nviewangle_lt_90[1]/(double)count1->nnottoosmall[1] << "\t\t"
  // 	  << count1->nviewangle_lt_90[0] << "\t" << count1->nviewangle_lt_90[1] << "\n";


  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Reality check:  EM and hadronic fractions both nonzero\t" << (double)count1->ngoodfracs[0]/(double)count1->nviewangle_lt_90[0] << "\t" <<
  //   (double)count1->ngoodfracs[1]/(double)count1->nviewangle_lt_90[1] << "\t\t"
  // 	  << count1->ngoodfracs[0] << "\t" << count1->ngoodfracs[1] << "\n";
  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "\tBoth EM and hadronic fractions are zero\t\t" << (double)count1->nbadfracs[0]/(double)count1->nviewangle_lt_90[0] << "\t" <<
  //   (double)count1->nbadfracs[1]/(double)count1->nviewangle_lt_90[1] << "\t\t" <<
  //   count1->nbadfracs[0] << "\t" << count1->nbadfracs[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "After finding neutrino direction,  \n\tchance of making through Earth\t\t\t" << (double)count_chanceofsurviving/(double)count1->ngoodfracs[0] << "\t\t\t";
  // icemc::report().foutput.precision(10);
  // icemc::report().foutput << count_chanceofsurviving << "\n";
  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Neutrino enters ice south of 60deg S latitude\t\t" << (double)count1->nentersice[0]/(double)count_chanceofsurviving << "\t" <<
  //   (double)count1->nentersice[1]/(double)count_chanceofsurviving <<
  //   "\t\t" <<
  //   count1->nentersice[0] << "\t" << count1->nentersice[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Neutrino reasonably likely to survive trip through Earth " << (double)count1->nabsorbed[0]/(double)count1->nentersice[0] << "\t" <<
  //   (double)count1->nabsorbed[1]/(double)count1->nentersice[1] << "\t\t"
  // 	  << count1->nabsorbed[0] << "\t" << count1->nabsorbed[1] << "\n";
  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Ray leaves the ice south of 60deg S latitude\t\t" << (double)count1->nraywithincontinent1[0]/(double)count1->nabsorbed[0] << "\t" <<
  //   (double)count1->nraywithincontinent1[1]/(double)count1->nabsorbed[1] << "\t" <<
  //   count1->nraywithincontinent1[0] << "\t" <<
  //   count1->nraywithincontinent1[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "After 1/r,  best guess ice attenuation,  \n\tmaximum signal is detectable\t\t\t" << (double)count_chanceinhell0/(double)count1->nraywithincontinent1[0] << "\t\t\t";
  // icemc::report().foutput.precision(10);
  // icemc::report().foutput <<count_chanceinhell0 << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Ray is not totally internally reflected\t\t\t" << (double)count1->nnottir[0]/(double)count_chanceinhell0 <<  "\t" <<
  //   (double)count1->nnottir[1]/(double)count_chanceinhell0 <<  "\t\t" <<
  //   count1->nnottir[0] << "\t" << count1->nnottir[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "From surface to balloon,  ray not intersected by earth \t" << (double)count1->nraypointsup2[0]/(double)count1->nnottir[0] << "\t" <<
  //   (double)count1->nraypointsup2[1]/(double)count1->nnottir[1] <<
  //   "\t\t"
  // 	  << count1->nraypointsup2[0] << "\t" << count1->nraypointsup2[1] << "\n";
  // icemc::report().foutput.precision(4);

  // icemc::report().foutput << "Ray leaves the ice south of 60deg S latitude\t\t" << (double)count1->nraywithincontinent2[0]/(double)count1->nraypointsup2[0] <<"\t" <<
  //   (double)count1->nraywithincontinent2[0]/(double)count1->nraypointsup2[1] <<
  //   "\t\t" << count1->nraywithincontinent2[0] << "\t" <<
  //   count1->nraywithincontinent2[1]     << "\n";
  // icemc::report().foutput.precision(4);

  // icemc::report().foutput << "Ray leaves where there is ice\t\t\t\t" << (double)count1->nacceptablerf[0]/(double)count1->nraywithincontinent2[0] << "\t" <<
  //   (double)count1->nacceptablerf[1]/(double)count1->nraywithincontinent2[1] << "\t\t"
  // 	  << count1->nacceptablerf[0] << "\t" << count1->nacceptablerf[1] << "\n";
  // icemc::report().foutput.precision(4);

  // icemc::report().foutput << "Ray tracing converges to within 10 m\t\t\t" << (double)count1->nconverges[0]/(double)count1->nacceptablerf[0] << "\t" <<
  //   (double)count1->nconverges[1]/(double)count1->nacceptablerf[1] <<
  //   "\t\t" << count1->nconverges[0] << "\t" << count1->nconverges[1] << "\n";
  // icemc::report().foutput.precision(4);

  // icemc::report().foutput << "After fresnel coefficient,  \n\tmaximum signal is detectable\t\t\t" << (double)count1->nchanceinhell_fresnel[0]/(double)count1->nconverges[0] << "\t" << (double)count1->nchanceinhell_fresnel[1]/(double)count1->nconverges[1] <<
  //   "\t\t" <<count1->nchanceinhell_fresnel[0] << "\t" << count1->nchanceinhell_fresnel[1] << "\n";
  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "After 1/r,  \n\tmaximum signal is detectable\t\t\t" << (double)count1->nchanceinhell_1overr[0]/(double)count1->nchanceinhell_fresnel[0] << "\t" <<  (double)count1->nchanceinhell_1overr[1]/(double)count1->nchanceinhell_fresnel[1] << "\t\t" <<count1->nchanceinhell_1overr[0] << "\t" << count1->nchanceinhell_1overr[1] << "\n";
  // icemc::report().foutput.precision(4);

  // icemc::report().foutput << "After ice attenuation,  \n\tmaximum signal is detectable\t\t\t" << (double)count1->nchanceinhell[0]/(double)count1->nchanceinhell_1overr[0] << "\t" <<
  //   (double)count1->nchanceinhell[1]/(double)count1->nchanceinhell_1overr[1] << "\t" <<count1->nchanceinhell[0] << "\t" << count1->nchanceinhell[1] << "\n";
  // icemc::report().foutput.precision(4);

  // icemc::report().foutput << "After viewing angle cut, \t\t\t\t" << (double)count1->nviewanglecut[0]/(double)count1->nchanceinhell[0] << "\t" << (double)count1->nviewanglecut[1]/(double)count1->nchanceinhell[1] << "\t\t" << count1->nviewanglecut[0] << " " << count1->nviewanglecut[1] << "\n";

  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "After factoring in off-Cerenkov cone tapering, \n\tMaximum signal is detectable\t\t\t" << (double)count1->nchanceinhell2[0]/(double)count1->nviewanglecut[0] << "\t" << (double)count1->nchanceinhell2[1]/(double)count1->nviewanglecut[1] << "\t\t" << count1->nchanceinhell2[0] << " " << count1->nchanceinhell2[1] << "\n";

  // icemc::report().foutput << "Survive dead time \t\t\t\t\t" << (double)count1->ndeadtime[0]/(double)count1->nchanceinhell2[0] << "\t" << (double)count1->ndeadtime[1]/(double)count1->nchanceinhell2[1] << "\t\t" << (double)count1->ndeadtime[0] << " " << count1->ndeadtime[1] << "\n";
  
  // icemc::report().foutput << "Passes trigger\t\t\t\t\t\t" << (double)count1->npassestrigger[0]/(double)count1->ndeadtime[0] << "\t" << (double)count1->npassestrigger[1]/(double)count1->ndeadtime[1] << "\t\t" << count1->npassestrigger[0] << "\t" << count1->npassestrigger[1] << "\n";
  // icemc::report().foutput << "Number of l1 triggers\t\t\t\t\t\t" << (double)count1->nl1triggers[0][0] << "\t" << (double)count1->nl1triggers[1][0] << "\n";

  // icemc::report().foutput << "Chord is good length\t\t\t\t\t" << (double)count_chordgoodlength/(double)count1->npassestrigger[0] << "\t\t\t";
  // icemc::report().foutput.precision(10);
  // icemc::report().foutput <<count_chordgoodlength << "\n";
  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Neutrino's path in ice more than 1m \t\t\t" << (double)count_d2goodlength/(double)count_chordgoodlength << "\t\t\t";
  // icemc::report().foutput.precision(10);
  // icemc::report().foutput << count_d2goodlength << "\n";
  // icemc::report().foutput.precision(4);
  // icemc::report().foutput << "Events that pass all cuts\t\t\t\t" << (double)count1->npass[0]/(double)count_d2goodlength << "\t" << (double)count1->npass[1]/(double)count_d2goodlength << "\t\t";
  // icemc::report().foutput.precision(10);
  // icemc::report().foutput <<count1->npass[0] << "\t" << count1->npass[1] << "\n";

  // std::cout << "Events that pass all cuts\t\t\t\t" << (double)count1->npass[0]/(double)count_d2goodlength << "\t" << (double)count1->npass[1]/(double)count_d2goodlength << "\t\t";
  // std::cout <<count1->npass[0] << "\t" << count1->npass[1] << "\n";

  //  if (EXPONENT<=10||EXPONENT>100) {
  if ( nuSpectra->IsSpectrum() ) {
    double sum_events=0.;
    double thisenergy=0.;
    double thislen_int_kgm2=0.;
    // double *energy=nuSpectra->GetEnergyArray();
    // double *EdNdEdAdt=nuSpectra->GetEdNdEdAdt();

    // for models which don't have even spaced energy bin,
    double even_E;
    int N_even_E = 12;
    double integral=0;
    even_E = ( nuSpectra->Getenergy()[nuSpectra->GetE_bin() - 1] - nuSpectra->Getenergy()[0] ) / ( (double) N_even_E );
    for (int i=0;i<N_even_E;i++) {
      thisenergy=pow(10., (nuSpectra->Getenergy())[0]+((double)i)*even_E);
      primary1->GetSigma(thisenergy, sigma, thislen_int_kgm2, settings1, xsecParam_nutype, xsecParam_nuint);

      // EdNdEdAdt is in #/cm^2/s
      // need to be multiplied by 1e4 to change 1/cm2 to 1/m^2
      // can also be written dN/d(lnE)dAdt
      // = dN*log(10)/d(icemc::report() E)dAdt
      // the bin spacing is 0.5
      // so # events ~ dN*log(10)*0.5/d(icemc::report() E)dAdt
      sum_events+=even_E*log(10.)*( nuSpectra->GetEdNdEdAdt(log10(thisenergy))*1e4 )/(thislen_int_kgm2/askFreqGen->RHOH20);
      integral+=even_E*log(10.)*( nuSpectra->GetEdNdEdAdt(log10(thisenergy)) );
      std::cout << "thisenergy,  EdNdEdAdt is " << thisenergy << " " <<  nuSpectra->GetEdNdEdAdt(log10(thisenergy)) << "\n";
      //foutput << "interaction length is " << thislen_int_kgm2/RHOH20 << "\n";
    }//end for N_even_E
    std::cout << "SUM EVENTS IS " << sum_events << std::endl;
    std::cout << "INTEGRAL : " << integral << std::endl;
    sum_events*=volume*anita1->LIVETIME*askFreqGen->RHOMEDIUM/askFreqGen->RHOH20*nevents/(double)NNU*constants::sr;
    // sum_events*=anita1->LIVETIME*km3sr*1e9;
    icemc::report().foutput << "volume,  LIVETIME,  askFreqGen->RHOMEDIUM,  RHOH20,  nevents,  NNU,  sr are " << volume << " " << anita1->LIVETIME << " " << askFreqGen->RHOMEDIUM << " " << askFreqGen->RHOH20 << " " << nevents << " " << NNU << " " << constants::sr << "\n";
    icemc::report().foutput << "Total events observed is " << sum_events << "\n";
  } //end if IsSpectrum

}
//end Summarize()







TVector3 icemc::EventGenerator::GetPolarization(const TVector3 &nnu, const TVector3 &nrf2_iceside, int inu) const {
  // Want to find a unit vector in the same plane as nnu and n_refr,
  // but perpendicular to n_refr,  pointing away from nnu.

  // cross nnu with n_refr to get the direction of the B field.
  TVector3 n_bfield = nnu.Cross(nrf2_iceside);
  // cross b-field with nrf2_iceside to get the polarization vector.
  TVector3 n_pol = n_bfield.Cross(nrf2_iceside);
  n_pol = n_pol.Unit();
  // check and make sure E-field is pointing in the right direction.
  if (nnu.Dot(nrf2_iceside)>0 && n_pol.Dot(nnu)>0){
    icemc::report() << severity::error << "In GetPolarization.  inu = " << inu << std::endl;
  }
  return n_pol;
}
//end GetPolarization()



double icemc::EventGenerator::IsItDoubleBang(double exitlength, double plepton) const {
  double gamma = plepton/constants::MTAU;
  return 1-exp(-1*exitlength/(constants::TAUDECAY_TIME*constants::CLIGHT*gamma));
}
//end IsItDoubleBang()


int icemc::EventGenerator::WhereIsSecondBang(const Geoid::Position &posnu, const TVector3 &nnu, double nuexitlength, double pnu, Antarctica *antarctica1, const Geoid::Position &r_bn, Geoid::Position &posnu2, Geoid::Position &rfexit_db, TVector3 &n_exit2bn_db) const {
  double rnd1=0;
  double rnd2=2;
  double gamma=pnu/constants::MTAU;

  if (exp(-1*nuexitlength/(constants::TAUDECAY_TIME*constants::CLIGHT*gamma))>0.999){
    rnd1=gRandom->Rndm()*nuexitlength;
  }
  else {
    while (rnd2>1-exp(-1*rnd1/(constants::TAUDECAY_TIME*constants::CLIGHT*gamma))) {
      rnd1=gRandom->Rndm()*nuexitlength;
      rnd2=gRandom->Rndm();
    } //while
  } //else
  posnu2 = posnu + rnd1*nnu;
  rfexit_db = antarctica1->WorldModel::Surface(posnu2)*posnu2.Unit();

  // unit vector pointing to antenna from exit point.
  n_exit2bn_db = (r_bn - rfexit_db).Unit();

  double cosangle=(n_exit2bn_db.Dot(posnu2)) / posnu2.Mag();
  if (cosangle<0){
    return 0;
  }
  return 1;
}
//end WhereIsSecondBang()


//the following is  a new function only for reflected case.
void icemc::EventGenerator::Attenuate_down(Antarctica *antarctica1, const Settings *settings1, double& vmmhz_max, const Geoid::Position &rfexit2, const Geoid::Position &posnu, const Geoid::Position &posnu_down) const {
  double ATTENLENGTH=700;
  if(!settings1->VARIABLE_ATTEN){
    ATTENLENGTH=antarctica1->EffectiveAttenuationLength(settings1, posnu, 1);
  }

  int position_in_iceshelves=antarctica1->IceOnWater(posnu);
  int position_in_rossexcept=antarctica1->RossExcept(posnu);
  // int position_in_ross = antarctica->RossIceShelf(posnu);
  // int position_in_ronne = antarctica->RonneIceShelf(posnu);
  double dtemp=posnu_down.Distance(rfexit2)/ATTENLENGTH;
  double scalefactor_attenuation = 0;
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


void icemc::EventGenerator::Attenuate(const Antarctica *antarctica1, const Settings *settings1, double& vmmhz_max,  double rflength, const Geoid::Position &posnu) const {
  double ATTENLENGTH=700;  // constant attenuation length for now.
  if (!settings1->VARIABLE_ATTEN){
    ATTENLENGTH = antarctica1->EffectiveAttenuationLength(settings1, posnu, 0);
  }

  double dtemp=(rflength/ATTENLENGTH);
  double scalefactor_attenuation = 0;
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


void icemc::EventGenerator::IsAbsorbed(double chord_kgm2, double len_int_kgm2, double &weight1) const {
  // see if neutrino is absorbed
  //  weighting works,  but not to much purpose since nu's always
  //   interact at these energies.
  double rtemp;

  rtemp=chord_kgm2/len_int_kgm2;
  if (rtemp<=20){
    weight1=exp(-rtemp);
  }
  else{
    weight1=0;
  }
}
//end IsAbsorbed()


void icemc::EventGenerator::GetSmearedIncidentAngle(TVector3 &specular, TVector3 &nrf_iceside, TVector3 &n_exit2bn, double SMEARINCIDENTANGLE) const{
  // Smear the incident angle for roughness studies
  specular+=nrf_iceside; // specular is the ray that we got from Snell's law
  TVector3 parallel_to_surface; // find vector parallel to surface to rotate the vector around
  parallel_to_surface+=n_exit2bn; // want to cross specular with n_exit2bn
  parallel_to_surface.Cross(specular);
  nrf_iceside.Rotate(SMEARINCIDENTANGLE*(2*gRandom->Rndm()-1.), parallel_to_surface); // smear the incident ray
  //   theta_inc_smeared=acos(nrf_iceside.Dot(nsurf_rfexit));
}
//end GetSmearedIncidentAngle()


int icemc::EventGenerator::GetRayIceSide(const TVector3 &n_exit2rx,  const TVector3 &nsurf_rfexit, double nexit,  double nenter,  TVector3 &nrf2_iceside) const {
  // this function performs snell's law in three dimensions
  double costh=0;
  double NRATIO=nexit/nenter;
  costh=(n_exit2rx.Dot(nsurf_rfexit))/(n_exit2rx.Mag() * nsurf_rfexit.Mag()); // cos(theta) of the transmission angle

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




int icemc::EventGenerator::GetDirection(const Settings *settings1, Interaction *interaction1, const TVector3 &refr,
					double deltheta_em,  double deltheta_had, const ShowerProperties& showerProps,
					double vmmhz1m_max,  double r_fromballoon,  RayTracer *ray1,
					const AskaryanFreqsGenerator *askFreqGen,  Geoid::Position posnu,  Anita *anita1,
					Balloon *bn1, TVector3 &nnu,  double& costhetanu,  double& theta_threshold) { 

  // In the specular (settings1->ROUGHNESS = 0) this function sets the neutrino direction according to a selection routine based on viewing within the Cerenkov cone
  // In the roughness case we just want to pick a random allowable direction,
  // so let's keep the original sampled neutrino direction from back in Antarctica::PickUnbiased() inside Ray::PickRoughnessInteractionPoint()

  int dont_count = 0;
  double theta_test = 0;
  double vmmhz1m_test = 0;
  double costhetanu1 = 0;
  double costhetanu2 = 0;

  if (settings1->SKIPCUTS || !settings1->USEDIRECTIONWEIGHTS) { // this is a setting that allows all neutrino angles,  no restriction.  Makes the code slower.
    costhetanu2=1.;
    costhetanu1=-1.;
    theta_threshold=1;
  }
  else {


    if (showerProps.emFrac<=1.E-10 && deltheta_had >1.E-10) {
      if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(showerProps.hadFrac*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)*sin(askFreqGen->GetChangle())>1)
	//if (Tools::dMax(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(hadfrac*vmmhz1m_max*heff_max*bw/1.E6)*sin(askFreqGen->GetChangle())>1)
	theta_threshold = -1;
      else {
	theta_threshold = sqrt(-1*deltheta_had*deltheta_had*log(anita1->VNOISE[0]/10.*anita1->maxthreshold/(showerProps.hadFrac*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)*sin(askFreqGen->GetChangle()))/constants::ALOG2);
	averaging_thetas1+=theta_threshold;
      } //else
      count_inthisloop1++;
    }



    if (showerProps.emFrac>1.E-10 && deltheta_had <=1.E-10) {
      dont_count++;
      if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(showerProps.emFrac*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)*sin(askFreqGen->GetChangle())>1)
	//if (Tools::dMax(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(showerProps.emFrac*vmmhz1m_max*heff_max*bw/1.E6)*sin(askFreqGen->GetChangle())>1)
	theta_threshold=-1;
      else {
	theta_threshold=sqrt(-1*deltheta_em*deltheta_em*log(anita1->VNOISE[0]/10.*anita1->maxthreshold/(showerProps.emFrac*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)*sin(askFreqGen->GetChangle()))/0.5);
	averaging_thetas2+=theta_threshold;
      } //else
      count_inthisloop2++;
    }




    if (showerProps.emFrac>1.E-10 && deltheta_had>1.E-10) {
      // if the electromagnetic and hadronic components of the shower are both non-negligible
      // then theta_threshold cannot be determined analytically so we step away from the cerenkov angle in steps equal to 1/2 * deltheta_em
      if (anita1->VNOISE[0]/10.*anita1->maxthreshold/((showerProps.sumFrac())*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.) {
	//if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/((showerProps.sumFrac())*vmmhz1m_max*heff_max*bw/1.E6)>1.) {
	theta_threshold=-1.; // if it's not detectable at all
      }
      else { // otherwise,  start stepping.
	theta_test=deltheta_em; // this is the angle we start stepping at
	vmmhz1m_test=vmmhz1m_max; // this will be the magnitude of the signal at theta_test away from the cerenkov cone.
	// find the magnitude of the signal at theta_test away from the cerenkov cone.
	askFreqGen->TaperVmMHz(askFreqGen->GetChangle()+theta_test, deltheta_em, deltheta_had, showerProps, vmmhz1m_test, djunk);
	//  if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(vmmhz1m_test*heff_max*bw/1.E6)>1.) { // is this electric field already too low to have a chance of passing the trigger threshold?
	if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(vmmhz1m_test/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.) { // is this electric field already too low to have a chance of passing the trigger threshold?
	  theta_threshold=theta_test; // then that is the maximum angular deviation
	}
	else { // otherwise increment by the step size and check again.
	  theta_test=1.5*deltheta_em;
	  vmmhz1m_test=vmmhz1m_max;
	  askFreqGen->TaperVmMHz(askFreqGen->GetChangle()+theta_test, deltheta_em, deltheta_had, showerProps, vmmhz1m_test, djunk);

	  if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(vmmhz1m_test/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.) {
	    theta_threshold=theta_test;
	  }
	  else { // otherwise increment by the step size and check again.
	    theta_test=2*deltheta_em;
	    vmmhz1m_test=vmmhz1m_max;
	    askFreqGen->TaperVmMHz(askFreqGen->GetChangle()+theta_test, deltheta_em, deltheta_had, showerProps, vmmhz1m_test, djunk);
	    //if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(vmmhz1m_test*heff_max*bw/1.E6)>1.)
	    if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(vmmhz1m_test/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.){
	      theta_threshold=theta_test;
	    }
	    else { // otherwise increment by the step size and check again.
	      theta_test=3*deltheta_em;
	      vmmhz1m_test=vmmhz1m_max;
	      askFreqGen->TaperVmMHz(askFreqGen->GetChangle()+theta_test, deltheta_em, deltheta_had, showerProps, vmmhz1m_test, djunk);
	      //if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(vmmhz1m_test*heff_max*bw/1.E6)>1.)

	      if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(vmmhz1m_test/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.){
		theta_threshold=theta_test;
	      }
	      else { // otherwise,  set is the the width of the hadronic component (much wider than the electromagnetic component)
		theta_test=deltheta_had;
		vmmhz1m_test=vmmhz1m_max;
		askFreqGen->TaperVmMHz(askFreqGen->GetChangle()+theta_test, deltheta_em, deltheta_had, showerProps, vmmhz1m_test, djunk);
		// if at the hadronic width,  you're below the threshold
		if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(vmmhz1m_test/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.)
		  //if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(vmmhz1m_test*heff_max*bw/1.E6)>1.) // if at the hadronic width,  you're below the threshold
		  theta_threshold=theta_test; // set theta_threshold
		else { // otherwise,  find theta_threshold considering the hadronic component alone.  This is conservative-- an electromagnetic component would only make it narrower.
		  theta_threshold=sqrt(-1*deltheta_had*deltheta_had*log(anita1->VNOISE[0]/10.*anita1->maxthreshold/(showerProps.hadFrac*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)*sin(askFreqGen->GetChangle()))/0.5);
		} // else: not below threshold at deltheta_had
	      } // else: not below threshold at 3*deltheta_em

	    } // else: not below threshold at 2.5*deltheta_em
	  } // else: not below threshold at 2.0*deltheta_em

	} // else: not below threshold at 1.5*deltheta_em

      } // not below threshold at 1.0*deltheta_em
      count_inthisloop3++;
      averaging_thetas3+=theta_threshold;

    }

    theta_threshold*=settings1->THETA_TH_FACTOR; // multiply theta_threshold by scale factor if requested,  for testing purposes.
    if (theta_threshold>0) { // we only pick the angle between 0 and pi so set the upper and lower limits accordingly.
      if (askFreqGen->GetChangle()-theta_threshold<0 && askFreqGen->GetChangle()+theta_threshold> constants::PI) {
	costhetanu2=1.;
	costhetanu1=-1.;
      } //if
      else if (askFreqGen->GetChangle()-theta_threshold>0 && askFreqGen->GetChangle()+theta_threshold> constants::PI) {
	costhetanu2=cos(askFreqGen->GetChangle()-theta_threshold);
	costhetanu1=-1.;
      } //else if
      else if (askFreqGen->GetChangle()-theta_threshold<0 && askFreqGen->GetChangle()+theta_threshold< constants::PI) {
	costhetanu2=1.;
	costhetanu1=cos(askFreqGen->GetChangle()+theta_threshold);
      } //else if
      else if (askFreqGen->GetChangle()-theta_threshold>0 && askFreqGen->GetChangle()+theta_threshold< constants::PI) {
	costhetanu2=cos(askFreqGen->GetChangle()-theta_threshold);
	costhetanu1=cos(askFreqGen->GetChangle()+theta_threshold);
      } //else if
    } // end if theta_threshold>0
  } // if SKIP_CUTS !=0

  if (theta_threshold>0) {
    // pick the neutrino direction,  in a coordinate system where the z axis lies along the cerenkov cone.
    costhetanu = costhetanu1+gRandom->Rndm()*(costhetanu2-costhetanu1);

    double phinu=constants::TWOPI*gRandom->Rndm(); // pick the phi of the neutrino direction,  in the same coordinate system.
    double sinthetanu=sqrt(1-costhetanu*costhetanu);
    // 3-vector of neutrino direction,  at that same coordinate system.
    nnu = TVector3(sinthetanu*cos(phinu), sinthetanu*sin(phinu), costhetanu);

    ///@todo FIX CHANGECOORD
    // nnu = nnu.ChangeCoord(refr); // rotate so it's in our normal coordinate system.
    // now the ray is aligned along the cerenkov cone and
    // the neutrino is rotated by that same angle

    //dtryingdirection+=4*PI/(2.*theta_threshold*sin(askFreqGen->GetChangle())*2*PI);
    interaction1->dtryingdirection=1/((costhetanu2-costhetanu1)/2.);
    return 1;
  } //end if theta_threshold
  else if (theta_threshold==-1.) {
    std::cout << "theta_threshold is " << theta_threshold << "\n";
    return 0;
  }
  else if (showerProps.emFrac<=1.E-10 && deltheta_had <= 1.E-10) {
    std::cout << "Error:  emfrac, hadfrac are (1st place)" << showerProps.emFrac << " " << showerProps.hadFrac << " " << "\n";
    return 0;
  } //else if

  return 0;
}
//end GetDirection()





double icemc::EventGenerator::GetThisAirColumn(const Settings* settings1,  Geoid::Position r_in, TVector3 nnu, Geoid::Position posnu,  double *col1,  double& cosalpha, double& mytheta, double& cosbeta0, double& mybeta) const {
  double myair=0; // this is the output
  // it is the column of air in kg/m^2
  cosalpha=(r_in.Dot(nnu)) / r_in.Mag(); // cosangle that the neutrino enters the earth wrt surface normal at its entrry point
  mytheta=(double)(acos(cosalpha)*constants::DEGRAD)-90.; // turn this into an angle

  //------------------added on Dec 8------------------------
  if (settings1->ATMOSPHERE) {
    int index11=int(mytheta*10.); // which index this theta corresponds to
    int index12=index11+1;
    // find column of air at this theta
    myair=(col1[index11]+(col1[index12]-col1[index11])*(mytheta*10.-double(index11)))*10.;//unit is kg/m^2
  }
  else{
    myair=0.;//don't include effect of atmosphere
  }

  //cout<<"mytheta="<<mytheta<<"; myair="<<myair<<std::endl;
  //------------------added on Dec 8------------------------
  cosbeta0= (posnu.Dot(nnu)) / posnu.Mag(); // cos angle of neutrino wrt person standing over the interaction point
  mybeta=(double)(acos(cosbeta0)*constants::DEGRAD)-90.; // turn that into a theta
  return myair;
}
//end GetThisAirColumn()


void icemc::EventGenerator::GetAir(double *col1) const {
  double nothing;
  std::string icemcSrcDir = icemc::EnvironmentVariable::ICEMC_SRC_DIR();
  
  std::ifstream air1(icemcSrcDir+"/data/atmosphere.dat"); // length of chord in air vs. theta (deg)
  //where theta is respect to "up"
  // binned in 0.1 degrees
  for(int iii=0;iii<900;iii++) {
    air1>>nothing>>col1[iii];
  } // read in chord lengths
}
//end GetAir()


int icemc::EventGenerator::TIR(const TVector3 &n_surf, const TVector3 &nrf2_iceside,  double N_IN, double N_OUT) const {
  double test=sin(nrf2_iceside.Angle(n_surf))*N_IN/N_OUT;
  if(test>=1){
    return 1;
  }
  else{
    return 0;
  }
}
//end TIR()


double icemc::EventGenerator::GetViewAngle(const TVector3 &nrf2_iceside, const TVector3 &nnu) const {
  // get viewing angle of shower
  double dtemp=nrf2_iceside.Dot(nnu);
  if (dtemp>=1 && dtemp<1.02)
    dtemp=0.999999;
  if (dtemp<=-1 && dtemp>-1.02)
    dtemp=-0.9999999;

  return acos(dtemp);
}
//end ())etViewAngle






void icemc::EventGenerator::GetFresnel(Roughness *rough1, int ROUGHNESS_SETTING,
				       const TVector3 &surface_normal, 
				       const TVector3 &air_rf, 
				       TVector3 &pol, 
				       const TVector3 &firn_rf, 
				       double efield, 
				       // const ShowerProperties& sp,
				       double deltheta_em_max, double deltheta_had_max, 
				       double &t_coeff_pokey, double &t_coeff_slappy,
				       double &fresnel, double &mag) const {

  // find angle of incidence and angle of transmission
  double incident_angle = surface_normal.Angle(firn_rf);
  double transmitted_angle = surface_normal.Angle(air_rf);

  //  double t_coeff_pokey, t_coeff_slappy;

  // this is perp the surface normal and transmitted ray,  parallel to surface
  TVector3 perp = air_rf.Cross(surface_normal).Unit();
  // this is in the bending plane
  TVector3 air_parallel = perp.Cross(air_rf).Unit();
  // this is in the bending plane
  TVector3 firn_parallel = perp.Cross(firn_rf).Unit();

  // component of polarization (in the air) perp to surface normal
  double pol_perp_firn = pol.Dot(perp); // this is the slappy component in the firn
  double pol_parallel_firn = pol.Dot(firn_parallel); // this is the pokey component in the firn
  double pol_perp_air=0, pol_parallel_air=0;

  double r_coeff_pokey =  tan(incident_angle - transmitted_angle) / tan(incident_angle + transmitted_angle);
  t_coeff_pokey = sqrt((1. - r_coeff_pokey*r_coeff_pokey));
  pol_parallel_air = t_coeff_pokey * pol_parallel_firn; // find pokey component in the air

  double r_coeff_slappy = sin(incident_angle - transmitted_angle) / sin(incident_angle + transmitted_angle);
  t_coeff_slappy = sqrt((1. - r_coeff_slappy*r_coeff_slappy));
  pol_perp_air = t_coeff_slappy * pol_perp_firn; // find slappy component in the firn

  mag=sqrt( tan(incident_angle) / tan(transmitted_angle) );

  fresnel = sqrt( pow(efield * pol_perp_air, 2) + pow(efield * pol_parallel_air, 2)) / efield;

  pol = (pol_perp_air * perp + pol_parallel_air * air_parallel).Unit();
  //cerr<<"(spec): inc "<<incident_angle*180./PI<<" : trans "<<transmitted_angle*180./PI<<" : tpokey "<<t_coeff_pokey<<" : tslappy "<<t_coeff_slappy<< std::endl;
}
//end GetFresnel()






void icemc::EventGenerator::WriteNeutrinoInfo(const int& inu, const Geoid::Position &posnu,  const TVector3 &nnu,  const Geoid::Position &r_bn,  double altitude_int,  Neutrino::Flavor nuflavor,  Neutrino::Current current,  double elast_y,  std::ofstream &nu_out) const {
  nu_out << "\n" << inu << "\t" << posnu[0] << " " << posnu[1] << " " << posnu[2] << "\t" << altitude_int << "\t" << nnu[0] << " " << nnu[1] << " " << nnu[2] << "\t" << r_bn[0] << " " << r_bn[1] << " " << r_bn[2] << "\t" << nuflavor << "\t" << current << "\t" << elast_y << "\n\n";
}




void icemc::EventGenerator::generateNeutrinos(const Settings& settings1, const CommandLineOpts& clOpts){

#ifdef ICEMC_FEEXCEPT
  feenableexcept(FE_INVALID | FE_DIVBYZERO); 
#endif

  icemc::report() << severity::info << "Seed is " << settings1.SEED << std::endl;

  // TRandom *rsave = gRandom;
  // TRandom3 *Rand3 = new TRandom3(settings1.SEED);//for generating random numbers
  // gRandom=Rand3;
  
  Secondaries sec1;
  Primaries* primary1 = new Primaries();
  AskaryanFreqsGenerator askFreqGen;
  // Antarctica* antarctica = new Antarctica(settings1.ICE_MODEL + settings1.NOFZ*10,
  // 					  settings1.CONSTANTICETHICKNESS * 1000 + settings1.CONSTANTCRUST * 100 + settings1.FIXEDELEVATION * 10 + 0,
  // 					  settings1.WEIGHTABSORPTION);
  Antarctica* antarctica = new Antarctica();
  icemc::report() << "Area of the earth's surface covered by antarctic ice is " << antarctica->ice_area << std::endl;
  
  Taumodel* taus1 = new Taumodel();
  Screen* panel1 = new Screen(0);  // create new instance of the screen class

  icemc::RootOutput ro(this, &settings1, clOpts.outputdir.c_str(), clOpts.run_no);

  if(!fDetector){
    fDetector = new ANITA(&settings1, &ro);
  }

  // input parameters
  settings1.ApplyInputs(fDetector,  &sec1,  &askFreqGen);
  askFreqGen.Initialize();

  NNU = settings1.NNU;
  RANDOMISEPOL = settings1.RANDOMISEPOL;
  
  gRandom->SetSeed(settings1.SEED); // settings seed is now updated with run_no to uniquify it elsewhere

  // @todo do this somewhere where it makes sense
  if (clOpts.trig_thresh!=0){
    fDetector->powerthreshold[4]=clOpts.trig_thresh;
  }

  Source::Spectra* nuSpectra = new Source::Spectra(&settings1);
  if(!interaction1){
    interaction1 = new Interaction(primary1, &settings1); //, 0, count1);
  }
  Roughness* rough1 = new Roughness(&settings1); // create new instance of the roughness class
  rough1->SetRoughScale(settings1.ROUGHSIZE);

  if(nuSpectra->IsSpectrum()){ ///@todo move to spectra constructor
    icemc::report() <<" Lowest energy for spectrum is 10^18 eV! \n";
  }
  time_t raw_start_time = time(NULL);
  struct tm*  start_time = localtime(&raw_start_time);
  icemc::report() << "Date and time at start of run are: " << asctime (start_time) << "\n";

  // for attenuation of neutrino in atmosphere
  // only important for black hole studies
  double col1[900];
  GetAir(col1);
  // double myair;//air column density, kg/m^2


  // what is this?
  costheta_nutraject2 = interaction1->costheta_nutraject;
  phi_nutraject2 = interaction1->phi_nutraject;
  altitude_int2 = interaction1->altitude_int;
  d12 = interaction1->d1;
  d22 = interaction1->d2;
  dtryingdirection2 = interaction1->dtryingdirection;
  logchord2 = interaction1->logchord;
  r_fromballoon2 = interaction1->r_fromballoon[0];
  chord_kgm2_bestcase2 = interaction1->chord_kgm2_bestcase;
  chord_kgm2_ice2 = interaction1->chord_kgm2_ice;
  weight_bestcase2 = interaction1->weight_bestcase;
  r_exit2bn2 = interaction1->r_exit2bn;
  r_exit2bn_measured2 = interaction1->r_exit2bn_measured;

  //added djg ////////////////////////////////////////////////////////
  icemc::report().al_voltages_direct<<"antenna #"<<"   "<<"volts chan 1"<<"   "<<"volts chan 2"<<"    "<<"volts chan 3"<<"    "<<"volts chan 4"<<"    "<<"noise chan 1"<<"    "<<"noise chan 2"<<"    "<<"noise chan 3"<<"   "<<"noise chan 4"<<"  "<<"weight"<<std::endl;
  ////////////////////////////////////////////////////////////////////


  // zeroing
  interaction1->dnutries=0;
  eventsfound=0.; // sums weights for events that pass

  // Tools::Zero(count1->npass, 2); // sums events that pass,  without weights
  Tools::Zero(sum, 3);
  eventsfound_db=0;
  eventsfound_nfb=0;

  Tools::Zero(eventsfound_binned, NBINS);
  Tools::Zero(eventsfound_binned_e, NBINS);
  Tools::Zero(eventsfound_binned_mu, NBINS);
  Tools::Zero(eventsfound_binned_tau, NBINS);
  Tools::Zero(eventsfound_nfb_binned, NBINS);

  fNeutrinoPath = new NeutrinoPath(); // init here for branch setting


  // fills arrays according to antenna specs
  fDetector->GetBeamWidths(&settings1); // this is used if GAINS set to 0
  // Antenna measured gain vs. frequency
  fDetector->ReadGains(); // this is used if GAINS set to 1
  fDetector->Set_gain_angle(&settings1, askFreqGen.NMEDIUM_RECEIVER);
  if(settings1.WHICH == Payload::Anita1Simple ||
     settings1.WHICH == Payload::Anita1){
    fDetector->SetDiffraction(); // for the upper ring
  }

  fDetector->saveGainsPlot(clOpts.outputdir+"/gains.eps");
  

  // sets position of balloon and related quantities
  fDetector->SetDefaultBalloonPosition();
  fDetector->SetNoise(&settings1, fDetector, antarctica);

  // find the maximum distance the interaction could be from the balloon and still be within the horizon.
  // antarctica->GetMAXHORIZON(bn1);
  antarctica->GetMAXHORIZON(fDetector);  

  // calculate the volume of ice seen by the balloon for all balloon positions
  // antarctica->CreateHorizons(&settings1, bn1, fDetector->theta_bn, fDetector->phi_bn, fDetector->altitude_bn, icemc::report().foutput);
  // antarctica->CreateHorizons(&settings1, fDetector);
  // antarctica->WorldModel::CreateHorizons(fDetector, 800e3);
  icemc::report() << "Done with CreateHorizons.\n";


  // sets neutrino energy
  if ( nuSpectra->IsMonoenergetic() ){
    pnu=pow(10., settings1.EXPONENT);
    primary1->GetSigma(pnu, sigma, len_int_kgm2, &settings1, xsecParam_nutype, xsecParam_nuint);    // get cross section and interaction length.
    icemc::report() << "pnu,  sigma,  len_int_kgm2 are " << pnu << " " << sigma << " " << len_int_kgm2 << "\n";
  }

  if (settings1.WRITEPOSFILE==1) {
    icemc::report().nu_out << "Neutrinos with energy " << pnu << "\n\n"; //Write header to file of neutrino positions
    icemc::report().nu_out << "nu #,  position of nu interaction,  depth of int.,  Direction of nu momentum,  Geoid::Position of balloon,  nu flavour,  current type,  elasticity\n\n\n\n";
  }



  // get positions of the anita payload during the slac test
  if (settings1.SLAC){
    icemc::report() << severity::error << "SLAC is no longer supported!" << std::endl;
    exit(1);
  }

  nuSpectra->savePlots2(clOpts.outputdir + "/GetG_test1.pdf");
  //if using energy spectrum
  if (nuSpectra->IsSpectrum()){
    nuSpectra->savePlots("Temp");
  }

  NeutrinoFactory nuFactory(&settings1);  

  time_t raw_loop_start_time = time(NULL);
  icemc::report() << "Starting loop over events. Time required for setup is "
       <<(int)((raw_loop_start_time - raw_start_time)/60) << ":"
       << ((raw_loop_start_time - raw_start_time)%60) << std::endl;


  double startTime = fDetector->getStartTime();
  double endTime = fDetector->getEndTime();

  icemc::report() << severity::info << "Start time is " << startTime << ", end time is " << endTime << ", time range is " << endTime - startTime << std::endl;
  // UInt_t simulationEndTime = fDetector->startTime();

  Source::DiffuseFlux directionModel;
  
  /**
   * Main loop over generated neutrinos
   */
  for (int inu = clOpts.startNu; inu < NNU && !ABORT_EARLY; inu++) {
    
    if (NNU >= 100 && ((inu % (NNU / 100)) == 0)){
      std::cout << inu << " neutrinos. " << (double(inu)/double(NNU)) * 100 << "% complete.\n";
    }
    else if(NNU < 100){
      std::cout << inu << " / " << NNU << "neutrinos" << std::endl;
    }
    
    eventNumber=(UInt_t)(clOpts.run_no)*NNU+inu;
    RNG::newSeeds(clOpts.run_no, eventNumber); // updates all RNG seeds

    // @todo Set seed of all random number generators to be dependent on eventNumber
    gRandom->SetSeed(eventNumber+6e7);

    double eventTime = fDetector->pickEventTime();

    // first pick a balloon position...
    fDetector->PickBalloonPosition(eventTime, &settings1,  fDetector);

    Geoid::Position interactionPos = antarctica->PickInteractionLocation(fDetector->position());

    RayTracer rayTracer(antarctica, fDetector->position());
    TVector3 rfDirFromInteraction = rayTracer.findPathToDetector(interactionPos);
    if(rfDirFromInteraction.Mag()==0){
      // no solution was found withing the fitter tolerance
      ///@todo Fill some tree indicating this!
      continue;
    }
    
    // if we get here, then there's a ray tracing solution
    // to get from our interaction point to the payload
    // now we have that, we can calculate the neutrino path

    TVector3 v = directionModel.pickNeutrinoDirection(interactionPos, rfDirFromInteraction); //, antarctica);

    // make a neutrino, we've picked energy, flavor, interaction current, 
    Neutrino nu = nuFactory.makeNeutrino();

    // ShowerProperties 
    
    // askFreqs = askFreqGen.generateAskaryanFreqs(vmmhz_max, vmmhz1m_max, pnu, fDetector->NFREQ, fDetector->freq, fDetector->NOTCH_MIN, fDetector->NOTCH_MAX, &showerProps);
    
    /// pick neutrino energy


    // propagate to ANITA

    
    
  }

  return;
  // for(int inu=0; inu < 1; inu++){

  //   // then pick an interation position

    

    

    

  //   // generate a new one for each loop...
  //   // is there a more elegant way to do this?
  //   fNeutrinoPath->reset();

  //   if(fGenNu) {delete fGenNu;}
  //   fGenNu = new GeneratedNeutrino(inu);

  //   if(fPassNu) {
  //     delete fPassNu;
  //     fPassNu = NULL; // will be initialized in the case that a neutrino triggers the detector
  //   }

  //   if (NNU >= 100 && (inu % (NNU / 100) == 0)){
  //     std::cout << inu << " neutrinos. " << (double(inu)/double(NNU)) * 100 << "% complete.\n";
  //   }
  //   else{
  //     std::cout << inu << " neutrinos.  " << (double(inu) / double(NNU)) * 100 << "% complete.\n";
  //   }

  //   // if (settings1.NOISEFROMFLIGHTDIGITIZER || settings1.NOISEFROMFLIGHTTRIGGER) {
  //   //   fDetector->fRand->SetSeed(eventNumber+8e9);
  //   // }
    
  //   //reset screen parameters (even for no roughness) for the new event
  //   panel1->ResetParameters();
  //   // fDetector->inu=inu;

  //   int whichray = 0;
  //   // minray = direct,  maxray = reflected
  //   // for (int whichray = settings1.MINRAY; whichray <= settings1.MAXRAY; whichray++) {

  //   fDetector->passglobtrig[0] = 0;
  //   fDetector->passglobtrig[1] = 0;
  //   passes_thisevent = 0;

  //   if (nuSpectra->IsSpectrum()){//if using energy spectrum

  //     if(settings1.USEDARTBOARD){
  // 	pnu = nuSpectra->GetNuEnergy();
  //     }
  //     else {
  // 	pnu = nuSpectra->GetCDFEnergy();
  //     }

  //     ierr = primary1->GetSigma(pnu, sigma, len_int_kgm2, &settings1, xsecParam_nutype, xsecParam_nuint);  // given neutrino momentum,  cross section and interaction length of neutrino.
  //     // ierr=0 if the energy is too low for the parameterization
  //     // ierr=1 otherwise
  //     fNeutrinoPath->len_int=1.0/(sigma*askFreqGen.RHOH20*(1./constants::M_NUCL)*1000); // in km (why interaction length in water?) //EH
  //   }// end IsSpectrum

  //   // count_pass=0;
  //   passestrigger=0;
  //   chanceinhell2=0;
  //   sec1.secondbang=false;
  //   count_total++;
  //   // initializing the voltage seen by each polarization of each antenna
  //   interaction1->dtryingdirection=0;
  //   // fDetector->dtryingposition=0;

  //   if (settings1.HIST && !settings1.ONLYFINAL && ro.prob_eachilon_bn.GetEntries() < settings1.HIST_MAX_ENTRIES) {
  //     ro.prob_eachilon_bn.Fill(fDetector->position().Longitude());
  //   }

  //   // pick random point in ice.
  //   // also get initial guess shower exit position
  //   // unit vector from shower exit to balloon.
  //   // get both positions of interaction point and its mirror point.
  //   //-------------------------------------------------------
  //   beyondhorizon = 0;

  //   if (interaction1){
  //     delete interaction1;
  //     interaction1 = NULL;
  //   }
  //   interaction1 = new Interaction(primary1,  &settings1); //,  whichray,  count1);

    // moved from deep inside interaction
    // if (settings1.MINRAY==whichray) {
    //   // only increment neutrino flavor on first ray so you don't count twice
    //   if (interaction1->nuflavor==Neutrino::Flavor::e){
    // 	count1->nnu_e++;
    //   }
    //   else if (interaction1->nuflavor==Neutrino::Flavor::mu){
    // 	count1->nnu_mu++;
    //   }
    //   else if (interaction1->nuflavor==Neutrino::Flavor::tau){
    // 	count1->nnu_tau++;
    //   }
    // }

    // if(taus1){
    //   delete taus1;
    //   taus1 = NULL;
    // }
    // taus1 = new Taumodel();

    // int taumodes = settings1.taumodes;
    // tauweighttrigger=0;
    // interaction1->weight_nu=0;
    // interaction1->weight_nu_prob=0;
    // taus1->weight_tau_prob=0;
      
    // if (taumodes==1 && interaction1->nuflavor=="nutau" && interaction1->current=="cc"){
    // if (taumodes==1 && interaction1->nuflavor==Neutrino::Flavor::tau && interaction1->current==Neutrino::Current::Charged){
    //   tautrigger = 1;///< tau trigger sets the chance to create tau particle
    // }
    // else{
    //   tautrigger=0;
    // }

    // fDetector->PickDownwardInteractionPoint(interaction1,  anita1,  &settings1,  antarctica,  ray1,  beyondhorizon);
    // beyondhorizon = interaction1->PickDownwardInteractionPoint(fDetector->position(), &settings1,  antarctica,  &rayTracer);
    // fDetector->PickDownwardInteractionPoint(interaction1,  &settings1,  antarctica,  &rayTracer,  beyondhorizon);      

    // if (interaction1->noway){
    //   fGenNu->passCutNoWay = 0;
    //   std::cerr << "All tree filled!" << ro.allTree.GetEntries() << "\n";
    //   ro.allTree.Fill();
    //   continue;
    // }
    // else{
    //   fGenNu->passCutNoWay = 1;
    // }
    // // count1->noway[whichray]++;

    // if (interaction1->wheredoesitleave_err){	
    //   continue;
    // }
    // // count1->wheredoesitleave_err[whichray]++;

    // if (interaction1->neverseesice){
    //   continue;
    // }
    // // count1->neverseesice[whichray]++;

    // if (interaction1->wheredoesitenterice_err){
    //   continue;
    // }
    // // count1->wheredoesitenterice_err[whichray]++;

    // if (interaction1->toohigh){
    //   continue;
    // }
    // // count1->toohigh[whichray]++;

    // if (interaction1->toolow){
    //   continue;
    // }
    // // count1->toolow[whichray]++;

    // if (!interaction1->iceinteraction){
    //   continue;
    // }
    // // count1->iceinteraction[whichray]++;

    // if (beyondhorizon) {
    //   fGenNu->passCutWithinHorizon = 0;
    //   ro.allTree.Fill();
    //   continue;
    // }
    // else{
    //   fGenNu->passCutWithinHorizon = 1;
    // }
    // // count1->inhorizon[whichray]++;

    // // cerenkov angle depends on depth because index of refraction depends on depth.

    // if (settings1.FIRN) {
    //   askFreqGen.SetNDepth(antarctica->GetN(interaction1->altitude_int));
    //   changle_deg=askFreqGen.GetChangle()*constants::DEGRAD;
    // }

    // rayTracer.GetSurfaceNormal(&settings1, antarctica, interaction1->posnu, slopeyangle, 0);

    // *** NOTE **** for Snell's law,  I call the ray on the air-side
    // the incident angle and the ice-side ray the refracted

    // ray's angle of incidence (in the air) onto ice
    // costheta_inc=rayTracer.n_exit2bn[0].Dot(rayTracer.nsurf_rfexit);    // just for plotting

    // // just for plotting
    // costheta_exit=cos(rayTracer.rfexit[0].Theta()); // just for plotting

    // if (!ray1->TraceRay(&settings1, anita1, 1, askFreqGen.N_DEPTH)) {
    // if (!rayTracer.TraceRay(&settings1, fDetector, 1, askFreqGen.N_DEPTH)) {
    //   continue;
    // }

    // // use snell's law to get the first guess at the
    // // direction of the rf as it leaves ice surface.
    // // 0th guess was simply radially outward from interaction position
    // // this now takes into account balloon position and surface normal.
    // // ray1->GetRFExit(&settings1, anita1, whichray, interaction1->posnu, interaction1->posnu_down, fDetector->r_bn, fDetector->r_boresights, 1, antarctica); // fills ray1->n_exit2bn[1]
    // const int inuDebug = -1; //1014; //1481; //1705; //34; //532; //543; //517; //523; //522;
    // bool debugRay = inu == inuDebug ? true : false;
    // rayTracer.setDebug(debugRay);
    // TVector3 rfDirFromInteraction = rayTracer.findPathToDetector(interaction1->posnu, fDetector->position());

    // if(rfDirFromInteraction.Mag()==0){
    //   // no solution was found withing the fitter tolerance
    //   ///@todo Fill some tree indicating this!
    //   continue;
    // }

    // rayTracer.GetRFExit(&settings1, fDetector, whichray, interaction1->posnu, interaction1->posnu_down, fDetector->position(), fDetector->r_boresights, 1, antarctica);

     // rayTracer.GetSurfaceNormal(&settings1, antarctica, interaction1->posnu, slopeyangle, 1);

    // std::cout << "inu = " << inu << std::endl;
    // if(debugRay){
    //   exit(1);
    // }
    // if (!ray1->TraceRay(&settings1, anita1, 2, askFreqGen.N_DEPTH)) {; // trace ray,  2nd iteration.
    // if (!rayTracer.TraceRay(&settings1, fDetector, 2, askFreqGen.N_DEPTH)) {; // trace ray,  2nd iteration.
    //   continue;
    // }

    // fills ray1->n_exit2bn[2] ?
    // ray1->GetRFExit(&settings1, anita1, whichray, interaction1->posnu, interaction1->posnu_down, fDetector->r_bn, fDetector->r_boresights, 2, antarctica);
    // rayTracer.GetRFExit(&settings1, fDetector, whichray, interaction1->posnu, interaction1->posnu_down, fDetector->position(), fDetector->r_boresights, 2, antarctica);

    // rayTracer.GetSurfaceNormal(&settings1, antarctica, interaction1->posnu, slopeyangle, 2);

    // intermediate counter
    // count1->nraypointsup1[whichray]++;

    // double emfrac_db = 0, hadfrac_db = 0;
    // sec1.GetTauDecay(interaction1->nuflavor, interaction1->current, taudecay,  emfrac_db,  hadfrac_db);

    // // pick elasticity
    // elast_y = primary1->pickY(&settings1, pnu, 0, Neutrino::Current::Charged);
    // if (settings1.CONSTANTY==1) { // if we ask to make y a constant=0.2
    //   elast_y = 0.2;
    //   interaction1->nuflavor = Neutrino::Flavor::e;//"nue";
    //   interaction1->current = Neutrino::Current::Charged; //"cc";
    // }

    // if (ro.ytree.GetEntries()<settings1.HIST_MAX_ENTRIES && !settings1.ONLYFINAL && settings1.HIST==1){
    //   ro.ytree.Fill();
    // }

    // //TAU STUFF. Pick whether it will stay as a neutrino or create tau
    // if(tautrigger==1){
    //   if (( !settings1.UNBIASED_SELECTION) && !settings1.SLAC ) {

    // 	  /**
    // 	   * @todo @warning bug discovered during =nicemc= refactor?
    // 	   * The original, unrefactored version of this code passed emfrac and hadfrac to GetDirection (which uses those numbers)
    // 	   * These values were not set in this iteration of the event loop, and so were probably the values from the last loop.
    // 	   * It seems likely that emfrac_db and hadfrac_db were meant to be used here based on how emfrac & hadfrac are set to these
    // 	   * values in the tau case later.
    // 	   * I'll put these in now, but since I didn't write the original code I'm not 100% sure of the author's intent.
    // 	   */
    // 	  ShowerProperties showerPropsTemp;
    // 	  showerPropsTemp.emFrac = emfrac_db;
    // 	  showerPropsTemp.hadFrac = hadfrac_db;
    //     err = GetDirection(&settings1, interaction1, rayTracer.nrf_iceside[4], deltheta_em_max, deltheta_had_max, showerPropsTemp, vmmhz1m_max*bestcase_atten, interaction1->r_fromballoon[whichray], &rayTracer, &askFreqGen, interaction1->posnu, fDetector, fDetector, interaction1->nnu, costhetanu, theta_threshold);	  
    //   }
    //   if(err==0){
    //     continue;//bad stuff has happened.
    // 	}
    //   interaction1->r_in = antarctica->WhereDoesItEnter(interaction1->posnu, interaction1->nnu);

    //   taus1->GetTauWeight(primary1,  &settings1,  antarctica,  interaction1,  pnu,  1,  ptauf, crust_entered);

    //   antarctica->Getchord(&settings1,
    // 			     len_int_kgm2,
    // 			     interaction1->r_in,
    // 			     interaction1->r_enterice,
    // 			     interaction1->nuexitice,
    // 			     interaction1->posnu, inu,
    // 			     interaction1->chord,
    // 			     interaction1->weight_nu_prob,
    // 			     interaction1->weight_nu,
    // 			     fNeutrinoPath->nearthlayers,
    // 			     myair,
    // 			     total_kgm2,
    // 			     crust_entered,
    // 			     mantle_entered,
    // 			     core_entered);

    //   nutauweight = interaction1->weight_nu_prob;
    //   tauweight = taus1->weight_tau_prob;

    //   nutauweightsum += nutauweight;
    //   tauweightsum += tauweight;
    //   double xrndm = gRandom->Rndm();

    //   if(xrndm <=taus1->weight_tau_prob/(taus1->weight_tau_prob+interaction1->weight_nu_prob)){
    //     pnu = ptauf;//set the energy we are looking at to the final energy of the tau. cuts out alot of if-else statements
    //     tauweighttrigger=1; //From now on, signifies a tau particle
    //   }
    //   if(fTauPtr){
    //     delete fTauPtr;
    // 	  fTauPtr = NULL;
    // 	}

    //   fTauPtr = new Taumodel();
    //   fTauPtr->inu = inu;
    //   fTauPtr->ptauf = ptauf;
    //   fTauPtr->weight_nu_prob = interaction1->weight_nu_prob;
    //   fTauPtr->weight_tau_prob = taus1->weight_tau_prob;
    //   ro.mytaus_tree.Fill();

    // }//end tautrigger ==1
    // /////////////////////// end of tau stuff

      
    // get fraction of shower that is electromagnetic.
    // pi^0's are counted as hadronic.
    // sec1->GetEMFrac(&settings1, interaction1->nuflavor, interaction1->current, taudecay, elast_y, &ro.hy, pnu, inu,emfrac, hadfrac, n_interactions, tauweighttrigger);

    // ShowerProperties showerProps = sec1.GetEMFrac(&settings1, interaction1->nuflavor, interaction1->current, taudecay, elast_y, &ro.hy, pnu, inu, tauweighttrigger);

    // // for double bangs, surely this should be in Secondaries?
    // if(sec1.secondbang && sec1.interestedintaus) {
    //   ptau=(1-elast_y)*pnu;
    //   showerProps.emFrac=emfrac_db;
    //   showerProps.hadFrac=hadfrac_db;
    // }

    // /**
    //  * Find the highest possible electric field emitted at this energy
    //  * (this corresponds to the electric field at the highest frequency detectable by the antennas.
    //  * Also find the maximum width of Cerenkov cone (which is at lowest frequency)
    //  * These are used to find the maximum angular deviation from Cerenkov cone where signal is still detectable.
    //  */

    // //@todo this could be better!
    // const double testTaperFreqHz = fDetector->freq[0] > 0 ? fDetector->freq[0] : fDetector->freq[1];
      
    // if(sec1.secondbang && sec1.interestedintaus) {
    //   vmmhz1m_max = askFreqGen.GetVmMHz1m(ptau, fDetector->FREQ_HIGH);
    //   askFreqGen.GetSpread(ptau, showerProps, testTaperFreqHz, deltheta_em_max, deltheta_had_max);
    // }
    // else {
    //   // get peak signal at highest edge of frequency band because that is where it is highest
    //   vmmhz1m_max = askFreqGen.GetVmMHz1m(pnu, fDetector->FREQ_HIGH);
    //   askFreqGen.GetSpread(pnu, showerProps, testTaperFreqHz, deltheta_em_max, deltheta_had_max);
    // } //end else (not secondbang or not interested in taus)

    // static bool firstTimeSpread = true;
    // TFile* fTestSpread = NULL;
    // if(firstTimeSpread){
    //   fTestSpread = new TFile("fTestSpread.root","recreate");	
    //   TGraph gr;
    //   for(int j=0; j < Anita::NFREQ; j++){
    // 	askFreqGen.GetSpread(pnu, showerProps, fDetector->freq[j], deltheta_em_max, deltheta_had_max);
    // 	gr.SetPoint(j, fDetector->freq[j], deltheta_em_max);
    // 	std::cout << j << "\t" <<  fDetector->freq[j] << "\t" << deltheta_em_max << "\t" << deltheta_had_max << std::endl;
    //   }

    //   gr.SetName("grTestSpread");
    //   std::cout << "grTestSpread!!!!!!!!!!!!!" << std::endl;
    //   std::cout << "grTestSpread!!!!!!!!!!!!!" << std::endl;
    //   std::cout << "grTestSpread!!!!!!!!!!!!!" << std::endl;
    //   std::cout << "grTestSpread!!!!!!!!!!!!!" << std::endl;
    //   gr.Write();
    //   fTestSpread->Close();
    //   firstTimeSpread = false;
    // }
      
    // /*
    //  * Using highest possible signal and minimum noise, accounting for distance from balloon,
    //  * and best case attenutation... reject if signal is undetectable.
    //  */
    // if (whichray==direct){ // if we are looking at direct rays	
    //   // the attenuation is obtained from the altitude of the interaction (shortest path is if the signal went straight up)	
    //   bestcase_atten = exp(interaction1->altitude_int/MAX_ATTENLENGTH);
    // }
    // if (whichray==downgoing){ // if we are looking at reflected rays
    //   //use the real path which seems from the mirror point.
    //   bestcase_atten=exp(interaction1->altitude_int_mirror/MAX_ATTENLENGTH);
    // }

    // // let's keep this even in the roughness case, since it still represents a ceiling value
    // if (fDetector->VNOISE[0]/10.*fDetector->maxthreshold/((showerProps.sumFrac())*vmmhz1m_max*bestcase_atten/interaction1->r_fromballoon[whichray]*heff_max*fDetector->bwmin/1.E6)>settings1.CHANCEINHELL_FACTOR
    // 	&& !settings1.SKIPCUTS) {
    //   // by comparing highest possible signal to the lowest possible noise, reject if there is just no way we could detect this event.
    //   continue;

    //   // vmmhz1m_max=signal at highest frequency
    //   // bestcase_atten=best case attenuation
    //   // r_fromballoon=distance from interaction to balloon
    //   //heff_max=maximum effective height over the frequency range
    // }

    // intermediate counter
    // count1->nnottoosmall[whichray]++;

    // pick neutrino direction.
    // This GetDirection() picks a neutrino direction such that its cerenkov cone
    // is close enough to the balloon line of sight that you have a chance in hell of seeing the signal.
    // if (whichray==downgoing) {
    //   chengji=rayTracer.nrf_iceside[4].Dot(rayTracer.nrf_iceside[0]);//the projection of nrf_iceside[2] on the direction of radius direction of interaction1->posnu
    //   //original nrf_iceside[4] is the upgoing direction of signals after being reflected.
    //   //now I get the corresponding downward direction of real signals in my case.
    //   //The two vectors are symmetric to the tangent plane of the Earth at interaction point
    //   rayTracer.nrf_iceside[4] = rayTracer.nrf_iceside[4] - 2*chengji*rayTracer.nrf_iceside[0];
    // } //if whichray==1

    // if(tautrigger==0){//did this for cc- taus already,  do again for all other particles
    //   if (( !settings1.UNBIASED_SELECTION) && !settings1.SLAC ){
    // 	err = GetDirection(&settings1, interaction1, rayTracer.nrf_iceside[4],
    // 			   deltheta_em_max, deltheta_had_max,
    // 			   showerProps,
    // 			   vmmhz1m_max*bestcase_atten, interaction1->r_fromballoon[whichray],
    // 			   &rayTracer,
    // 			   &askFreqGen,
    // 			   interaction1->posnu,
    // 			   // anita1, bn1,
    // 			   fDetector, fDetector,
    // 			   interaction1->nnu,
    // 			   costhetanu,
    // 			   theta_threshold);
    //   }
    // }//end tau trigger ==0

    // // gets angle between ray and neutrino direction
    // viewangle = GetViewAngle(rayTracer.nrf_iceside[4], interaction1->nnu);
    // if(viewangle>1.57 && !settings1.SKIPCUTS) { //discard the event if viewangle is greater than 90 degrees	
    //   continue;
    // }
    // // count1->nviewangle_lt_90[whichray]++; // add to counter

    // if (!RayTracer::WhereDoesItLeave(interaction1->posnu, interaction1->nnu, antarctica, interaction1->nuexit)){
    //   std::cout << "You fucked up the ray tracer!?!?!?!" << std::endl;
    //   continue; // doesn't give a real value from quadratic formula
    // }
      
    // nuexitlength=interaction1->posnu.Distance(interaction1->nuexit);
    // // probability a tau would decay within this length at this
    // // energy
    // nuexitice=interaction1->posnu.Distance(interaction1->nuexitice);
    // theta_threshold_deg=theta_threshold*constants::DEGRAD;

    // // neutrino direction in frame where balloon is up,  0=east, 1=north, 2=up
    // // n_nutraject_ontheground = Vector(fDetector->n_east.Dot(interaction1->nnu),  fDetector->n_north.Dot(interaction1->nnu),  fDetector->n_bn.Dot(interaction1->nnu));

    // cosviewangle=cos(viewangle); // cosine angle
    // viewangle_deg=viewangle*constants::DEGRAD; // same angle but in degrees
    // dviewangle_deg=(askFreqGen.GetChangle()-viewangle)*constants::DEGRAD; // deviation from cerenkov angle

    // if (ro.viewangletree.GetEntries()<settings1.HIST_MAX_ENTRIES && !settings1.ONLYFINAL && settings1.HIST==1){
    //   ro.viewangletree.Fill(); // fills variables related to viewing angle
    // }

    // if (whichray==downgoing) {
    //   //return it to the upgoing direction that is after being reflected
    //   rayTracer.nrf_iceside[4] = rayTracer.nrf_iceside[4] + 2*chengji*rayTracer.nrf_iceside[0];
    // }

    // if (err==0) {
    //   // count1->nbadfracs[whichray]++;
    //   std::cout<<"err==0,  so leaving.\n";
    //   continue;
    // }
    // // count1->ngoodfracs[whichray]++;

    // // these variables are just for plotting
    // nsigma_em_threshold=theta_threshold/deltheta_em_max;
    // nsigma_had_threshold=theta_threshold/deltheta_had_max;

    // For each neutrino,  multiply the number of tries
    // necessary to generate the appropriate direction,
    // and the number of tries necessary to generate
    // an appropriate position,  and assume they are
    // independent.

    
    // double dtryingposition = 1;///@todo dtryingposition IS AN IMPORTANT WEIGHT FACTOR THAT MUST BE RESTORED
    // interaction1->dnutries=interaction1->dtryingdirection*dtryingposition;

    // // where the neutrino enters the earth
    // if (tautrigger==0){//did for cc-taus already,  do for all other particles
    //   interaction1->r_in = antarctica->WhereDoesItEnter(interaction1->posnu, interaction1->nnu);
    // }

    // // total chord
    // double chord_kgm2_test=interaction1->posnu.Distance(interaction1->r_in)*askFreqGen.RHOMEDIUM;

    // double weight_test=0;  // weight if the whole chord from interaction to earth entrance is ice.
    // // take best case scenario chord length and find corresponding weight

    // IsAbsorbed(chord_kgm2_test, len_int_kgm2, weight_test);
    // // if the probably the neutrino gets absorbed is almost 1,  throw it out.

    // if (!settings1.SKIPCUTS) {
    //   if (weight_test<CUTONWEIGHTS) {
    // 	continue;
    //   }
    // }
    // count_chanceofsurviving++;

    // // theta of nu entrance point,  in earth frame
    // // and latitude
    // fNeutrinoPath->theta_in = interaction1->r_in.Theta();
    // fNeutrinoPath->lat_in = -90+fNeutrinoPath->theta_in*constants::DEGRAD;

    // // find quantities relevent for studying impact of atmosphere
    // // for black hole studies
    // // costheta and mytheta: theta of neutrino wrt surface normal where neutrino enters earth
    // // cosbeta0, mybeta: theta of neutrino wrt surface normal for a person standing above the interaction point
    // myair = GetThisAirColumn(&settings1,  interaction1->r_in, interaction1->nnu, interaction1->posnu, col1, cosalpha, mytheta, cosbeta0, mybeta);

    // // where the neutrino enters the ice
    // // reject if it enters beyond the borders of the continent.
    // // step size is 1/10 of interaction length
    // if (!settings1.UNBIASED_SELECTION) {
    //   if (!antarctica->WhereDoesItEnterIce(interaction1->posnu, interaction1->nnu, len_int_kgm2/askFreqGen.RHOMEDIUM/10., interaction1->r_enterice)) {
    // 	//r_enterice.Print();
    // 	if (antarctica->OutsideAntarctica(interaction1->r_enterice)) {
    // 	  std::cout<<"Warning!  Neutrino enters beyond continent,  program is rejecting neutrino! inu = "<<inu<<std::endl;
    // 	  continue;
    // 	}// end outside antarctica
    //   }// end wheredoesitenterice
    // }// end if unbiased
    // // intermediate counter
    // // count1->nentersice[whichray]++;

    // // d1=earth entrance to rock-ice interface
    // // d2=rock-ice interface to position of neutrino interaction
    // interaction1->d1=interaction1->r_enterice.Distance(interaction1->r_in);
    // interaction1->d2=interaction1->r_enterice.Distance(interaction1->posnu);


    // // get a lower limit on the chord that the neutrino traverses,
    // // so that later we can see if the signal is detectable in
    // // the best case scenario.
    // if(sec1.secondbang && sec1.interestedintaus) {
    //   std::cout << "Need to bring back GetFirstBang before you can simulate taus.\n";
    //   std::cout << "I removed it because it required Earth and I wanted Secondaries to be a stand-alone class to use in the embedded simulation.\n";
    //   icethickness=interaction1->r_enterice.Distance(interaction1->nuexit);
    //   interaction1->chord_kgm2_bestcase=nuentrancelength*TMath::MinElement(3, icemc::densities);
    // }
    // else {
    //   // finds minimum chord (in kg/m^2) traversed by neutrino
    //   // only keeping events with weight > 10^-3
    //   // periodically need to make sure this is still valid
    //   // chord_kgm2_bestcase=(d1+d2)*askFreqGen->RHOMEDIUM;
    //   interaction1->chord_kgm2_bestcase=(interaction1->d1+interaction1->d2)*TMath::MinElement(3, icemc::densities);
    // }

    // // chord just through ice.
    // interaction1->chord_kgm2_ice=interaction1->d2*askFreqGen.RHOMEDIUM;

    // // take best case scenario chord length and find corresponding weight
    // IsAbsorbed(interaction1->chord_kgm2_bestcase, len_int_kgm2, interaction1->weight_bestcase);

    // // if the probability that the neutrino gets absorbed is almost 1,  throw it out.
    // if (interaction1->weight_bestcase<CUTONWEIGHTS && !settings1.SKIPCUTS) {
    //   //
    //   continue;
    // }
    // //intermediate counter
    // // count1->nabsorbed[whichray]++;

    // // intermediate counter
    // // count1->nraywithincontinent1[whichray]++;

    // // now we have second guess for rf exit point,  which is
    // // pretty close to the right answer.
    // // now get our best case attenuation again,
    // // and see if we can reject the event.
    // if (whichray==direct){
    //   bestcase_atten = exp(-1*rayTracer.rfexit[1].Distance(interaction1->posnu)/MAX_ATTENLENGTH);
    // }
    // if (whichray==downgoing){
    //   bestcase_atten = exp(-1*rayTracer.rfexit[1].Distance(interaction1->posnu_down)/MAX_ATTENLENGTH);//use the real distance
    // }
    // if (fDetector->VNOISE[0]/10.*fDetector->maxthreshold/((showerProps.sumFrac())*vmmhz1m_max*bestcase_atten/interaction1->r_fromballoon[whichray]*heff_max*fDetector->bwmin/1.E6)>settings1.CHANCEINHELL_FACTOR && !settings1.SKIPCUTS) {
    //   //
    //   continue;
    // }
    // count_chanceinhell0++;

    // // intermediate counting
    // // count1->nraypointsup2[whichray]++;
     
    // double nbelowsurface = 0;
    // // reject if it is totally internally reflected at the surface AND NOT CONSIDERING ROUGHNESS
    // if (settings1.FIRN){
    //   nbelowsurface=constants::NFIRN;
    // }
    // else{
    //   nbelowsurface=askFreqGen.NICE;
    // }
    // // this is purely a sanity check.
    // // if everything is working,  events should pass with 100% efficiency
    // if (!settings1.ROUGHNESS && TIR(rayTracer.nsurf_rfexit, rayTracer.nrf_iceside[3], nbelowsurface, askFreqGen.N_AIR)) {
    //   continue;
    // }
    // // count1->nnottir[whichray]++;

    // // this sets n_exit2bn[2] to the ray from the exit point to the balloon,
    // // last iteration.  Now we're ready to do some calculations!!!!
    // // rayTracer.GetRFExit(&settings1, fDetector, whichray, interaction1->posnu, interaction1->posnu_down, fDetector->position(), fDetector->r_boresights, 2, antarctica);

    // // count1->nraywithincontinent2[whichray]++;

    // // for plotting- cos(theta) of neutrino direction standing on earth below balloon.
    // interaction1->costheta_nutraject=(interaction1->nnu.Dot(fDetector->position()))/sqrt(fDetector->position().Dot(fDetector->position()));

    // theta_rf_atbn = rayTracer.n_exit2bn[2].Angle(fDetector->position()); // polar angle of the rf signal as seen at the balloon.
    // // measured theta of the rf,  which is actual smeared by SIGMA_THETA,  whose default is 0.5 degrees.
    // theta_rf_atbn_measured = theta_rf_atbn+gRandom->Gaus()*fDetector->SIGMA_THETA;
    // interaction1->r_exit2bn=fDetector->position().Distance(rayTracer.rfexit[2]);
    // interaction1->r_exit2bn_measured=fDetector->altitude_bn/cos(theta_rf_atbn_measured);

    // // reject if the rf leaves the ice where there is water,  for example.
    // if (!antarctica->AcceptableRfexit(rayTracer.nsurf_rfexit, rayTracer.rfexit[2], rayTracer.n_exit2bn[2])){
    //   continue;
    // }//end if acceptableRFexit

    // // intermediate counting
    // // count1->nacceptablerf[whichray]++;

    // // difference between exit points of 2nd and 3rd iterations.
    // diff_3tries=rayTracer.rfexit[1].Distance(rayTracer.rfexit[2]);

    // // reject if 2nd and 3rd tries
    // // don't converge within 10m.
    // if (diff_3tries>10) {
    //   continue;
    // }
    // // count1->nconverges[whichray]++;

    // // Time to start assembling signal information...      

    // // Get Polarization vector.  See Jackson,  Cherenkov section.
    // TVector3 n_pol = GetPolarization(interaction1->nnu, rayTracer.nrf_iceside[4], inu);

    // if (settings1.FIRN){
    //   // now rotate that polarization vector according to ray paths in firn and air.
    //   // fresnel factor at ice-firn interface

    //   GetFresnel(rough1, settings1.ROUGHNESS, rayTracer.nsurf_rfexit, rayTracer.nrf_iceside[3], n_pol, rayTracer.nrf_iceside[4], vmmhz1m_max, deltheta_em_max, deltheta_had_max, t_coeff_pokey, t_coeff_slappy, fresnel1, mag1);
    //   //The gradual transition in the firn means that there is no fresnel factor, only magnification
    //   // and the magnification factor is upside down compared to what it is
    //   // for the firn-air interface
    //   vmmhz1m_fresneledonce = vmmhz1m_max/mag1;

    //   //  get fresnel factor at firn-air interface
    //   GetFresnel(rough1, settings1.ROUGHNESS, rayTracer.nsurf_rfexit, rayTracer.n_exit2bn[2], n_pol, rayTracer.nrf_iceside[3], vmmhz1m_fresneledonce, deltheta_em_max, deltheta_had_max, t_coeff_pokey, t_coeff_slappy, fresnel2, mag2);	
    //   // use both fresnel and magnification factors at firn-air interface.  Notice that magnification factor is upside-down compared to what it is in the firn.
    //   vmmhz1m_fresneledtwice = vmmhz1m_fresneledonce*fresnel2*mag2;

    // }//end if firn
    // else {
    //   askFreqGen.GetSpread(pnu, showerProps, (fDetector->bwslice_min[2]+fDetector->bwslice_max[2])/2., deltheta_em_mid2, deltheta_had_mid2);

    //   GetFresnel(rough1, settings1.ROUGHNESS, rayTracer.nsurf_rfexit, rayTracer.n_exit2bn[2], n_pol, rayTracer.nrf_iceside[4], vmmhz1m_max, deltheta_em_mid2, deltheta_had_mid2, t_coeff_pokey, t_coeff_slappy,  fresnel1, mag1);

    //   vmmhz1m_fresneledtwice = vmmhz1m_max*fresnel1*mag1;  //  only the ice-air interface

    // }//end else firn
    // //cerr<<inu<<" -- here"<<std::endl;      //}

    // if(settings1.ROUGHNESS){
    //   // applyRoughness(settings1, inu, interaction1, ray1, panel1, antarctica, bn1, &askFreqGen, anita1, showerProps);
    //   applyRoughness(settings1, inu, interaction1, &rayTracer, panel1, antarctica, fDetector, &askFreqGen, fDetector, showerProps);	
    // }

    // if( settings1.ROUGHNESS && !panel1->GetNvalidPoints() ){
    //   continue;
    // }
      
    // count1->nchanceinhell_fresnel[whichray]++;
      
    // // for plotting
    // diffexit = rayTracer.rfexit[0].Distance(rayTracer.rfexit[1]);
    // diffnorm = acos(rayTracer.nsurf_rfexit[0]*rayTracer.nsurf_rfexit[1]);
    // diffrefr = acos(rayTracer.nrf_iceside[4].Dot(rayTracer.nrf_iceside[0]));
      
    // scale by 1/r once you've found the 3rd iteration exit point
    // ALREADY DEALT WITH IN CASE OF ROUGHNESS
  // //   if (!settings1.ROUGHNESS) {
  // //     // this is what was in ScaleVmMHz
  // //     const Geoid::Position&posnu = whichray == direct ? interaction1->posnu : interaction1->posnu_down;
  // //     double r_meters = fDetector->position().Distance(rayTracer.rfexit[2]) + rayTracer.rfexit[2].Distance(posnu);
  // //     vmmhz_max = vmmhz1m_fresneledtwice/r_meters;
  // //   }



      
      

  // //   // reject if the event is undetectable.
  // //   if (!settings1.ROUGHNESS){
  // //     if (settings1.CHANCEINHELL_FACTOR*vmmhz_max*heff_max*0.5*(fDetector->bwmin/1.E6)<fDetector->maxthreshold*fDetector->VNOISE[0]/10. && !settings1.SKIPCUTS) {
  // // 	//
  // // 	continue;
  // //     } //if
  // //   }
  // //   // count1->nchanceinhell_1overr[whichray]++;

  // //   // distance ray travels through ice.
  // //   if (!settings1.ROUGHNESS) {
  // //     if (whichray==direct) {
  // // 	rflength=interaction1->posnu.Distance(rayTracer.rfexit[2]);
  // //     }
  // //     if (whichray==downgoing) {
  // // 	rflength=interaction1->posnu_down.Distance(rayTracer.rfexit[2]);//use the real distance that singals pass
  // //     }
  // //   }


  // //   // applying ice attenuation factor
  // //   if (!settings1.ROUGHNESS) {
  // //     if (whichray==direct){
  // // 	Attenuate(antarctica, &settings1, vmmhz_max,  rflength,  interaction1->posnu);
  // //     }
  // //     if (whichray==downgoing){
  // // 	Attenuate_down(antarctica, &settings1, vmmhz_max,  rayTracer.rfexit[2],  interaction1->posnu, interaction1->posnu_down);
  // //     }
  // //   }

  // //   // intermediate counting
  // //   count_dbexitsice++;

  // //   // reject if the event is undetectable.
  // //   if (!settings1.ROUGHNESS){
  // //     if (settings1.CHANCEINHELL_FACTOR*vmmhz_max*heff_max*0.5*(fDetector->bwmin/1.E6)<fDetector->maxthreshold*fDetector->VNOISE[0]/10. && !settings1.SKIPCUTS) {
  // // 	continue;
  // //     }
  // //   }

  // //   // count1->nchanceinhell[whichray]++;
      
  // //   // index for each antenna so you can use it to fill arrays
  // //   // count_rx=0;
  // //   // keeps track of maximum voltage seen on either polarization of any antenna
  // //   // voltsRX.max=0;

  // //   // Make a vector of V/m/MHz scaled by 1/r and attenuated.
  // //   // Calculates Jaime's V/m/MHz at 1 m for each frequency
  // //   // then multiplies by scale factor vmmhz_max/vmmhz1m_max
  // //   // this will need to be improved once frequency-dependent
  // //   // attenuation length is included.
  // //   AskaryanFreqs askFreqs;
  // //   if (!settings1.ROUGHNESS){

  // //     // the 1/r through ice is controlled by vmmhz_max,  which has been scaled in Attenuate (Attenuate_down)
  // //     askFreqs = askFreqGen.generateAskaryanFreqs(vmmhz_max, vmmhz1m_max, pnu, fDetector->NFREQ, fDetector->freq, fDetector->NOTCH_MIN, fDetector->NOTCH_MAX, &showerProps);

  // //     // here we get the array vmmhz by taking vmmhz1m_max (signal at lowest frequency bin) and
  // //     // vmmhz_max (signal at lowest frequency after applying 1/r factor and attenuation factor)
  // //     // and making an array across frequency bins by putting in frequency dependence.
  // //   }

  // //   // For each frequency,  get the width of Cerenkov cone
  // //   // and size of signal once position of viewing angle is taken into account

  // //   // these variables are for energy reconstruction studies

  // //   if (!settings1.ROUGHNESS){
  // //     // don't loop over frequencies if the viewing angle is too far off
  // //     double rtemp = TMath::Min((viewangle-askFreqGen.GetChangle())/(deltheta_em_max), (viewangle-askFreqGen.GetChangle())/(deltheta_had_max));
  // //     if (rtemp > AskaryanFreqsGenerator::VIEWANGLE_CUT && !settings1.SKIPCUTS) {
  // // 	continue;
  // //     }
  // //     // count1->nviewanglecut[whichray]++;

  // //     askFreqs.taperAmplitudesForOffConeViewing(viewangle);

  // //     // store low frequency post-tapering
  // //     // vmmhz_lowfreq=askFreqs[0]; // for plotting,  vmmhz at the lowest frequency
  // //     pdgcode = interaction1->getPdgCode();

  // //     // reject if it is undetectable now that we have accounted for viewing angle
  // //     if (settings1.CHANCEINHELL_FACTOR*askFreqs.maxElement()*heff_max*0.5*(fDetector->bwmin/1.E6)<fDetector->maxthreshold*fDetector->VNOISE[0]/10. && !settings1.SKIPCUTS) {
  // // 	continue;
  // //     }
  // //   }//end if roughness==0 before the Anita::NFREQ k loop, this isolates the TaperVmMHz()
      

  // //   // just for plotting
  // //   if(!settings1.ROUGHNESS){
  // //     vmmhz_max=askFreqs.maxElement();
  // //     vmmhz_min=askFreqs.minElement();
  // //   }
  // //   // intermediate counting
  // //   // count1->nchanceinhell2[whichray]++;
  // //   chanceinhell2=1;

  // //   // Dead time
  // //   if (settings1.USEDEADTIME){
  // //     if ( (fDetector->deadTime>0.9) || (gRandom->Uniform(1)<fDetector->deadTime) ){
  // // 	continue;
  // //     }
  // //   }
	    
  // //   // count1->ndeadtime[whichray]++;

  // //   //if no-roughness case, add its parameters to the saved screen parameters so specular and roughness simulations use the same code in the waveform construction
  // //   if(!settings1.ROUGHNESS){
  // //     panel1->SetNvalidPoints(1);
  // //     for (int k=0;k<Anita::NFREQ;k++) {
  // // 	panel1->AddVmmhz_freq(askFreqs[k]);
  // //     }
  // //     panel1->AddVmmhz0(askFreqs[0]);
  // //     panel1->AddVec2bln(rayTracer.n_exit2bn[2]); //
  // //     panel1->AddPol(n_pol);
  // //     panel1->AddDelay( 0. );
  // //     panel1->AddImpactPt(rayTracer.rfexit[2]); // impact point is the surface position?

  // //     // std::cout << inu << " geom vectors:" << "\n" << rayTracer.rfexit[2] << "\n" << rayTracer.n_exit2bn[2] << "\n" << std::endl;
	
  // //     panel1->AddViewangle(viewangle);
  // //     panel1->AddIncidenceAngle(rayTracer.nsurf_rfexit.Angle(rayTracer.nrf_iceside[3]));
  // //     panel1->AddTransmissionAngle(rayTracer.nsurf_rfexit.Angle(rayTracer.n_exit2bn[2]));
  // //     panel1->AddWeight( 1. );
  // //     panel1->SetWeightNorm( 1. );
  // //     panel1->AddFacetLength( 1. );
  // //     panel1->AddTparallel_polParallel(t_coeff_pokey);
  // //     panel1->AddTperpendicular_polPerpendicular(t_coeff_slappy);

  // //     panel1->AddTparallel_polPerpendicular(0.);
  // //     panel1->AddTperpendicular_polParallel(0.);

  // //     // for (int k=0;k<Anita::NFREQ;k++) {
  // //     //   if (fDetector->whichPath()==FlightPath::PeterEvent){
  // //     //     // IntegrateBands(anita1, k, panel1, fDetector->freq, vmmhz1m_max/(vmmhz_max*1.E6), sumsignal_aftertaper);
  // //     //     IntegrateBands(fDetector, k, panel1, fDetector->freq, vmmhz1m_max/(vmmhz_max*1.E6), sumsignal_aftertaper);	    
  // //     //   }
  // //     // }
  // //   }

  // //   //just added this temporarily - will make it run slower
  // //   //this gets the weight due to stopping in earth
  // //   //returns 0 if chord<1m
  // //   if (!antarctica->Getchord(&settings1, len_int_kgm2, interaction1->r_in, interaction1->r_enterice, interaction1->nuexitice, interaction1->posnu, inu, interaction1->chord, interaction1->weight_nu_prob, interaction1->weight_nu, fNeutrinoPath->nearthlayers, myair, total_kgm2, crust_entered,  mantle_entered, core_entered)){
  // //     interaction1->weight_nu_prob = -1.;
  // //   }

  // //   if(tauweighttrigger==1){
  // //     fNeutrinoPath->weight1=interaction1->weight_nu_prob + taus1->weight_tau_prob;
  // //   }
  // //   else{
  // //     fNeutrinoPath->weight1=interaction1->weight_nu_prob;
  // //   }

  // //   fNeutrinoPath->weight = fNeutrinoPath->weight1 / interaction1->dnutries * settings1.SIGMA_FACTOR;  // total weight is the earth absorption factor
  // //   // divided by the factor accounting for the fact that we only chose our interaction point within the horizon of the balloon
  // //   // then multiply by the cross section multiplier,  to account for the fact that we get more interactions when the cross section is higher
  // //   if (fNeutrinoPath->weight<CUTONWEIGHTS) {
  // //     continue;
  // //   }
  // //   eventsfound_beforetrigger += fNeutrinoPath->weight;



  // //   // intermediate counter
  // //   if(sec1.secondbang && sec1.interestedintaus){
  // //     count_asktrigger_nfb++;  // just for taus
  // //   }
  // //   else{
  // //     count_asktrigger++;
  // //   }

  // //   panel1->PropagateSignalsToDetector(&settings1, fDetector, inu);

  // //   bool eventPassedTrigger = fDetector->applyTrigger(inu);
  // //   if(eventPassedTrigger){

	
  // //     std::cout << "inu " <<  inu << ", eventNumber = " << eventNumber << "... It passed!" << std::endl;
  // //     TString canName = TString::Format("rayTracer_%d.C", inu);
  // //     rayTracer.makePlot(canName, antarctica, interaction1->posnu, interaction1->nnu, fDetector->position());

  // //     // the neutrino has passed the trigger...
  // //     fPassNu = new PassingNeutrino(*fGenNu, askFreqs, showerProps); // forced to be NULL at loop start

  // //     // fDetector->passglobtrig[0]=thispasses[0];
  // //     // fDetector->passglobtrig[1]=thispasses[1];

  // //     //calculate the phi angle wrt +x axis of the ray from exit to balloon
  // //     n_exit_phi = Tools::AbbyPhiCalc(rayTracer.n_exit2bn[2][0], rayTracer.n_exit2bn[2][1]);

  // //     // keep track of events passing trigger
  // //     // count1->npassestrigger[whichray]++;
  // //     // tags this event as passing
  // //     passestrigger=1;

  // //     // for taus
  // //     if(sec1.secondbang && sec1.interestedintaus){
  // // 	count_passestrigger_nfb++;
  // //     }
  // //     crust_entered=0; //These are switches that let us tell how far a given neutrino penetrated.  Clear them before entering Getchord.
  // //     mantle_entered=0;
  // //     core_entered=0;

  // //     // this gets the weight due to stopping in earth
  // //     // returns 0 if chord<1m
  // //     if (tautrigger==1 || antarctica->Getchord(&settings1,
  // // 						len_int_kgm2,
  // // 						interaction1->r_in,
  // // 						interaction1->r_enterice,
  // // 						interaction1->nuexitice,
  // // 						interaction1->posnu,
  // // 						inu,
  // // 						interaction1->chord,
  // // 						interaction1->weight_nu_prob,
  // // 						interaction1->weight_nu,
  // // 						fNeutrinoPath->nearthlayers,
  // // 						myair,
  // // 						total_kgm2,
  // // 						crust_entered,
  // // 						mantle_entered,
  // // 						core_entered)) {
  // // 	//cout << "passes chord.\n";
  // // 	if (ro.nupathtree.GetEntries()<settings1.HIST_MAX_ENTRIES && !settings1.ONLYFINAL && settings1.HIST==1){
  // // 	  ro.nupathtree.Fill();
  // // 	}
  // // 	// counts how many have a good chord length
  // // 	count_chordgoodlength++;

  // // 	// divide phase space factor into weight1
  // // 	if(tauweighttrigger==1){
  // // 	  fNeutrinoPath->weight_prob=interaction1->weight_nu_prob + taus1->weight_tau_prob;
  // // 	}
  // // 	else{
  // // 	  fNeutrinoPath->weight_prob=interaction1->weight_nu_prob;
  // // 	}
  // // 	fNeutrinoPath->weight1=interaction1->weight_nu;
  // // 	fNeutrinoPath->weight=fNeutrinoPath->weight1/interaction1->dnutries*settings1.SIGMA_FACTOR;
  // // 	fNeutrinoPath->weight_prob=fNeutrinoPath->weight_prob/interaction1->dnutries*settings1.SIGMA_FACTOR;

  // // 	fNeutrinoPath->pieceofkm2sr=fNeutrinoPath->weight*antarctica->volume*pow(1.E-3, 3)*askFreqGen.RHOMEDIUM/askFreqGen.RHOH20*constants::sr/(double)NNU/fNeutrinoPath->len_int;
  // // 	// if (ro.h10.GetEntries()<settings1.HIST_MAX_ENTRIES && !settings1.ONLYFINAL && settings1.HIST){
  // // 	//   ro.h10.Fill(hitangle_e_all[0], fNeutrinoPath->weight);
  // // 	// }
  // // 	fNeutrinoPath->logweight=log10(fNeutrinoPath->weight);
  // // 	interaction1->logchord=log10(interaction1->chord);
	  
  // // 	// if neutrino travels more than one meter in ice
  // // 	if (interaction1->d2>1) {
  // // 	  // intermediate counter
  // // 	  count_d2goodlength++;

  // // 	  // for taus
  // // 	  // add to tally of neutrinos found,  weighted.
  // // 	  if(sec1.secondbang && sec1.interestedintaus) {
  // // 	    eventsfound_nfb+=fNeutrinoPath->weight;
  // // 	    index_weights=(int)(((fNeutrinoPath->logweight-MIN_LOGWEIGHT)/(MAX_LOGWEIGHT-MIN_LOGWEIGHT))*(double)NBINS);
  // // 	    eventsfound_nfb_binned[index_weights]++;
  // // 	  }//end if secondbang & interestedintaus
  // // 	  else {
  // // 	    allcuts[whichray]++;
  // // 	    allcuts_weighted[whichray]+=fNeutrinoPath->weight;
  // // 	    // if (thispasses[0] && thispasses[1]) {
  // // 	    // 	allcuts_weighted_polarization[2]+=fNeutrinoPath->weight;
  // // 	    // } else if (thispasses[0]){
  // // 	    // 	allcuts_weighted_polarization[0]+=fNeutrinoPath->weight;
  // // 	    // } else if (thispasses[1]){
  // // 	    // 	allcuts_weighted_polarization[1]+=fNeutrinoPath->weight;
  // // 	    // }
  // // 	    fDetector->weight_inanita=fNeutrinoPath->weight;

  // // 	    if (ro.h1mybeta.GetEntries()<settings1.HIST_MAX_ENTRIES && !settings1.ONLYFINAL && settings1.HIST==1){
  // // 	      ro.h1mybeta.Fill(mybeta, fNeutrinoPath->weight); //get the angle distribution of mybeta
  // // 	    }
  // // 	    eventsfound+=fNeutrinoPath->weight; // counting events that pass,  weighted.
  // // 	    eventsfound_prob+=fNeutrinoPath->weight_prob; // counting events that pass,  probabilities.
  // // 	    if (cosalpha>0){
  // // 	      eventsfound_belowhorizon+=fNeutrinoPath->weight;
  // // 	    }
  // // 	    // count1->npass[whichray]++;  // counting events that pass,  unweighted.
  // // 	    // for calculating errors on sensitivity
  // // 	    // need to find how many events as a function of weight
  // // 	    // here,  we find how to index weight
  // // 	    if (fNeutrinoPath->logweight<MIN_LOGWEIGHT){  // underflows,  set to 0th bin
  // // 	      index_weights=0;
  // // 	    }
  // // 	    else if (fNeutrinoPath->logweight>MAX_LOGWEIGHT){ // overflows,  set to last bin
  // // 	      index_weights=NBINS-1;
  // // 	    }
  // // 	    else{ // which index weight corresponds to.
  // // 	      index_weights=(int)(((fNeutrinoPath->logweight-MIN_LOGWEIGHT)/(MAX_LOGWEIGHT-MIN_LOGWEIGHT))*(double)NBINS);
  // // 	    }
  // // 	    // count number of events that pass,  binned in weight
  // // 	    if (index_weights<NBINS){
  // // 	      eventsfound_binned[index_weights]++;
  // // 	    }
  // // 	    // number of events in a ring at distance from balloon
  // // 	    // if (index_distance<NBINS_DISTANCE){
  // // 	    // 	eventsfound_binned_distance[index_distance]+= fNeutrinoPath->weight;
  // // 	    // }
  // // 	    // // same,  now binned in weight,  for calculating errors
  //     // if (index_distance<NBINS_DISTANCE && index_weights<NBINS){
  // 	    // 	eventsfound_binned_distance_forerror[index_distance][index_weights]++;
  // 	    // }
  // 	    // for debugging
  // 	    if (fNeutrinoPath->logweight>-3){
  // 	      eventsfound_weightgt01+=fNeutrinoPath->weight;
  // 	    }
  // 	    // how many events just pass through crust,  for same purpose.
  // 	    if (fNeutrinoPath->nearthlayers==1){
  // 	      eventsfound_crust+=fNeutrinoPath->weight;
  // 	    }
  // 	    if (ro.h1mybeta.GetEntries()<settings1.HIST_MAX_ENTRIES && !settings1.ONLYFINAL && settings1.HIST==1) {
  // 	      ro.h1mybeta.Fill(mybeta, fNeutrinoPath->weight);
  // 	      ro.h1mytheta.Fill(mytheta, fNeutrinoPath->weight);//fill mytheta
  // 	    }
  // 	  }//end else secondbang & interestedintaus

  // 	  //for plotting events distribution map only
  // 	  if(fNeutrinoPath->weight>0.0001){
  // 	    double int_lon, int_lat;
  // 	    int event_e_coord=0, event_n_coord=0;
  // 	    float event_e, event_n;
  // 	    //here are the longitude and altitude which Amy defined
  // 	    int_lon = interaction1->posnu.Longitude(); // what latitude,  longitude does interaction occur at
  // 	    int_lat = interaction1->posnu.Latitude();
  // 	    antarctica->IceLonLattoEN(int_lon, int_lat, event_e_coord, event_n_coord);
  // 	    event_e=float(antarctica->xLowerLeft_ice+event_e_coord*antarctica->cellSize)/1000.;
  // 	    event_n=float(-1*(antarctica->yLowerLeft_ice+(antarctica->cellSize*event_n_coord)))/1000.;
  // 	    if(whichray==direct){//direct
  // 	      ro.dir_int_coord.Fill(event_e, event_n);
  // 	    }
  // 	    else if(whichray==downgoing){
  // 	      ro.ref_int_coord.Fill(event_e, event_n);
  // 	    }
  // 	  }

  // 	  // just for plotting.
  // 	  offaxis=(double)fabs(viewangle-askFreqGen.GetChangle());
  // 	  nsigma_offaxis=offaxis/deltheta_had_max;
	    
  // 	  // ro.hundogaintoheight_e.Fill(undogaintoheight_e, fNeutrinoPath->weight);
  // 	  // ro.hundogaintoheight_h.Fill(undogaintoheight_h, fNeutrinoPath->weight);
  // 	  // ro.rec_diff.Fill((rec_efield-true_efield)/true_efield, fNeutrinoPath->weight);
  // 	  // ro.rec_diff0.Fill((rec_efield_array[0]-true_efield_array[0])/true_efield_array[0], fNeutrinoPath->weight);
  // 	  // ro.rec_diff1.Fill((rec_efield_array[1]-true_efield_array[1])/true_efield_array[1], fNeutrinoPath->weight);
  // 	  // ro.rec_diff2.Fill((rec_efield_array[2]-true_efield_array[2])/true_efield_array[2], fNeutrinoPath->weight);
  // 	  // ro.rec_diff3.Fill((rec_efield_array[3]-true_efield_array[3])/true_efield_array[3], fNeutrinoPath->weight);
  // 	  // ro.recsum_diff.Fill((rec_efield_array[0]+rec_efield_array[1]+rec_efield_array[2]+rec_efield_array[3]-true_efield)/true_efield, fNeutrinoPath->weight);

  // 	  sourceLon = rayTracer.rfexit[2].Longitude();
  // 	  sourceLat = rayTracer.rfexit[2].Latitude();
  // 	  sourceAlt = antarctica->SurfaceAboveGeoid(sourceLon, sourceLat);

  // 	  //Now put data in Vectors and Positions into arrays for output to the ROOT file.
  // 	  if (settings1.HIST && ro.finaltree.GetEntries()<settings1.HIST_MAX_ENTRIES) {
  // 	    TVector3 n_bn = fDetector->position().Unit();
  // 	    for (int i=0;i<3;i++) {
  // 	      nnu_array[i] = interaction1->nnu[i];
  // 	      r_in_array[i] = interaction1->r_in[i];
  // 	      r_bn_array[i] = fDetector->position()[i];
  // 	      n_bn_array[i] = n_bn[i];
  // 	      posnu_array[i] = interaction1->posnu[i];
  // 	      // ant_max_normal0_array[i] = ant_max_normal0[i];
  // 	      // ant_max_normal1_array[i] = ant_max_normal1[i];
  // 	      // ant_max_normal2_array[i] = ant_max_normal2[i];
  // 	      // n_pol_array[i] = n_pol[i];
  // 	      r_enterice_array[i] = interaction1->r_enterice[i];
  // 	      nsurf_rfexit_array[i] = rayTracer.nsurf_rfexit[i];
  // 	      nsurf_rfexit_db_array[i] = rayTracer.nsurf_rfexit_db[i];
  // 	    } //end for (fill arrays)
  // 	    for (int j=0;j<5;j++) {
  // 	      for (int i=0;i<3;i++) {
  // 		nrf_iceside_array[j][i] = rayTracer.nrf_iceside[j][i];
  // 		nrf_iceside_db_array[j][i] = nrf_iceside_db[j][i];
  // 		n_exit2bn_array[j][i] = rayTracer.n_exit2bn[j][i];
  // 		n_exit2bn_db_array[j][i] = n_exit2bn_db[j][i];
  // 		rfexit_array[j][i] = rayTracer.rfexit[j][i];
  // 		rfexit_db_array[j][i] = rayTracer.rfexit_db[j][i];
  // 	      } //end for
  // 	    } //end for

  // 	    costheta_nutraject2	= interaction1->costheta_nutraject;
  // 	    phi_nutraject2		= interaction1->phi_nutraject;
  // 	    altitude_int2		= interaction1->altitude_int;
  // 	    d12			= interaction1->d1;
  // 	    d22			= interaction1->d2;
  // 	    dtryingdirection2		= interaction1->dtryingdirection;
  // 	    logchord2			= interaction1->logchord;
  // 	    r_fromballoon2		= interaction1->r_fromballoon[0];
  // 	    chord_kgm2_bestcase2	= interaction1->chord_kgm2_bestcase;
  // 	    chord_kgm2_ice2		= interaction1->chord_kgm2_ice;
  // 	    weight_bestcase2		= interaction1->weight_bestcase;
  // 	    r_exit2bn2		= interaction1->r_exit2bn;
  // 	    r_exit2bn_measured2	= interaction1->r_exit2bn_measured;

  // 	    sourceMag = rayTracer.rfexit[2].Mag();

  // 	    ro.finaltree.Fill();
  // 	    // count1->IncrementWeights_r_in(interaction1->r_in, fNeutrinoPath->weight);
  // 	  } //end if HIST & HISTMAXENTRIES

  // 	    // Adds an entry to header, event, gps and truth trees
  // 	    // ro.fillRootifiedAnitaDataTrees(this, settings1, &rayTracer, panel1);

  // 	  sum_weights += fNeutrinoPath->weight;
  // 	  neutrinos_passing_all_cuts++;
  // 	  times_crust_entered_det += crust_entered;  //Increment counter for neutrino numbers in each earth layer - passing neutrinos
  // 	  times_mantle_entered_det += mantle_entered;
  // 	  times_core_entered_det += core_entered;

  // 	  if (settings1.WRITEPOSFILE==1){
  // 	    WriteNeutrinoInfo(inu, interaction1->posnu, interaction1->nnu, fDetector->position(), interaction1->altitude_int, interaction1->nuflavor, interaction1->current, elast_y, icemc::report().nu_out);
  // 	  }

  // 	  // sample first 1000 events that pass to see the distribution of weights
  // 	  if (settings1.HIST && !settings1.ONLYFINAL && ro.sampleweights.GetEntries()<settings1.HIST_MAX_ENTRIES) {
  // 	    if (fNeutrinoPath->weight>1.E-6){
  // 	      ro.sampleweights.Fill(log10(fNeutrinoPath->weight));
  // 	    }
  // 	    else{
  // 	      ro.sampleweights.Fill(-6.);
  // 	    }

  // 	    // on the 1000th one,  see how low you should make the cut so that you catch 99% of the events (weighted)
  // 	    if (ro.sampleweights.GetEntries()==1000) {
  // 	      double sum_sampleintegral=0.;
  // 	      double sum_sample=0.;
  // 	      // first calculate total integral of all the weights
  // 	      for (int k=ro.sampleweights.GetNbinsX();k>=1;k--) {
  // 		sum_sampleintegral+=ro.sampleweights.GetBinContent(k)*pow(10., ro.sampleweights.GetBinLowEdge(k));
  // 	      }
  // 	      // treat the underflow bin specially
  // 	      sum_sampleintegral+=ro.sampleweights.GetBinContent(0)*pow(10., ro.sampleweights.GetBinLowEdge(1));
  // 	      // now sum until you reach 99% of the integral.
  // 	      for (int k=ro.sampleweights.GetNbinsX();k>=1;k--) {
  // 		sum_sample+=ro.sampleweights.GetBinContent(k)*pow(10., ro.sampleweights.GetBinLowEdge(k));
  // 		if (sum_sample>0.99*sum_sampleintegral) {
  // 		  // reset the cut value.
  // 		  CUTONWEIGHTS=pow(10., ro.sampleweights.GetBinLowEdge(k));
  // 		  std::cout << "CUTONWEIGHTS is " << CUTONWEIGHTS << "\n";
  // 		  k=0;
  // 		}
  // 	      }
  // 	    }
  // 	  }//end if HIST & ONLYFINAL & sampleweights HISTMAXENTRIES

  // 	  // outputs to text file variables relevant to sky map.
	    
  // 	  // icemc::report().forbrian << interaction1->costheta_nutraject << " " << n_nutraject_ontheground.Phi() << " " << fDetector->phi_bn << " " << fNeutrinoPath->logweight << "\n";
	    
  // 	  // incrementing by flavor	    
  // 	  // also bin in weight for error calculation.
  // 	  if (interaction1->nuflavor==Neutrino::Flavor::e) {
  // 	    sum[0]+=fNeutrinoPath->weight;
  // 	    eventsfound_binned_e[index_weights]++;
  // 	  } //if
  // 	  if (interaction1->nuflavor==Neutrino::Flavor::mu) {
  // 	    sum[1]+=fNeutrinoPath->weight;
  // 	    eventsfound_binned_mu[index_weights]++;
  // 	  } //if
  // 	  if(!sec1.secondbang || !sec1.interestedintaus) {
  // 	    if (interaction1->nuflavor==Neutrino::Flavor::tau) {
  // 	      sum[2]+=fNeutrinoPath->weight;
  // 	      eventsfound_binned_tau[index_weights]++;
  // 	    } //if
  // 	  } //if

  // 	} //end if interaction1->d2>1

  //     } //end if tautrigger || GetChord
  //     else {
  // 	std::cout << "Chord is less than 1m.\n";
  //     } //end else GetChord

  //     if (settings1.HIST==1 && !settings1.ONLYFINAL && fDetector->tglob->GetEntries()<settings1.HIST_MAX_ENTRIES) {// all events
  // 	// std::cout << "Filling global trigger tree.  inu is " << inu << "\n";
  // 	fDetector->tglob->Fill();
  //     }

  //     passes_thisevent = 1; // flag this event as passing
  //     fDetector->tdata->Fill();
  //     fDetector->tgaryanderic->Fill();
  //   } // end if passing global trigger conditions
  //   else {
  //     passes_thisevent = 0; // flag this event as not passing
  //   }// end else event does not pass trigger
      
  //   // if it passes the trigger,  then go ahead and
  //   // calculate the chord length,  etc.

      

	
  //   // keeping track of intermediate counters,  incrementing by weight1.
  //   // weight1 was not yet determined when integer counters were incremented.
  //   if (chanceinhell2){
  //     count_chanceinhell2_w += fNeutrinoPath->weight;
  //   }
  //   if (passestrigger){
  //     count_passestrigger_w += fNeutrinoPath->weight;
  //   }
      
  //   // volume_thishorizon=antarctica->volume_inhorizon[fDetector->Getibnposition()]/1.E9;


  //   if(fPassNu){
  //     ro.passTree.Fill();
  //   }
  //   ro.allTree.Fill(); // should be called at each "continue" inu loop, but we made it though.    
    
  //   if (ABORT_EARLY){
  //     icemc::report() << "\n***********************************************************";
  //     icemc::report() << "\n* SIGINT received,  aborting loop over events early.";
  //     icemc::report() << "\n* Stopped after event " << inu << " instead of " << NNU;
  //     icemc::report() << "\n* Any output which relied on NNU should be corrected for.";
  //     icemc::report() << "\n***********************************************************\n";
  //     break;
  //   }
  // }//end NNU neutrino loop

  // gRandom = rsave;
  // delete Rand3;

  // icemc::report() << "about to close tsignals tree.\n";
  // fDetector->fsignals=fDetector->tsignals->GetCurrentFile();
  // fDetector->fdata=fDetector->tdata->GetCurrentFile();
  // fDetector->fdata=fDetector->tgaryanderic->GetCurrentFile();
  // fDetector->fsignals->Write();
  // fDetector->fsignals->Close();

  // fDetector->fdata=fDetector->tglob->GetCurrentFile();
  // fDetector->fdata->Write();
  // fDetector->fdata->Close();

  // if (settings1.EVENTSMAP){
  //   //draw the S80-degree-latitude circle
  //   TH2F *lat80deg=new TH2F("lat80deg", "", 600, -3000, 3000, 500, -2500, 2500);
  //   lat80deg->SetMarkerColor(kRed);
  //   int tmp_e_coord=0, tmp_n_coord=0;
  //   float tmp_e, tmp_n=0;
  //   for(double lon=0;lon<360.;lon+=0.5){
  //     double lat=10.;
  //     antarctica->IceLonLattoEN(lon, lat, tmp_e_coord, tmp_n_coord);
  //     tmp_e=float(antarctica->xLowerLeft_ice+tmp_e_coord*antarctica->cellSize)/1000.;
  //     tmp_n=float(-1*(antarctica->yLowerLeft_ice+tmp_n_coord*antarctica->cellSize))/1000.;
  //     lat80deg->Fill(tmp_e, tmp_n);
  //   }//end for lon loop
  // }// end if EVENTSMAP

  // fDetector->rms_rfcm[0] = sqrt(fDetector->rms_rfcm[0] / (double)fDetector->count_getnoisewaveforms)*1000.;
  // fDetector->rms_rfcm[1] = sqrt(fDetector->rms_rfcm[1] / (double)fDetector->count_getnoisewaveforms)*1000.;
  // fDetector->rms_lab[0] = sqrt(fDetector->rms_lab[0] / (double)fDetector->count_getnoisewaveforms)*1000.;
  // fDetector->rms_lab[1] = sqrt(fDetector->rms_lab[1] / (double)fDetector->count_getnoisewaveforms)*1000.;

  // icemc::report() << "RMS noise in rfcm e-pol is " << fDetector->rms_rfcm[0] << " mV.\n";
  // icemc::report() << "RMS noise in rfcm h-pol is " << fDetector->rms_rfcm[1] << " mV.\n";
  // icemc::report() << "RMS noise in lab e-pol is " << fDetector->rms_lab[0] << "mV.\n";
  // icemc::report() << "RMS noise in lab h-pol is " << fDetector->rms_lab[1] << "mV.\n";
  // for (int i=0;i<Anita::NFREQ;i++) {
  //   fDetector->avgfreq_rfcm[i]/=(double)fDetector->count_getnoisewaveforms;
  //   fDetector->avgfreq_rfcm_lab[i]/=(double)fDetector->count_getnoisewaveforms;
  // }

  // rms_rfcm_e=fDetector->rms_rfcm[0];
  // rms_rfcm_h=fDetector->rms_rfcm[1];
  // rms_lab_e=fDetector->rms_lab[0];
  // rms_lab_h=fDetector->rms_lab[1];
  // for (int i=0;i<Anita::NFREQ;i++) {
  //   avgfreq_rfcm[i]=fDetector->avgfreq_rfcm[i];
  //   avgfreq_rfcm_lab[i]=fDetector->avgfreq_rfcm_lab[i];
  //   freq[i]=fDetector->freq[i];
  // }
  // icemc::report() << "Filling summarytree.  rms_rfcm_e is " << rms_rfcm_e << "\n";
  // ro.summarytree.Fill();


  // // maks the output file
  // // Summarize(&settings1, anita1, count1, nuSpectra, &askFreqGen, primary1, pnu, eventsfound, eventsfound_db, eventsfound_nfb,
  // // 	    sigma, sum, antarctica->volume, antarctica->ice_area, km3sr, km3sr_e, km3sr_mu, km3sr_tau, clOpts.outputdir);
  // Summarize(&settings1, fDetector, nuSpectra, &askFreqGen, primary1, pnu, eventsfound, eventsfound_db, eventsfound_nfb,
  // 	    sigma, sum, antarctica->volume, antarctica->ice_area, km3sr, km3sr_e, km3sr_mu, km3sr_tau, clOpts.outputdir);

  // icemc::report().veff_out << settings1.EXPONENT << "\t" << km3sr << "\t" << km3sr_e << "\t" << km3sr_mu << "\t" << km3sr_tau << "\t" << settings1.SIGMA_FACTOR << std::endl;//this is for my convenience

  // // for each neutrino flavor,  fraction each contributes to sensitivity.
  // sum_frac[0]=sum[0]/eventsfound;
  // sum_frac[1]=sum[1]/eventsfound;
  // sum_frac[2]=sum[2]/eventsfound;

  // // for taus.
  // sum_frac_db[0]=sum[0]/(eventsfound+eventsfound_db+eventsfound_nfb);
  // sum_frac_db[1]=sum[1]/(eventsfound+eventsfound_db+eventsfound_nfb);
  // sum_frac_db[2]=(sum[2]+eventsfound_db+eventsfound_nfb)/(eventsfound+eventsfound_db+eventsfound_nfb);
  //if (tree17->GetEntries()<settings1.HIST_MAX_ENTRIES && !settings1.ONLYFINAL && HIST==1)
  //tree17->Fill();

  // std::cout << "closing file.\n";

  time_t raw_end_time = time(NULL);
  struct tm * end_time = localtime(&raw_end_time);
  icemc::report() << "Date and time at end of run are: " << asctime (end_time) << "\n";
  icemc::report() << "\nTotal time elapsed in run is " <<(int)((raw_end_time - raw_start_time)/60)<<":"<< ((raw_end_time - raw_start_time)%60)<<std::endl;

  // heap allocated non-members need deleting
  if(primary1)       delete primary1;
  // if(ray1)        delete ray1;
  // if(count1)         delete count1;
  // if(globaltrig1) delete globaltrig1;
  if(taus1)          delete taus1;
  if(rough1)         delete rough1;
  if(panel1)         delete panel1;
  if(antarctica)     delete antarctica;

  return;
}





// void icemc::EventGenerator::applyRoughness(const Settings& settings1, const int& inu, Interaction* interaction1,
// 					   RayTracer* ray1, Screen* panel1, Antarctica* antarctica,
// 					   Balloon* bn1, const AskaryanFreqsGenerator* askFreqGen, Anita* anita1, const ShowerProperties& showerProps){
  
//   //(vector) ray1->nsurf_rfexit:  surface normal at RFexit position
//   //(pos)        ->rfexit[2]:     final iterated position of RF exit
//   //(vector)     ->n_exit2bn[2]:  vector from RF exit position TO balloon
//   //(pos)    fDetector->r_bn:           position of balloon
//   //(vector) n_pol:               polarization vector
//   //(pos)    posnu:               position of neutrino interaction

//   int num_validscreenpoints = 0;
//   Geoid::Position pos_current;
//   TVector3 vec_pos_current_to_balloon;

//   Geoid::Position pos_projectedImpactPoint;
//   TVector3 vec_localnormal;         //normalized, normal vector at projected ground point
//   TVector3 vec_nnu_to_impactPoint;  //normalized
//   TVector3 vec_inc_perp;            //normalized, vector perp. to incident and surface normal (out-of-inc place)
//   TVector3 vec_inc_parl;            //normalized, vector parl. to incident and surface normal (in-inc plane)
//   double pol_perp_inc, pol_parl_inc;  //component of incident polarization
//   TVector3 vec_local_grnd_perp;     //normalized, vector perp. to transmitted and surface normal (out-of-trans place)
//   TVector3 vec_local_grnd_parl;     //normalized, vector parl. to transmitted and surface normal (in-trans plane)
//   double pol_perp_trans, pol_parl_trans;  //component of transmitted polarization
//   TVector3 vec_grndcomp2bln;
//   TVector3 vec_grndcomp2IP;

//   double time_reference_specular, time_reference_local;
//   double pathlength_local;        // set for each screen point
//   double viewangle_local;
//   double azimuth_local;           // azimuthal angle between local surface normal and vector to balloon [radians]
//   double theta_local;             // polar angle between local surface normal and vector to balloon [radians]
//   double theta_0_local;                 //angle between local surface normal and incident direction [radians]
//   double tcoeff_perp_polperp, tcoeff_parl_polperp;  //for perpendicular polarization (in-ground comp)
//   double tcoeff_perp_polparl, tcoeff_parl_polparl;  //for parallel polarization
//   double power_perp_polperp, power_parl_polperp;
//   double power_perp_polparl, power_parl_polparl;
//   power_perp_polperp = power_parl_polperp = power_perp_polparl = power_parl_polparl = 0.;
//   double fresnel_r, mag_r;

//   TVector3 npol_local_inc, npol_local_trans;
//   TVector3 temp_a;

//   double Emag_local;
//   double taperfactor;
//   //cerr<<inu<<": "<<vmmhz1m_max<<std::endl;

//   //double pathlength_specular = interaction1->posnu.Distance(ray1->rfexit[2]) + ray1->rfexit[2].Distance(fDetector->r_bn);
//   if (settings1.FIRN){
//     time_reference_specular = (interaction1->posnu.Distance(ray1->rfexit[2])*constants::NFIRN / constants::CLIGHT) + (ray1->rfexit[2].Distance(fDetector->position())/constants::CLIGHT);
//   }
//   else{
//     time_reference_specular = (interaction1->posnu.Distance(ray1->rfexit[2])*constants::NICE / constants::CLIGHT) + (ray1->rfexit[2].Distance(fDetector->position())/constants::CLIGHT);
//   }

//   double slopeyx, slopeyy, slopeyz, rtemp;
//   TVector3 ntemp2;
//   TVector3 xaxis = TVector3(1.,0.,0.);
//   TVector3 yaxis = TVector3(0.,1.,0.);
//   TVector3 zaxis = TVector3(0.,0.,1.);

//   double basescreenedgelength = settings1.SCREENEDGELENGTH;
//   double grd_stepsize = settings1.SCREENSTEPSIZE;
//   int grd_nsteps;
//   if(settings1.ROUGHSIZE>0)
//     grd_nsteps = int(basescreenedgelength/2. / grd_stepsize);
//   else
//     grd_nsteps = 0;

//   //#########
//   //iterate points on the screen, get their position and project back to find ground impact
//   //calculate incident and transmitted angles, look up power fraction, and add to running total

//   //reset
//   panel1->ResetParameters();

//   panel1->SetNsamples( grd_nsteps );
//   panel1->SetEdgeLength( basescreenedgelength );

//   panel1->SetCentralPoint( ray1->rfexit[2] );
//   vec_localnormal = antarctica->GetSurfaceNormal(ray1->rfexit[2]).Unit();
//   panel1->SetNormal( vec_localnormal );
//   panel1->SetCosineProjectionFactor( 1. );

//   panel1->SetUnitY( (vec_localnormal.Cross(ray1->n_exit2bn[2])).Unit() );
//   panel1->SetUnitX( (panel1->GetUnitY().Cross(vec_localnormal)).Unit() );
//   //cerr<<panel1->GetCentralPoint()<<"  "<<  fDetector->r_bn<<std::endl;
//   // loop over grid point on ground and see if it's valid
//   for (int ii= -2*panel1->GetNsamples(); ii< 2*panel1->GetNsamples()+1; ii++){
//     for (int jj= -2*panel1->GetNsamples(); jj< 2*panel1->GetNsamples()+1; jj++){
//       //cerr<<"+++++++++++++"<<std::endl;
//       //cerr<<inu<<": "<<ii<<"  "<<jj<<std::endl;
//       //cerr<<"+ seed point: "<<jj<<" / "<<panel1->GetNsamples()*panel1->GetNsamples()<<std::endl;
//       Emag_local = vmmhz1m_max;
//       taperfactor = fresnel_r = mag_r =  1.;
//       tcoeff_perp_polparl = tcoeff_parl_polparl = 0.;
//       tcoeff_perp_polperp = tcoeff_parl_polperp = 0.;
//       pos_projectedImpactPoint = panel1->GetPosition(ii, jj);        // this gets the new screen position
//       vec_pos_current_to_balloon = fDetector->position() - pos_projectedImpactPoint;
//       // local angles of transmission and incidence in their respective planes
//       vec_localnormal = antarctica->GetSurfaceNormal(pos_projectedImpactPoint).Unit();
//       if (settings1.SLOPEY) {
// 	slopeyx=ray1->slopeyx;
// 	slopeyy=ray1->slopeyy;
// 	slopeyz=ray1->slopeyz;
// 	ntemp2 = vec_localnormal + slopeyx*xaxis + slopeyy*yaxis + slopeyz*zaxis;
// 	ntemp2 = ntemp2.Unit();
// 	rtemp= ntemp2.Dot(vec_localnormal);
// 	if (rtemp<=1) {
// 	  vec_localnormal = ntemp2;
// 	}//if
//       }//end local slopeyness
//       //cerr<<inu<<"  "<<pos_projectedImpactPoint<<std::endl;
//       vec_nnu_to_impactPoint =  TVector3( pos_projectedImpactPoint[0]-interaction1->posnu[0], pos_projectedImpactPoint[1]-interaction1->posnu[1], pos_projectedImpactPoint[2]-interaction1->posnu[2] ).Unit();

//       vec_grndcomp2IP = (vec_nnu_to_impactPoint - (vec_nnu_to_impactPoint.Dot(vec_localnormal)*vec_localnormal)).Unit();
//       vec_grndcomp2bln = (vec_pos_current_to_balloon - (vec_pos_current_to_balloon.Dot(vec_localnormal)*vec_localnormal)).Unit();
//       temp_a = vec_localnormal.Cross(vec_pos_current_to_balloon).Unit();
//       azimuth_local = vec_grndcomp2IP.Angle(vec_grndcomp2bln); //[rad]
//       if( temp_a.Dot(vec_nnu_to_impactPoint) < 0 )
// 	azimuth_local *= -1.;
//       if( panel1->GetCentralPoint().Distance(pos_projectedImpactPoint)<0.75*grd_stepsize ){
// 	azimuth_local = 0.;
//       }
//       //cerr<<inu<<":  "<<jj<<"  "<<vec_grndcomp2IP<<" : "<<vec_grndcomp2bln<<" : "<<azimuth_local*180./PI<<std::endl;
//       theta_local = vec_localnormal.Angle( (const TVector3)vec_pos_current_to_balloon ); //[rad]
//       theta_0_local = vec_localnormal.Angle(vec_nnu_to_impactPoint); //[rad]
//       //cerr<<inu<<"  "<<ii<<"  "<<jj<<";  "<<panel1->GetCentralPoint()<<" : "<<  pos_projectedImpactPoint<<" : "<<theta_local*180./PI<<"  "<<theta_0_local*180./PI<<"  "<< azimuth_local*180./PI<< std::endl;
//       //cerr<< panel1->GetCentralPoint() - pos_projectedImpactPoint<<std::endl;
//       if( isnan(theta_local) | isnan(theta_0_local) | isnan(azimuth_local) ){
// 	continue;
//       }
//       viewangle_local = GetViewAngle(vec_nnu_to_impactPoint, interaction1->nnu);

//       // at this point, only figure out if taper will kill the geometry, but don't actually apply the factor
//       deltheta_em[0] = deltheta_em_max*fDetector->FREQ_LOW/fDetector->freq[0];
//       deltheta_had[0] = deltheta_had_max*fDetector->FREQ_LOW/fDetector->freq[0];

//       // askFreqGen->TaperVmMHz(viewangle_local, deltheta_em[0], deltheta_had[0], showerProps, taperfactor, vmmhz_em[0]);// this applies the angular dependence.
//       double dummy_vmmhz_em_0 = 0;
//       askFreqGen->TaperVmMHz(viewangle_local, deltheta_em[0], deltheta_had[0], showerProps, taperfactor, dummy_vmmhz_em_0);// this applies the angular dependence.
//       if(taperfactor==0){
// 	continue;
//       }
//       //cerr<<inu<< ": past E=0"<<std::endl;
//       /////
//       // Field Magnitude
// #ifdef USE_HEALPIX
//       if (settings1.FIRN)
// 	rough1->InterpolatePowerValue(power_perp_polperp, power_parl_polperp, power_perp_polparl, power_parl_polparl, theta_0_local*180./PI, theta_local*180./PI, azimuth_local *180./PI);
//       else
// 	rough1->InterpolatePowerValue(power_perp_polperp, power_parl_polperp, power_perp_polparl, power_parl_polparl, theta_0_local*180./PI, theta_local*180./PI, azimuth_local *180./PI);
// #endif
//       //cerr<<"P: "<<power_perp<<"  "<<power_parl<<std::endl;
//       if( (power_perp_polperp==0.)&(power_parl_polperp==0.)&(power_perp_polparl==0.)&(power_parl_polparl==0.) ){
// 	//continue;
//       }
//       //cerr<<"survived power cut"<<std::endl;
//       if (settings1.FIRN){
// 	tcoeff_perp_polparl = sqrt(power_perp_polparl);//*NFIRN*cos(theta_0_local)*cos(theta_local));
// 	tcoeff_parl_polparl = sqrt(power_parl_polparl);//*NFIRN*cos(theta_0_local)*cos(theta_local));
// 	tcoeff_perp_polperp = sqrt(power_perp_polperp);//*NFIRN*cos(theta_0_local)*cos(theta_local));
// 	tcoeff_parl_polperp = sqrt(power_parl_polperp);//*NFIRN*cos(theta_0_local)*cos(theta_local));
//       }
//       else{
// 	tcoeff_perp_polparl = sqrt(power_perp_polparl);//*NFIRN*cos(theta_0_local)*cos(theta_local));
// 	tcoeff_parl_polparl = sqrt(power_parl_polparl);//*NFIRN*cos(theta_0_local)*cos(theta_local));
// 	tcoeff_perp_polperp = sqrt(power_perp_polperp);//*NFIRN*cos(theta_0_local)*cos(theta_local));
// 	tcoeff_parl_polperp = sqrt(power_parl_polperp);//*NFIRN*cos(theta_0_local)*cos(theta_local));
//       }
//       //
//       //cerr<<"T: "<<tcoeff_perp<<"  "<<tcoeff_parl<<std::endl;
//       //cerr<<"T (spec): "<<t_coeff_slappy<<"  "<<t_coeff_pokey<<std::endl;
//       //Emag_local *= sqrt((tcoeff_perp*tcoeff_perp + tcoeff_parl*tcoeff_parl)) * mag_r;// * (antennalength*antennalength/(vec_pos_current_to_balloon.Mag()*vec_pos_current_to_balloon.Mag()))/HP_64_binarea);
//       //cerr<<"E: "<<Emag_local<<std::endl;
//       // account for 1/r for 1)interaction point to impact point and 2)impact point to balloon, and attenuation in ice
//       pathlength_local = interaction1->posnu.Distance(pos_projectedImpactPoint) + pos_projectedImpactPoint.Distance(fDetector->position());
//       //cerr<<"P: "<<pathlength_local<<std::endl;
//       Emag_local /= pathlength_local ;
//       //cerr<<"E: "<<Emag_local<<std::endl;
//       Attenuate(antarctica, &settings1, Emag_local,  interaction1->posnu.Distance(pos_projectedImpactPoint),  interaction1->posnu);
//       //cerr<<"E: "<<Emag_local<<std::endl;
//       /////
//       // Incident and Transmitted Polarizations
//       // set incident polarization
//       npol_local_inc = GetPolarization(interaction1->nnu, vec_nnu_to_impactPoint, inu).Unit();
//       vec_inc_perp = (vec_localnormal.Cross(vec_nnu_to_impactPoint)).Unit();
//       vec_inc_parl = (vec_nnu_to_impactPoint.Cross(vec_inc_perp)).Unit();
//       pol_perp_inc = npol_local_inc.Dot(vec_inc_perp);
//       pol_parl_inc = npol_local_inc.Dot(vec_inc_parl);
//       //
//       pol_perp_trans = pol_perp_inc * tcoeff_perp_polperp + pol_parl_inc * tcoeff_perp_polparl;
//       pol_parl_trans = pol_parl_inc * tcoeff_parl_polparl + pol_perp_inc * tcoeff_parl_polperp;
//       //
//       vec_local_grnd_perp = (vec_localnormal.Cross(vec_pos_current_to_balloon)).Unit();
//       vec_local_grnd_parl = (vec_pos_current_to_balloon.Cross(vec_local_grnd_perp)).Unit();
//       //
//       // set transmitted polarization
//       npol_local_trans= (pol_perp_trans*vec_local_grnd_perp + pol_parl_trans*vec_local_grnd_parl).Unit();
//       //cerr<<inu<<":  v_nu "<<interaction1->nnu<<" : 2IP "<<vec_nnu_to_impactPoint<<" : npol "<<npol_local_trans<< std::endl;
//       // check if transmitted polarization is undefined
//       if( isnan(npol_local_trans[0]) ){
// 	continue;
//       }
//       //cerr<<"past pol cut"<<std::endl;
//       //
//       fresnel_r = sqrt( pow(vmmhz1m_max*pol_perp_trans,2) + pow(vmmhz1m_max*pol_parl_trans,2) ) / vmmhz1m_max;
//       mag_r = sqrt( tan(theta_0_local) / tan(theta_local) );
//       Emag_local *= fresnel_r * mag_r;
//       //cerr<<"E: "<<Emag_local<<std::endl;
//       if (settings1.FIRN){
// 	time_reference_local = (interaction1->posnu.Distance(pos_projectedImpactPoint)*constants::NFIRN / constants::CLIGHT) + (pos_projectedImpactPoint.Distance(fDetector->position())/constants::CLIGHT);
//       }
//       else{
// 	time_reference_local = (interaction1->posnu.Distance(pos_projectedImpactPoint)*constants::NICE / constants::CLIGHT) + (pos_projectedImpactPoint.Distance(fDetector->position())/constants::CLIGHT);
//       // increment counter so we can track the size of the screen's vector arrays
//       }
//       num_validscreenpoints++;

//       //add the contribution to the running total
//       panel1->AddVmmhz0(Emag_local);                // pre-taper Efield
//       panel1->AddVec2bln(vec_pos_current_to_balloon);
//       panel1->AddPol(npol_local_trans);
//       panel1->AddDelay( time_reference_specular - time_reference_local );
//       panel1->AddImpactPt(pos_projectedImpactPoint);
//       panel1->AddViewangle(viewangle_local);
//       panel1->AddIncidenceAngle(theta_0_local);
//       panel1->AddTransmissionAngle(theta_local);
//       panel1->AddWeight( (panel1->GetEdgeLength() / panel1->GetNsamples()) * (panel1->GetEdgeLength() / panel1->GetNsamples()) );
//       panel1->AddFacetLength(panel1->GetEdgeLength() / panel1->GetNsamples());
//       panel1->AddTparallel_polParallel(tcoeff_parl_polparl);
//       panel1->AddTperpendicular_polParallel(tcoeff_perp_polparl);
//       panel1->AddTparallel_polPerpendicular(tcoeff_parl_polperp);
//       panel1->AddTperpendicular_polPerpendicular(tcoeff_perp_polperp);
//     }
//   }

//   panel1->SetNvalidPoints(num_validscreenpoints);
//   //now construct the Screen's vmmhz array for all points, so it gets passed to the trigger object later to make the waveforms
//   // here we get the array vmmhz by taking vmmhz1m_max (signal at lowest frequency bin) and vmmhz_max (signal at lowest frequency after applying 1/r factor and attenuation factor) and making an array across frequency bins by putting in frequency dependence.
//   double validScreenSummedArea = 0.;
//   double vmmhz_local_array[Anita::NFREQ];
//   // AskaryanFreqs askFreqsLocal;
//   for (int jj=0; jj<panel1->GetNvalidPoints(); jj++){
//     // fill the frequency array vmmhz_local_array
//     askFreqGen->GetVmMHz(panel1->GetVmmhz0(jj), vmmhz1m_max, pnu, fDetector->freq, fDetector->NOTCH_MIN, fDetector->NOTCH_MAX, vmmhz_local_array, Anita::NFREQ);
//     // apply the off-angle tapering
//     for (int k=0;k<Anita::NFREQ;k++) {
//       deltheta_em[k]=deltheta_em_max*fDetector->FREQ_LOW/fDetector->freq[k];
//       deltheta_had[k]=deltheta_had_max*fDetector->FREQ_LOW/fDetector->freq[k];
//       // askFreqGen->TaperVmMHz(panel1->GetViewangle(jj), deltheta_em[k], deltheta_had[k], emfrac, hadfrac, vmmhz_local_array[k], vmmhz_em[k]);// this applies the angular dependence.
//       double dummy_vmmhz_em_k = 0;
//       askFreqGen->TaperVmMHz(panel1->GetViewangle(jj), deltheta_em[k], deltheta_had[k], showerProps.emFrac, showerProps.hadFrac, vmmhz_local_array[k], dummy_vmmhz_em_k);// this applies the angular dependence.
//       panel1->AddVmmhz_freq(vmmhz_local_array[k]);
//     }

//     validScreenSummedArea += panel1->GetWeight(jj);
//   }//end jj over panel Nvalid points
//   panel1->SetWeightNorm(validScreenSummedArea);
//   vmmhz_max = 0.;
//   for(int jj=0; jj<panel1->GetNvalidPoints(); jj++){
//     vmmhz_max = Tools::dMax(vmmhz_max, panel1->GetVmmhz_freq(jj*Anita::NFREQ));
//   }
// }
