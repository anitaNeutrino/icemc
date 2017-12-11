#include <array>
#include "vector.hh"
#include "position.hh"

#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "earthmodel.hh"
#include "icemodel.hh"
#include "TF1.h"

#include "rx.hpp"
#include "anita.hh"

#include "ChanTrigger.h"
#include "balloon.hh"
#include <iostream>
#include <cmath>
#include <vector>

#include <string>

#include "Tools.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include "TMath.h"

#include "EnvironmentVariable.h"

#ifdef ANITA_UTIL_EXISTS
#include "FFTtools.h"
#include "AnitaEventCalibrator.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#endif

const std::string ICEMC_SRC_DIR=EnvironmentVariable::ICEMC_SRC_DIR();
const std::string ICEMC_DATA_DIR=ICEMC_SRC_DIR+"/data/";


using std::cout;
using std::cin;
using std::endl;
using std::vector;

Anita::Anita() {

  stemp = "";

  for (int irx=0;irx<Anita::NLAYERS_MAX*Anita::NPHI_MAX;irx++) {
    arrival_times[0][irx]=arrival_times[1][irx]=0.;
  }
  NCH_PASS=4;
  rx_minarrivaltime=0.;

  Tools::Zero(VNOISE_ANITALITE,16);
  INCLINE_TOPTHREE = 10.; // cant angle of top three layers of antennas
  INCLINE_NADIR = 10.; // cant angle of nadir (bottom) layer of antennas
  SIGMA_THETA=0.5*RADDEG; // resolution on the polar angle of the signal
    
  FREQ_LOW=0.;//200.E6;
  FREQ_HIGH=1300.E6; //1200.E6;
    
  antennatosurf[0]=2;
  antennatosurf[1]=4;
  antennatosurf[2]=6;
  antennatosurf[3]=8;
  antennatosurf[4]=2;
  antennatosurf[5]=4;
  antennatosurf[6]=6;
  antennatosurf[7]=8;
  // layer 1
  antennatosurf[8]=1;
  antennatosurf[9]=3;
  antennatosurf[10]=5;
  antennatosurf[11]=7;
  antennatosurf[12]=1;
  antennatosurf[13]=3;
  antennatosurf[14]=5;
  antennatosurf[15]=7;
  // layer 2
  antennatosurf[16]=1;
  antennatosurf[17]=2;
  antennatosurf[18]=3;
  antennatosurf[19]=4;
  antennatosurf[20]=5;
  antennatosurf[21]=6;
  antennatosurf[22]=7;
  antennatosurf[23]=8;
  antennatosurf[24]=1;
  antennatosurf[25]=2;
  antennatosurf[26]=3;
  antennatosurf[27]=4;
  antennatosurf[28]=5;
  antennatosurf[29]=6;
  antennatosurf[30]=7;
  antennatosurf[31]=8; // layer 3
      
  maxthreshold=0.;
  Tools::Zero(bwslice_thresholds,5); // thresholds for each band -- this is just an initialization- this is set in the input file
  for (int i=0;i<5;i++) {
    bwslice_allowed[i]=1; // these bands are allowed to contribute to the trigger sum -- this is set in the input file
  }
  bwslice_required[0]=0;
  bwslice_required[1]=0;
  bwslice_required[2]=0;
  bwslice_required[3]=0;
  bwslice_required[4]=1; // these bands are required to be among the channels that pass -- this is set in the input file
  pol_allowed[0]=1;// which polarisations are allowed to have channels that fire (V,H)
  pol_allowed[1]=1;// which polarisations are allowed to have channels that fire (V,H)
  pol_required[0]=1;// which polarisations are required to have channels that fire (V,H)
  pol_required[1]=0;// which polarisations are required to have channels that fire (V,H)
    
    
  bwslice_center[0]=265.e6;
  bwslice_center[1]=435.e6;
  bwslice_center[2]=650.e6;
  bwslice_center[3]=980.e6;
  bwslice_center[4]=682.5E6; // center frequencies
    
  bwslice_width[0]=130.e6;
  bwslice_width[1]=160.e6;
  bwslice_width[2]=250.e6;
  bwslice_width[3]=370.e6;
  bwslice_width[4]=965.E6; // 3 dB bandwidths, without overlap
    
  for (int i=0;i<4;i++) {
    bwslice_min[i]=bwslice_center[i]-bwslice_width[i]/2.;
    bwslice_max[i]=bwslice_center[i]+bwslice_width[i]/2.;
  }
    
  bwslice_min[4]= bwslice_center[0]-bwslice_width[0]/2.; //minimum of each bandwidth slice
  bwslice_max[4]=bwslice_center[3]+bwslice_width[3]/2.; //minimum of each bandwidth slice
    
  bwmin=0.; // minimum width of any allowed bandwidth slice
    
  /*** Used for the Coherent Sum Trigger ***/
  summed_power_trigger_digitizer_zero_random = new TRandom3();
  // Prepare the file and tree for the coherent sum trigger's data tree
  coherent_datafile = new TFile("outputs/coherent_sum_data_file.root","RECREATE");
  coherent_waveform_sum_tree = new TTree("coherent_waveform_sum_tree", "Coherent Waveform Sum");
  coherent_waveform_sum_tree->Branch("event_number", &cwst_event_number);
  coherent_waveform_sum_tree->Branch("center_phi_sector", &cwst_center_phi_sector);
  coherent_waveform_sum_tree->Branch("rms_noise", &cwst_rms_noise);
  coherent_waveform_sum_tree->Branch("actual_rms", &cwst_actual_rms);
  coherent_waveform_sum_tree->Branch("threshold", &cwst_threshold);
  coherent_waveform_sum_tree->Branch("window_start", &cwst_window_start);
  coherent_waveform_sum_tree->Branch("window_end", &cwst_window_end);
    
  coherent_waveform_sum_tree->Branch("deg_theta", &cwst_deg_theta);
  coherent_waveform_sum_tree->Branch("deg_phi", &cwst_deg_phi);
    
  coherent_waveform_sum_tree->Branch("actual_deg_theta", &cwst_actual_deg_theta);
  coherent_waveform_sum_tree->Branch("actual_deg_phi", &cwst_actual_deg_phi);
    
  coherent_waveform_sum_tree->Branch("timesteps", cwst_timesteps);
    
  for (unsigned i = 0; i < 48; ++i) {
    cwst_RXs[i].waveform = new vector <double>(HALFNFOUR, 0.);
    cwst_RXs[i].digitized = new vector <double>(HALFNFOUR, 0.);
    coherent_waveform_sum_tree->Branch(Form("rx%u", i), &(cwst_RXs[i]));
  }
    
  for (unsigned int i = 0; i < 9; i++) {
    cwst_aligned_wfms[i].digitized = new vector <double>(HALFNFOUR, 0.);
    coherent_waveform_sum_tree->Branch(Form("aligned_wfms%u", i), &(cwst_aligned_wfms[i]));
    /*
    //coherent_waveform_sum_tree->Branch(Form("whole_wfms%u", i), &(cwst_whole_wfms[i]));
    cwst_whole_wfms[i] = new vector <double>(HALFNFOUR, 0.);
    cwst_wfms[i] = new vector <double>(HALFNFOUR, 0.);
    cwst_aligned_wfms[i] = new vector <double>(HALFNFOUR, 0.);
    coherent_waveform_sum_tree->Branch(Form("whole_wfms%u", i), &(cwst_whole_wfms[i]));
    coherent_waveform_sum_tree->Branch(Form("wfms%u", i), &(cwst_wfms[i]));
    coherent_waveform_sum_tree->Branch(Form("aligned_wfms%u", i), &(cwst_aligned_wfms[i]));
    */
  }
	
  coherent_waveform_sum_tree->Branch("summed_wfm", &cwst_summed_wfm);
  coherent_waveform_sum_tree->Branch("power_of_summed_wfm", &cwst_power_of_summed_wfm);
  coherent_waveform_sum_tree->Branch("power", &cwst_power);
  // End of preparition for the coherent sum trigger's data tree
}

Anita::~Anita(){
  coherent_datafile->cd();
    
  coherent_waveform_sum_tree->Write();
  coherent_datafile->Write();
    
  coherent_waveform_sum_tree->Delete();
  coherent_datafile->Close();
  coherent_datafile->Delete();
    
  delete summed_power_trigger_digitizer_zero_random;
	
  return;
}

int Anita::Match(int ilayer,int ifold,int rx_minarrivaltime) {
    
  if (ilayer==GetLayer(rx_minarrivaltime) && ifold==GetIfold(rx_minarrivaltime))
    return 1;
  else
    return 0;
    
}
int Anita::GetRx(int ilayer, int ifold) { // get antenna number based on which layer and position it is
    
  int irx=0;
  for (int i=0;i<ilayer;i++) {
    irx+=NRX_PHI[i];
  }
  irx+=ifold;
  return irx;
    
}

int Anita::GetRxTriggerNumbering(int ilayer, int ifold) { // get antenna number based on which layer and position it is
  // make the top trigger layer count 1-16 left to right
  if (ilayer==0)
    //cout << "ilayer, ifold, getrx are " << ilayer << "\t" << ifold << "\t" << 2*ifold+ilayer << "\n";
    return 2*ifold;
  else if(ilayer==1) {
    return 2*ifold+1;
  }
  else {
    //cout << "ilayer, ifold, getrx are " << ilayer << "\t" << ifold << "\t" << GetRx(ilayer,ifold) << "\n";
    return GetRx(ilayer,ifold);
  }
}

void Anita::SetNoise(Settings *settings1,Balloon *bn1,IceModel *antarctica) {
    
  // these should only be used for the frequency domain trigger.
  if (settings1->WHICH==2 || settings1->WHICH==6) { //this is for anita 1
    for (int il=0;il<NLAYERS_MAX;il++) {
      bwslice_vnoise[il][0]=5.5125E-6;
      bwslice_vnoise[il][1]=6.1157E-6;
      bwslice_vnoise[il][2]=7.64446E-6;
      bwslice_vnoise[il][3]=9.29989E-6;
      bwslice_vnoise[il][4]=0.;
      for (int i=0;i<4;i++) {
	bwslice_vnoise[il][4]+=bwslice_vnoise[il][i]*bwslice_vnoise[il][i];
      }
      bwslice_vnoise[il][4]=sqrt(bwslice_vnoise[il][4]);
    }
  }
  else { // if not anita 1
    for (int il=0;il<NLAYERS_MAX;il++) {
      for (int ibw=0;ibw<5;ibw++) {
	bwslice_vnoise[il][ibw]=ChanTrigger::GetNoise(settings1,bn1->altitude_bn,antarctica->SurfaceAboveGeoid(bn1->latitude,bn1->longitude),THETA_ZENITH[il],bwslice_max[ibw]-bwslice_min[ibw],0.);
      }
    }
  }
    
    
    
}
void Anita::Initialize(Settings *settings1,ofstream &foutput,int thisInu, TString outputdir)
{
    
    
  count_getnoisewaveforms=0;
  rms_lab[0]=rms_lab[1]=0.;
  rms_rfcm[0]=rms_rfcm[1]=0.;
    
  NBANDS=4; // subbands (not counting full band)
  inu=thisInu;
    
  PERCENTBW=10; // subbands (not counting full band)

  TIMESTEP=(1./2.6)*1.E-9; // time step between samples

  for (int i=0;i<HALFNFOUR;i++)   fTimes[i] = i * TIMESTEP * 1.0E9; 
 
  for (int i=0;i<NFREQ;i++) {
    freq[i]=FREQ_LOW+(FREQ_HIGH-FREQ_LOW)*(double)i/(double)NFREQ; // freq. of each bin.
    avgfreq_rfcm[i]=0.;
    avgfreq_rfcm_lab[i]=0.;
  } //for

  initializeFixedPowerThresholds(foutput);
  
  /// TEMP HACK FOR ANITA-4 !!!!
  if (settings1->WHICH==10) powerthreshold[4] /= TMath::Sqrt(2.);
  
  if (settings1->TRIGGERSCHEME==5)
    l1window=3.75E-9;
  else
    l1window=11.19E-9; // l1 coincidence window
    
  minsignalstrength=0.1;
    
  impedence=50.;
  phase=90.; // phase for positive frequencies
  // assuming v(t) is real, then phase(neg)=-phase(pos)
  INTEGRATIONTIME=3.5E-9; // integration time of trigger diode

  DEADTIME=10.E-9; // dead time after a trigger
  energythreshold=3.;  // power threshold

    
  // for ANITA 3 trigger
  TRIG_TIMESTEP=3.75E-9;
  N_STEPS_PHI = 100;
  N_STEPS_THETA = 100;
  MIN_PHI_HYPOTHESIS=-22.5;
  MAX_PHI_HYPOTHESIS=22.5;
  MIN_THETA_HYPOTHESIS=-50.;
  MAX_THETA_HYPOTHESIS=20.;
  //MIN_THETA_HYPOTHESIS=-55.;
  //MAX_THETA_HYPOTHESIS=20.;

  double freqstep=1./(double)(NFOUR/2)/TIMESTEP;
  // double freqstep_long=1./(double)(NFOUR)/TIMESTEP; // freqstep_long is actually shorter so that time domain waveform is longer
    
  for (int i=0;i<HALFNFOUR;i++) {
    freq_forfft[2*i]=(double)i*freqstep;
    freq_forfft[2*i+1]=(double)i*freqstep;
      
    if (i<HALFNFOUR/2) freq_forplotting[i]=freq_forfft[2*i];
		
  }


  readVariableThresholds(settings1);

  phiTrigMask=0;
  phiTrigMaskH=0;
  l1TrigMask=0;
  l1TrigMaskH=0;


  readAmplification();    
    
  // get lab attenuation data
  getLabAttn(NPOINTS_LAB,freqlab,labattn);
    

  getDiodeDataAndAttenuation(settings1, outputdir);

  if (settings1->PULSER==1) {
    getPulserData();
  } 
    
  // for antenna gains
  reference_angle[0]=0.;
  reference_angle[1]=5.;
  reference_angle[2]=10.;
  reference_angle[3]=20.;
  reference_angle[4]=30.;
  reference_angle[5]=45.;
  reference_angle[6]=90.; // reference angles for finding gains of antennas

  
  THERMALNOISE_FACTOR=settings1->THERMALNOISE_FACTOR;
  for (int j=0;j<settings1->NLAYERS;j++) {
    // noise depends on cant angle of antenna
    //     //   VNOISE[j]=ChanTrigger::GetNoise(altitude_bn,surface_under_balloon,THETA_ZENITH[j],BW_SEAVEYS,0.);
    VNOISE[j]=1.52889E-5; // this comes from V^2/R=kT*bw -> V=sqrt(kT*bw*R)
    VNOISE[j]*=THERMALNOISE_FACTOR;
  }//for

#ifdef ANITA_UTIL_EXISTS
  if (settings1->NOISEFROMFLIGHTDIGITIZER || settings1->NOISEFROMFLIGHTTRIGGER){
    readNoiseFromFlight(settings1);
  }
  if (settings1->APPLYIMPULSERESPONSEDIGITIZER){
    readImpulseResponseDigitizer(settings1);
    if(settings1->TUFFSON){
      readTuffResponseDigitizer(settings1);
      readTuffResponseTrigger(settings1);
    }
  }
  if (settings1->APPLYIMPULSERESPONSETRIGGER){
    readImpulseResponseTrigger(settings1);
  }
  if (settings1->TRIGGEREFFSCAN){
    readTriggerEfficiencyScanPulser(settings1);
  }
#endif

  setDiodeRMS(settings1, outputdir);

  
  // Setting up output files
  
  string stemp=string(outputdir.Data())+"/signals.root";
  fsignals=new TFile(stemp.c_str(),"RECREATE");
  tsignals=new TTree("tsignals","tsignals");
    
  stemp=string(outputdir.Data())+"/data.root";
  fdata=new TFile(stemp.c_str(),"RECREATE");
  tdata=new TTree("tdata","tdata");
    
    
  tsignals->Branch("signal_vpol_inanita",&signal_vpol_inanita,"signal_vpol_inanita[5][512]/D");
  tsignals->Branch("timedomainnoise_rfcm_banding",&timedomainnoise_rfcm_banding,"timedomainnoise_rfcm_banding[2][5][512]/D");
  tsignals->Branch("total_vpol_inanita",&total_vpol_inanita,"total_vpol_inanita[5][512]/D");
  tsignals->Branch("total_diodeinput_1_inanita",&total_diodeinput_1_inanita,"total_diodeinput_1_inanita[5][512]/D"); // this is the waveform that is input to the tunnel diode in the first (LCP or vertical) polarization
  tsignals->Branch("total_diodeinput_2_inanita",&total_diodeinput_2_inanita,"total_diodeinput_2_inanita[5][512]/D"); // this is the waveform that is input to the tunnel diode in the first (RCP or horizontal) polarization
  tsignals->Branch("timedomain_output_corrected_forplotting",&timedomain_output_corrected_forplotting,"timedomain_output_1_corrected_forplotting[2][6][512]/D"); 
  tsignals->Branch("timedomain_output_inanita",&timedomain_output_inanita,"timedomain_output_inanita[2][5][512]/D");
    
    
  tsignals->Branch("peak_rx_rfcm_lab",&peak_rx_rfcm_lab,"peak_rx_rfcm_lab[2]/D");
  tsignals->Branch("inu",&inu,"inu/I");
  tsignals->Branch("dangle",&dangle_inanita,"dangle/D");
  tsignals->Branch("emfrac",&emfrac_inanita,"emfrac/D");
  tsignals->Branch("hadfrac",&hadfrac_inanita,"hadfrac/D");
  tsignals->Branch("ston",&ston,"ston[5]/D");

  tsignals->Branch("peak",                     &peak_v_banding_rfcm,      "peak_v_banding_rfcm[2][5]/D"   );
  tsignals->Branch("peak_rx",                  &peak_rx_signalonly,       "peak_rx[2]/D"                  );
  tsignals->Branch("peak_rx_rfcm",             &peak_rx_rfcm,             "peak_rx_rfcm[2]/D"             );
  tsignals->Branch("peak_rx_rfcm_signalonly",  &peak_rx_rfcm_signalonly,  "peak_rx_rfcm_signalonly[2]/D"  );
  tsignals->Branch("peak_rx_rfcm_lab",         &peak_rx_rfcm_lab,         "peak_rx_rfcm_lab[2]/D"         );
  tsignals->Branch("bwslice_vrms",&bwslice_vrms,"bwslice_vrms[5]/D");
  tsignals->Branch("iminbin",&iminbin,"iminbin[5]/I");
  tsignals->Branch("imaxbin",&imaxbin,"imaxbin[5]/I");
  tsignals->Branch("maxbin_fortotal",&maxbin_fortotal,"maxbin_fortotal[5]/I");
  tsignals->Branch("channels_passing",&channels_passing,"channels_passing[2][5]/I");
  tsignals->Branch("bwslice_rmsdiode",&bwslice_rmsdiode,"bwslice_rmsdiode[5]/D");
  tsignals->Branch("l1_passing",&l1_passing,"l1_passing/I");
  tsignals->Branch("integral_vmmhz",&integral_vmmhz_foranita,"integral_vmmhz/D");
  //tsignals->Branch("dnutries",&dnutries,"dnutries/D");
  //tsignals->Branch("weight_test",&weight_test,"weight_test/D");
  tsignals->Branch("flag_e_inanita",&flag_e_inanita,"flag_e_inanita[5][512]/I");
  tsignals->Branch("flag_h_inanita",&flag_h_inanita,"flag_h_inanita[5][512]/I");
    
    
  tdata=new TTree("tdata","tdata");
  tdata->Branch("total_diodeinput_1_allantennas",&total_diodeinput_1_allantennas,"total_diodeinput_1_allantennas[48][512]/D"); // this is the waveform that is input to the tunnel diode in the first (LCP or vertical) polarization
  tdata->Branch("total_diodeinput_2_allantennas",&total_diodeinput_2_allantennas,"total_diodeinput_2_allantennas[48][512]/D"); // this is the waveform that is input to the tunnel diode in the first (LCP or vertical) polarization
  tdata->Branch("timedomain_output_allantennas",&timedomain_output_allantennas,"timedomain_output_allantennas[2][48][512]/D"); // this is the waveform that is output to the tunnel diode in the first (LCP or vertical) polarization
  tdata->Branch("arrival_times",&arrival_times,"arrival_times[2][48]/D");
  tdata->Branch("inu",&inu,"inu/I");
  tdata->Branch("powerthreshold",&powerthreshold,"powerthreshold[5]/D");
  tdata->Branch("bwslice_rmsdiode",&bwslice_rmsdiode,"bwslice_rmsdiode[5]/D");

  //std::array< std::array< std::array< std::array<std::vector<int>,5>, 2>, 16>, 3>  arrayofhits_inanita; 

  tdata->Branch("arrayofhits_inanita",&arrayofhits_inanita,"arrayofhits_inanita[3][16][2][512]/I");

  tdata->Branch("l1trig_anita3and4_inanita",&l1trig_anita3and4_inanita,"l1trig_anita3and4_inanita[2][16][512]/I");

  tdata->Branch("l1trig_anita4lr_inanita",&l1trig_anita4lr_inanita,"l1trig_anita4lr_inanita[3][16][512]/I");


  tdata->Branch("l2trig_anita4lr_inanita",&l2trig_anita4lr_inanita,"l2trig_anita4lr_inanita[16][3][512]/I");

  tdata->Branch("l3type0trig_anita4lr_inanita",&l3type0trig_anita4lr_inanita,"l3type0trig_anita4lr_inanita[16][512]/I");
  tdata->Branch("l3trig_anita4lr_inanita",&l3trig_anita4lr_inanita,"l3trig_anita4lr_inanita[16][512]/I");


  //tdata->Branch("arrayofhits_inanita",&arrayofhits_inanita,"std::array< std::array< std::array< std::array<std::vector<int>,5>, 2>, 16>, 3>");
  tdata->Branch("passglobtrig",&passglobtrig,"passglobtrig[2]/I");
    

  tgaryanderic=new TTree("tgaryanderic","tgaryanderic");
  tgaryanderic->Branch("arrayofhits",&arrayofhits_forgaryanderic,"arrayofhits_forgaryanderic[3][16][2][512]/I");
  tgaryanderic->Branch("l1trig",&l1trig_anita4lr_forgaryanderic,"l1trig_anita4lr_forgaryanderic[3][16][512]/I");

  tgaryanderic->Branch("l2trig",&l2trig_anita4lr_forgaryanderic,"l2trig_anita4lr_forgaryanderic[16][512]/I"); 
  tgaryanderic->Branch("l3type0trig",&l3type0trig_anita4lr_forgaryanderic,"l3type0trig_anita4lr_forgaryanderic[16][512]/I"); 
  tgaryanderic->Branch("l3type1trig",&l3type1trig_anita4lr_forgaryanderic,"l3type1trig_anita4lr_forgaryanderic[16][512]/I"); 
  tgaryanderic->Branch("passglobtrig",&passglobtrig,"passglobtrig[2]/I");
  tgaryanderic->Branch("weight",&weight_inanita,"weight_inanita/D");
  tgaryanderic->Branch("time",&time_trig,"time_trig[512]/D");


  tglob=new TTree("tglob","tglob");
  tglob->Branch("inu",&inu,"inu/I");
  tglob->Branch("passglobtrig",&passglobtrig,"passglobtrig[2]/I");
  tglob->Branch("l1_passing_allantennas",&l1_passing_allantennas,"l1_passing_allantennas[48]/I");
  
  
    
}

void Anita::initializeFixedPowerThresholds(ofstream &foutput){
    
  if (BANDING==2) { //anita 2
		
    // I derive these thresholds from Ryan's plot
    //of measured rates vs. threshold
    // that appears on p. 9 of his talk at the Anita meeting 19th Feb 2008
    //  Using 14 MHz, 8 MHz, 8MHz and 1 MHz for L, M, H and full bands
    powerthreshold[0]=-1.87; // low band
    powerthreshold[1]=-2.53; // middle band
    powerthreshold[2]=-2.15; // high band
    powerthreshold[3]=-1.; // not used for Anita 2
    powerthreshold[4]=-4.41; // full band - use this threshold when other bands are used
    //powerthreshold[4]=-5.3; // full band - use this threshold for full band only trigger
		
		
    powerthreshold_nadir[0]=-6.7; // low band
    powerthreshold_nadir[1]=-6.5; // middle band
    powerthreshold_nadir[2]=-5.1; // high band
    powerthreshold_nadir[3]=-1.; // not used for Anita 2
    powerthreshold_nadir[4]=-6.7; // full band - use this threshold when other bands are used
		
		
    foutput << "Thresholds are (in p/<p>):  " <<
      powerthreshold[0] << " (L)\t" <<
      powerthreshold[1] << " (M)\t" <<
      powerthreshold[2] << " (H)\t" <<
      powerthreshold[4] << " (F)\n";
  }
  else if (BANDING==0 || BANDING==1) { // anita 1 or set your own
		
    powerthreshold[0]=-3.27;
    powerthreshold[1]=-3.24;
    powerthreshold[2]=-2.48;
    powerthreshold[3]=-2.56;
    powerthreshold[4]=-3.;
		
    foutput << "Thresholds are (in p/<p>):  " <<
      powerthreshold[0] << " (L)\t" <<
      powerthreshold[1] << " (M1)\t" <<
      powerthreshold[2] << " (M2)\t" <<
      powerthreshold[3] << " (H)\t" <<
      powerthreshold[4] << " \n";
		
  } else if (BANDING==4 || BANDING==5){ // anita-3
    powerthreshold[0]=-1; // not used 
    powerthreshold[1]=-1; // not used 
    powerthreshold[2]=-1; // not used 
    powerthreshold[3]=-1; // not used
    powerthreshold[4]=-5.40247; // Average Anita-3 scaler is 450kHz, which corresponds to this threshold as seen in
    // p. 9 of Ryan's talk at the Anita meeting 19th Feb 2008
		
    foutput << "Thresholds are (in p/<p>):  " <<
      powerthreshold[0] << " (L)\t" <<
      powerthreshold[1] << " (M1)\t" <<
      powerthreshold[2] << " (M2)\t" <<
      powerthreshold[3] << " (H)\t" <<
      powerthreshold[4] << " \n";
      
  }
}

void Anita::readVariableThresholds(Settings *settings1){

  if (settings1->WHICH==8) { // ANITA-2
    fturf=new TFile((ICEMC_DATA_DIR+"/turfrate_icemc.root").c_str());
    turfratechain=(TTree*)fturf->Get("turfrate_icemc");
    turfratechain->SetMakeClass(1);
    turfratechain->SetBranchAddress("phiTrigMask",&phiTrigMask);
    turfratechain->SetBranchAddress("realTime",&realTime_turfrate);
    turfratechain->BuildIndex("realTime");
    turfratechain->GetEvent(0);
    realTime_tr_min=realTime_turfrate; // realTime of first event in the file
    turfratechain->GetEvent(turfratechain->GetEntries()-1);
    realTime_tr_max=realTime_turfrate; // realTime of last event in file

    
  }else if (settings1->WHICH==9 || settings1->WHICH==10){ // ANITA-3 and 4
    
    string turfFile="";
    string surfFile="";
    if (settings1->WHICH==9){
      turfFile+=ICEMC_DATA_DIR+"/SampleTurf_icemc_anita3.root";
      surfFile+=ICEMC_DATA_DIR+"/SampleSurf_icemc_anita3.root";
    }else{
      turfFile+=ICEMC_DATA_DIR+"/SampleTurf_run42to367_anita4.root";
      surfFile+=ICEMC_DATA_DIR+"/SampleSurf_run42to367_anita4.root";
    }

    // Reading in masking every 60 seconds
    fturf=new TFile(turfFile.c_str());
    turfratechain=(TTree*)fturf->Get("turfrate_icemc");
    turfratechain->SetMakeClass(1);
    turfratechain->SetBranchAddress("phiTrigMask",&phiTrigMask);
    turfratechain->SetBranchAddress("phiTrigMaskH",&phiTrigMaskH);
    turfratechain->SetBranchAddress("l1TrigMask",&l1TrigMask);
    turfratechain->SetBranchAddress("l1TrigMaskH",&l1TrigMaskH);

    /////////// DEAD TIME ONLY DEFINED FOR ANITA-3 !!!!!!!!!!!!
    if (settings1->WHICH==9) turfratechain->SetBranchAddress("deadTime",&deadTime);
    turfratechain->SetBranchAddress("realTime",&realTime_turfrate);
    turfratechain->BuildIndex("realTime");
    turfratechain->GetEvent(0);
    realTime_tr_min=realTime_turfrate; // realTime of first event in the file
    turfratechain->GetEvent(turfratechain->GetEntries()-1);
    realTime_tr_max=realTime_turfrate; // realTime of last event in file

    // Reading in thresholds/scalers every 60 seconds
    fsurf=new TFile(surfFile.c_str());
    surfchain=(TTree*)fsurf->Get("surf_icemc");
    surfchain->SetMakeClass(1);
    surfchain->SetBranchAddress("thresholds",       &thresholds       );
    surfchain->SetBranchAddress("scalers",          &scalers          );
    ///////////// FAKE SCALERS AND THRESHOLDS ONLY FOR ANITA 3 !!!!!!!!!!
    if (settings1->WHICH==9){
      surfchain->SetBranchAddress("fakeThreshold",    &fakeThresholds   );
      surfchain->SetBranchAddress("fakeThreshold2",   &fakeThresholds2  );
      surfchain->SetBranchAddress("fakeScaler",       &fakeScalers      );
    }
    surfchain->SetBranchAddress("realTime",         &realTime_surf    );
    surfchain->BuildIndex("realTime");
    surfchain->GetEvent(0);
    realTime_surf_min=realTime_surf; // realTime of first event in the file
    surfchain->GetEvent(surfchain->GetEntries()-1);
    realTime_surf_max=realTime_surf; // realTime of last event in file
  }
}


void Anita::readAmplification(){
    
  // get rfcm amplification data
  // read from tree with amplification data
  // this tree contains a different event for each antenna and polarization
  TFile *f2=new TFile((ICEMC_DATA_DIR+"/gains.root").c_str());
  TTree *tgain=(TTree*)f2->Get("tree1");
    
  float freq_ampl_eachantenna[NPOINTS_AMPL];
  float ampl_eachantenna[NPOINTS_AMPL];
  float noisetemp_eachantenna[NPOINTS_AMPL];
    
  tgain->SetBranchAddress("freq",freq_ampl_eachantenna);
  tgain->SetBranchAddress("ampl",ampl_eachantenna);
  tgain->SetBranchAddress("noisetemp",noisetemp_eachantenna);

  for (int iant=0;iant<48;iant++) {
    tgain->GetEvent(iant);
    for (int j=0;j<NPOINTS_AMPL;j++) {
	
      freq_ampl[iant][j]=(double)freq_ampl_eachantenna[j];
	
      ampl[iant][j]=(double)ampl_eachantenna[j];
	
      ampl[iant][j]+=32.; // add 32 dB to correct for attenuation that was used during the test
      ampl_notdb[iant][j]=pow(10.,ampl[iant][j]/10.); // convert to regular fraction
	
      noisetemp[iant][j]=(double)noisetemp_eachantenna[j]; // so far we don't use this for anything
    }
  }
  f2->Close();
}




void Anita::getDiodeDataAndAttenuation(Settings *settings1, TString outputdir){

  // get vnoise data
  string sdiode;
  if (BANDING==0)
    sdiode=ICEMC_DATA_DIR+"/diode_anita1.root";
  else if (BANDING==1) 
    sdiode=ICEMC_DATA_DIR+"diode_nobanding.root";
  else if (BANDING==2)
    sdiode=ICEMC_DATA_DIR+"/diode_anita2.root";
  else if (BANDING==4 || BANDING==5) // Linda
    sdiode=ICEMC_DATA_DIR+"/diode_anita3.root";
    
  fnoise=new TFile(sdiode.c_str());
  tdiode=(TTree*)fnoise->Get("diodeoutputtree");
    
    
    
  string sbands;
  if (BANDING==0)
    sbands=ICEMC_DATA_DIR+"/bands_anita1.root";
  else if (BANDING==1)   
    sbands=ICEMC_DATA_DIR+"/bands_nobanding.root";
  else if (BANDING==2)
    sbands=ICEMC_DATA_DIR+"/bands_anita2.root";
  else if (BANDING==4 || BANDING==5) // Linda 
    sbands=ICEMC_DATA_DIR+"/bands_anita2.root";
    
  TFile *fbands=new TFile(sbands.c_str());
  TTree *tbands=(TTree*)fbands->Get("bandstree");
    
    
  for (int i=0;i<HALFNFOUR;i++) {
    time[i]=(double)i*TIMESTEP;
    time_long[i]=time[i];
    //cout << "time is " << time[i] << "\n";
    time_centered[i]=time[i]-(double)HALFNFOUR/2*TIMESTEP;
  }
  for (int i=HALFNFOUR;i<NFOUR;i++) {
    time_long[i]=(double)i*TIMESTEP;
  }
    
    
  // get diode model
  getDiodeModel();
    
  int m=(int)(maxt_diode/TIMESTEP);
  for (int j=0;j<5;j++) {
    for (int i=0;i<m;i++) {
      fdiode_real[j][i]=diode_real[j][i];
    }
    for (int i=m;i<NFOUR;i++) {
      fdiode_real[j][i]=0.;   // now fdiode_real is NFOUR array which is double sized then the signal we have. This is for zero padding for later convolution.
    }
      
    Tools::realft(fdiode_real[j],1,NFOUR);  // now fdiode_real is in freq domain
  }
  // try applying an exponential to the frequency domain
    
    
  TCanvas *cdiode=new TCanvas("cdiode","cdiode",880,800);
  cdiode->Divide(1,2);
  TGraph *gdiode=new TGraph(NFOUR/2,time,diode_real[4]);
  cdiode->cd(1);
  gdiode->Draw("al");
  gdiode=new TGraph(NFOUR/2,freq_forfft,fdiode_real[4]);
  cdiode->cd(2);
  gdiode->Draw("al");
  
  stemp=string(outputdir.Data())+"/diode.eps";
  cdiode->Print((TString)stemp);
    
  tdiode->SetBranchAddress("avgfreqdomain_lab",&(avgfreqdomain_lab[0]));
  tdiode->SetBranchAddress("freqdomain_amp_icemc",&(freqdomain_rfcm[0]));
  tdiode->SetBranchAddress("freqdomain_rfcm_banding",&(freqdomain_rfcm_banding[0][0]));
    

    
  tbands->SetBranchAddress("freq_bands",freq_bands);
  tbands->SetBranchAddress("bandsattn",bandsattn);
  //tbands->SetBranchAddress("correl",&(correl[0][0]));
  tbands->SetBranchAddress("correl_banding",&(correl_banding[0][0]));
  tbands->SetBranchAddress("correl_lab",&(correl_lab[0]));
  tbands->GetEntry(0);
    
    
    
    
    
  //  cout << "WARNING!! Altering bandsattn.\n";
  for (int j=0;j<5;j++) {
    //BoxAverage(bandsattn[j],NPOINTS_BANDS,10);
    for (int i=0;i<NPOINTS_BANDS;i++) {
			
      //    bandsattn[j][i]*=4.;
      if (bandsattn[j][i]>1.)
	bandsattn[j][i]=1.;
			
      //    if (BANDING==0 && j==4) // make this the same as the full band in anita 2
      //bandsattn[j][i]=1.;
    }
  }
    
    
  TGraph *gbandsattn[5];
  TGraph *gcorr[5];
  TH2F *hbandsattn=new TH2F("hbandsattn","hbandsattn",100,0.,2.E9,100,0.,1.);
  TH2F *hcorr=new TH2F("hcorr","hcorr",100,0.,2.E9,100,0.,2.);
  for (int i=0;i<5;i++) {
    gbandsattn[i]=new TGraph(NPOINTS_BANDS,freq_bands[i],bandsattn[i]);
    gbandsattn[i]->SetLineColor(2+i);
    gcorr[i]=new TGraph(NPOINTS_BANDS,freq_bands[i],correl_banding[i]);
    gcorr[i]->SetLineColor(2+i);
  }
  TCanvas *cbands=new TCanvas("cbands","cbands",880,800);
  cbands->Divide(1,2);
  cbands->cd(1);
  hbandsattn->Draw();
  for (int i=0;i<5;i++) {
    gbandsattn[i]->Draw("l");
  }
  cbands->cd(2);
  hcorr->Draw();
  for (int i=0;i<5;i++) {
    gcorr[i]->Draw("l");
  }
  stemp=string(outputdir.Data())+"/bands.eps";
  cbands->Print((TString)stemp);
    
    
    
  //   if (BANDING==0)
  //     sbands=ICEMC_DATA_DIR+"/bands_anita2.root";
    
    
  //   TFile *fbands_temp=new TFile(sbands.c_str());
  //   TTree *tbands_temp=(TTree*)fbands_temp->Get("bandstree");
  //   double bandsattn_temp[5][NPOINTS_BANDS];
  //   tbands_temp->SetBranchAddress("bandsattn",bandsattn_temp);
  //   tbands_temp->GetEvent(0);
  //   for (int i=0;i<NPOINTS_BANDS;i++) {
  //     bandsattn[4][i]=bandsattn_temp[4][i];
  //   }
    
}



void Anita::setDiodeRMS(Settings *settings1, TString outputdir){

  double mindiodeconvl[5];
  double onediodeconvl[5];
    
  double power_noise_eachband[5][NFOUR];
  double timedomain_output[5][NFOUR];
    
  // average result of diode integrator during quiet time
    
  for (int i=0;i<5;i++) {
    bwslice_enoise[i]=0.;
    bwslice_rmsdiode[i]=0.;
    bwslice_vrms[i]=0.;
  }

  
  nnoiseevents=tdiode->GetEntries();
  noiseeventcounter=0;
    
    
  int nhnoisebins=100;
  TH1F *hnoise[5];
    
  char histname[150];
  for (int j=0;j<3;j++) {
		
    sprintf(histname,"hnoise_%d",j);
    if (BANDING==2)
      hnoise[j]=new TH1F(histname,histname,nhnoisebins,-40.,20.);
    else
      hnoise[j]=new TH1F(histname,histname,nhnoisebins,-80.,40.);
  }
  if (BANDING==2) {
    sprintf(histname,"hnoise_3");
    hnoise[3]=new TH1F(histname,histname,nhnoisebins,-60.0,40.0);
    sprintf(histname,"hnoise_4");
    hnoise[4]=new TH1F(histname,histname,nhnoisebins,-60.0,40.0);
  }
  else {
    sprintf(histname,"hnoise_3");
    hnoise[3]=new TH1F(histname,histname,nhnoisebins,-120.0,80.0);
    sprintf(histname,"hnoise_4");
    hnoise[4]=new TH1F(histname,histname,nhnoisebins,-120.0,80.0);
		
  }
    
  // just take the average noise arrays from the tdiode tree
  tdiode->GetEvent(0);
    
  cout << "after getting event, freqdomain_rfcm_banding is " << freqdomain_rfcm_banding[0][NFOUR/8-1] << "\n";
    
  // TGraph *gfreqdomain_rfcm=new TGraph(NFOUR/4,freq_forplotting,freqdomain_rfcm);
  // TGraph *gavgfreqdomain_lab=new TGraph(NFOUR/4,freq_forplotting,avgfreqdomain_lab);
  // TGraph *gfreqdomain_rfcm_banding[5];

  // for (int i=0;i<5;i++) {
  //   gfreqdomain_rfcm_banding[i]=new TGraph(NFOUR/4,freq_forplotting,freqdomain_rfcm_banding[i]);
  // }
    
    
    
  // TCanvas *cfreq=new TCanvas("cfreq","cfreq",880,800);
  // cfreq->Divide(1,3);
  // cfreq->cd(1);
  // gfreqdomain_rfcm->Draw("al");
  // cfreq->cd(2);
  // gavgfreqdomain_lab->Draw("al");
  // cfreq->cd(3);
  // gfreqdomain_rfcm_banding[4]->Draw("al");
  // for (int i=0;i<5;i++) {
  //   gfreqdomain_rfcm_banding[i]->Draw("l");
  // }
  // stemp=string(outputdir.Data())+"/freqdomainplots.eps";
  // cfreq->Print((TString)stemp);
    
    
    
  // do the box banding for BANDING==1
  if (BANDING==1) {
    for (int j=0;j<5;j++) {
	
      for (int k=0;k<NFOUR/4;k++) {
	  
	if (bwslice_min[j]>freq_forplotting[k] || bwslice_max[j]<freq_forplotting[k]) {
	  freqdomain_rfcm_banding[j][k]=0.;
	}
      }
    }// end loop over bands
  }
    
  
  
  double power=0.;
  for (int j=0;j<5;j++) {
    for (int k=0;k<NFOUR/4;k++) {
      power+=freqdomain_rfcm_banding[j][k]/((double)NFOUR/4); // in V^2
    }
    //  cout << "power is " << power << "\n";
  }
    
  int ngeneratedevents=1000;
  int passes[5]={0,0,0,0,0};
  double averageoutput[5]={0.,0.,0.,0.,0.};

  if (!settings1->NOISEFROMFLIGHTTRIGGER){
    for (int i=0;i<ngeneratedevents;i++) {
      // put phases to freqdomain_rfcm_banding_rfcm array
      // assign phases (w/ correlations) to freqdomain_rfcm_banding_rfmc_banding array
      GetNoiseWaveforms();
      for (int j=0;j<5;j++) {
	//       // timedomainnoise_rfcm_banding_e[j] contain noise waveforms for
	//       // each band.  They come from taking the waveforms measured in
	//       // the signal stream, removing the rfcm's and lab attn.,
	//       // then reinserting rfcm's and band attn.
	//       // so they should be narrow band
	//       // myconvlv performs a convolution on that time domain waveform
	//       // based on the function that you define in getDiodeModel
	    
	//       cout << "timedomainnoise_rfcm_banding_e is " << timedomainnoise_rfcm_banding_e[j][NFOUR/4-1] << "\n";
	myconvlv(timedomainnoise_rfcm_banding[0][j],NFOUR,fdiode_real[j],mindiodeconvl[j],onediodeconvl[j],power_noise_eachband[j],timedomain_output[j]);
			
	hnoise[j]->Fill(onediodeconvl[j]*1.E15);
			
	bwslice_enoise[j]+=onediodeconvl[j];
			
			
	for (int m=(int)(maxt_diode/TIMESTEP);m<NFOUR/2;m++) {
	  bwslice_meandiode[j]+=timedomain_output[j][m]/((double)ngeneratedevents*((double)NFOUR/2-maxt_diode/TIMESTEP));
	  bwslice_vrms[j]+=timedomainnoise_rfcm_banding[0][j][m]*timedomainnoise_rfcm_banding[0][j][m]/((double)ngeneratedevents*((double)NFOUR/2-maxt_diode/TIMESTEP)); // this is the rms of the diode input voltage
	  averageoutput[j]+=timedomain_output[j][m]*timedomain_output[j][m]/((double)ngeneratedevents*((double)NFOUR/2-maxt_diode/TIMESTEP));
	} // end loop over samples where diode function is fully contained
			
      } // end loop over bands
		
    } // end loop over generated events
    
    
    // TCanvas *ctest=new TCanvas("ctest","ctest",880,800);
    // ctest->Divide(1,5);
    // TGraph *gtest[5];
    // for (int i=0;i<5;i++) {
    //   ctest->cd(i+1);
    //   gtest[i]=new TGraph(NFOUR,time_long,timedomain_output[i]);
    //   gtest[i]->Draw("al");
    // }
    // stemp = string(outputdir.Data())+"/test.eps";
    // ctest->Print((TString)stemp);
    
    for (int j=0;j<5;j++) {
      
      bwslice_vrms[j]=sqrt(bwslice_vrms[j]); // this is the rms input voltage
    }
    
    
    for (int j=0;j<5;j++) {
      bwslice_enoise[j]=hnoise[j]->GetMean()/1.E15; // mean diode output with no correlations
      
      bwslice_fwhmnoise[j]=Tools::GetFWHM(hnoise[j]);
      
    }
    
    
    for (int i=0;i<ngeneratedevents;i++) {// now we need to get the rms
      GetNoiseWaveforms();
      for (int j=0;j<5;j++) {
	
	myconvlv(timedomainnoise_rfcm_banding[0][j],NFOUR,fdiode_real[j],mindiodeconvl[j],onediodeconvl[j],power_noise_eachband[j],timedomain_output[j]);
	
	for (int m=(int)(maxt_diode/TIMESTEP);m<NFOUR/2;m++) {	  
	  bwslice_rmsdiode[j]+=(timedomain_output[j][m]-bwslice_meandiode[j])*(timedomain_output[j][m]-bwslice_meandiode[j])/((double)ngeneratedevents*((double)NFOUR/2-maxt_diode/TIMESTEP));
	}
      
      }
    }
    
    for (int j=0;j<5;j++) {
      bwslice_rmsdiode[j]=sqrt(bwslice_rmsdiode[j]);
      cout << "mean, rms are " << bwslice_meandiode[j] << " " << bwslice_rmsdiode[j] << "\n";
    }
  
    double thresh_begin=-1.;
    double thresh_end=-11.;
    double thresh_step=1.;
  
    double rate[5][100];
    for (double testthresh=thresh_begin;testthresh>=thresh_end;testthresh-=thresh_step) {
      for (int j=0;j<5;j++) {
	passes[j]=0;
      }
      for (int i=0;i<ngeneratedevents;i++) {
  	
	GetNoiseWaveforms();
	for (int j=0;j<5;j++) {

	  myconvlv(timedomainnoise_rfcm_banding[0][j],NFOUR,fdiode_real[j],mindiodeconvl[j],onediodeconvl[j],power_noise_eachband[j],timedomain_output[j]);
  	
	  for (int m=(int)(maxt_diode/TIMESTEP);m<NFOUR/2;m++) {
  	    
	    if (timedomain_output[j][m+1]<bwslice_rmsdiode[j]*testthresh) {
	      passes[j]++;
	      m+=(int)(DEADTIME/TIMESTEP);
	    }
  	    	    
	  } // end loop over samples where diode function is fully contained
	}
      }
      // brian m.
      for (int j=0;j<5;j++) {
	int ibin=(int)fabs((testthresh-thresh_begin)/thresh_step);
	// relthresh[j][ibin]=fabs(testthresh);
	rate[j][ibin]=passes[j]/((double)ngeneratedevents*((double)NFOUR/2*(TIMESTEP)-maxt_diode));
	if (rate[j][ibin]!=0)
	  rate[j][ibin]=log10(rate[j][ibin]);
	//cout << "Threshold " << testthresh << " For the " << j << "th band, rate is " << rate[j][ibin] << "\n";
      }
    }
    
  } else { // IF WE HAVE NOISE FROM FLIGHT


#ifdef ANITA_UTIL_EXISTS
    double quickNoise[HALFNFOUR];

    memset(bwslice_diodemean_fullband_allchan,  0, sizeof(bwslice_diodemean_fullband_allchan)  );
    memset(bwslice_dioderms_fullband_allchan,   0, sizeof(bwslice_dioderms_fullband_allchan)   );

    static double tempdiodeoutput[1000][NFOUR];

    for (int ipol=0; ipol<2; ipol++){
      for (int iant=0; iant<48; iant++){
	
	memset(tempdiodeoutput, 0, sizeof(tempdiodeoutput) );

	for (int i=0;i<ngeneratedevents;i++) {
	  
	  getQuickTrigNoiseFromFlight(quickNoise, ipol, iant);
	  
	  myconvlv(quickNoise,NFOUR,fdiode_real[4],mindiodeconvl[4],onediodeconvl[4],power_noise_eachband[4],tempdiodeoutput[i]);

	  // First calculate the mean
	  for (int m=(int)(maxt_diode/TIMESTEP);m<NFOUR/2;m++) {
	    bwslice_diodemean_fullband_allchan[ipol][iant]+=tempdiodeoutput[i][m]/((double)ngeneratedevents*((double)NFOUR/2-maxt_diode/TIMESTEP));
	    //	  cout << m << " " << timedomain_output[j][m] << " " << ((double)ngeneratedevents*((double)NFOUR/2-maxt_diode/TIMESTEP)) << endl;
	  
	  }
	}

	// Then get the RMS
	for (int i=0;i<ngeneratedevents;i++) {
	  
	  for (int m=(int)(maxt_diode/TIMESTEP);m<NFOUR/2;m++) {
	    bwslice_dioderms_fullband_allchan[ipol][iant]+=(tempdiodeoutput[i][m]-bwslice_diodemean_fullband_allchan[ipol][iant])*(tempdiodeoutput[i][m]-bwslice_diodemean_fullband_allchan[ipol][iant])/((double)ngeneratedevents*((double)NFOUR/2-maxt_diode/TIMESTEP));
	  }

	}
	
	bwslice_dioderms_fullband_allchan[ipol][iant]=sqrt(bwslice_dioderms_fullband_allchan[ipol][iant]);
	//cout << "EACH CHAN MEAN, RMS " <<  ipol << " " << iant << " " << bwslice_diodemean_fullband_allchan[ipol][iant] << " , " << bwslice_dioderms_fullband_allchan[ipol][iant] << endl;  
	
      }
	
    }
    
    
    // double thresh_begin=-1.;
    // double thresh_end=-11.;
    // double thresh_step=1.;
  
    // double rate[5][100];
    // for (double testthresh=thresh_begin;testthresh>=thresh_end;testthresh-=thresh_step) {
    //   for (int j=0;j<5;j++) {
    // 	passes[j]=0;
    //   }
    //   for (int i=0;i<ngeneratedevents;i++) {
  	
    // 	for (int j=4;j<5;j++) {

    // 	  getQuickTrigNoiseFromFlight(quickNoise, 0, 0);

    // 	  myconvlv(quickNoise,NFOUR,fdiode_real[j],mindiodeconvl[j],onediodeconvl[j],power_noise_eachband[j],timedomain_output[j]);
  	
    // 	  for (int m=(int)(maxt_diode/TIMESTEP);m<NFOUR/2;m++) {
  	    
    // 	    if (timedomain_output[j][m+1]<bwslice_rmsdiode[j]*testthresh) {
    // 	      passes[j]++;
    // 	      m+=(int)(DEADTIME/TIMESTEP);
    // 	    }
  	    	    
    // 	  } // end loop over samples where diode function is fully contained
    // 	}
    //   }
    //   // brian m.
    //   for (int j=0;j<5;j++) {
    // 	int ibin=(int)fabs((testthresh-thresh_begin)/thresh_step);
    // 	// relthresh[j][ibin]=fabs(testthresh);
    // 	rate[j][ibin]=passes[j]/((double)ngeneratedevents*((double)NFOUR/2*(TIMESTEP)-maxt_diode));
    // 	if (rate[j][ibin]!=0)
    // 	  rate[j][ibin]=log10(rate[j][ibin]);
    // 	//cout << "Threshold " << testthresh << " For the " << j << "th band, rate is " << rate[j][ibin] << "\n";
    //   }
    // }
   
#endif 
  }
    
  TF1 *frice[5];
    
  int jplot;
  TCanvas *c4=new TCanvas("c4","c4",880,800);
  c4->Divide(1,5);
  for (int j=0;j<5;j++) {
	
	
    //    cout << "cd to " << j+1 << "\n";
    frice[j]=new TF1("frice","[2]*(-1.*x-[1])/[0]^2*exp(-1.*(-1.*x-[1])^2/2/[0]^2)",-50.,50.);
    frice[j]->SetParameter(0,4.);
    //    frice[j]->SetParameter(1,-5.);
    frice[j]->SetParameter(1,5.);
    frice[j]->SetParameter(2,400.);
    //frice[j]->SetParameter(2,frice[j]->GetParameter(2)*hnoise[j]->GetEntries()/frice[j]->Integral(-500.,500.));
    jplot=j;
    //    if (BANDING==2 && j==4)
    //jplot=j-1;
	
    if ((BANDING==2 && j!=3) ||
	//	(BANDING!=2 && j!=4)) {
	(BANDING!=2 && BANDING!=4) ||
	((BANDING==4||BANDING==5) && j==4) ){
      c4->cd(jplot+1);
	    
      if (j==0 && BANDING==2)
	hnoise[j]->Fit(frice[j],"Q","",-20.,1.);
      if (j==1 && BANDING==2)
	hnoise[j]->Fit(frice[j],"Q","",-20.,2.);
      if (j==2 && BANDING==2)
	hnoise[j]->Fit(frice[j],"Q","",-20.,10.);
      if (j==4 && (BANDING==2 || BANDING==4 || BANDING==5))
	hnoise[j]->Fit(frice[j],"Q","",-15.,15.);
	    
      if (j==0 && BANDING==0)
	hnoise[j]->Fit(frice[j],"Q","",-20.,5.);
      if (j==1 && BANDING==0)
	hnoise[j]->Fit(frice[j],"Q","",-20.,5.);
      if (j==2 && BANDING==0)
	hnoise[j]->Fit(frice[j],"Q","",-20.,15.);
      if (j==3 && BANDING==0)
	hnoise[j]->Fit(frice[j],"Q","",-20.,20.);
	    
      frice[j]->GetHistogram()->SetLineWidth(3);
	    
      hnoise[j]->SetLineWidth(3);
      hnoise[j]->SetXTitle("Diode Output (10^{-15} Joules)");
      hnoise[j]->SetYTitle("Number of Noise Waveforms");
      hnoise[j]->SetTitleSize(0.05,"X");
      hnoise[j]->SetTitleSize(0.05,"Y");
	    
      hnoise[j]->Draw();
      frice[j]->Draw("same");
      //    cout << "bwslice_enoise is " << bwslice_enoise[j] << "\n";
    }
	
  }
  stemp=string(outputdir.Data())+"/hnoise.eps";
  c4->Print((TString)stemp);

  
  
}



#ifdef ANITA_UTIL_EXISTS
void Anita::getQuickTrigNoiseFromFlight(double justNoise[HALFNFOUR], int ipol, int iant){

  FFTWComplex *phasorsTrig = new FFTWComplex[numFreqs];
  phasorsTrig[0].setMagPhase(0,0);
  Double_t sigma, realPart, imPart, norm;
  int iring=2;
  if (iant<16) iring=0;
  else if (iant<32) iring=1;

  int iphi = iant - (iring*16);

  for(int i=1;i<numFreqs;i++) {
    norm           = fRatioTriggerDigitizerFreqDomain[ipol][iring][iphi][i];
    sigma          = RayleighFits[ipol][iant]->Eval(freqs[i])*4./TMath::Sqrt(numFreqs);
    sigma*=norm;
    realPart       = fRand->Gaus(0,sigma);
    imPart         = fRand->Gaus(0,sigma);

    //    cout << " " << i << " " << trig << " " << dig << " " << trig/dig << " " << norm << endl;
    phasorsTrig[i] = FFTWComplex(realPart, imPart);
  }

  RFSignal *rfNoiseTrig   = new RFSignal(numFreqs,freqs,phasorsTrig,1);  
  Double_t *justNoiseTemp = rfNoiseTrig->GetY();
  
  for (int i=0; i<HALFNFOUR; i++){
    justNoise[i]  = justNoiseTemp[i]*THERMALNOISE_FACTOR;
  }
  
  delete rfNoiseTrig;
  delete[] phasorsTrig;
  
}
#endif

void Anita::getPulserData(){
      
  TFile *fpulser=new TFile((ICEMC_DATA_DIR+"/pulser.root").c_str());
	
  TGraph *gpulser=(TGraph*)fpulser->Get("pulser");
  TGraph *gphases=(TGraph*)fpulser->Get("phases");
  TGraph *gnoise=(TGraph*)fpulser->Get("noise");
	
	
  double *temp1=gpulser->GetX();
  for (int i=0;i<NFOUR/4;i++) {
    f_pulser[i]=temp1[i];
  }
  double *temp2=gphases->GetX();
  for (int i=0;i<NFOUR/4;i++) {
    f_phases[i]=temp2[i];
  }
  double *temp3=gpulser->GetY();
  double *temp4=gnoise->GetY();
	
	
  for (int i=0;i<NFOUR/4;i++) {
    v_pulser[i]=sqrt(temp3[i]); // is this right
    v_noise[i]=sqrt(temp4[i]);
    if (f_phases[i]<75.E6)
      v_noise[i]=0.;
	    
  }
	
  TGraph *gpulser_eachband;
  TGraph *gnoise_eachband;
  TCanvas *cpulser=new TCanvas("cpulser","cpulser",880,800);
  cpulser->Divide(1,5);
  int iplot;
  for (int i=0;i<5;i++) {
    iplot=i;
			
    if (i!=3) {
      cpulser->cd(iplot+1);
				
      gpulser_eachband=new TGraph(NFOUR/4,f_pulser,v_pulser);
				
      gpulser_eachband->Draw("al");
    }
  }
		
  cpulser->Print("pulser.eps");
		
  TCanvas *cnoise=new TCanvas("cnoise","cnoise",880,800);
  cnoise->Divide(1,5);
  for (int i=0;i<5;i++) {
    iplot=i;
			
    if (i!=3) {
      cnoise->cd(iplot+1);
      gnoise_eachband=new TGraph(NFOUR/4,f_pulser,v_noise);
				
      gnoise_eachband->Draw("al");
    }
  }
		
  cnoise->Print("noise.eps");
	       
		
  //	gpulser_eachband=new TGraph(NFOUR/4,f_pulser,v_pulser[iplot]);
  //gpulser_eachband->Draw("al");
		
		
		
  double sumpulserpower=0.;
  double sumnoisepower=0.;
		
  for (int i=0;i<NFOUR/4;i++) {
    sumpulserpower+=v_pulser[i]*v_pulser[i];
    sumnoisepower+=v_noise[i]*v_noise[i];
  }
		
  double *temp5=gphases->GetY();
  for (int i=0;i<NFOUR/4;i++) {
    v_phases[i]=temp5[i];
  }
		
  gpulser->Delete();
  gphases->Delete();
  gnoise->Delete();
		
  fpulser->Close();
}








void Anita::ReadGains(void) {
  // gains from university of hawaii measurements.
  double sfrequency;
  int iii;
  ifstream gainsfile;
  gainsfile.open((ICEMC_DATA_DIR+"/hh_0").c_str()); // gains for horizontal polarization
  if(gainsfile.fail()) {
    cout << "can't open `$ICEMC_DATA_DIR`/hh_0\n";
    exit(1);
  }
  for(iii = 0; iii < NPOINTS_GAIN; iii++)
    gainsfile >> frequency_forgain_measured[iii] >> gainh_measured[iii];
  gainsfile.close();
    
  gainsfile.open((ICEMC_DATA_DIR+"/vv_0").c_str()); // gains for vertical polarization
  if(gainsfile.fail()) {
    cout << "can't open `$ICEMC_DATA_DIR`/vv_0\n";
    exit(1);
  }
  for(iii = 0; iii < NPOINTS_GAIN; iii++) {
    gainsfile >> sfrequency >> gainv_measured[iii];
    if(sfrequency != frequency_forgain_measured[iii])
      cout << "warning: sfrequency = " << sfrequency << ", frequency_forgain_measured[iii] = " << frequency_forgain_measured[iii] << endl;
  }
  gainsfile.close();
    
  gainsfile.open((ICEMC_DATA_DIR+"/hv_0").c_str()); // gains for h-->v cross polarization
  if(gainsfile.fail()) {
    cout << "can't open `$ICEMC_DATA_DIR`/hv_0\n";
    exit(1);
  }
  for(iii = 0; iii < NPOINTS_GAIN; iii++) {
    gainsfile >> sfrequency >> gainhv_measured[iii];
    if(sfrequency != frequency_forgain_measured[iii])
      cout << "warning: sfrequency = " << sfrequency << ", frequency_forgain_measured[iii] = " << frequency_forgain_measured[iii] << endl;
  }
  gainsfile.close();
    
  gainsfile.open((ICEMC_DATA_DIR+"/vh_0").c_str()); // gains for v-->h cross polarization
  if(gainsfile.fail()) {
    cout << "can't open `$ICEMC_DATA_DIR`/vh_0\n";
    exit(1);
  }
  for(iii = 0; iii < NPOINTS_GAIN; iii++) {
    gainsfile >> sfrequency >> gainvh_measured[iii];
    if(sfrequency != frequency_forgain_measured[iii])
      cout << "warning: sfrequency = " << sfrequency << ", frequency_forgain_measured[iii] = " << frequency_forgain_measured[iii] << endl;
  }
  gainsfile.close();
} //ReadGains




void Anita::AntennaGain(Settings *settings1,double hitangle_e,double hitangle_h,double e_component,double h_component,int k,double &vsignalarray_e,double &vsignalarray_h) {
    
  if (freq[k]>=settings1->FREQ_LOW_SEAVEYS && freq[k]<=settings1->FREQ_HIGH_SEAVEYS) {
		
    double relativegains[4]; // fill this for each frequency bin for each antenna.  It's the gain of the antenna given the angle that the signal hits the balloon, for vv, vh, hv, hh, relative to the gain at boresight
		
    for (int pols=0;pols<2;pols++) {// loop over vv, hv
      if (fabs(hitangle_e)<PI/2)
	relativegains[pols]=Get_gain_angle(pols,k,hitangle_e);//change here to be constant if need be (ABBY).  Add a setting to the code for use constant gain or dish or something.  Make this a member function.
      else
	relativegains[pols]=0.;
      //if (fabs(hitangle_e)<PI/12)
      //cout << "relative gains is " << relativegains[pols] << "\n";
    }
		
		
    //      if (fabs(hitangle_e)<PI/12)
    //cout << "vsignalarray_e before is " << vsignalarray_e << "\n";
    vsignalarray_e = vsignalarray_e * 0.5 * sqrt(vvGaintoHeight[k] * vvGaintoHeight[k] * e_component * e_component * relativegains[0] + hvGaintoHeight[k] * hvGaintoHeight[k] * h_component * h_component * relativegains[1]); // 0.5 is for voltage dividing
		
    //cout << "In AntennaGain, vsignalarray_e is " << vsignalarray_e << "\n";
		
    //if (fabs(hitangle_e)<PI/12)
    //cout << "vsignalarray_e after is " << vsignalarray_e << "\n";
		
    for (int pols=2;pols<4;pols++) { // loop over hh, vh
      if (fabs(hitangle_h)<PI/2)
	relativegains[pols]=Get_gain_angle(pols,k,hitangle_h);
      else
	relativegains[pols]=0.;
    }
		
    // V/MHz
		
    //if (fabs(hitangle_h)<PI/12)
    //cout << "vsignalarray_h before is " << vsignalarray_h << "\n";
		
    vsignalarray_h=vsignalarray_h*0.5*
      sqrt(hhGaintoHeight[k]*hhGaintoHeight[k]*h_component*h_component*relativegains[2] + vhGaintoHeight[k]*vhGaintoHeight[k]*e_component*e_component*relativegains[3]);
		
    //      if (fabs(hitangle_h)<PI/12)
    //cout << "vsignalarray_h after is " << vsignalarray_h << "\n";
		
  }
}


// The deck affects signals reaching the upper ring. These equations are for Fresnel diffraction around a half plane. The edge of the half plane passes just over and in front of the lower ring antenna that's facing the radio pulse. This assumption should be okay for the three upper ring antenns in the phi sector facing the pulse and should underestimate the effect of diffraction for the other upper ring antennas.
// This only corrects for signal strength, not for changes in polarization and signal direction.
// page and figure references are from Principles of Optics, by Born & Wolf, 7th edition
void Anita::SetDiffraction() {
  const double cc = 299792458.0;
  const double height[2] = {3.21571, 2.25625}; // height of the phase centers of antennas 1 and 9 above the deck
  const double radius[2] = {1.40185, 1.5677}; // radius of the deck minus how far the phase centers of antennas 1 and 9 are from the payload's axis.
  const int ntau = 200;
  double ww, xx, ss, theta, sintheta, squigglyC, squigglyS, tau, dtau;
  for(int jj = 0; jj < 2; jj++) // loop through antenna layers
    for(int kk = 0; kk < 89; kk++) { // loop through pulse directions
      sintheta = 0.37 + 0.005*double(kk);
      theta = asin(sintheta); // angle between the direction from the center of the Earth to the payload and the direction of the pulse. This should be delta in figure 8.38
      ss = height[jj]/cos(theta); // page 382 fig 8.5 . This should be the distance between the phase center and wherever the pulse passes the deck.
      xx = height[jj]*tan(theta)-radius[jj]; // page 429 fig 8.38
      for(int ll = 0; ll < NFREQ; ll++) { // loop through frequencies
	ww = sqrt(2./cc*freq[ll]/ss)*xx*sintheta; // page 433 eq 23
	dtau = ww/double(ntau);
	squigglyS = 0.0;
	squigglyC = 0.0;
	for(int ii = 0; ii < ntau; ii++) { // page 430 eq 12
	  tau = dtau*(double(ii)+0.5);
	  squigglyC += cos(M_PI*0.5*tau*tau);
	  squigglyS += sin(M_PI*0.5*tau*tau);
	}
	squigglyC *= dtau;
	squigglyS *= dtau;
	diffraction[jj][kk][ll] = sqrt(0.5*((0.5+squigglyC)*(0.5+squigglyC)+(0.5+squigglyS)*(0.5+squigglyS))); // page 434 eq 28
	if(jj == 0 && kk > 66) diffraction[jj][kk][ll] = 1.0; // in the top layer, diffraction[][][] has values greater than 1 in the low frequencies for sintheta > 0.73
      }
    }
  return;
} // void SetDiffraction()

double Anita::GetDiffraction(int ilayer, double zenith_angle, int ifreq) {
  double reference1 = (sin(zenith_angle)-0.37)*200.0;
  int whichbin = int(reference1);
  double factor = reference1 - double(whichbin);
  return factor*diffraction[ilayer][whichbin+1][ifreq] + (1.-factor)*diffraction[ilayer][whichbin][ifreq];
}

int Anita::SurfChanneltoBand(int ichan) {
    
  // takes surf channel (0-31) and finds what band it is
  if (ichan<0 || ichan>31) {
    cout << "surf channel out of range!\n";
    exit(1);
  }
    
  int iband= (ichan%8-ichan%2)/2;
    
    
  return iband;
    
}

//  int Anita::GetScalarNumber(int lr,int ibw, int ilayer,int ifold) {

//   int rx = GetAntennaNumber(ilayer,ifold);
//   int rx_on_surf = (rx-1)%4 + 1; // surf holds four antennas.  This number is 1-4.
//   int scalar= (rx_on_surf-1)*8 + ibw*2 + lr;

//   // returns the scalar number
//   return scalar;
// }
int Anita::GetAntennaNumber(int ilayer,int ifold) {
    
  return 8*(ilayer)+ifold+1; // antenna number 1-32
    
}
int Anita::GetLayer(int rx) {
  if (rx>=0 && rx<8)
    return 0;
  else if (rx>=8 && rx<16)
    return 1;
  else if (rx>=16 && rx<32)
    return 2;
    
  return -1;
}
int Anita::GetIfold(int rx) {
  if (rx>=0 && rx<16)
    return rx%8;
  else if (rx<32)
    return rx%16;
    
  return -1;
}
int Anita::AntennaWaveformtoSurf(int ilayer,int ifold) {
    
  int antenna=GetAntennaNumber(ilayer,ifold); // antenna number 1-32
    
  return antennatosurf[antenna-1]; // returns the surf number
    
}
int Anita::AntennaNumbertoSurfNumber(int ilayer,int ifold) {
  int antenna=GetAntennaNumber(ilayer,ifold); // antenna number 1-32
    
  return (antenna-1-(antenna-1)%4)/4+1; // returns the surf number 1-9
    
    
}
int Anita::GetSurfChannel(int antenna,int ibw,int ipol) { // 1 to 32
    
  return ((antenna-1)%4)*8+WhichBand(ibw,ipol);//which scalar channel, numbered 1 to 32
    
}
int Anita::WhichBand(int ibw,int ipol) {
    
  return 2*ibw+ipol+1; // from 1 to 8, which scalar channel on an antenna this band corresponds to, in order as they are on the surfs
}
void Anita::Banding(int j,double *freq_noise,double *vmmhz,int NPOINTS_NOISE) {
    
  if (j<0 || j>4) {
    cout << "Band out of range.\n";
    exit(1);
  }
  if (BANDING!=1) {
    int iband; // which index of freq_bands
    for (int i=0;i<NPOINTS_NOISE;i++) {
			
      iband=Tools::findIndex(freq_bands[j],freq_noise[i],NPOINTS_BANDS,freq_bands[j][0],freq_bands[j][NPOINTS_BANDS-1]);
      //      cout << "freq_bands, freq_noise, NPOINTS_BANDS, freq_bands are " << freq_noise[i] << " " << NPOINTS_BANDS << " " << freq_bands[j][0] << " " << iband << "\n";
      vmmhz[i]=vmmhz[i]*sqrt(bandsattn[j][iband]); // sqrt since it's voltage and not power
			
			
    }
  }
  else if (BANDING==1) {
		
    for (int i=0;i<NPOINTS_NOISE;i++) {
      if (freq_noise[i]<bwslice_min[j] || freq_noise[i]>bwslice_max[j])
	vmmhz[i]=0.;
    } // loop over frequency bins
		
  } // if it's choose-your-own banding
}
void Anita::RFCMs(int ilayer,int ifold,double *vmmhz) {
    
  int irx=Anita::GetAntennaNumber(ilayer,ifold)-1; // want this to be 0-31
    
  int iampl;
  for (int i=0;i<NFREQ;i++) {
    // first find the index of the amplification array
    iampl=Tools::findIndex(freq_ampl[irx],freq[i],NPOINTS_AMPL,freq_ampl[irx][0],freq_ampl[irx][NPOINTS_AMPL-1]);
		
    vmmhz[i]=vmmhz[i]*sqrt(ampl_notdb[irx][iampl]);
		
  }
    
}


// reads in the effect of a signal not hitting the antenna straight on
// also reads in gainxx_measured and sets xxgaintoheight
void Anita::Set_gain_angle(Settings *settings1,double nmedium_receiver) {
  string gain_null1, gain_null2;
  double sfrequency;
  int iii, jjj;
  ifstream anglefile;
  for(jjj = 0; jjj < 4; jjj++)
    for(iii = 0; iii < 131; iii++)
      gain_angle[jjj][iii][0] = 1.;
    
  anglefile.open((ICEMC_DATA_DIR+"/vv_az").c_str()); // v polarization, a angle
  if(anglefile.fail()) {
    cout << "can't open `$ICEMC_DATA_DIR`/vv_az\n";
    exit(1);
  }
  for(jjj = 1; jjj < 7; jjj++)
    for(iii = 0; iii < 131; iii++) {
      anglefile >> sfrequency >> gain_angle[0][iii][jjj];
      if(sfrequency != frequency_forgain_measured[iii]) cout << "check frequencies for vv_az\n";
    }
  anglefile.close();
    
    
    
  anglefile.open((ICEMC_DATA_DIR+"/hh_az").c_str()); // h polarization, a angle
  if(anglefile.fail()) {
    cout << "can't open `$ICEMC_DATA_DIR`/hh_az\n";
    exit(1);
  }
  for(jjj = 1; jjj < 7; jjj++)
    for(iii = 0; iii < 131; iii++) {
      anglefile >> sfrequency >> gain_angle[1][iii][jjj];
      if(sfrequency != frequency_forgain_measured[iii]) cout << "check frequencies for hh_az\n";
    }
  anglefile.close();
    
  anglefile.open((ICEMC_DATA_DIR+"/hh_el").c_str()); // h polarization, e angle
  if(anglefile.fail()) {
    cout << "can't open `$ICEMC_DATA_DIR`/hh_el\n";
    exit(1);
  }
  for(jjj = 1; jjj < 7; jjj++)
    for(iii = 0; iii < 131; iii++) {
      anglefile >> sfrequency >> gain_angle[2][iii][jjj];
      if(sfrequency != frequency_forgain_measured[iii]) cout << "check frequencies for hh_el\n";
    }
  anglefile.close();
    
  anglefile.open((ICEMC_DATA_DIR+"/vv_el").c_str()); // v polarization, e angle
  if(anglefile.fail()) {
    cout << "can't open `$ICEMC_DATA_DIR`/vv_el\n";
    exit(1);
  }
  for(jjj = 1; jjj < 7; jjj++)
    for(iii = 0; iii < 131; iii++) {
      anglefile >> sfrequency >> gain_angle[3][iii][jjj];
      if(sfrequency != frequency_forgain_measured[iii]) cout << "check frequencies for vv_el\n";
    }
  anglefile.close();
  for(jjj = 0; jjj < 6; jjj++) inv_angle_bin_size[jjj] = 1. /
				 (reference_angle[jjj+1] - reference_angle[jjj]); // this is used for interpolating gains at angles between reference angles
    
  double gainhv, gainhh, gainvh, gainvv;
  double gain_step = frequency_forgain_measured[1]-frequency_forgain_measured[0]; // how wide are the steps in frequency;
    
  cout << "GAINS is " << GAINS << "\n";
  for (int k = 0; k < NFREQ; ++k) {
    whichbin[k] = int((freq[k] - frequency_forgain_measured[0]) / gain_step); // finds the gains that were measured for the frequencies closest to the frequency being considered here
    if((whichbin[k] >= NPOINTS_GAIN || whichbin[k] < 0) && !settings1->FORSECKEL) {
      cout << "Set_gain_angle out of range, freq = " << freq[k] << endl;
      exit(1);
    }
		
    //now a linear interpolation for the frequency
    scalef2[k] = (freq[k] - frequency_forgain_measured[whichbin[k]]) / gain_step;
    // how far from the lower frequency
    scalef1[k] = 1. - scalef2[k]; // how far from the higher frequency
		
    // convert the gain at 0 degrees to the effective antenna height at 0 degrees for every frequency
    if(whichbin[k] == NPOINTS_GAIN - 1) { // if the frequency is 1.5e9 or goes a little over due to rounding
      gainhv = gainhv_measured[whichbin[k]];
      gainhh = gainh_measured[whichbin[k]];
      gainvh = gainvh_measured[whichbin[k]];
      gainvv = gainv_measured[whichbin[k]];
    } else {
      // These gains should be dimensionless numbers, not in dBi
      gainhv = scalef1[k] * gainhv_measured[whichbin[k]] + scalef2[k] * gainhv_measured[whichbin[k] + 1];
      gainhh = scalef1[k] * gainh_measured[whichbin[k]] + scalef2[k] * gainh_measured[whichbin[k] + 1];
      gainvh = scalef1[k] * gainvh_measured[whichbin[k]] + scalef2[k] * gainvh_measured[whichbin[k] + 1];
      gainvv = scalef1[k] * gainv_measured[whichbin[k]] + scalef2[k] * gainv_measured[whichbin[k] + 1];
    }
    if (GAINS==0) {
			
      gainhv=0.;
      gainhh=gain[1][k];
      gainvh=0.;
      gainvv=gain[0][k];
			
    }
		
		
    hvGaintoHeight[k] = GaintoHeight(gainhv,freq[k],nmedium_receiver);
    hhGaintoHeight[k] = GaintoHeight(gainhh,freq[k],nmedium_receiver);
    vhGaintoHeight[k] = GaintoHeight(gainvh,freq[k],nmedium_receiver);
    vvGaintoHeight[k] = GaintoHeight(gainvv,freq[k],nmedium_receiver);
		
  } // loop through frequency bins
    
    
    
    
    
    
    
    
    
} // void Set_gain_angle()

// determines the effect of a signal not hitting the antenna straight on
// for the gain type, 0 means V polarization and el angle, 1 means H polarization and el angle, 2 means H polarization and az angle, 3 means V polarization and az angle
// gain_angle[gain_type][][] lists the decrease in gain for different frequencies and angles. This subroutine finds the two angles closest to hitangle and the two frequencies closest to freq. It then returns a linear interpolation in 2 dimensions.
int Anita::GetBeamWidths(Settings *settings1) {
    
  // first component is frequency
  // second component is which plane and which polarization
  // it goes  e-plane: vp/hp, h-plane: vp/hp
    
  // these number were read from antenna specs
    
  // these are beam widths
    
  double freq_specs[5];
  int NFREQ_FORGAINS; // Number of
    
  if (settings1->WHICH==7) { // EeVa
    NFREQ_FORGAINS=4;
    freq_specs[0]=265.E6; // lower edge of
    freq_specs[1]=435.E6;
    freq_specs[2]=650.E6;
    freq_specs[3]=992.5E6;
		
		
  }
  else if (settings1->WHICH==11) { // Satellite
    NFREQ_FORGAINS=4;
    freq_specs[0]=265.E6; // lower edge of
    freq_specs[1]=435.E6;
    freq_specs[2]=650.E6;
    freq_specs[3]=992.5E6;
		
  }
    
  else { // everything else
    NFREQ_FORGAINS=5;
    freq_specs[0]=300.E6;
    freq_specs[1]=600.E6;
    freq_specs[2]=900.E6;
    freq_specs[3]=1200.E6;
    freq_specs[4]=1500.E6;
		
  }
    
    
    
  double specs[5][4];
    
  if (settings1->WHICH==7) { // EeVA
		
		
    // these are all for 200 MHz
    specs[0][0]=2.5; // eplane, vp
    specs[0][1]=2.5;// eplane, hp
    specs[0][2]=1.5;// hplane, vp
    specs[0][3]=1.5;// hplane, hp
		
    // these are for 355 MHz
    specs[1][0]=2.2;
    specs[1][1]=2.2;
    specs[1][2]=1.2;
    specs[1][3]=1.2;
		
    // 515 MHz
    specs[2][0]=2.2;
    specs[2][1]=2.2;
    specs[2][2]=1.0;
    specs[2][3]=1.0;
		
    // 785 MHz
    specs[3][0]=2.2;
    specs[3][1]=2.2;
    specs[3][2]=0.9;
    specs[3][3]=0.9;
		
  }
  else { // everything but EeVA
    // for now use this for satellite even
    // these are beam widths in degrees
    // these are all for 300 MHz
    specs[0][0]=57.5; // eplane, vp
    specs[0][1]=58.5;// eplane, hp
    specs[0][2]=66;// hplane, vp
    specs[0][3]=57;// hplane, hp
		
    // these are for 600 MHz
    specs[1][0]=33.5;
    specs[1][1]=34.5;
    specs[1][2]=36.5;
    specs[1][3]=38;
		
    // 900 MHz
    specs[2][0]=50.5;
    specs[2][1]=53;
    specs[2][2]=33;
    specs[2][3]=32;
		
    // 1200 MHz
    specs[3][0]=43.5;
    specs[3][1]=43;
    specs[3][2]=39;
    specs[3][3]=41.5;
		
    //1500 MHz
    specs[4][0]=36.5;
    specs[4][1]=46.5;
    specs[4][2]=32;
    specs[4][3]=31;
  }
    
  // These are gains, which are not used because we use Ped's numbers instead
  double specs2[5][2];
    
  if (settings1->WHICH==7) {// EeVA
    // this is in dBi
    // 200 MHz
    specs2[0][0]=26.85; // vp
    specs2[0][1]=26.85; // hp
		
    // 355 MHz
    specs2[1][0]=35.8;
    specs2[1][1]=35.8;
		
    // 515 MHz
    specs2[2][0]=25.;
    specs2[2][1]=25.;
		
    // 785 MHz
    specs2[3][0]=10.5;
    specs2[3][1]=10.5;
		
		
  }
  else if (settings1->WHICH==11) { // satellite
		
    // this is in dBi
    // 40 MHz
    specs2[0][0]=7.5; // vp
    specs2[0][1]=7.5; // hp
    //specs2[0][0]=0.0; // vp
    //specs2[0][1]=0.0; // hp
		
    //  MHz
    specs2[1][0]=6.5;
    specs2[1][1]=6.5;
    //specs2[1][0]=0.0;
    //specs2[1][1]=0.0;
		
    //  MHz
    specs2[2][0]=6.5;
    specs2[2][1]=6.5;
    //specs2[2][0]=0.0;
    //specs2[2][1]=0.0;
		
    //  MHz
    specs2[3][0]=5.5;
    specs2[3][1]=5.5;
    //specs2[3][0]=0.0;
    //specs2[3][1]=0.0;
		
		
		
		
  }
  else { // everything else
    // these are dimensionless gains
    // turn to dBi below
    // 300 MHz
    specs2[0][0]=8.5; // vp
    specs2[0][1]=8.8; // hp
		
    // 600 MHz
    specs2[1][0]=11.0;
    specs2[1][1]=9.2;
		
    // 900 MHz
    specs2[2][0]=9.3;
    specs2[2][1]=9.6;
		
    // 1200 MHz
    specs2[3][0]=10.1;
    specs2[3][1]=11.5;
		
    // 1500 MHz
    specs2[4][0]=8.9;
    specs2[4][1]=9.0;
		
  }
    
  // now convert to dimensionless units
  for (int i=0;i<4;i++) {
    for (int j=0;j<2;j++) {
      specs2[i][j]=pow(10.,specs2[i][j]/10.);
    }
  }
    
    
    
    
    
  double scale=0;
  // loop through frequency bins and fill the flare array
  for (int k=0;k<NFREQ;k++) {
    // if the frequency is below the lowest from the specs, just use the 300 MHz value
    if (freq[k]<freq_specs[0]) {
      for (int j=0;j<4;j++) {
	flare[j][k]=specs[0][j]*RADDEG;
      } //for
    } //if
    // if the frequency is higher than the highest frequency from the specs, just use the 1500 MHz value
    else if (freq[k]>=freq_specs[NFREQ_FORGAINS-1]) {
      for (int j=0;j<4;j++) {
	flare[j][k]=specs[NFREQ_FORGAINS-1][j]*RADDEG;
				
      } //for
    } //else if
    // if it is in between, interpolate
    else {
      for (int i=0;i<NFREQ_FORGAINS;i++) {
	if (freq[k]>=freq_specs[i] && freq[k]<freq_specs[i+1]) {
	  scale = (freq[k]-freq_specs[i])/(freq_specs[i+1]-freq_specs[i]);
					
	  for (int j=0;j<4;j++) {
	    flare[j][k]=(specs[i][j]+scale*(specs[i+1][j]-specs[i][j]))*RADDEG;
						
	  } //for
	  i=NFREQ_FORGAINS;
	} //if
      } //for
    } //else
  } //for (frequencies)
    
    // loop through frequencies to fill the gain array
  for (int k=0;k<NFREQ;k++) {
    // if below 300 MHz, use the 300 MHz value
    if (freq[k]<freq_specs[0]) {
      for (int j=0;j<2;j++) {
	gain[j][k]=specs2[0][j];
      } //for
    } //if
    // if higher than 1500 MHz, use the 1500 MHz value
    else if (freq[k]>=freq_specs[NFREQ_FORGAINS-1]) {
      for (int j=0;j<2;j++) {
	gain[j][k]=specs2[NFREQ_FORGAINS-1][j];
      } //for
    } //else if
    // if in between, interpolate
    else {
      for (int i=0;i<NFREQ_FORGAINS;i++) {
	if (freq[k]>=freq_specs[i] && freq[k]<freq_specs[i+1]) {
	  scale = (freq[k]-freq_specs[i])/(freq_specs[i+1]-freq_specs[i]);
					
	  for (int j=0;j<2;j++) {
	    gain[j][k]=specs2[i][j]+scale*(specs2[i+1][j]-specs2[i][j]);
						
	  } //for
	  i=NFREQ_FORGAINS;
	} //if
      } //for
    } //else
  } //for (frequencies)
    
  return 1;
} //GetBeamWidths
double Anita::Get_gain_angle(int gain_type, int k, double hitangle) {
  double scaleh1, scaleh2;
  if(gain_type < 0 || gain_type > 3) {
    cout << "gain_type out of range\n";
    exit(1);
  }
    
  hitangle = fabs(hitangle)*DEGRAD;
  flare[gain_type][k]=flare[gain_type][k]*DEGRAD;
  if(hitangle > 90.00001) {
    cout << "hitangle out of range\n";
    exit(1);
  }
  if(hitangle >= 90.) hitangle = 89.99999;
    
  // int pol;
  if (GAINS==0) {
    // for this simple model just treat the cross pols the same with regard to effect of being hit off-axis
    // the absolute gain will of course be different for these
    // pol=(int)((double)gain_type/2.);
		
    //    cout << "gain_type, k, flare, gain are " << gain_type << " " << k << " " << flare[gain_type][k] << " " << gain[pol][k] << "\n";
    //cout << "hitangle, flare, factor are " << hitangle << " " << flare[gain_type][k] << " " << exp(-1.*(hitangle*hitangle)/(2*flare[gain_type][k]*flare[gain_type][k])) << "\n";
    return exp(-1.*(hitangle*hitangle)/(2*flare[gain_type][k]*flare[gain_type][k]));
  } else {
    for(int iii = 1; iii < 7; iii++) { // linear interpolation for the angle
      if(hitangle <= reference_angle[iii]) {
	scaleh2 = (hitangle - reference_angle[iii-1]) *
	  inv_angle_bin_size[iii-1]; // how far from the smaller angle
	scaleh1 = 1. - scaleh2; // how far from the larger angle
				
	if(whichbin[k] == NPOINTS_GAIN - 1) // if the frequency is 1.5e9 or goes a little over due to rounding
	  return scaleh1 * gain_angle[gain_type][whichbin[k]][iii-1] +
	    scaleh2 * gain_angle[gain_type][whichbin[k]][iii];
				
	return scaleh1 * scalef1[k] * gain_angle[gain_type][whichbin[k]][iii-1] +
	  scaleh1 * scalef2[k] * gain_angle[gain_type][whichbin[k]+1][iii-1] +
	  scaleh2 * scalef1[k] * gain_angle[gain_type][whichbin[k]][iii] +
	  scaleh2 * scalef2[k] * gain_angle[gain_type][whichbin[k]+1][iii];
				
      }
    }
  }
    
    
  cout << "Get_gain_angle should have returned a value\n";
  exit(1);
} // double Get_gain_angle(int gain_type, double freq, double hitangle)


void Anita::getDiodeModel( ) {
    
    
  //  this is our homegrown diode response function which is a downgoing gaussian followed by an upward step function
  TF1 *fdown1=new TF1("fl_down1","[3]+[0]*exp(-1.*(x-[1])*(x-[1])/(2*[2]*[2]))",-300.E-9,300.E-9);
  fdown1->SetParameter(0,-0.8);
  //  fdown1->SetParameter(1,15.E-9);
  fdown1->SetParameter(1,15.E-9);
  fdown1->SetParameter(2,2.3E-9);
  //fdown1->SetParameter(2,0.5E-9);
  fdown1->SetParameter(3,0.);
    
  TF1 *fdown2=new TF1("fl_down2","[3]+[0]*exp(-1.*(x-[1])*(x-[1])/(2*[2]*[2]))",-300.E-9,300.E-9);
  fdown2->SetParameter(0,-0.2);
  //  fdown2->SetParameter(1,15.E-9);
  fdown2->SetParameter(1,15.E-9);
  fdown2->SetParameter(2,4.0E-9);
  //fdown2->SetParameter(2,0.5E-9);
  fdown2->SetParameter(3,0.);
    
    
  maxt_diode=70.E-9;
  idelaybeforepeak[0]=(int)(5.E-9/TIMESTEP);
  iwindow[0]=(int)(20.E-9/TIMESTEP);
  idelaybeforepeak[1]=(int)(5.E-9/TIMESTEP);
  iwindow[1]=(int)(20.E-9/TIMESTEP);
  idelaybeforepeak[2]=(int)(5.E-9/TIMESTEP);
  iwindow[2]=(int)(20.E-9/TIMESTEP);
  idelaybeforepeak[3]=(int)(5.E-9/TIMESTEP);
  iwindow[3]=(int)(20.E-9/TIMESTEP);
  idelaybeforepeak[4]=(int)(13.E-9/TIMESTEP);
  iwindow[4]=(int)(4.E-9/TIMESTEP);

    
  fdown1->Copy(fdiode);
    
  TF1 *f_up=new TF1("f_up","[0]*([3]*(x-[1]))^2*exp(-(x-[1])/[2])",-200.E-9,100.E-9);
    
  f_up->SetParameter(2,7.0E-9);
  f_up->SetParameter(0,1.);
  f_up->SetParameter(1,18.E-9);
  f_up->SetParameter(3,1.E9);
    
  f_up->SetParameter(0,-1.*sqrt(2.*PI)*(fdiode.GetParameter(0)*fdiode.GetParameter(2)+fdown2->GetParameter(0)*fdown2->GetParameter(2))/(2.*pow(f_up->GetParameter(2),3.)*1.E18));
    
  double sum=0.;
  for (int j=0;j<5;j++) {
		
    if (j==0)
      f_up->SetParameter(0,0.00275);
    else if (j==1)
      f_up->SetParameter(0,0.004544);
    else
      f_up->SetParameter(0,-1.*sqrt(2.*PI)*(fdiode.GetParameter(0)*fdiode.GetParameter(2)+fdown2->GetParameter(0)*fdown2->GetParameter(2))/(2.*pow(f_up->GetParameter(2),3.)*1.E18));
		
    for (int i=0;i<NFOUR/2;i++) {
			
      if (time[i]>0. && time[i]<maxt_diode) {
		    
	diode_real[j][i]=fdiode.Eval(time[i])+fdown2->Eval(time[i]);
	if (time[i]>f_up->GetParameter(1))
	  diode_real[j][i]+=f_up->Eval(time[i]);
		    
	sum+=diode_real[j][i];
      }
      else {
	diode_real[j][i]=0.;
	//cout << "diode: " << time[i] << " " << diode_real[i] << "\n";
      }
    }
  }
    
    
  // diode_real is the time domain response of the diode
}

void Anita::myconvlv(double *data,const int NFOUR,double *fdiode,double &mindiodeconvl,double &onediodeconvl,double *power_noise,double *diodeconv) {
    
  // NFOUR is the size array of fdiode_real, while data is NFOUR/2 sized array.
  // so we are going to make data array to NFOUR sized array with zero padding (see Numerical Recipes 3rd ED, 643 page)
  // We are actually using NFOUR/2 array for zero padding which might be too long in terms of efficiency but it's a simple way.
  // Returned diodeconv array is NFOUR sized array but actually initial NFOUR/2 is the information what we need
  // Last half NFOUR/2 array of diodeconv is result from zero padding and some spoiled information.
    
  const int length=NFOUR;
  // double data_copy[length];
  //double fdiode_real[length];
  double power_noise_copy[length];
    
    
  for (int i=0;i<NFOUR/2;i++) {
    power_noise_copy[i]=(data[i]*data[i])/impedence*TIMESTEP;
  }
  for (int i=NFOUR/2;i<NFOUR;i++) {
    power_noise_copy[i]=0.;
  }

  
  // TCanvas *c = new TCanvas("c");
  // string name = "oldnoise";
  // TGraph *gV = new TGraph(HALFNFOUR, fTimes, power_noise_copy);
  // gV->Draw("Al");
  // c->Print(Form("%s_%d_diodePower.png", name.c_str(), inu));
  // delete gV;
  // delete c;
    
  Tools::realft(power_noise_copy,1,length);
  //  realft(fdiode_real,1,length);
    
  double ans_copy[length];
    
    
    
  // the +, - sign in the equation might look different then Numerical Recipes 3rd ED 649 page.
  // our equation (below) is the situation when power_noise_copy's 0th bin starts to meet fdiode's 0th bin
  // while code in Numerical Recipes 649 page is the situation when power_noise_copy's final bin starts to meet fdiode's 0th bin
  // so, in our case, below equation is correct.
  for (int j=1;j<length/2;j++) {
    ans_copy[2*j]=(power_noise_copy[2*j]*fdiode[2*j]-power_noise_copy[2*j+1]*fdiode[2*j+1])/((double)length/2);
    ans_copy[2*j+1]=(power_noise_copy[2*j+1]*fdiode[2*j]+power_noise_copy[2*j]*fdiode[2*j+1])/((double)length/2);
  }
  ans_copy[0]=power_noise_copy[0]*fdiode[0]/((double)length/2);
  ans_copy[1]=power_noise_copy[1]*fdiode[1]/((double)length/2);
   
    
    
  Tools::realft(ans_copy,-1,length);
    
    
  // even if there are NFOUR array in diodeconv, only first NFOUR/2 is meaningful for us.
  for (int i=0;i<NFOUR;i++) {
    power_noise[i]=power_noise_copy[i];
    diodeconv[i]=ans_copy[i];
  }
    
    
    
  int iminsamp,imaxsamp; // find min and max samples such that
  // the diode response is fully overlapping with the noise waveform
  iminsamp=(int)(maxt_diode/TIMESTEP);
  // the noise waveform is NFOUR/2 long
  // then for a time maxt_diode/TIMESTEP at the end of that, the
  // diode response function is only partially overlappying with the
  // waveform in the convolution
  imaxsamp=NFOUR/2;
    
  //  cout << "iminsamp, imaxsamp are " << iminsamp << " " << imaxsamp << "\n";
  if (imaxsamp<iminsamp) {
    cout << "Noise waveform is not long enough for this diode response.\n";
    exit(1);
  }
    
  int ibin=((iminsamp+imaxsamp)-(iminsamp+imaxsamp)%2)/2;
  //cout << "ibin is " << ibin << "\n";
    
  // return the 50th sample, right in the middle
  onediodeconvl=diodeconv[ibin];
    
    
  mindiodeconvl=0.;
    
  for (int i=iminsamp;i<imaxsamp;i++) {
		
    if (diodeconv[i]<mindiodeconvl)
      mindiodeconvl=diodeconv[i];
		
		
  }
    
    
    
    
    
    
}







// get the frequency bin for a given frequency, the lower and upper limits and the
// number of bins

void Anita::GetPhases() {

  int iband;
  double corr,uncorr;
  double phase_corr,phase_uncorr;
  double phasor_x,phasor_y;
    
  for (int k=0;k<NFOUR/4;k++) { // loop through samples
    iband=Tools::findIndex(freq_bands[0],freq_forfft[2*k],NPOINTS_BANDS,freq_bands[0][0],freq_bands[0][NPOINTS_BANDS-1]);
		
		
    phases_rfcm[0][k]=2*PI*gRandom->Rndm(); // set phases at output of rfcm randoml
    phases_rfcm[1][k]=2*PI*gRandom->Rndm(); // set phases at output of rfcm randomly
		
		
    // now set phases at the lab chip
		
    if(iband < 0) corr = 0;
    else corr=correl_lab[iband];
    uncorr=1-corr;
    phase_corr=phases_rfcm[0][k];
    phase_uncorr=2*PI*gRandom->Rndm();
    phasor_x=corr*cos(phase_corr)+uncorr*cos(phase_uncorr);
    phasor_y=corr*sin(phase_corr)+uncorr*sin(phase_uncorr);
    phases_lab[0][k]=TMath::ATan2(phasor_y,phasor_x);
		
		
		
    phase_corr=phases_rfcm[1][k];
    phase_uncorr=2*PI*gRandom->Rndm();
    phasor_x=corr*cos(phase_corr)+uncorr*cos(phase_uncorr);
    phasor_y=corr*sin(phase_corr)+uncorr*sin(phase_uncorr);
    phases_lab[1][k]=TMath::ATan2(phasor_y,phasor_x);
		
		
    // do the same thing for the bands
    for (int j=0;j<5;j++) {
			
      if(iband < 0) corr = 0;
      else corr=correl_banding[j][iband];
      uncorr=1-corr;
      phase_corr=phases_rfcm[0][k];
      phase_uncorr=2*PI*gRandom->Rndm();
      phasor_x=corr*cos(phase_corr)+uncorr*cos(phase_uncorr);
      phasor_y=corr*sin(phase_corr)+uncorr*sin(phase_uncorr);
      phases_rfcm_banding[0][j][k]=TMath::ATan2(phasor_y,phasor_x);
      phase_corr=phases_rfcm[1][k];
      phase_uncorr=2*PI*gRandom->Rndm();
      phasor_x=corr*cos(phase_corr)+uncorr*cos(phase_uncorr);
      phasor_y=corr*sin(phase_corr)+uncorr*sin(phase_uncorr);
      phases_rfcm_banding[1][j][k]=TMath::ATan2(phasor_y,phasor_x);

    }
  }
    
    

    
}

void Anita::normalize_for_nsamples(double *spectrum, double nsamples, double nsamp){
  for (int k = 0; k < NFOUR / 4; k++){
    spectrum[2 * k] *= sqrt((double) nsamples / (double) nsamp);
    spectrum[2 * k + 1] *= sqrt((double) nsamples / (double) nsamp);
  }
  return;
}

void Anita::convert_power_spectrum_to_voltage_spectrum_for_fft(int length,double *spectrum, double domain[], double phase[]){
  double current_amplitude, current_phase;
  //    for (int k = 0; k < NFOUR / 4; k++){
  for (int k = 0; k < length/2 ; k++){
    current_amplitude = domain[k];
    current_phase = phase[k];
      
    // Uncomment the following line of code to allow the the noise generated to have a Rician/Rayleigh distribution, which should increase the number of weighted neutrinos that pass *slightly*.
    Tools::get_random_rician(0., 0., sqrt(2. / M_PI) * domain[k], current_amplitude, current_phase);
    // The above function uses 0. as the mean of the rician distribution, and uses sqrt(2. / M_PI) * domain[k] as the standard deviation for the distribution, as discussed in Goodman's "Statistical Optics", equation 2.9-13 on page 49
      
    spectrum[2 * k] = sqrt(current_amplitude) * cos(phase[k]) / (( (double) NFOUR / 2) / 2);
    spectrum[2 * k + 1] = sqrt(current_amplitude) * sin(phase[k]) / (( (double) NFOUR / 2) / 2);
  }
  return;
}

void Anita::GetNoiseWaveforms() {
  GetPhases();
  int nsamples = NFOUR / 2;
  // int nsamples_long = NFOUR;
  double sumfreqdomain = 0.;
  double sumtimedomain = 0.;

  count_getnoisewaveforms++;


  // This is done in a stupid way for the moment to provide the same order of gRandom calls
  // So that I have the exact same results as masters
  // This should be done in 1 loop once I merge the trigger branch with master
  // LC, 16/02/17

  
  for (int ipol=0;ipol<2;ipol++)
    convert_power_spectrum_to_voltage_spectrum_for_fft(nsamples, timedomainnoise_rfcm[ipol], freqdomain_rfcm,   phases_rfcm[ipol] );
  for (int ipol=0;ipol<2;ipol++)
    convert_power_spectrum_to_voltage_spectrum_for_fft(nsamples, timedomainnoise_lab[ipol],  avgfreqdomain_lab, phases_lab[ipol]  );

  //     convert_power_spectrum_to_voltage_spectrum_for_fft(nsamples_long, timedomainnoise_rfcm_long[ipol], freqdomain_rfcm_long,   phases_rfcm_long[ipol] );
  //     convert_power_spectrum_to_voltage_spectrum_for_fft(nsamples_long, timedomainnoise_lab_long[ipol],  avgfreqdomain_lab_long, phases_lab_long[ipol]  );
    
  // want to restrict it to NFOUR/2 samples -# samples that equal twice
  // maxt_diode
  // freqdomain_rfcm_banding was made by averaging many nsamp=100-sample noise waveforms in the
  // frequency domain
  // the noise waveform that we create from it has NFOUR/2 samples
  // I don't know how to restrict it to nsamp=100 samples.
  // So in order to restrict the power to nsamp=100 samples we zero
  // everything but the middle nsamp samples and scale the waveforms
  // by sqrt((NFOUR/2)/nsamp)
    
  for (int ipol=0;ipol<2;ipol++){
    normalize_for_nsamples(timedomainnoise_rfcm[ipol], (double) nsamples, (double) nsamp);
    normalize_for_nsamples(timedomainnoise_lab[ipol],  (double) nsamples, (double) nsamp);

    //     normalize_for_nsamples(timedomainnoise_rfcm_long[ipol], (double) nsamples_long, (double) nsamp);
    //     normalize_for_nsamples(timedomainnoise_lab_long[ipol],  (double) nsamples_long, (double) nsamp);
    
    for (int k = 0; k < NFOUR / 4; k++){
      sumfreqdomain += avgfreqdomain_lab[k];
    }
    
    Tools::realft(timedomainnoise_rfcm[ipol], -1, NFOUR / 2);
    Tools::realft(timedomainnoise_lab[ipol],  -1, NFOUR / 2);

    //    Tools::realft(timedomainnoise_rfcm_long[ipol], -1, NFOUR );
    //     Tools::realft(timedomainnoise_lab_long[ipol],-1, NFOUR );
    
    for (int k = 0; k < NFOUR / 2; k++) {
      timedomainnoise_lab[ipol][k] *=THERMALNOISE_FACTOR;
      timedomainnoise_rfcm[ipol][k]*=THERMALNOISE_FACTOR;
      
      if (ipol==0){
	sumtimedomain += timedomainnoise_lab[ipol][k] * timedomainnoise_lab[ipol][k];
	rms_rfcm_e_single_event += timedomainnoise_rfcm[ipol][k] * timedomainnoise_rfcm[ipol][k];
      }
      rms_rfcm[ipol] += timedomainnoise_rfcm[ipol][k] * timedomainnoise_rfcm[ipol][k] / ((double) NFOUR / 2);
      rms_lab[ipol]  += timedomainnoise_lab[ipol][k] * timedomainnoise_lab[ipol][k] / ((double) NFOUR / 2);
            
    }
    
  }

  for (int iband=0; iband<5; iband++) {
    for (int ipol=0;ipol<2;ipol++)
      convert_power_spectrum_to_voltage_spectrum_for_fft(NFOUR/2,timedomainnoise_rfcm_banding[ipol][iband], freqdomain_rfcm_banding[iband], phases_rfcm_banding[ipol][iband]);
    for (int ipol=0;ipol<2;ipol++)
      normalize_for_nsamples(timedomainnoise_rfcm_banding[ipol][iband], (double) nsamples, (double) nsamp);
    for (int ipol=0;ipol<2;ipol++)
      Tools::realft(timedomainnoise_rfcm_banding[ipol][iband], -1, NFOUR / 2);

    //      convert_power_spectrum_to_voltage_spectrum_for_fft(NFOUR,timedomainnoise_rfcm_banding_long[ipol][iband], freqdomain_rfcm_banding_long[iband], phases_rfcm_banding_e_long[iband]);
    //       normalize_for_nsamples(timedomainnoise_rfcm_banding_long[ipol][iband], (double) nsamples_long, (double) nsamp);
    //       Tools::realft(timedomainnoise_rfcm_banding_long[ipol][iband], -1, NFOUR);

    
    
  }
}

void Anita::GetArrayFromFFT(double *tmp_fftvhz, double *vhz_rx){
  
  int firstNonZero = Tools::Getifreq(freq[0],freq_forfft[0],freq_forfft[NFOUR/2-1],NFOUR/4);
  int lastNonZero  = Tools::Getifreq(freq[NFREQ-1],freq_forfft[0],freq_forfft[NFOUR/2-1],NFOUR/4);
  double norm=TMath::Sqrt(double(lastNonZero-firstNonZero)/double(NFREQ));
  //  cout << firstNonZero << " " << lastNonZero << " " << lastNonZero-firstNonZero << " " << norm << endl;

  for (int ifreq=0;ifreq<NFOUR/4;ifreq++){
    tmp_fftvhz[ifreq]=TMath::Sqrt(tmp_fftvhz[2*ifreq]*tmp_fftvhz[2*ifreq] + tmp_fftvhz[2*ifreq+1]*tmp_fftvhz[2*ifreq+1]);
  }
    
  for (int ifreq=0; ifreq<NFREQ; ifreq++){
      
    int ifour=Tools::Getifreq(freq[ifreq],freq_forfft[0],freq_forfft[NFOUR/2-1],NFOUR/4);
    vhz_rx[ifreq]=tmp_fftvhz[ifour]*norm;
    //cout << ifour << " " << freq[ifreq] << " " << vhz_rx[ifreq] << " " << endl;
  }

}

void Anita::MakeArraysforFFT(double *vsignalarray_e,double *vsignalarray_h,double *vsignal_e_forfft,double *vsignal_h_forfft, double phasedelay, bool useconstantdelay) {
    
  Tools::Zero(vsignal_e_forfft,NFOUR/2);
  Tools::Zero(vsignal_h_forfft,NFOUR/2);
    
  double previous_value_e_even=0.;
  double previous_value_e_odd=0.;
  double previous_value_h_even=0.;
  double previous_value_h_odd=0.;
  int count_nonzero=0;
  int iprevious=0;
  int ifirstnonzero=-1;
  int ilastnonzero=2000;
  for (int i=0;i<NFREQ;i++) {
    // freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
    // but there are only NFOUR/4 different values
    // it's the index among the NFOUR/4 that we're interested in
    int ifour=Tools::Getifreq(freq[i],freq_forfft[0],freq_forfft[NFOUR/2-1],NFOUR/4);
      
    if (ifour!=-1 && 2*ifour+1<NFOUR/2) {
      count_nonzero++;
      if (ifirstnonzero==-1){
	ifirstnonzero=ifour;
      }

      vsignal_e_forfft[2*ifour]=vsignalarray_e[i]*2/((double)NFOUR/2); // phases is 90 deg.

      //      cout << "ifour, vsignal is " << ifour << " " << vsignal_e_forfft[2*ifour] << "\n";

      vsignal_e_forfft[2*ifour+1]=vsignalarray_e[i]*2/((double)NFOUR/2); // phase is 90 deg.
      // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting

      vsignal_h_forfft[2*ifour]=vsignalarray_h[i]*2/((double)NFOUR/2); // what to do about phases
      vsignal_h_forfft[2*ifour+1]=vsignalarray_h[i]*2/((double)NFOUR/2); // what to do about phases
      // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting
      // how about we interpolate instead of doing a box average

      for (int j=iprevious+1;j<ifour;j++) {
        vsignal_e_forfft[2*j]=previous_value_e_even+(vsignal_e_forfft[2*ifour]-previous_value_e_even)*(double)(j-iprevious)/(double)(ifour-iprevious);
        //	cout << "j, vsignal is " << j << " " << vsignal_e_forfft[2*j] << "\n";
        vsignal_h_forfft[2*j]=previous_value_h_even+(vsignal_h_forfft[2*ifour]-previous_value_h_even)*(double)(j-iprevious)/(double)(ifour-iprevious);

        vsignal_e_forfft[2*j+1]=previous_value_e_odd+(vsignal_e_forfft[2*ifour+1]-previous_value_e_odd)*(double)(j-iprevious)/(double)(ifour-iprevious);
        vsignal_h_forfft[2*j+1]=previous_value_h_odd+(vsignal_h_forfft[2*ifour+1]-previous_value_h_odd)*(double)(j-iprevious)/(double)(ifour-iprevious);
      }

      ilastnonzero=ifour;
      iprevious=ifour;
      previous_value_e_even=vsignal_e_forfft[2*ifour];
      previous_value_e_odd=vsignal_e_forfft[2*ifour+1];
      previous_value_h_even=vsignal_h_forfft[2*ifour];
      previous_value_h_odd=vsignal_h_forfft[2*ifour+1];
    }
      
  } // end loop over nfreq
    
    // EH check
    // cout << "ifirstnonzero, ilastnonzero are " << ifirstnonzero << " " << ilastnonzero << "\n";
    // cout << "ratio is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
  for (int j=0;j<NFOUR/4;j++) {
    vsignal_e_forfft[2*j]*=sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
    vsignal_e_forfft[2*j+1]*=sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
    vsignal_h_forfft[2*j]*=sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
    vsignal_h_forfft[2*j+1]*=sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
  }
    
  //  Tools::InterpolateComplex(vsignal_e_forfft,NFOUR/4);
  //Tools::InterpolateComplex(vsignal_h_forfft,NFOUR/4);

  if (useconstantdelay){
    double cosphase=cos(phasedelay*PI/180.);
    double sinphase=sin(phasedelay*PI/180.);
    for (int ifour=0;ifour<NFOUR/4;ifour++) {      
      if (PULSER) {
	cosphase = cos(v_phases[ifour]*PI/180.);
	sinphase = sin(v_phases[ifour]*PI/180.);
      }
      vsignal_e_forfft[2*ifour]*=cosphase;
      vsignal_e_forfft[2*ifour+1]*=sinphase;
      vsignal_h_forfft[2*ifour]*=cosphase;
      vsignal_h_forfft[2*ifour+1]*=sinphase;	
    }
  }
}

void Anita::MakeArrayforFFT(double *vsignalarray_e,double *vsignal_e_forfft, double phasedelay, bool useconstantdelay) {
    
  Tools::Zero(vsignal_e_forfft,NFOUR/2);
    
  double previous_value_e_even=0.;
  double previous_value_e_odd=0.;
  int count_nonzero=0;
  int iprevious=0;
  int ifirstnonzero=-1;
  int ilastnonzero=2000;
  for (int i=0;i<NFREQ;i++) {
    // freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
    // but there are only NFOUR/4 different values
    // it's the index among the NFOUR/4 that we're interested in
    int ifour=Tools::Getifreq(freq[i],freq_forfft[0],freq_forfft[NFOUR/2-1],NFOUR/4);
      
    if (ifour!=-1 && 2*ifour+1<NFOUR/2) {
      count_nonzero++;
      if (ifirstnonzero==-1){
	ifirstnonzero=ifour;
      }

      vsignal_e_forfft[2*ifour]=vsignalarray_e[i]*2/((double)NFOUR/2); // phases is 90 deg.

      //      cout << "ifour, vsignal is " << ifour << " " << vsignal_e_forfft[2*ifour] << "\n";

      vsignal_e_forfft[2*ifour+1]=vsignalarray_e[i]*2/((double)NFOUR/2); // phase is 90 deg.
      // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting

      // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting
      // how about we interpolate instead of doing a box average

      for (int j=iprevious+1;j<ifour;j++) {
        vsignal_e_forfft[2*j]=previous_value_e_even+(vsignal_e_forfft[2*ifour]-previous_value_e_even)*(double)(j-iprevious)/(double)(ifour-iprevious);
        //	cout << "j, vsignal is " << j << " " << vsignal_e_forfft[2*j] << "\n";

        vsignal_e_forfft[2*j+1]=previous_value_e_odd+(vsignal_e_forfft[2*ifour+1]-previous_value_e_odd)*(double)(j-iprevious)/(double)(ifour-iprevious);
      }

      ilastnonzero=ifour;
      iprevious=ifour;
      previous_value_e_even=vsignal_e_forfft[2*ifour];
      previous_value_e_odd=vsignal_e_forfft[2*ifour+1];
    }
      
  } // end loop over nfreq
    
    // EH check
    // cout << "ifirstnonzero, ilastnonzero are " << ifirstnonzero << " " << ilastnonzero << "\n";
    // cout << "ratio is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
  for (int j=0;j<NFOUR/4;j++) {
    vsignal_e_forfft[2*j]*=sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
    vsignal_e_forfft[2*j+1]*=sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
  }
    
  //  Tools::InterpolateComplex(vsignal_e_forfft,NFOUR/4);

  if (useconstantdelay){
    double cosphase=cos(phasedelay*PI/180.);
    double sinphase=sin(phasedelay*PI/180.);
    for (int ifour=0;ifour<NFOUR/4;ifour++) {      
      if (PULSER) {
	cosphase = cos(v_phases[ifour]*PI/180.);
	sinphase = sin(v_phases[ifour]*PI/180.);
      }
      vsignal_e_forfft[2*ifour]*=cosphase;
      vsignal_e_forfft[2*ifour+1]*=sinphase;
    }
  }
}


void Anita::BoxAverageComplex(double *array, const int n,int navg) {
  // to get rid of the zero bins
  double array_temp[2*n];
  for (int i=0;i<n;i++) {
		
    array_temp[2*i]=0.;
    array_temp[2*i+1]=0.;
		
		
    for (int k=i;k<i+navg;k++) {
      if (k<n) {
	array_temp[2*i]+=array[2*k];
	array_temp[2*i+1]+=array[2*k+1];
      }
    }
		
    array[2*i]=array_temp[2*i]/(double)navg;
    array[2*i+1]=array_temp[2*i+1]/(double)navg;
		
    array_temp[2*i]=array[2*i];
    array_temp[2*i+1]=array[2*i+1];
		
  }
}
void Anita::BoxAverage(double *array, const int n,int navg) {
  // to get rid of the zero bins
  double array_temp[n];
  for (int i=0;i<n;i++) {
		
    array_temp[i]=0.;
		
    for (int k=i;k<i+navg;k++) {
      if (k<n) {
	array_temp[i]+=array[k];
      }
    }
		
    array[i]=array_temp[i]/(double)navg;
		
		
    array_temp[i]=array[i];
		
  }
}


void Anita::labAttn(double *vhz) {
  for (int i=0;i<NFREQ;i++) {
    // next find the index of the lab attenuation array
    int ilab=Tools::findIndex(freqlab,freq[i],NPOINTS_LAB,freqlab[0],freqlab[NPOINTS_LAB-1]);
    if (ilab!=-1) {
      vhz[i]*=sqrt(labattn[ilab]);
      vhz[i]*=sqrt(labattn[ilab]);
    }
    else
      cout << "Lab attenuation outside of band.\n";
  }
}

int Anita::getLabAttn(int NPOINTS_LAB,double *freqlab,double *labattn) {
    
  ifstream flab((ICEMC_DATA_DIR+"/surfatten_run294_ch23v.dat").c_str());
  if (flab.fail()) {
    cout << "Cannot open lab data file.\n";
    exit(1);
  }
  int index=0;
    
  string line;
  getline(flab,line); // first two lines are junk
  getline(flab,line);
    
  float ffreqlab,flabattn;
    
  // read in the rest of the lab data from a file
  while (!flab.eof() && index<NPOINTS_LAB) {
		
    flab >> ffreqlab >> flabattn;
    labattn[index]=(double)pow(10.,flabattn/10.);
    //cout << "Mofifying labattn!\n";
    //labattn[index]*=2.;
		
    freqlab[index]=(double)ffreqlab;
		
    index++;
		
  }
    
    
  flab.close();
    
  return 1;
    
}


double Anita::GaintoHeight(double gain,double freq,double nmedium_receiver) {
    
  // from gain=4*pi*A_eff/lambda^2
  // and h_eff=2*sqrt(A_eff*Z_rx/Z_air)
  // gain is in dB
  return 2*sqrt(gain/4/PI*CLIGHT*CLIGHT/(freq*freq)*Zr/(Z0*nmedium_receiver));
} //GaintoHeight


void Anita::fill_coherent_waveform_sum_tree(unsigned event_number, unsigned center_phi_sector_index, Settings* settings1, double rms_noise, double actual_rms, unsigned window_start, unsigned window_end, double deg_theta, double deg_phi, double actual_deg_theta, double actual_deg_phi, vector <double>& summed_wfm, vector <double>& power_of_summed_wfm, double power){
    
  cwst_event_number = event_number;
  cwst_center_phi_sector = center_phi_sector_index;
  cwst_rms_noise = rms_noise;
  cwst_actual_rms = actual_rms;
  cwst_threshold = settings1->COHERENT_THRESHOLD;
  cwst_window_start = window_start;
  cwst_window_end = window_end;
    
  cwst_deg_theta = deg_theta;
  cwst_deg_phi = deg_phi;
    
  cwst_actual_deg_theta = actual_deg_theta;
  cwst_actual_deg_phi = actual_deg_phi;
    
  // Reassign the array associated with the timesteps branch
  for (int i = 0; i < HALFNFOUR; i++){
    cwst_timesteps[i] = i;
  }
    
  cwst_summed_wfm = summed_wfm;
  cwst_power_of_summed_wfm = power_of_summed_wfm;
  cwst_power = power;
    
  coherent_waveform_sum_tree->Fill();
    
  return;
}

void Anita::GetPayload(Settings* settings1, Balloon* bn1){
  // anita-lite payload
  // see comments next to variable definitions
  double temp_eachrx[Anita::NPHI_MAX]; // temperature of each antenna (for the anita-lite configuration)
    
  const double gps_offset = atan2(-0.7042,0.71), MINCH = 0.0254, phase_center = 0.17;
  // const double phase_center_anita2=0.17;
  const double phase_center_anita2_analysis=.2;
  //const double gps_offset_anita2=atan2(0.89,-0.29);
  const double gps_offset_anita2=atan2(-0.7085,0.7056); // from elog 473
  const double gps_offset_anita3= 45*RADDEG; // Linda: 45 degrees from EventReader
  const double phase_center_anita3=0.20; // Linda: phase-centers are around 20 cm inwards of antennas face-end


  if (settings1->WHICH==0) { // anita-lite
		
    settings1->NFOLD=3;
		
    settings1->CYLINDRICALSYMMETRY=0;
    PHI_EACHLAYER[0][0]=0.;
    PHI_EACHLAYER[0][1]=22.5*RADDEG;
    NRX_PHI[0]=2;
    PHI_OFFSET[0]=0;
    THETA_ZENITH[0]=PI/2+10.*RADDEG;
    LAYER_VPOSITION[0]=0.; // vertical separation between layers
    LAYER_HPOSITION[0]=0.;  // position of layer relative to center axis in horizontal plane
    LAYER_PHIPOSITION[0]=0.; // phi of position of layer
		
    for (int i=0;i<5;i++) {
      RRX[i]=3.006;
    }
		
    for (int i=0;i<NRX_PHI[0];i++) {
      VNOISE_ANITALITE[i]=ChanTrigger::GetNoise(settings1,bn1->altitude_bn,bn1->surface_under_balloon,THETA_ZENITH[i],settings1->BW_SEAVEYS,temp_eachrx[i]);
    }
  } //if (ANITA-lite)
    
    //Ross Payload
  else if (settings1->WHICH==1) {
    //settings1->NFOLD=8;
    settings1->CYLINDRICALSYMMETRY=1;
		
    NRX_PHI[0]=5;
    NRX_PHI[1]=5;
    NRX_PHI[2]=5;
    NRX_PHI[3]=5;
    NRX_PHI[4]=4;
		
    PHI_OFFSET[0]=0;
    PHI_OFFSET[1]=2*PI/(double)NRX_PHI[1]/2;
    PHI_OFFSET[2]=0;
    PHI_OFFSET[3]=2*PI/(double)NRX_PHI[3]/2;
    PHI_OFFSET[4]=0;
		
    THETA_ZENITH[0]=PI/2+10.*RADDEG;
    THETA_ZENITH[1]=PI/2+10.*RADDEG;
    THETA_ZENITH[2]=PI/2+10.*RADDEG;
    THETA_ZENITH[3]=PI/2+10.*RADDEG;
    THETA_ZENITH[4]=3*PI/4;
		
    // anita proposal "says that the separation between upper and lower
    // 2 layers of antennas is just under 4m.
    LAYER_VPOSITION[0]=3.5;
    LAYER_VPOSITION[1]=2.0;
    LAYER_VPOSITION[2]=-0.5;
    LAYER_VPOSITION[3]=-2.0;
    LAYER_VPOSITION[4]=-3.5;
		
    LAYER_HPOSITION[0]=0.;
    LAYER_HPOSITION[1]=0.;
    LAYER_HPOSITION[2]=0.;
    LAYER_HPOSITION[3]=0.;
    LAYER_HPOSITION[4]=0.;
		
    LAYER_PHIPOSITION[0]=0.;
    LAYER_PHIPOSITION[1]=0.;
    LAYER_PHIPOSITION[2]=0.;
    LAYER_PHIPOSITION[3]=0.;
    LAYER_PHIPOSITION[4]=0.;
		
    // position of layers in z relative to vertical center of the payload
    // radius that antennas sit at on the payload
    for (int i=0;i<5;i++) {
      RRX[i]=0.5;
    }
		
  } //else if (Ross payload)
    // Smex payload (full ANITA flown 2006-2007)
  else if (settings1->WHICH==2) {
    settings1->CYLINDRICALSYMMETRY=1;
		
    //these are physical layers
    NRX_PHI[0]=8;
    NRX_PHI[1]=8;
    NRX_PHI[2]=16;
    NRX_PHI[3]=8;
		
    PHITRIG[0]=16; // number of positions in phi in each *trigger* layer
    PHITRIG[1]=16;
    PHITRIG[2]=8;
		
    //these are physical layers again
    PHI_OFFSET[0]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    // it's rotated by 1/2 the azimuth that separates two antennas on the 0th layer
    PHI_OFFSET[1]=-2.*PI/(double)NRX_PHI[0]/2.;
    PHI_OFFSET[2]=-2.*PI/(double)NRX_PHI[0]/2.;
    PHI_OFFSET[3]=-2.*PI/(double)NRX_PHI[0]/4.;
		
    //double INCLINE_NADIR=55; // this is set in the input file now EWG: so why not remove it?
		
    // sets their declination
    THETA_ZENITH[0]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[1]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[2]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[3]=PI/2+INCLINE_NADIR*RADDEG;
		
    // radius from center axis of the payload
    RRX[0] = 0.9210802;
    RRX[1] = 0.7553198;
    RRX[2] = 2.0645374;
    RRX[3]=1.40; // this is wrong, but I put this version of icemc here just to show one possible strategy to use V polarization and nadirs so it doesn't matter
		
    // vertical separation between layers.
    LAYER_VPOSITION[0]=0;
    LAYER_VPOSITION[1] = -0.9670034;
    LAYER_VPOSITION[2] = -3.730752;
    LAYER_VPOSITION[3]=(-3.175-0.948-0.889); // this is wrong too, but nadirs are probably gone so doesn't matter.
		
    LAYER_HPOSITION[0]=0.;
    LAYER_HPOSITION[1] = 0.;
    LAYER_HPOSITION[2] = 0.;
    LAYER_HPOSITION[3]=0.; // this is wrong too, but nadirs are probably gone so doesn't matter.
		
		
  } //else if (SMEX payload - default)
  else if (settings1->WHICH==3) {
		
    cout << "Is this configuration cylindrically symmetric? Yes(1) or No(0)\n";
    cin >> settings1->CYLINDRICALSYMMETRY;
		
    cout << "How many layers?\n";
    cin >> settings1->NLAYERS;
		
    for (int i=0;i<settings1->NLAYERS;i++) {
			
      cout << "How many antennas in the " << i << "th layer?\n";
      cin >> NRX_PHI[i];
			
      cout << "What is the offset in phi for the " << i << "th layer?\n";
      cin >> PHI_OFFSET[i];
			
      cout << "What is the theta ascent for the " << i << "th layer (0 if pointed straight upwards, PI if pointed straight downwards)?\n";
      cin >> THETA_ZENITH[i];
			
      cout << "What is the vertical position of this layer relative to the vertical center of the payload?";
      cin >> LAYER_VPOSITION[i];
			
      cout << "What is the distance between of the vertical axis of the payload and the center of this layer?";
      cin >> LAYER_HPOSITION[i];
			
      cout << "What is the phi of this layer relative to the vertical center of the payload in the horizontal plane?";
      cin >> LAYER_PHIPOSITION[i];
			
			
			
      if (settings1->CYLINDRICALSYMMETRY==0) {
	for (int j=0;j<NRX_PHI[i];j++) {
	  cout << "What is the phi of the " << j << "th antenna is this layer?\n";
	  cin >> PHI_EACHLAYER[i][j];
	} //for (read antenna phi)
      }//if (not cylindrically symmetric)
    } //for (antenna layers)
		
    cout << "How many polarizations must pass a voltage threshold?\n";
    cin >> settings1->NFOLD;
		
    cout << "How many times the expected noise level should the voltage threshold be?\n";
    cin >> maxthreshold;
  } //else if (custom payload)
    
  else if (settings1->WHICH==4) {// anita hill
		
    //settings1->NFOLD=1; // how many channels must pass the trigger
    //maxthreshold=2.3;
		
    if (settings1->NLAYERS!=2)
      cout << "Warning!!! Did not enter the right number of layers in the input file.  For Anita Hill, it's 2.";
		
    settings1->CYLINDRICALSYMMETRY=1;
		
    NRX_PHI[0]=1; // this is how many antennas we have in phi on each "layer"
    NRX_PHI[1]=1; // for anita hill, we are calling each station a different "layer"
		
		
    PHI_OFFSET[0]=(bn1->BN_LONGITUDE+0.)*RADDEG; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    // it's rotated by 1/2 the azimuth that separates two antennas on the 0th layer
    PHI_OFFSET[1]=(bn1->BN_LONGITUDE+0.)*RADDEG;
		
    cout << "phi_offsets are " << PHI_OFFSET[0] << " " << PHI_OFFSET[1] << "\n";
		
    //double INCLINE_NADIR=55; SET IN INPUT NOW!!!!
		
    // sets their declination
    THETA_ZENITH[0]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[1]=PI/2+INCLINE_TOPTHREE*RADDEG;
		
    // radius from center axis of the payload
    RRX[0] = 1.;
    RRX[1] = 1.;
		
    // vertical separation between layers.
    LAYER_VPOSITION[0]=0.;
    LAYER_VPOSITION[1] = 0.;
		
    LAYER_HPOSITION[0]=0.;
    LAYER_HPOSITION[1] = 100.; // in meters
		
    LAYER_PHIPOSITION[0]=0.;
    LAYER_PHIPOSITION[1] = 30.*RADDEG;// in radians
  }
    
  else if(settings1->WHICH==6) { // Kurt's measurements for the first flight in elog 345
    //settings1->NFOLD=8;
    settings1->CYLINDRICALSYMMETRY=0;
		
    NRX_PHI[0]=8;
    NRX_PHI[1]=8;
    NRX_PHI[2]=16;
		
    PHITRIG[0]=16; // number of positions in phi in each trigger layer
    PHITRIG[1]=16;
    PHITRIG[2]=8;
		
		
    THETA_ZENITH[0]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[1]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[2]=PI/2+INCLINE_TOPTHREE*RADDEG;
		
    PHI_OFFSET[0]=0.;
    PHI_OFFSET[1]=0.;
    PHI_OFFSET[2]=0.;
		
    ANTENNA_POSITION_START[0][0][0] = MINCH * Vector(40.957,-39.29,125.38).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][1] = MINCH * Vector(57.608,0.606,124.97).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][2] = MINCH * Vector(41.049,40.605,124.896).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][3] = MINCH * Vector(1.032,57.128,124.93).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][4] = MINCH * Vector(-38.832,40.508,125.607).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][5] = MINCH * Vector(-55.545,0.549,125.851).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][6] = MINCH * Vector(-38.793,-39.423,126.105).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][0][7] = MINCH * Vector(1.113,-55.918,125.731).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][0] = MINCH * Vector(19.841,-45.739,87.738).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][1] = MINCH * Vector(46.959,-18.601,87.364).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][2] = MINCH * Vector(46.983,19.646,87.208).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][3] = MINCH * Vector(19.823,46.633,87.194).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][4] = MINCH * Vector(-18.429,46.496,87.486).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][5] = MINCH * Vector(-45.439,19.34,88.0).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][6] = MINCH * Vector(-45.446,-18.769,88.183).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][1][7] = MINCH * Vector(-18.297,-45.857,88.066).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][0] = MINCH * Vector(38.622,-94.184,-20.615).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][1] = MINCH * Vector(71.662,-72.036,-20.966).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][2] = MINCH * Vector(93.857,-39.130,-21.512).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][3] = MINCH * Vector(101.664,-0.202,-21.942).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][4] = MINCH * Vector(93.993,38.733,-22.351).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][5] = MINCH * Vector(71.92,71.815,-22.386).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][6] = MINCH * Vector(38.944,93.885,-22.274).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][7] = MINCH * Vector(0.036,101.619,-21.813).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][8] = MINCH * Vector(-38.991,93.809,-21.281).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][9] = MINCH * Vector(-71.899,71.704,-20.777).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][10] = MINCH * Vector(-94.121,38.754,-20.825).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][11] = MINCH * Vector(-101.986,-0.133,-20.71).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][12] = MINCH * Vector(-94.175,-39.024,-20.671).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][13] = MINCH * Vector(-72.138,-72.132,-20.295).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][14] = MINCH * Vector(-39.111,-94.25,-20.23).RotateZ(-gps_offset);
    ANTENNA_POSITION_START[0][2][15] = MINCH * Vector(-0.163,-101.975,-20.229).RotateZ(-gps_offset);
    PHI_EACHLAYER[0][0] = -45.103 * RADDEG - gps_offset;
    PHI_EACHLAYER[0][1] = -0.14 * RADDEG - gps_offset;
    PHI_EACHLAYER[0][2] = 44.559 * RADDEG - gps_offset;
    PHI_EACHLAYER[0][3] = 89.959 * RADDEG - gps_offset;
    PHI_EACHLAYER[0][4] = 135.555 * RADDEG - gps_offset;
    PHI_EACHLAYER[0][5] = 179.651 * RADDEG - gps_offset;
    PHI_EACHLAYER[0][6] = -135.14 * RADDEG - gps_offset;
    PHI_EACHLAYER[0][7] = -90.18 * RADDEG - gps_offset;
    PHI_EACHLAYER[1][0] = -67.283 * RADDEG - gps_offset;
    PHI_EACHLAYER[1][1] = -23.004 * RADDEG - gps_offset;
    PHI_EACHLAYER[1][2] = 22.72 * RADDEG - gps_offset;
    PHI_EACHLAYER[1][3] = 67.82 * RADDEG - gps_offset;
    PHI_EACHLAYER[1][4] = 112.698 * RADDEG - gps_offset;
    PHI_EACHLAYER[1][5] = 157.565 * RADDEG - gps_offset;
    PHI_EACHLAYER[1][6] = -157.376 * RADDEG - gps_offset;
    PHI_EACHLAYER[1][7] = -112.449 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][0] = -67.26 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][1] = -45.284 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][2] = -22.457 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][3] = 0.227 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][4] = 22.318 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][5] = 45.008 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][6] = 67.751 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][7] = 89.913 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][8] = 113.016 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][9] = 135.608 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][10] = 157.487 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][11] = 179.709 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][12] = -157.569 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][13] = -135.021 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][14] = -112.773 * RADDEG - gps_offset;
    PHI_EACHLAYER[2][15] = -89.959 * RADDEG - gps_offset;
    ANTENNA_DOWN[0][0] = 10.422 * RADDEG;
    ANTENNA_DOWN[0][1] = 10.207 * RADDEG;
    ANTENNA_DOWN[0][2] = 10.714 * RADDEG;
    ANTENNA_DOWN[0][3] = 10.381 * RADDEG;
    ANTENNA_DOWN[0][4] = 10.026 * RADDEG;
    ANTENNA_DOWN[0][5] = 9.515 * RADDEG;
    ANTENNA_DOWN[0][6] = 9.677 * RADDEG;
    ANTENNA_DOWN[0][7] = 9.544 * RADDEG;
    ANTENNA_DOWN[1][0] = 10.183 * RADDEG;
    ANTENNA_DOWN[1][1] = 10.44 * RADDEG;
    ANTENNA_DOWN[1][2] = 10.562 * RADDEG;
    ANTENNA_DOWN[1][3] = 10.655 * RADDEG;
    ANTENNA_DOWN[1][4] = 10.265 * RADDEG;
    ANTENNA_DOWN[1][5] = 9.77 * RADDEG;
    ANTENNA_DOWN[1][6] = 9.422 * RADDEG;
    ANTENNA_DOWN[1][7] = 9.526 * RADDEG;
    ANTENNA_DOWN[2][0] = 9.364 * RADDEG;
    ANTENNA_DOWN[2][1] = 9.712 * RADDEG;
    ANTENNA_DOWN[2][2] = 9.892 * RADDEG;
    ANTENNA_DOWN[2][3] = 10.253 * RADDEG;
    ANTENNA_DOWN[2][4] = 10.574 * RADDEG;
    ANTENNA_DOWN[2][5] = 10.62 * RADDEG;
    ANTENNA_DOWN[2][6] = 10.416 * RADDEG;
    ANTENNA_DOWN[2][7] = 10.189 * RADDEG;
    ANTENNA_DOWN[2][8] = 9.776 * RADDEG;
    ANTENNA_DOWN[2][9] = 9.596 * RADDEG;
    ANTENNA_DOWN[2][10] = 9.561 * RADDEG;
    ANTENNA_DOWN[2][11] = 9.695 * RADDEG;
    ANTENNA_DOWN[2][12] = 9.445 * RADDEG;
    ANTENNA_DOWN[2][13] = 9.387 * RADDEG;
    ANTENNA_DOWN[2][14] = 9.398 * RADDEG;
    ANTENNA_DOWN[2][15] = 9.288 * RADDEG;
    for(int iii = 0; iii < 3; iii++) // move from the square centers to the phase centers
      for(int jjj = 0; jjj < NRX_PHI[iii]; jjj++)
	ANTENNA_POSITION_START[0][iii][jjj] = ANTENNA_POSITION_START[0][iii][jjj] - phase_center * Vector(cos(PHI_EACHLAYER[iii][jjj])*sin(90.*RADDEG+ANTENNA_DOWN[iii][jjj]), sin(PHI_EACHLAYER[iii][jjj])*sin(90.*RADDEG+ANTENNA_DOWN[iii][jjj]), cos(90.*RADDEG+ANTENNA_DOWN[iii][jjj]));
  }
  else if (settings1->WHICH==7) {
		
    // Just one layer of antennas around balloon
    // EeVEX
		
    //   settings1->NFOLD=1; //
    //maxthreshold=2.3;
		
    settings1->CYLINDRICALSYMMETRY=1;
		
    NRX_PHI[0]=360;
    NRX_PHI[1]=360;
    NRX_PHI[2]=360;
    NRX_PHI[3]=360;
    NRX_PHI[4]=360;
		
    PHITRIG[0]=360;
    PHITRIG[1]=360;
    PHITRIG[2]=360;
    PHITRIG[3]=360;
    PHITRIG[4]=360;
		
		
    PHI_OFFSET[0]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    PHI_OFFSET[1]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    PHI_OFFSET[2]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    PHI_OFFSET[3]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    PHI_OFFSET[4]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    // it's rotated by 1/2 the azimuth that separates two antennas on the 0th layer
		
    // sets their declination
    THETA_ZENITH[0]=PI/2+5.5*RADDEG;
    THETA_ZENITH[1]=PI/2+7.5*RADDEG;
    THETA_ZENITH[2]=PI/2+9.5*RADDEG;
    THETA_ZENITH[3]=PI/2+11.5*RADDEG;
    THETA_ZENITH[4]=PI/2+13.5*RADDEG;
		
		
    // radius from center axis of the payload
    RRX[0] = 0.9210802;
    RRX[1] = 0.9210802;
    RRX[2] = 0.9210802;
    RRX[3] = 0.9210802;
    RRX[4] = 0.9210802;
		
    // vertical separation between layers.
    LAYER_VPOSITION[0]=0;
    LAYER_VPOSITION[1]=10.;
    LAYER_VPOSITION[2]=20.;
    LAYER_VPOSITION[3]=30.;
    LAYER_VPOSITION[4]=40.;
		
    LAYER_HPOSITION[0]=0.;
    LAYER_HPOSITION[1]=0.;
    LAYER_HPOSITION[2]=0.;
    LAYER_HPOSITION[3]=0.;
    LAYER_HPOSITION[4]=0.;
		
  } //else if (EeVEX)
    //anitaII or satellite
  else if (settings1->WHICH==8) {
    cout<<"initializing and using anitaII payload geometry"<<endl;
		
    // layer 0 is antennas 1-8 on the payload
    // layer 1 is antennas 9-15
    // layer 2 is antennas 16-32
		
    //settings1->NFOLD=8;
    //maxthreshold=2.3;
		
    settings1->CYLINDRICALSYMMETRY=0;
		
    NRX_PHI[0]=8;
    NRX_PHI[1]=8;
    NRX_PHI[2]=16;
    NRX_PHI[3]=8;
		
    PHITRIG[0]=16; // number of positions in phi in each trigger layer
    PHITRIG[1]=16;
    PHITRIG[2]=8;
		
    PHI_OFFSET[0]=0.;
    PHI_OFFSET[1]=0;
    PHI_OFFSET[2]=0;
    PHI_OFFSET[3]=0;
		
    // sets their declination
    THETA_ZENITH[0]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[1]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[2]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[3]=PI/2+INCLINE_TOPTHREE*RADDEG;
		
    ANTENNA_POSITION_START[0][0][0] = MINCH * Vector(40.438,-36.958,147.227);
    ANTENNA_POSITION_START[0][0][1] = MINCH * Vector(57.134,3.109,146.476);
    ANTENNA_POSITION_START[0][0][2] = MINCH * Vector(40.549,43.106,145.871);
    ANTENNA_POSITION_START[0][0][3] = MINCH * Vector(0.624,59.688,145.361);
    ANTENNA_POSITION_START[0][0][4] = MINCH * Vector(-39.455,43.147,145.928);
    ANTENNA_POSITION_START[0][0][5] = MINCH * Vector(-56.096,3.177,146.894);
    ANTENNA_POSITION_START[0][0][6] = MINCH * Vector(-39.355,-36.753,147.757);
    ANTENNA_POSITION_START[0][0][7] = MINCH * Vector(0.645,-53.539,147.876);
    ANTENNA_POSITION_START[0][1][0] = MINCH * Vector(19.554,-43.890,109.531);
    ANTENNA_POSITION_START[0][1][1] = MINCH * Vector(46.600,-16.625,108.889);
    ANTENNA_POSITION_START[0][1][2] = MINCH * Vector(46.587,21.659,108.220);
    ANTENNA_POSITION_START[0][1][3] = MINCH * Vector(19.476,48.539,107.671);
    ANTENNA_POSITION_START[0][1][4] = MINCH * Vector(-18.798,48.502,107.852);
    ANTENNA_POSITION_START[0][1][5] = MINCH * Vector(-45.899,21.424,108.516);
    ANTENNA_POSITION_START[0][1][6] = MINCH * Vector(-45.895,-16.821,109.354);
    ANTENNA_POSITION_START[0][1][7] = MINCH * Vector(-18.691,-43.864,109.843);
    ANTENNA_POSITION_START[0][2][0] = MINCH * Vector(38.636,-93.988,2.636);
    ANTENNA_POSITION_START[0][2][1] = MINCH * Vector(71.690,-72.108,1.953);
    ANTENNA_POSITION_START[0][2][2] = MINCH * Vector(93.897,-39.211,0.498);
    ANTENNA_POSITION_START[0][2][3] = MINCH * Vector(101.790,-0.212,-0.661);
    ANTENNA_POSITION_START[0][2][4] = MINCH * Vector(94.047,38.773,-1.788);
    ANTENNA_POSITION_START[0][2][5] = MINCH * Vector(72.080,71.816,-2.223);
    ANTENNA_POSITION_START[0][2][6] = MINCH * Vector(39.065,93.999,-2.561);
    ANTENNA_POSITION_START[0][2][7] = MINCH * Vector(0.121,101.815,-2.314);
    ANTENNA_POSITION_START[0][2][8] = MINCH * Vector(-38.815,94.002,-2.034);
    ANTENNA_POSITION_START[0][2][9] = MINCH * Vector(-71.809,71.912,-1.102);
    ANTENNA_POSITION_START[0][2][10] = MINCH * Vector(-93.886,39.000,-0.673);
    ANTENNA_POSITION_START[0][2][11] = MINCH * Vector(-101.885,0.048,0.102);
    ANTENNA_POSITION_START[0][2][12] = MINCH * Vector(-94.017,-38.841,0.865);
    ANTENNA_POSITION_START[0][2][13] = MINCH * Vector(-72.079,-71.902,1.864);
    ANTENNA_POSITION_START[0][2][14] = MINCH * Vector(-39.152,-93.935,2.464);
    ANTENNA_POSITION_START[0][2][15] = MINCH * Vector(-0.290,-101.771,2.991);
    ANTENNA_POSITION_START[0][3][0] = MINCH * Vector(32.625,-82.045,-71.140);
    ANTENNA_POSITION_START[0][3][1] = MINCH * Vector(79.071,-35.639,-72.809);
    ANTENNA_POSITION_START[0][3][2] = MINCH * Vector(79.172,30.988,-74.893);
    ANTENNA_POSITION_START[0][3][3] = MINCH * Vector(32.608,77.414,-75.342);
    ANTENNA_POSITION_START[0][3][4] = MINCH * Vector(-33.398,78.088,-74.957);
    ANTENNA_POSITION_START[0][3][5] = MINCH * Vector(-79.367,31.568,-73.922);
    ANTENNA_POSITION_START[0][3][6] = MINCH * Vector(-78.900,-34.192,-72.645);
    ANTENNA_POSITION_START[0][3][7] = MINCH * Vector(-33.046,-81.696,-70.907);
    PHI_EACHLAYER[0][0] = -45.012 * RADDEG ;//ant 7
    PHI_EACHLAYER[0][1] = -0.588 * RADDEG ;//ant 0
    PHI_EACHLAYER[0][2] = 45.694 * RADDEG ;//ant 1
    PHI_EACHLAYER[0][3] = 90.310 * RADDEG ;//ant 2
    PHI_EACHLAYER[0][4] = 135.161 * RADDEG ;//ant3
    PHI_EACHLAYER[0][5] = 179.861 * RADDEG ;//ant4
    PHI_EACHLAYER[0][6] = -134.930 * RADDEG ;//ant5
    PHI_EACHLAYER[0][7] = -90.638 * RADDEG ;//ant 6
    PHI_EACHLAYER[1][0] = -67.412 * RADDEG ;//ant 15
    PHI_EACHLAYER[1][1] = -23.005 * RADDEG ;//ant 8
    PHI_EACHLAYER[1][2] = 22.503 * RADDEG ;//ant 9
    PHI_EACHLAYER[1][3] = 67.722 * RADDEG ;//ant 10
    PHI_EACHLAYER[1][4] = 112.614 * RADDEG ;//ant 11
    PHI_EACHLAYER[1][5] = 157.685 * RADDEG ;//ant 12
    PHI_EACHLAYER[1][6] = -156.639 * RADDEG ;//ant 13
    PHI_EACHLAYER[1][7] = -112.587 * RADDEG ;//ant 14
    PHI_EACHLAYER[2][0] = -67.365 * RADDEG ;//ant 29 
    PHI_EACHLAYER[2][1] = -45.135 * RADDEG ;//ant 30
    PHI_EACHLAYER[2][2] = -23.002 * RADDEG ;//ant 31
    PHI_EACHLAYER[2][3] = -1.013 * RADDEG ;//ant 16
    PHI_EACHLAYER[2][4] = 21.934 * RADDEG ;//ant 17
    PHI_EACHLAYER[2][5] = 44.467 * RADDEG ;//ant 18
    PHI_EACHLAYER[2][6] = 67.288 * RADDEG ;//ant 19
    PHI_EACHLAYER[2][7] = 89.971 * RADDEG ;//ant 20
    PHI_EACHLAYER[2][8] = 112.390 * RADDEG ;//ant 21
    PHI_EACHLAYER[2][9] = 134.988 * RADDEG ;//ant 22
    PHI_EACHLAYER[2][10] = 157.387 * RADDEG ;//ant 23
    PHI_EACHLAYER[2][11] = 179.843 * RADDEG ;//ant 24
    PHI_EACHLAYER[2][12] = -157.444 * RADDEG ;//ant 25
    PHI_EACHLAYER[2][13] = -134.877 * RADDEG ;//ant 26
    PHI_EACHLAYER[2][14] = -112.406 * RADDEG ;//ant 27
    PHI_EACHLAYER[2][15] = -90.081 * RADDEG ;//ant 28
    PHI_EACHLAYER[3][0] = -67.997 * RADDEG ;//ant 
    PHI_EACHLAYER[3][1] = -22.948 * RADDEG ;//ant 
    PHI_EACHLAYER[3][2] = 22.382 * RADDEG ;//ant
    PHI_EACHLAYER[3][3] = 67.583 * RADDEG ;//ant
    PHI_EACHLAYER[3][4] = 112.844 * RADDEG ;//ant
    PHI_EACHLAYER[3][5] = 157.761 * RADDEG ;//ant 
    PHI_EACHLAYER[3][6] = -157.896 * RADDEG ;//ant
    PHI_EACHLAYER[3][7] = -112.791 * RADDEG ;//ant
    ANTENNA_DOWN[0][0] = 9.637 * RADDEG;
    ANTENNA_DOWN[0][1] = 10.108 * RADDEG;
    ANTENNA_DOWN[0][2] = 11.245 * RADDEG;
    ANTENNA_DOWN[0][3] = 11.291 * RADDEG;
    ANTENNA_DOWN[0][4] = 10.988 * RADDEG;
    ANTENNA_DOWN[0][5] = 9.491 * RADDEG;
    ANTENNA_DOWN[0][6] = 9.027 * RADDEG;
    ANTENNA_DOWN[0][7] = 8.743 * RADDEG;
    ANTENNA_DOWN[1][0] = 9.445 * RADDEG;
    ANTENNA_DOWN[1][1] = 10.061 * RADDEG;
    ANTENNA_DOWN[1][2] = 10.772 * RADDEG;
    ANTENNA_DOWN[1][3] = 11.484 * RADDEG;
    ANTENNA_DOWN[1][4] = 11.122 * RADDEG;
    ANTENNA_DOWN[1][5] = 10.376 * RADDEG;
    ANTENNA_DOWN[1][6] = 9.410 * RADDEG;
    ANTENNA_DOWN[1][7] = 9.039 * RADDEG;
    ANTENNA_DOWN[2][0] = 8.233 * RADDEG;
    ANTENNA_DOWN[2][1] = 8.807 * RADDEG;
    ANTENNA_DOWN[2][2] = 9.120 * RADDEG;
    ANTENNA_DOWN[2][3] = 10.352 * RADDEG;
    ANTENNA_DOWN[2][4] = 10.889 * RADDEG;
    ANTENNA_DOWN[2][5] = 11.315 * RADDEG;
    ANTENNA_DOWN[2][6] = 11.402 * RADDEG;
    ANTENNA_DOWN[2][7] = 11.379 * RADDEG;
    ANTENNA_DOWN[2][8] = 10.842 * RADDEG;
    ANTENNA_DOWN[2][9] = 10.725 * RADDEG;
    ANTENNA_DOWN[2][10] = 10.143 * RADDEG;
    ANTENNA_DOWN[2][11] = 10.067 * RADDEG;
    ANTENNA_DOWN[2][12] = 9.503 * RADDEG;
    ANTENNA_DOWN[2][13] = 9.021 * RADDEG;
    ANTENNA_DOWN[2][14] = 8.453 * RADDEG;
    ANTENNA_DOWN[2][15] = 8.268 * RADDEG;
    ANTENNA_DOWN[3][0] = 8.007 * RADDEG;
    ANTENNA_DOWN[3][1] = 9.817 * RADDEG;
    ANTENNA_DOWN[3][2] = 10.259 * RADDEG;
    ANTENNA_DOWN[3][3] = 11.648 * RADDEG;
    ANTENNA_DOWN[3][4] = 10.271 * RADDEG;
    ANTENNA_DOWN[3][5] = 10.015 * RADDEG;
    ANTENNA_DOWN[3][6] = 10.889 * RADDEG;
    ANTENNA_DOWN[3][7] = 7.314 * RADDEG;

    SIMON_DELTA_R[0][0] = -0.0384839;
    SIMON_DELTA_R[0][1] = 0.00634697;
    SIMON_DELTA_R[0][2] = -0.0861167;
    SIMON_DELTA_R[0][3] = 0.0461873;
    SIMON_DELTA_R[0][4] = 0.0153388;
    SIMON_DELTA_R[0][5] = -0.00927728;
    SIMON_DELTA_R[0][6] = 0.0239867;
    SIMON_DELTA_R[0][7] = 0.0125282;
    SIMON_DELTA_R[1][0] = -0.0111636;
    SIMON_DELTA_R[1][1] = -0.0959452;
    SIMON_DELTA_R[1][2] = -0.0330808;
    SIMON_DELTA_R[1][3] = -0.0475617;
    SIMON_DELTA_R[1][4] = 0.0196292;
    SIMON_DELTA_R[1][5] = -0.0190837;
    SIMON_DELTA_R[1][6] = -0.00922367;
    SIMON_DELTA_R[1][7] = -0.0294811;
    SIMON_DELTA_R[2][0] = 0.0140245;
    SIMON_DELTA_R[2][1] = -0.0621836;
    SIMON_DELTA_R[2][2] = -0.0379325;
    SIMON_DELTA_R[2][3] = -0.0108062;
    SIMON_DELTA_R[2][4] = -0.0601935;
    SIMON_DELTA_R[2][5] = -0.0968276;
    SIMON_DELTA_R[2][6] = -0.0348523;
    SIMON_DELTA_R[2][7] = 0.0121726;
    SIMON_DELTA_R[2][8] = 0.0405193;
    SIMON_DELTA_R[2][9] = 0.0239992;
    SIMON_DELTA_R[2][10] = -0.0405203;
    SIMON_DELTA_R[2][11] = -0.00401756;
    SIMON_DELTA_R[2][12] = -0.0362955;
    SIMON_DELTA_R[2][13] = -0.00587152;
    SIMON_DELTA_R[2][14] = -0.00611182;
    SIMON_DELTA_R[2][15] = -0.00321244;
    SIMON_DELTA_R[3][0] = -0.0437687;
    SIMON_DELTA_R[3][1] = -0.0643475;
    SIMON_DELTA_R[3][2] = -0.0804245;
    SIMON_DELTA_R[3][3] = -0.0112675;
    SIMON_DELTA_R[3][4] = 0.0337428;
    SIMON_DELTA_R[3][5] = -0.0525977;
    SIMON_DELTA_R[3][6] = -0.101587;
    SIMON_DELTA_R[3][7] = -0.0401037;

    SIMON_DELTA_PHI[0][0] = -0.0100608;
    SIMON_DELTA_PHI[0][1] = -0.00313443;
    SIMON_DELTA_PHI[0][2] = -0.015312;
    SIMON_DELTA_PHI[0][3] = 0.00206827;
    SIMON_DELTA_PHI[0][4] = -0.0227948;
    SIMON_DELTA_PHI[0][5] = 0.00750385;
    SIMON_DELTA_PHI[0][6] = 0.00388065;
    SIMON_DELTA_PHI[0][7] = -0.00131021;
    SIMON_DELTA_PHI[1][0] = -0.0299233;
    SIMON_DELTA_PHI[1][1] = -0.00165365;
    SIMON_DELTA_PHI[1][2] = -0.0107407;
    SIMON_DELTA_PHI[1][3] = 0.0145914;
    SIMON_DELTA_PHI[1][4] = -0.0150373;
    SIMON_DELTA_PHI[1][5] = -0.0121967;
    SIMON_DELTA_PHI[1][6] = -0.0038106;
    SIMON_DELTA_PHI[1][7] = 0.0106842;
    SIMON_DELTA_PHI[2][0] = -0.0087849;
    SIMON_DELTA_PHI[2][1] = 0.000682206;
    SIMON_DELTA_PHI[2][2] = -0.00516052;
    SIMON_DELTA_PHI[2][3] = -0.00770935;
    SIMON_DELTA_PHI[2][4] = -0.00862535;
    SIMON_DELTA_PHI[2][5] = -0.00920648;
    SIMON_DELTA_PHI[2][6] = 0.00037431;
    SIMON_DELTA_PHI[2][7] = 0.00310935;
    SIMON_DELTA_PHI[2][8] = -0.00546085;
    SIMON_DELTA_PHI[2][9] = -0.00901249;
    SIMON_DELTA_PHI[2][10] = -0.0145529;
    SIMON_DELTA_PHI[2][11] = -0.00666063;
    SIMON_DELTA_PHI[2][12] = -0.00372999;
    SIMON_DELTA_PHI[2][13] = 0.00197442;
    SIMON_DELTA_PHI[2][14] = -0.000789595;
    SIMON_DELTA_PHI[2][15] = 0.000188257;
    SIMON_DELTA_PHI[3][0] = -0.00289577;
    SIMON_DELTA_PHI[3][1] = -0.0203117;
    SIMON_DELTA_PHI[3][2] = -0.00503387;
    SIMON_DELTA_PHI[3][3] = -0.000220575;
    SIMON_DELTA_PHI[3][4] = -0.00416114;
    SIMON_DELTA_PHI[3][5] = -0.0223176;
    SIMON_DELTA_PHI[3][6] = 0.0058874;
    SIMON_DELTA_PHI[3][7] = 0.00899651;
		
    for(int iii = 0; iii < 4; iii++){ // move from the square centers to the phase centers
      for(int jjj = 0; jjj < NRX_PHI[iii]; jjj++){
			 
	//ANTENNA_DOWN is measured from horiztonal. Put negatives in correct places. Verified with analysis code 
	ANTENNA_POSITION_START[0][iii][jjj] = ANTENNA_POSITION_START[0][iii][jjj] - phase_center_anita2_analysis * Vector(cos(PHI_EACHLAYER[iii][jjj])*cos(-1*ANTENNA_DOWN[iii][jjj]), sin(PHI_EACHLAYER[iii][jjj])*cos(-1*ANTENNA_DOWN[iii][jjj]), sin(-1*ANTENNA_DOWN[iii][jjj]));
      }//jjj
    }//iii
    double r;
    double phi;
		
    double x;
    double y;
    double z;

	
    for(int iii = 0; iii < 4; iii++){ // move from the square centers to the phase centers
      for(int jjj = 0; jjj < NRX_PHI[iii]; jjj++){
	x = ANTENNA_POSITION_START[0][iii][jjj][0];
	y = ANTENNA_POSITION_START[0][iii][jjj][1];
	z = ANTENNA_POSITION_START[0][iii][jjj][2];

	r = sqrt(pow(x,2)+pow(y,2));
	phi = atan2(y,x);
		   
	ANTENNA_POSITION_START[0][iii][jjj]= Vector((r+SIMON_DELTA_R[iii][jjj])*cos(phi+SIMON_DELTA_PHI[iii][jjj]),(r+SIMON_DELTA_R[iii][jjj])*sin(phi+SIMON_DELTA_PHI[iii][jjj]),z);
		   
	ANTENNA_POSITION_START[1][iii][jjj]=ANTENNA_POSITION_START[0][iii][jjj]=ANTENNA_POSITION_START[0][iii][jjj].RotateZ(-gps_offset_anita2);
	PHI_EACHLAYER[iii][jjj]=atan2(ANTENNA_POSITION_START[0][iii][jjj][1],ANTENNA_POSITION_START[0][iii][jjj][0]);//set phi of each antennas to correct starting position

	//cout<<"Antenna pos is "<<ANTENNA_POSITION_START[0][iii][jjj]<<" PHI is "<<PHI_EACHLAYER[iii][jjj]<<"\n";
      }
    }

		
  } 
  else if (settings1->WHICH==9 || settings1->WHICH==10) { // ANITA-3 and ANITA-4
    cout<<"initializing and using ANITA-III payload geometry"<<endl;
    // layer 0 is antennas 1-8 on the payload
    // layer 1 is antennas 9-15
    // layer 2 is antennas 16-32
    // layer 3 is antennas 32-48
        
      
    settings1->CYLINDRICALSYMMETRY=0;
      
    //these are physical layers
    NRX_PHI[0]=8;
    NRX_PHI[1]=8;
    NRX_PHI[2]=16;
    NRX_PHI[3]=16;
      
    PHITRIG[0]=16; // number of positions in phi in each *trigger* layer
    PHITRIG[1]=16;
    PHITRIG[2]=16;
      
    //these are physical layers again
    PHI_OFFSET[0]=0.; 
    PHI_OFFSET[1]=0.; // 2.*PI/(double)NRX_PHI[0]/2.; // Linda: changed this offset to 0 as  it shouldn't be needed
    PHI_OFFSET[2]=0.;
    PHI_OFFSET[3]=0.;
      
    //double INCLINE_NADIR=55; // this is set in the input file now So should be removed
      
    // sets their declination
    THETA_ZENITH[0]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[1]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[2]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[3]=PI/2+INCLINE_TOPTHREE*RADDEG;
      
    // Read photogrammetry positions
    string whichANITAroman="";
    if (settings1->WHICH==9) whichANITAroman+="III";
    else whichANITAroman+="IV";
    string photoFile;
#ifdef ANITA_UTIL_EXISTS
    photoFile += ( (string)getenv("ANITA_UTIL_INSTALL_DIR") +"/share/anitaCalib/anita"+whichANITAroman+"Photogrammetry.csv");
#else
    photoFile += (ICEMC_DATA_DIR+"/anita"+whichANITAroman+"Photogrammetry.csv");
#endif
    
    std::ifstream Anita3PhotoFile(photoFile.c_str());
      if (!Anita3PhotoFile){
      std::cerr << "Couldn't open photogrammetry!" << std::endl;
      return;
    }

    //First up are the antenna positions
    TString line;
    for(int i=0;i<2;i++) {
      line.ReadLine(Anita3PhotoFile);
      //std::cout << line.Data() << "\n";
    }

    //Array with photogrammetry values
    Double_t xAntPhoto[48]; //inch
    Double_t yAntPhoto[48]; //inch
    Double_t zAntPhoto[48]; //inch
    Double_t rAntPhoto[48]; //inch
    Double_t azCentrePhoto[48]; //deg
    Double_t apertureAzPhoto[48]; //deg
    Double_t apertureElPhoto[48]; //deg  

    for(int ant=0;ant<48;ant++) {
      line.ReadLine(Anita3PhotoFile);
      //std::cout << "Seavey:\t" << line.Data() << "\n";
      TObjArray *tokens = line.Tokenize(",");
      for(int j=0;j<8;j++) {
	const TString subString = ((TObjString*)tokens->At(j))->GetString();
	//	TString *subString = (TString*) tokens->At(j);
	//	std::cout << j << "\t" << subString.Data() << "\n";
	switch(j) {
	case 0:
	  if (ant+1 != subString.Atoi()) {
	    std::cerr << "Antenna number mismatch\n";
	  }
	  break;
	case 1:	   
	  xAntPhoto[ant]=subString.Atof(); //inch
	  break;
	case 2:	   
	  yAntPhoto[ant]=subString.Atof(); //inch
	  break;
	case 3:	   
	  zAntPhoto[ant]=subString.Atof(); //inch
	  break;
	case 4:	   
	  rAntPhoto[ant]=subString.Atof(); //inch
	  break;
	case 5:	   
	  azCentrePhoto[ant]=subString.Atof(); //deg
	  break;
	case 6:	   
	  apertureAzPhoto[ant]=subString.Atof(); //deg
	  break;
	case 7:	   
	  apertureElPhoto[ant]=subString.Atof()*(-1); //deg // photogrammetry elevation defined as negative, here positive
	  break;
	default:	   
	  break;
	}
	  
      }
      tokens->Delete();
	
    }
    Anita3PhotoFile.close();

    // Fill photogrammetry position for top rings
    for (int iant=0; iant<8;iant++){
      ANTENNA_POSITION_START[0][0][iant] = MINCH * Vector(xAntPhoto[iant*2], yAntPhoto[iant*2], zAntPhoto[iant*2]).RotateZ(-gps_offset_anita3);	    // top ring top antennas
      ANTENNA_POSITION_START[0][1][iant] = MINCH * Vector(xAntPhoto[iant*2+1], yAntPhoto[iant*2+1], zAntPhoto[iant*2+1]).RotateZ(-gps_offset_anita3);  // top ring bottom antennas

      PHI_EACHLAYER[0][iant] = azCentrePhoto[iant*2] * RADDEG - gps_offset_anita3;
      PHI_EACHLAYER[1][iant] = azCentrePhoto[iant*2+1] * RADDEG - gps_offset_anita3;
      ANTENNA_DOWN[0][iant] = apertureElPhoto[iant*2] * RADDEG; 
      ANTENNA_DOWN[1][iant] = apertureElPhoto[iant*2+1] * RADDEG; 
 
    }

    // Fill photogrammetry position for middle and bottom rings
    for (int iant=0; iant<16;iant++){
      ANTENNA_POSITION_START[0][2][iant] = MINCH * Vector(xAntPhoto[iant+16], yAntPhoto[iant+16], zAntPhoto[iant+16]).RotateZ(-gps_offset_anita3);	    // middle ring antennas
      ANTENNA_POSITION_START[0][3][iant] = MINCH * Vector(xAntPhoto[iant+32], yAntPhoto[iant+32], zAntPhoto[iant+32]).RotateZ(-gps_offset_anita3);  // bottom ring antennas

      PHI_EACHLAYER[2][iant] = azCentrePhoto[iant+16] * RADDEG - gps_offset_anita3;
      ANTENNA_DOWN[2][iant] = apertureElPhoto[iant+16] * RADDEG; 

      PHI_EACHLAYER[3][iant] = azCentrePhoto[iant+32] * RADDEG - gps_offset_anita3;
      ANTENNA_DOWN[3][iant] = apertureElPhoto[iant+32] * RADDEG; 

    }
      
    // HERE HPOL IS 0 AND VPOL IS 1
    string whichANITAcard="";
    if (settings1->WHICH==9) whichANITAcard+="3";
    else whichANITAcard+="4";
    string phaseCenterName;
#ifdef ANITA_UTIL_EXISTS
    phaseCenterName += ( (string)getenv("ANITA_UTIL_INSTALL_DIR") +"/share/anitaCalib/phaseCenterPositionsRelativeToPhotogrammetryAnita"+whichANITAcard+".dat");
#else
    phaseCenterName += (ICEMC_DATA_DIR+"/phaseCenterPositionsRelativeToPhotogrammetryAnita"+whichANITAcard+".dat");
#endif
    
    std::ifstream PhaseCenterFile(phaseCenterName.c_str());
    Int_t antNum, tpol, pol;
    Double_t deltaR,deltaPhi,deltaZ;
    char firstLine[180];
    Double_t deltaRPhaseCentre[2][4][16]; //Relative to photogrammetry + ring offset
    Double_t deltaZPhaseCentre[2][4][16]; //Relative to photogrammetry + ring offset
    Double_t deltaPhiPhaseCentre[2][4][16]; //Relative to photogrammetry + ring offset
      
    PhaseCenterFile.getline(firstLine,179);
    // HERE HPOL IS 0 AND VPOL IS 1 that's why we invert pol here
    while(PhaseCenterFile >> antNum >> tpol >> deltaR >> deltaPhi >> deltaZ) {
      int ilayer = (antNum<16)*((antNum%2==0)*0 + (antNum%2==1)*1)+ (antNum>15)*(antNum<32)*2+(antNum>31)*3;
      int ifold = (ilayer<2)*((antNum-ilayer)/2)+(ilayer>1)*(antNum%16);

      if (tpol==1) pol=0;
      else if (tpol==0) pol=1;
	
      deltaRPhaseCentre[pol][ilayer][ifold]=deltaR;
      deltaPhiPhaseCentre[pol][ilayer][ifold]=deltaPhi*TMath::DegToRad();
      deltaZPhaseCentre[pol][ilayer][ifold]=deltaZ;
    } 
    PhaseCenterFile.close();

    std::ifstream relativePhaseCenterToAmpaDelaysFile((ICEMC_DATA_DIR+"/relativePhaseCenterToAmpaDelaysAnita"+whichANITAcard+".dat").c_str());

#ifdef ANITA_UTIL_EXISTS
    AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
    AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance(stoi(whichANITAcard));
    int tempAnt, intTempPol;
    AnitaPol::AnitaPol_t tempPol;
    for(Int_t surf=0; surf<NUM_SURF; surf++){
      for(Int_t chan=0; chan<NUM_CHAN; chan++){
    	fGeomTool->getAntPolFromSurfChan(surf,chan, tempAnt, tempPol);
    	if (tempAnt!=-1){
	  // in EventReaderRoot 0: HPOL, 1: VPOL, in icemc it's the opposite
	  intTempPol = (tempPol==0) ? 1 : 0;
	  extraCableDelays[intTempPol][tempAnt] = cal->relativePhaseCenterToAmpaDelays[surf][chan]*1e-9 ; 
	}
      }
    }

#else
    for(int ipol=0; ipol<2; ipol++){
      for (int iant=0; iant<48; iant++){
	extraCableDelays[ipol][iant] = 0;
      }
    }
    
#endif

    
    double x, y, z, r, phi;
    for (int ipol = 0; ipol < 2; ipol++){
      for(int ilayer = 0; ilayer < 4; ilayer++){ 
	for(int ifold = 0; ifold < NRX_PHI[ilayer]; ifold++){

	  // First attempt to define phase-centers is done by coming 20 cm
	  // inwards from the antenna face end (that was only defined for VPOL)
	  // that's why the initial ipol=0
	  ANTENNA_POSITION_START[ipol][ilayer][ifold] = ANTENNA_POSITION_START[0][ilayer][ifold] - phase_center_anita3 * Vector(cos(PHI_EACHLAYER[ilayer][ifold])*cos(ANTENNA_DOWN[ilayer][ifold]), sin(PHI_EACHLAYER[ilayer][ifold])*cos(ANTENNA_DOWN[ilayer][ifold]), sin(ANTENNA_DOWN[ilayer][ifold]));
	  x = ANTENNA_POSITION_START[ipol][ilayer][ifold].GetX();
	  y = ANTENNA_POSITION_START[ipol][ilayer][ifold].GetY();
	  
	  r = sqrt(pow(x,2)+pow(y,2)) + deltaRPhaseCentre[ipol][ilayer][ifold];
	  phi = atan2(y,x) + deltaPhiPhaseCentre[ipol][ilayer][ifold];
	  z = ANTENNA_POSITION_START[ipol][ilayer][ifold].GetZ() + deltaZPhaseCentre[ipol][ilayer][ifold];	  

	  if(phi<0) phi+=TMath::TwoPi();
	  if(phi>TMath::TwoPi()) phi-=TMath::TwoPi();
	  
	  ANTENNA_POSITION_START[ipol][ilayer][ifold]= Vector(r*cos(phi),r*sin(phi),z);

	  PHI_EACHLAYER[ilayer][ifold]=phi;
	  
	}
      }
    }
  }
    
    
  else if (settings1->WHICH==11) { // satellite
      
		
    // layer 0 is antennas 1-8 on the payload
    // layer 1 is antennas 9-15
    settings1->CYLINDRICALSYMMETRY=1;
		
    //these are physical layers
    NRX_PHI[0]=8;
    NRX_PHI[1]=8;
		
    PHITRIG[0]=8; // number of positions in phi in each *trigger* layer
    PHITRIG[1]=8;
		
    //these are physical layers again
    PHI_OFFSET[0]=0.; // antenna 1 on 0th layer is rotated in phi wrt antenna 9 and antenna 17
    // it's rotated by 1/2 the azimuth that separates two antennas on the 0th layer
    PHI_OFFSET[1]=-2.*PI/(double)NRX_PHI[0]/2.;
		
    // sets their declination
    THETA_ZENITH[0]=PI/2+INCLINE_TOPTHREE*RADDEG;
    THETA_ZENITH[1]=PI/2+INCLINE_TOPTHREE*RADDEG;
		
    // radius from center axis of the payload
    RRX[0] = 0.9210802;
    RRX[1] = 0.7553198;   
		
    // vertical separation between layers.
    LAYER_VPOSITION[0]=0;
    LAYER_VPOSITION[1] = -7.5;
		
    LAYER_HPOSITION[0]=0.;
    LAYER_HPOSITION[1] = 0.;
		
		
  } //else if (satellite)
    
  settings1->NANTENNAS=0;
  for (int i=0;i<settings1->NLAYERS;i++)
    settings1->NANTENNAS+=NRX_PHI[i];
    
  cout << "nantennas is " << settings1->NANTENNAS << "\n";
  number_all_antennas=settings1->NANTENNAS;
    
  // gets noise (vrms) for each bandwidth slice and antenna layer according to antenna theta
    
  if (settings1->WHICH==0) {
    temp_eachrx[0]=919.4;
    temp_eachrx[1]=1051.8;
    for (int j=2;j<16;j++) {
      temp_eachrx[j]=0.; // noise temperature of each antenna?
    }
  }
  else {
    for (int j=0;j<16;j++) {
	
      temp_eachrx[j]=(240.+200.); // temp of each antenna in kelvin
    }
  }
    
  for (int i=0;i<NRX_PHI[0];i++) {
    VNOISE_ANITALITE[i]=ChanTrigger::GetNoise(settings1,bn1->altitude_bn,bn1->surface_under_balloon,THETA_ZENITH[0],settings1->BW_SEAVEYS,temp_eachrx[i]);
  }
}//GetPayload

void Anita::calculate_all_offsets(void) {
    
  double angle_phi, angle_theta;
    
  double hypothesis_offset[3][3];
    
  vector<double> angles_tmp;
  vector< vector <double> > angles_diffthetas_tmp;

  double step_phi=(MAX_PHI_HYPOTHESIS-MIN_PHI_HYPOTHESIS)/(double)N_STEPS_PHI;
  double step_theta=(MAX_THETA_HYPOTHESIS-MIN_THETA_HYPOTHESIS)/(double)N_STEPS_THETA;

  cout << "step_theta is " << step_theta << "\n";

  for (unsigned center_phi_sector_index = 0; center_phi_sector_index < 1; ++center_phi_sector_index) {
    angles_tmp.clear();
    for (unsigned index_phi = 0; index_phi < N_STEPS_PHI; ++index_phi) {
      //angle_phi = (center_phi_sector_index + 3)%16 * 22.5 - 22.5 + index_phi * 3.;
      // angle of this phi sector (=index_phi) rel to center phi sector 
      angle_phi = (double)center_phi_sector_index * 22.5 + MIN_PHI_HYPOTHESIS + (double)index_phi * step_phi;
      // Note that the zeroth phi sector is actually not the one with the phi = 0 direction, it's actually the third.
      // The above statement is probably false because of the gps_offset_anita2 rotation! [0][0] should be right after all...
      angles_diffthetas_tmp.clear();
      for (unsigned index_theta = 0; index_theta < N_STEPS_THETA; ++index_theta) {
	angle_theta = MIN_THETA_HYPOTHESIS + (double)index_theta*step_theta;

	//	  if (angle_phi==-11. && angle_theta==-20.)
	//if (fabs(angle_phi)<1. && angle_theta>15 && angle_theta<17)
	//cout << "angle phi, angle_theta are " << angle_phi << "\t" << angle_theta << "\n";
	angles_tmp.clear();
	angles_tmp.push_back(angle_phi);
	angles_tmp.push_back(angle_theta);
	  
	angles_diffthetas_tmp.push_back(angles_tmp);
	  
	calculate_single_offset(center_phi_sector_index, angle_phi, angle_theta, hypothesis_offset);
	//cout << "index_phi is " << index_phi << " " << "index_theta is " << index_theta << " ";
	for (unsigned i_layer = 0; i_layer < N_SUMMED_LAYERS; ++i_layer) {
	  for (unsigned i_sector = 0; i_sector < N_SUMMED_PHI_SECTORS; ++i_sector) {
	    hypothesis_offsets[center_phi_sector_index][index_phi][index_theta][i_sector][i_layer] = int(Tools::round(hypothesis_offset[i_sector][i_layer] / TRIG_TIMESTEP));
	      
	  }
	}
	  
      }
      hypothesis_angles.push_back(angles_diffthetas_tmp);
    }


  }
  getDifferentOffsets();
  printDifferentOffsets();
  return;
}
void Anita::printDifferentOffsets() {
  ofstream ofile("outputs/offsets.txt");
  ofile << "number of offsets is " << vdifferent_offsets.size() << "\n";
  

  for (unsigned int i=0;i<vdifferent_offsets.size();i++) {
    for (int j=0;j<2;j++) {
      ofile << vdifferent_angles[i][j] << "\t";
    }
    for (unsigned int j=0;j<N_SUMMED_PHI_SECTORS;j++) {
      for (unsigned int k=0;k<N_SUMMED_LAYERS;k++) {
	//    for (int j=0;j<vdifferent_offsets[i].size();j++) {
	ofile << vdifferent_offsets[i][N_SUMMED_LAYERS*j+k] << " ";
	//}
      }
      ofile << "\t";
    }
    ofile << "\n";
  }
  ofile.close();
}
void Anita::getDifferentOffsets() {
 
  vector<int> vtmp;
  vector<double> vangles_tmp;

  for (int center_phi_sector_index=0;center_phi_sector_index<1;center_phi_sector_index++) {
    for (unsigned index_phi = 0; index_phi < N_STEPS_PHI; ++index_phi) {
      for (unsigned index_theta = 0; index_theta < N_STEPS_THETA; ++index_theta) {
	vtmp.clear();
	vangles_tmp.clear();
	vangles_tmp.push_back(hypothesis_angles[index_phi][index_theta][0]);
	vangles_tmp.push_back(hypothesis_angles[index_phi][index_theta][1]);
	
	for (unsigned i_sector = 0; i_sector < N_SUMMED_PHI_SECTORS; ++i_sector) {	
	  for (unsigned i_layer = 0; i_layer < N_SUMMED_LAYERS; ++i_layer) {	    
	    vtmp.push_back(hypothesis_offsets[center_phi_sector_index][index_phi][index_theta][i_sector][i_layer]);
	    //if (hypothesis_angles[index_phi][index_theta][0]<-21 && hypothesis_angles[index_phi][index_theta][1] >-29 && hypothesis_angles[index_phi][index_theta][1]<-27)
	    //cout << vtmp[vtmp.size()-1] << " ";
	  }  // end loop over layer
	} // end loop over sector

	
	// now loop over all previous offsets and see if it is new or a repeat
	int foundone=0;


	for (unsigned int i=0;i<vdifferent_offsets.size();i++) {
	  if (vtmp==vdifferent_offsets[i]) {
	    foundone++;
	    //if (hypothesis_angles[index_phi][index_theta][0]<-21 && hypothesis_angles[index_phi][index_theta][1] >-29 && hypothesis_angles[index_phi][index_theta][1]<-27)
	    //cout << "found one is " << foundone << "\n";
	  }
	}
	//if (hypothesis_angles[index_phi][index_theta][0]<-21 && hypothesis_angles[index_phi][index_theta][1] >-29 && hypothesis_angles[index_phi][index_theta][1]<-27)
	//cout << "just before if statement, found one is " << foundone << "\n";

	if (foundone==0) {
	  vdifferent_offsets.push_back(vtmp);
	  vdifferent_angles.push_back(vangles_tmp);
	}
	//if (hypothesis_angles[index_phi][index_theta][0]<-21 && hypothesis_angles[index_phi][index_theta][1] >-29 && hypothesis_angles[index_phi][index_theta][1]<-27)
	//cout << "size of vdifferent_offsets is " << vdifferent_offsets.size() << "\n";
      } // end loop over hypotheses in theta
    }  // end loop over hypotheses in phi
  } // end loop over center_phi_sector
}
void Anita::calculate_single_offset(const unsigned center_phi_sector_index, const double angle_phi, const double angle_theta, double hypothesis_offset[][3]) {
  double maximum_time = -2000E-9;
   
  double to_center_of_summed_phi_sectors=((double)N_SUMMED_PHI_SECTORS/2.)*22.5-11.25;
  //    cout << "to_center_of_summed_phi_sectors is " << to_center_of_summed_phi_sectors << "\n";
  Vector normal_vector = Vector(cos(angle_theta * RADDEG) * cos((angle_phi+to_center_of_summed_phi_sectors) * RADDEG), cos(angle_theta * RADDEG) * sin((angle_phi+to_center_of_summed_phi_sectors) * RADDEG), sin(angle_theta * RADDEG));
    
  //    cout << "normal vector is ";
  //normal_vector.Print();

  //    Vector first_antenna_pos = ANTENNA_POSITION_START[0][3][center_phi_sector_index];
  Vector one_antenna_pos = antenna_positions[0][2*16+center_phi_sector_index];
  //cout << "one_antenna_pos is ";
  //one_antenna_pos.Print();

  int phi_start=-1*(int)((double)N_SUMMED_PHI_SECTORS/2.)+1;
    
  //    cout << "new set of 6.\n";
  //cout << "phi_start, phi_end are " << phi_start << "\t" << phi_start+N_SUMMED_PHI_SECTORS << "\n";
  //for (signed int phi_sector_offset = phi_start; phi_sector_offset < phi_start+N_SUMMED_PHI_SECTORS; phi_sector_offset++) {
  for (signed int phi_sector_offset = phi_start; phi_sector_offset < (int)(phi_start+N_SUMMED_PHI_SECTORS); phi_sector_offset++) {
    unsigned i_sector = (phi_sector_offset + center_phi_sector_index + 16)%16; // Must map to {0, ..., 15}.
    //cout << "i_sector is " << i_sector << "\n";
    for (unsigned i_layer = 0; i_layer < N_SUMMED_LAYERS; ++i_layer) {
	

      Vector antenna_pos = antenna_positions[0][16*i_layer+i_sector];
      //cout << "i_layer, i_sector are " << i_layer << "\t" << i_sector << "\n";
      //cout << "antenna_pos is ";
      //antenna_pos.Print();
	
      double offset = (-1. / CLIGHT) * normal_vector * (antenna_pos - one_antenna_pos);

      //	cout << "offset is " << offset << "\n";
      if (offset >= maximum_time) {
	maximum_time = offset;
      }
	
      //unsigned i_layer_trig = ((i_layer > 1) ? (i_layer - 1) : (0));	// Set i_layer_triger according to the following map:
      // 0 -> 0 (if in even phi_sector)
      // 1 -> 0 (if in odd phi_sector)
      // 2 -> 1
      // 3 -> 2
      //	cout << "phi_sector_offset-phi_start is " << phi_sector_offset-phi_start << "\n";
      hypothesis_offset[phi_sector_offset-phi_start][i_layer] = offset; // Use phi_sector_offset + 1 so that it maps {-1, 0, 1} to {0, 1, 2}.
    }
  }
  //    cout << "maximum_time is " << maximum_time << "\n";
  for (unsigned i_layer = 0; i_layer < N_SUMMED_LAYERS; ++i_layer) {
    for (unsigned i_sector = 0; i_sector < N_SUMMED_PHI_SECTORS; ++i_sector) {
      hypothesis_offset[i_sector][i_layer] -= maximum_time;
      hypothesis_offset[i_sector][i_layer]*=-1.;
    }
  }
  return;
}


void Anita::GetArrivalTimes(const Vector& rf_direction,Balloon *bn1, Settings *settings1) {
  //cout << "inside getarrivaltimes.\n";
  

  for (int ipol=0; ipol<2; ipol++){
  
    for (int antenna_index = 0; antenna_index < (number_all_antennas); antenna_index++) { //loop over layers on the payload
      arrival_times[ipol][antenna_index] = (antenna_positions[ipol][antenna_index] * rf_direction) / CLIGHT;

      //arrival_times[antenna_index] += extraCableDelays[0][antenna_index];
    
      // cout << "index is " << antenna_index << "\n";
      // cout << "antenna_positions are " << antenna_positions[antenna_index] << "\n";
      // cout << "rf direction is " << rf_direction << "\n";
      // cout << "arrival_times is " << arrival_times[antenna_index] << "\n";
    } // for: loop over antenna layers
    
    if( settings1->WHICH==8 ){//ANITA-II offboresight delay
      
      Vector rf_tmp_dir = bn1->unRotatePayload(-1*rf_direction);
     
      double theta_deg =rf_tmp_dir.Theta() * DEGRAD;//
     
      double phi_deg = rf_tmp_dir.Phi() *DEGRAD;
      double totalAngledeg;
      double extra_delay;
     
      double phi_eachlayer;
      double theta_offset;
      int ant_ctr=0;

      theta_offset = 10;//boresight_vector[ant_ctr].Theta()*DEGRAD;
       
      theta_deg = theta_deg -90;
       
     
      theta_deg = theta_deg - theta_offset;
      for(int iii = 0; iii < 4; iii++){ 
	for(int jjj = 0; jjj < NRX_PHI[iii]; jjj++){
		   
	  phi_deg = rf_tmp_dir.Phi();
      
	  phi_eachlayer =atan2(ANTENNA_POSITION_START[ipol][iii][jjj][1],ANTENNA_POSITION_START[ipol][iii][jjj][0]);
       
	  phi_deg =phi_deg- phi_eachlayer;
	 
	  if(fabs(phi_deg) > fabs(phi_deg+2*PI)) phi_deg+=2*PI;
	  if(fabs(phi_deg) > fabs(phi_deg-2*PI)) phi_deg-=2*PI;
	  phi_deg =phi_deg*DEGRAD;
	  totalAngledeg = phi_deg*phi_deg + theta_deg*theta_deg;
	  if(totalAngledeg > 2500) totalAngledeg=2500;
	 
	  extra_delay  = (totalAngledeg*totalAngledeg)*1.45676e-8;//pulled from Abby analysis
	  extra_delay -= (totalAngledeg)*5.01452e-6;//pulled from Abby analysis
	
	  arrival_times[ipol][ant_ctr]+=extra_delay*1E-9;
	  ant_ctr++;
	}
      }
    }

  }
  
  //    double last_trigger_time=Tools::dMax(arrival_times,(number_all_antennas));
  //cout << "last_trigger_time is " << last_trigger_time << "\n";
  double minV = Tools::dMin(arrival_times[0],(number_all_antennas));
  double minH = Tools::dMin(arrival_times[1],(number_all_antennas));
  double first_trigger_time = Tools::dMin(minV, minH);
  for (int ipol=0; ipol<2; ipol++){
    for (int i=0;i<(number_all_antennas);i++){
      // cout << "antenna_positions is ";
      //  antenna_positions[i].Print();
      // cout << "diff is ";
      // (antenna_positions[i]-one_antenna_position).Print();

      arrival_times[ipol][i] -= first_trigger_time;
      // cout << "arrivaltimes is " << arrival_times[i] << "\n";
      // arrival_times[i] -= last_trigger_time;

      // if (arrival_times[i] == 0){
      // 	 first_phi_sector_hit = (int)((double)i/16.);
      //        }
    }
  }
  // cout << "end of GetArrivalTimes.\n";
} // GetArrivalTimes


void Anita::GetArrivalTimesBoresights(const Vector rf_direction[NLAYERS_MAX][NPHI_MAX],Balloon *bn1, Settings *settings1) {

  for (int ipol=0; ipol<2; ipol++){
    for (int antenna_index = 0; antenna_index < (number_all_antennas); antenna_index++) { //loop over layers on the payload
      int ilayer = (antenna_index<8)*0 + (antenna_index>7)*(antenna_index<16)*1+ (antenna_index>15)*(antenna_index<32)*2+(antenna_index>31)*3;
      int ifold = (ilayer<2)*(antenna_index%8)+(ilayer>1)*(antenna_index%16);
      arrival_times[ipol][antenna_index] = (antenna_positions[ipol][antenna_index] * rf_direction[ilayer][ifold]) / CLIGHT;

      // cout << antenna_index << " " << arrival_times[antenna_index] << " " << extraCableDelays[0][antenna_index] << " " ;
      //      arrival_times[ipol][antenna_index] += extraCableDelays[ipol][antenna_index];
      // cout << arrival_times[antenna_index] << endl;
    
      //  arrival_times[antenna_index]=0;
      

      if(settings1->WHICH==8 ){//ANITA-II offboresight delay
      
	Vector rf_tmp_dir = bn1->unRotatePayload(-1*rf_direction[ilayer][ifold]);
     
	double theta_deg =rf_tmp_dir.Theta() * DEGRAD;//
     
	double phi_deg = rf_tmp_dir.Phi() *DEGRAD;
	double totalAngledeg;
	double extra_delay;
     
	double phi_eachlayer;
	double theta_offset;
	int ant_ctr=0;

	theta_offset = 10;//boresight_vector[ant_ctr].Theta()*DEGRAD;
       
	theta_deg = theta_deg -90;
       
     
	theta_deg = theta_deg - theta_offset;
		   
	phi_deg = rf_tmp_dir.Phi();
      
	phi_eachlayer =atan2(ANTENNA_POSITION_START[ipol][ilayer][ifold][1],ANTENNA_POSITION_START[ipol][ilayer][ifold][0]);
       
	phi_deg =phi_deg- phi_eachlayer;
	 
	if(fabs(phi_deg) > fabs(phi_deg+2*PI)) phi_deg+=2*PI;
	if(fabs(phi_deg) > fabs(phi_deg-2*PI)) phi_deg-=2*PI;
	phi_deg =phi_deg*DEGRAD;
	totalAngledeg = phi_deg*phi_deg + theta_deg*theta_deg;
	if(totalAngledeg > 2500) totalAngledeg=2500;
	 
	extra_delay  = (totalAngledeg*totalAngledeg)*1.45676e-8;//pulled from Abby analysis
	extra_delay -= (totalAngledeg)*5.01452e-6;//pulled from Abby analysis
	
	arrival_times[ipol][ant_ctr]+=extra_delay*1E-9;
	ant_ctr++;
      
      }  
    }
  }
  
  double minV = Tools::dMin(arrival_times[0],(number_all_antennas));
  double minH = Tools::dMin(arrival_times[1],(number_all_antennas));
  double first_trigger_time = Tools::dMin(minV, minH);
  for (int ipol=0; ipol<2; ipol++){
    for (int i=0;i<(number_all_antennas);i++){
      
      arrival_times[ipol][i] -= first_trigger_time;
      //cout<<"arrival_times boresight["<<i<<"] is "<<arrival_times[i]<<"\n";
    }
  }
} // GetArrivalTimesBoresights


void Anita::GetArrivalTimesBoresights(const Vector rf_direction[NLAYERS_MAX][NPHI_MAX]) {

  for (int ipol=0; ipol<2; ipol++){
    for (int antenna_index = 0; antenna_index < (number_all_antennas); antenna_index++) { //loop over layers on the payload
      int ilayer = (antenna_index<8)*0 + (antenna_index>7)*(antenna_index<16)*1+ (antenna_index>15)*(antenna_index<32)*2+(antenna_index>31)*3;
      int ifold = (ilayer<2)*(antenna_index%8)+(ilayer>1)*(antenna_index%16);
      arrival_times[ipol][antenna_index] = (antenna_positions[ipol][antenna_index] * rf_direction[ilayer][ifold]) / CLIGHT;

      // cout << antenna_index << " " << arrival_times[antenna_index] << " " << extraCableDelays[0][antenna_index] << " " ;
      //      arrival_times[ipol][antenna_index] += extraCableDelays[ipol][antenna_index];
      // cout << arrival_times[antenna_index] << endl;
    
      //  arrival_times[antenna_index]=0;
    } // for: loop over antenna layers
  }

  
  double minV = Tools::dMin(arrival_times[0],(number_all_antennas));
  double minH = Tools::dMin(arrival_times[1],(number_all_antennas));
  double first_trigger_time = Tools::dMin(minV, minH);
  for (int ipol=0; ipol<2; ipol++){
    for (int i=0;i<(number_all_antennas);i++){
      
      arrival_times[ipol][i] -= first_trigger_time;
      //cout<<"arrival_times boresight["<<i<<"] is "<<arrival_times[i]<<"\n";
    }
  }
} // GetArrivalTimesBoresights




void Anita::setphiTrigMask(UInt_t realTime_flightdata) {

  if (realTime_flightdata<realTime_tr_min || realTime_flightdata>realTime_tr_max) {
    phiTrigMask=0; // if the realTime for this balloon position is out of range then just set mask to 0
    phiTrigMaskH=0;
    l1TrigMask=0;
    l1TrigMaskH=0;
    deadTime=0;
  }
  else { // if it's in range
		
    iturf=turfratechain->GetEntryNumberWithBestIndex(realTime_flightdata); // find entry in turfratechain that is closest to this realTime_flightdata
    if (iturf<0){ // if it didn't find one
      phiTrigMask=0; // set to zero
      phiTrigMaskH=0;
      l1TrigMask=0;
      l1TrigMaskH=0;
      deadTime=0;
    }else{
      turfratechain->GetEvent(iturf);
    }
  } // end if it's in range
  
}



void Anita::setTimeDependentThresholds(UInt_t realTime_flightdata){
  
  if (realTime_flightdata<realTime_surf_min || realTime_flightdata>realTime_surf_max) {
    for(int ipol=0;ipol<2;ipol++){
      for (int iant=0;iant<48;iant++){
	scalers[ipol][iant]=thresholds[ipol][iant]=0.;
      }
    }
  }
  else { // if it's in range
		
    isurf=surfchain->GetEntryNumberWithBestIndex(realTime_flightdata); // find entry in surfchain that is closest to this realTime_flightdata
    if (isurf<0){ // if it didn't find one
      for(int ipol=0;ipol<2;ipol++){
	for (int iant=0;iant<48;iant++){
	  scalers[ipol][iant]=thresholds[ipol][iant]=0.;
	}
      }
    }else{
      surfchain->GetEvent(isurf);
    }
  } // end if it's in range
  

}


#ifdef ANITA_UTIL_EXISTS
void Anita::readImpulseResponseDigitizer(Settings *settings1){
  
  // Set deltaT to be used in the convolution
  deltaT = 1/(2.6*16);
  string graphNames[2][3][16];
  string fileName;
  double norm=1;
 
  // For ANITA-2 we have 1 impulse response for VPOL and 1 for HPOL
  // For ANITA-3 we have 3 impulse responses (Top, Middle, Bottom ring) for VPOL and 3 for HPOL.
  // Set Graph names for ANITA-2 and ANITA-3
  if (settings1->WHICH==8){
    fileName = ICEMC_DATA_DIR+"/sumPicoImpulse.root";
    
    for (int iring=0;iring<3;iring++){
      for (int iphi=1;iphi<17;iphi++){
	graphNames[0][iring][iphi]="grImpRespV";
	graphNames[1][iring][iphi]="grImpRespH";
      }
    }
    //Now need to scale our impulse response from unit areas to the area of kronecker-delta (i.e dt)
    norm*=0.1;
  } else if(settings1->WHICH==9 || settings1->WHICH==10){

    fileName = ICEMC_DATA_DIR+"/Anita3_ImpulseResponseDigitizer.root";

    string spol[2] ={"V", "H"};
    string sring[3]={"T", "M", "B"};
    
    for (int ipol=0;ipol<2;ipol++){
      for (int iring=0;iring<3;iring++){
	for (int iphi=0;iphi<16;iphi++){
	  graphNames[ipol][iring][iphi] = Form("g%02d%s%s", iphi+1, sring[iring].c_str(), spol[ipol].c_str() ) ;
	}
      }
    }

    // Normalisation of digitizer impulse response might be off by 3dB
    // See LC's talk at icemc meeting of 2017 Nov 13
    norm *= TMath::Power(10., -3./20.);

    norm *= TMath::Power(10., +1./20.);
    
  }

  // Read in input file
  TFile fImpulse(fileName.c_str());
  
  if(!fImpulse.IsOpen()) {
    std::cerr << "Couldn't read siganl chain impulse response from " << fileName << "\n";
    exit(0);
  } else {

    for (int ipol=0;ipol<2;ipol++){
      for (int iring=0;iring<3;iring++){
	for (int iphi=0;iphi<16;iphi++){
	  // Read graph
	  TGraph *grTemp = (TGraph*) fImpulse.Get(graphNames[ipol][iring][iphi].c_str());
	  if(!grTemp) {
	    std::cerr << "Couldn't read signal chain impulse response" << graphNames[ipol][iring][iphi] << " from file " << fileName << "\n";
	    exit(0);
	  }
	  // Interpolate to high sampling rate that will be used for the convolution
	  TGraph *grInt = FFTtools::getInterpolatedGraph(grTemp,deltaT); 
	  Int_t nPoints  = grInt->GetN();
	  Double_t *newx = grInt->GetX();
	  Double_t *newy = grInt->GetY();
	  // Normalise
	  for (int i=0;i<nPoints;i++){
	    newy[i]=newy[i]*norm;
	    // change time axis from ns to s
	    newx[i]=newx[i]*1E-9;
	  }
	  // Pave to 0
	  int paveNum = 8533;
	  grTemp = new TGraph(nPoints,  newx, newy);
	  
	  fSignalChainResponseDigitizer[ipol][iring][iphi] = new RFSignal(FFTtools::padWaveToLength(grTemp, paveNum));
	  
	  delete grInt;
	  delete grTemp;

	  TGraph *gDig  = fSignalChainResponseDigitizer[ipol][iring][iphi]->getFreqMagGraph();
	  // Smooth out the high frequency 
	  double temparray[512];
	  for(int i=0;i<numFreqs;i++) {
	    temparray[i] =  gDig->Eval(freqs[i]*1e6);
	    // cout <<  i <<  " " << ipol << " " << iring << " " << iphi << " " << freqs[i] << " " << fSignalChainResponseDigitizerFreqDomain[ipol][iring][iphi][i]<< endl;
	  }
	  
	  // Smoothing magnitude response a bit to avoid trig/dig ratio explodes
	  for (int i=0; i<numFreqs;i++){
	    if (freqs[i]<900.){
	      fSignalChainResponseDigitizerFreqDomain[ipol][iring][iphi][i]  = temparray[i];
	    } else {
	      fSignalChainResponseDigitizerFreqDomain[ipol][iring][iphi][i]  = (temparray[i-2] + temparray[i-1] + temparray[i] + temparray[i+1] + temparray[i+2])/5.;
	    }
	  }
	  
	  delete gDig;
	  
	}
      }
    }
  }

}

void Anita::readTuffResponseDigitizer(Settings *settings1){
  // for loops to make the RFSignal array that can be used in applyImpulseResponseDigitizer of ChanTrigger.cc
  // ipol is the polarization "v" or "H" 
  // ring is the number 3 for tmb or bottom middle top
  // iphi is the antenna number
  // ituff is the notch directory
  TString filename;
  string snotch_dir[6]={"notches_260_0_0","notches_260_375_0","notches_260_0_460","notches_260_385_0","notches_260_365_0","notches_260_375_460"};
  string spol[2] = {"V","H"};
  string sring[3] = {"T","M","B"};
 // Set deltaT to be used in the convolution
  deltaT = 1/(2.6*16);
  for(int ipol=0; ipol<=1; ipol++) {
    for(int iring = 0; iring<=2; iring++){
      for(int iphi=0; iphi<=15; iphi++) {
	for(int ituff=0; ituff <=5; ituff++) {
	  if(iphi+1 < 10) {
	    filename = Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/%s/0%d%s%s.imp",getenv("ANITA_UTIL_INSTALL_DIR"), snotch_dir[ituff].c_str(), iphi+1, sring[iring].c_str(), spol[ipol].c_str());
	  } 
          else {
	    filename = Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/%s/%d%s%s.imp",getenv("ANITA_UTIL_INSTALL_DIR"), snotch_dir[ituff].c_str(), iphi+1, sring[iring].c_str(), spol[ipol].c_str());
	  }
	  TGraph *gtemp = new TGraph(filename);
	  // interpolate
	  TGraph *gint = Tools::getInterpolatedGraph(gtemp,deltaT); 
	  int paveNum=8533; // change for 0 to just signal back 
	  fSignalChainResponseDigitizerTuffs[ipol][iring][iphi][ituff] = new RFSignal(FFTtools::padWaveToLength(gint, paveNum)); 
	  delete gint;
	  delete gtemp;
	}// end for loop ituff
      } // end for loop iphi
    }// end for loop iring
  }// end for loop ipol
}

void Anita::readTuffResponseTrigger(Settings *settings1){
  // for loops to make the RFSignal array that can be used in applyImpulseResponseTrigger of ChanTrigger.cc Do we need one for each antenna???
  TString filename;
  string snotch_dir[6]={"trigconfigA.imp","trigconfigB.imp","trigconfigC.imp","trigconfigG.imp","trigconfigO.imp","trigconfigP.imp"};
  string spol[2] = {"V","H"};
  string sring[3] = {"T","M","B"};
 // Set deltaT to be used in the convolution
  deltaT = 1/(2.6*16);
  for(int ipol=0; ipol<=1; ipol++) {
    for(int iring = 0; iring<=2; iring++){
      for(int iphi=0; iphi<=15; iphi++) {
        for(int ituff=0; ituff <=5; ituff++) {
            filename = Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/%s",getenv("ANITA_UTIL_INSTALL_DIR"), snotch_dir[ituff].c_str());
            //debugging
            //cout << Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/%s",getenv("ANITA_UTIL_INSTALL_DIR"), snotch_dir[ituff].c_str()) << endl;
          TGraph *gtemp = new TGraph(filename);
          // interpolate
          TGraph *gint = Tools::getInterpolatedGraph(gtemp,deltaT); 
          int paveNum=8533; // change for 0 to just signal back 
          fSignalChainResponseTriggerTuffs[ipol][iring][iphi][ituff] = new RFSignal(FFTtools::padWaveToLength(gint, paveNum)); 
          delete gint;
          delete gtemp;
        }// end for loop ituff
      } // end for loop iphi
    }// end for loop iring
  }// end for loop ipol
}

void Anita::readNoiseFromFlight(Settings *settings1){
  
  TFile *fRayleighAnita3 = new TFile((ICEMC_DATA_DIR+"/RayleighAmplitudesAnita3_noSun_Interp.root").c_str(), "read");
  
  for (int iant=0;iant<48;iant++){
    RayleighFits[0][iant] = (TGraph*)fRayleighAnita3->Get(Form("grSigma%dV_interp", iant+1));
    RayleighFits[1][iant] = (TGraph*)fRayleighAnita3->Get(Form("grSigma%dH_interp", iant+1));
  }
  
  Double_t *timeVals = new Double_t [780];
  Double_t *voltVals = new Double_t [780];
  for(int i=0;i<780;i++){
    timeVals[i] = i*1./2.6;
    voltVals[i] = 0.1;
  }
  
  RFSignal *rfTemplate = new RFSignal(780,timeVals,voltVals,true);
  
  numFreqs=rfTemplate->getNumFreqs();
  freqs=rfTemplate->getFreqs();
  
  fRand = new TRandom3(settings1->SEED);
  
}



void Anita::readImpulseResponseTrigger(Settings *settings1){

  // So far only available for ANITA-3
  
  // Set deltaT to be used in the convolution
  deltaT = 1/(2.6*16);
  string graphNames[2][3][16];
  string fileName;
  double norm=1;

  if(settings1->WHICH==9 || settings1->WHICH==10){

    fileName = ICEMC_DATA_DIR+"/Anita3_ImpulseResponseTrigger.root";

    string spol[2] ={"V", "H"};
    string sring[3]={"T", "M", "B"};
    
    for (int ipol=0;ipol<2;ipol++){
      for (int iring=0;iring<3;iring++){
	for (int iphi=0;iphi<16;iphi++){
	  graphNames[ipol][iring][iphi]= Form("gTrigPath") ;
	}
      }
    }

    // Impulse response accounts for trigger/digitizer splitter
    // norm *= sqrt(2);
    // if we are using the trigger impulse response and the old noise
    // we need to add this 7dB attenuation to have sensible results
    // see LC's talk on 2017 Nov 13
    if (!settings1->NOISEFROMFLIGHTTRIGGER) norm *= TMath::Power(10, -7/20.);
  }

  // Read in input file
  TFile fImpulse(fileName.c_str());
  
  if(!fImpulse.IsOpen()) {
    std::cerr << "Couldn't read siganl chain impulse response from " << fileName << "\n";
    exit(0);
  } else {

    for (int ipol=0; ipol<2; ipol++){
      for (int iring=0; iring<3; iring++){
	for (int iphi=0; iphi<16; iphi++){
	  // Read graph
	  TGraph *grTemp = (TGraph*) fImpulse.Get(graphNames[ipol][iring][iphi].c_str());
	  if(!grTemp) {
	    std::cerr << "Couldn't read signal chain impulse response" << graphNames[ipol][iring][iphi] << " from file " << fileName << "\n";
	    exit(0);
	  }
	  // Interpolate to high sampling rate that will be used for the convolution
	  TGraph *grInt = FFTtools::getInterpolatedGraph(grTemp,deltaT); 
	  Int_t nPoints  = grInt->GetN();
	  Double_t *newx = grInt->GetX();
	  Double_t *newy = grInt->GetY();
	  // Normalise
	  for (int i=0;i<nPoints;i++){
	    newy[i]=newy[i]*norm;
	    // change time axis from ns to s
	    newx[i]=newx[i]*1E-9;
	  }
	  // Pave to 0
	  int paveNum = 8533;
	  grTemp = new TGraph(nPoints,  newx, newy);
	
	  fSignalChainResponseTrigger[ipol][iring][iphi] = new RFSignal(FFTtools::padWaveToLength(grTemp, paveNum));
	
	  delete grInt;
	  delete grTemp;

	  
	  TGraph *gTrig = fSignalChainResponseTrigger[ipol][iring][iphi]->getFreqMagGraph();
	  for(int i=0;i<numFreqs;i++) {
	    fSignalChainResponseTriggerFreqDomain[ipol][iring][iphi][i]    = gTrig->Eval(freqs[i]*1e6);
	    // cout <<  i <<  " " << ipol << " " << iring << " " << iphi << " " << freqs[i] << " " << fSignalChainResponseDigitizerFreqDomain[ipol][iring][iphi][i] << " " << fSignalChainResponseTriggerFreqDomain[ipol][iring][iphi][i] << endl;
	  }
	  delete gTrig;

	}
      }
    }
  }


  double dig, trig;
  // TFile *fout = new TFile("RatioTrigDigResponses.root", "recreate");
  
  for (int ipol=0; ipol<2; ipol++){
    for (int iring=0; iring<3; iring++){
      for (int iphi=0; iphi<16; iphi++){
	
	for(int i=0;i<numFreqs;i++) {
	  if (freqs[i]<160.) {
	    fRatioTriggerDigitizerFreqDomain[ipol][iring][iphi][i]=0.1;  
	  } else {
	    dig    = fSignalChainResponseDigitizerFreqDomain[ipol][iring][iphi][i];
	    trig   = fSignalChainResponseTriggerFreqDomain[ipol][iring][iphi][i];
	    fRatioTriggerDigitizerFreqDomain[ipol][iring][iphi][i]    = (trig/dig);
	  }
	}
	// TGraph *temp = new TGraph (numFreqs, freqs, fRatioTriggerDigitizerFreqDomain[ipol][iring][iphi]);
	// temp->Write(Form("gratio_%d_%d_%d", ipol, iring, iphi));
	// delete temp;
    	
      }
    }
  }
  // delete fout;


  


  
}


void Anita::readTriggerEfficiencyScanPulser(Settings *settings1){
  
  if(settings1->WHICH==9){
     
    string fileName = ICEMC_DATA_DIR+"/TriggerEfficiencyScanPulser_anita3.root";
    TFile *f = new TFile(fileName.c_str(), "read");

    // Get average pulse as measured by scope
    gPulseAtAmpa  = (TGraph*)f->Get("gAvgPulseAtAmpa");
    
    bool useDelayGenerator = false;

    double maxDelays =  (Tools::dMax(trigEffScanRingDelay, 3) + Tools::dMax(trigEffScanPhiDelay,5) );
    maxDelays       -=  (Tools::dMin(trigEffScanRingDelay, 3) + Tools::dMin(trigEffScanPhiDelay,5) );
    
    if (maxDelays!=0) useDelayGenerator=true;
    
    for (int i=0;i<gPulseAtAmpa->GetN();i++){
      // 7db fixed attenuation
      gPulseAtAmpa->GetY()[i]*=TMath::Power(10,-7./20.);

      // Variable attenuation of central phi sector
      gPulseAtAmpa->GetY()[i]*=TMath::Power(10, trigEffScanAtt[2]*1./20.);

      // Signal in a 12-way splitter 
      gPulseAtAmpa->GetY()[i]*=TMath::Power(10, -10.8/20.);

      // Attenutation due to delay generator
      if (useDelayGenerator){
	gPulseAtAmpa->GetY()[i]*=TMath::Power(10, -12/20.);
      }

      // Splitter between digitizer and trigger path
      gPulseAtAmpa->GetY()[i]*=TMath::Power(10,-3./20.);
       
    }

    // To get the correct interpolation we need to shift the waveform
    // We shift it by 8*(1/(2.6*16)) ns
    
    double *x = gPulseAtAmpa->GetX();
    double *y = gPulseAtAmpa->GetY();
    double xnew[10000], ynew[10000];
    int N = gPulseAtAmpa->GetN();
    
    int newn = N - 6;
    
    for (int j=6; j<N; j++){
      xnew[j-6] = x[j];
      ynew[j-6] = y[j];
    }
    
    TGraph *gtemp = new TGraph (newn, xnew, ynew);
    
    TGraph *gPulseAtAmpaInt = FFTtools::getInterpolatedGraph(gtemp, 1/2.6);
    
    for (int i=0;i<HALFNFOUR;i++){
      trigEffScanPulseAtAmpa[i] = gPulseAtAmpaInt->Eval(fTimes[i]);
    }

    gPulseAtAmpa  = FFTtools::translateGraph(gPulseAtAmpa, 77.5721);
    
    delete gPulseAtAmpaInt;
    delete gtemp;
     
    for (int isample=0;isample<250;isample++){
      // Get average waveform at SURF as measured by scope
      TGraph *gPulseAtSurf = (TGraph*)f->Get(Form("gSamplePulseAtSurf_%i", isample));
     
      TGraph *gPulseAtSurfInt = FFTtools::getInterpolatedGraph(gPulseAtSurf, 1/(2.6));
      double *y2 = gPulseAtSurfInt->GetY();
      // 20dB attenuation was applied at the scope
      for (int i=0;i<HALFNFOUR;i++){
	trigEffScanPulseAtSurf[isample][i]=y2[i]/10.;
      }
      
      delete gPulseAtSurfInt;
      delete gPulseAtSurf;
        
    }


     
    f->Close();
  }else{
    cout << "Impulse response on trigger path can only be used with ANITA-3" << endl;
    exit(1);
  }
 
}

#endif


void Anita::calculateDelaysForEfficiencyScan(){

  int irx, phiIndex;

  for (int iant=0; iant<48; iant++){
    
    phiIndex = trigEffScanPhi - (iant%16);
    if (phiIndex>8) phiIndex=phiIndex-16;

    if(TMath::Abs(phiIndex)<=2){

      irx = iant;
      if (iant<16) {
	if (iant%2==0) irx = iant/2;
	else          irx = 8 + iant/2;
      }

      // Add phi sector delay
      arrival_times[0][irx] += trigEffScanPhiDelay[phiIndex+2];

      // Check if we are adding the ring delay to this phi sector
      if (trigEffScanApplyRingDelay[phiIndex+2]>0){
	// Add ring delay (T-M, M-B, T-B)
	if (iant<16)       arrival_times[0][irx] += trigEffScanRingDelay[0] + trigEffScanRingDelay[2];
	else if (iant<32)  arrival_times[0][irx] += trigEffScanRingDelay[1];
      
      }
    }
    
  }

}
