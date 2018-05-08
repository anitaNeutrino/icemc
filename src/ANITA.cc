#include "ANITA.h"
#include "GlobalTrigger.h"
#include "ChanTrigger.h"
#include "Tools.h"
#include "counting.hh"
#include "Settings.h"
#include "IcemcLog.h"
#include "Constants.h"
#include "screen.hh"
#include "ray.hh"
#include "VoltsRX.h"


icemc::ANITA::ANITA(const Settings* settings, Counting* retardedClass, Ray* sillyRay, Screen* sillyPanel)
  : fSettingsPtrIDontOwn(settings), fCountingPtrIDontOwn(retardedClass), fRayPtrIDontOwn(sillyRay), fScreenPtrIDontOwn(sillyPanel)
{
  for(int i=0; i < getNumRX(); i++){
    Vector v;
    testVecNotRealYet.push_back(v);
  }
  
  for(int i=0; i < getNumRX(); i++){
    fWaveformsRX.push_back(TGraph());
  }
}
  
icemc::ANITA::~ANITA(){
  // does nothing
}






void icemc::ANITA::addSignalToRX(const icemc::AskaryanSignal& signal, int rx){

  TGraph& gr = fWaveformsRX.at(rx);
  // apply the antenna gain...

  double e_component=0;
  double h_component=0;
  double n_component=0;
  double e_component_kvector=0;
  double h_component_kvector=0;
  double n_component_kvector=0;
  double hitangle_e=0;
  double hitangle_h=0;


  Vector n_eplane;
  Vector n_hplane;
  Vector n_normal;


  // this->GetAntennaOrientation(fSettingsPtrIDontOwn,  this,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);  
  
  //this->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal,  panel1->GetVec2bln(jpt), e_component_kvector,  h_component_kvector,  n_component_kvector);
  this->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal,  signal.poynting, e_component_kvector,  h_component_kvector,  n_component_kvector);
  // this->GetEcompHcompEvector(fSettingsPtrIDontOwn,  n_eplane,  n_hplane,  panel1->GetPol(jpt),  e_component,  h_component,  n_component);
  this->GetEcompHcompEvector(fSettingsPtrIDontOwn,  n_eplane,  n_hplane,  signal.polarization,  e_component,  h_component,  n_component);  
  this->GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector, hitangle_e, hitangle_h);

  // this->AntennaGain(fSettingsPtrIDontOwn, hitangle_e, hitangle_h, e_component, h_component, k, tmp_vhz[0][k], tmp_vhz[1][k]);

}





// bool icemc::ANITA::applyTrigger(const std::vector<TGraph>& pureSignalVoltageTimeGraphs, const TVector& poyntingVector, const TVector& polarizationVector){
bool icemc::ANITA::applyTrigger(){
  
  //////////////////////////////////////
  //       EVALUATE GLOBAL TRIGGER    //
  //          FOR VPOL AND HPOL       //
  //////////////////////////////////////

  int thispasses[Anita::NPOL]={0,0};
  GlobalTrigger globalTrigger(fSettingsPtrIDontOwn, dynamic_cast<Anita*>(this));

  // make a global trigger object (but don't touch the electric fences)
  // globaltrig1 = new GlobalTrigger(fSettingsPtrIDontOwn, anita1);
  // globaltrig1 = new GlobalTrigger(fSettingsPtrIDontOwn, fDetector);

  Tools::Zero(this->arrival_times[0], Anita::NLAYERS_MAX*Anita::NPHI_MAX);
  Tools::Zero(this->arrival_times[1], Anita::NLAYERS_MAX*Anita::NPHI_MAX);
  if (!fSettingsPtrIDontOwn->TRIGGEREFFSCAN){
    if(fSettingsPtrIDontOwn->BORESIGHTS){
      this->GetArrivalTimesBoresights(fRayPtrIDontOwn->n_exit2bn_eachboresight[2]);
    }
    else{
      // this->GetArrivalTimes(fRayPtrIDontOwn->n_exit2bn[2],bn1,fSettingsPtrIDontOwn);
      this->GetArrivalTimes(fRayPtrIDontOwn->n_exit2bn[2],this,fSettingsPtrIDontOwn);
    }	
  }
  this->rx_minarrivaltime=Tools::WhichIsMin(this->arrival_times[0], fSettingsPtrIDontOwn->NANTENNAS);



  //For verification plots - added by Stephen
  int max_antenna0 = -1;  //antenna with the peak voltage,  top layer
  int max_antenna1 = -1;  //antenna with the peak voltage,  middle layer
  int max_antenna2 = -1;  //antenna with the peak voltage,  bottom layer
  double voltagearray[Anita::NLAYERS_MAX*Anita::NPHI_MAX]; //Records max voltages on each antenna for one neutrino
  int discones_passing;  // number of discones that pass

  double max_antenna_volts0 = 0; //Voltage on the antenna with maximum signal,  top layer
  double max_antenna_volts0_em = 0; //Component of voltage from em shower on the antenna with maximum signal,  top layer
  double max_antenna_volts1 = 0; //Voltage on the antenna with maximum signal,  middle layer
  double max_antenna_volts2 = 0; //Voltage on the antenna with maximum signal,  bottom layer

  double e_comp_max1=0;
  double h_comp_max1=0;
  double e_comp_max2=0;
  double h_comp_max2=0;
  double e_comp_max3=0;
  double h_comp_max3=0;

  double hitangle_e = 0; // angle the ray hits the antenna wrt e-plane
  double hitangle_h = 0; // angle the ray hits the antenna wrt h-plane
  double hitangle_e_all[Anita::NANTENNAS_MAX] = {0};         // hit angles rel. to e plane stored for each antenna
  double hitangle_h_all[Anita::NANTENNAS_MAX] = {0};         // hit angles rel. to h plane stored for each antenna

  // for comparing with peter
  double sumsignal[5]={0.};
  double sumsignal_aftertaper[5]={0.};
  

  int l3trig[Anita::NPOL];  // 16 bit number which says which phi sectors pass L3 V-POL
  // For each trigger layer,  which "clumps" pass L2.  16 bit,  16 bit and 8 bit for layers 1 & 2 and nadirs
  int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];
  //For each trigger layer,  which antennas pass L1.  16 bit,  16 bit and 8 bit and layers 1,  2 and nadirs
  int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];


  double justNoise_trig[2][48][512] = {{{0}}};
  double justSignal_trig[2][48][512] = {{{0}}};
  double justNoise_dig[2][48][512] = {{{0}}};
  double justSignal_dig[2][48][512] = {{{0}}};

  
  // these variables are for energy reconstruction studies
  double undogaintoheight_e=0;
  double undogaintoheight_h=0;

  double undogaintoheight_e_array[4];
  double undogaintoheight_h_array[4];
  double nbins_array[4];

  double rec_efield=0;
  double true_efield=0;

  double rec_efield_array[4];
  double true_efield_array[4];
  // end energy reconstruction variables

  int count_pass=0;  // how many total trigger channels pass (4 bandwidth slices*2 pol * nrx)


  
  undogaintoheight_e=0;
  undogaintoheight_h=0;

  for (int k=0;k<4;k++) {
    undogaintoheight_e_array[k]=0.;
    undogaintoheight_h_array[k]=0.;
    nbins_array[k]=0;
    true_efield_array[k]=0.;
    rec_efield_array[k]=0.;
  }

  rec_efield=0;
  true_efield=0;
  
  
  //Zeroing
  for (int i=0;i<fSettingsPtrIDontOwn->NANTENNAS;i++) {
    voltagearray[i]=0;
    discones_passing=0;
  } //Zero the trigger array

  max_antenna0=-1;
  max_antenna1=-1;
  max_antenna2=-1;
  max_antenna_volts0 =0;
  max_antenna_volts1 =0;
  max_antenna_volts2 =0;
  e_comp_max1 = 0;
  h_comp_max1 = 0;
  e_comp_max2 = 0;
  h_comp_max2 = 0;
  e_comp_max3 = 0;
  h_comp_max3 = 0;
  //End zeroing


  // start looping over antennnas.
  // ilayer loops through vertical layers

  if (fSettingsPtrIDontOwn->SLAC){
    Log().fslac_hitangles << this->sslacpositions[this->islacposition] << "\n";
  }
  // if (RANDOMISEPOL) {
  //   double rotateangle=gRandom->Gaus(RANDOMISEPOL*constants::RADDEG);
  //   n_pol=n_pol.Rotate(rotateangle, fRayPtrIDontOwn->n_exit2bn[2]);
  // }

  // if (this->WHICHPATH==4) {
  //   Tools::Zero(sumsignal, 5);
  //   for (int k=0;k<Anita::NFREQ;k++){
  //     // IntegrateBands(anita1, k, panel1, this->freq, this->r_bn.Distance(interaction1->posnu)/1.E6, sumsignal);
  //     IntegrateBands(this, k, panel1, this->freq, this->r_bn.Distance(interaction1->posnu)/1.E6, sumsignal);
  //   }
  // }//end if whichpath==4


  if (fSettingsPtrIDontOwn->CENTER){
    this->CenterPayload(hitangle_e);
  }

  // if (fSettingsPtrIDontOwn->MAKEVERTICAL) {
  //   n_pol=this->n_bn;
  //   // rotate n_exit2bn too
  //   // rotation axis n_bn crossed with n_exit2bn
  //   Vector rotationaxis=fRayPtrIDontOwn->n_exit2bn[2].Cross(this->n_bn);
  //   double rotateangle=constants::PI/2.-fRayPtrIDontOwn->n_exit2bn[2].Dot(this->n_bn);
  //   fRayPtrIDontOwn->n_exit2bn[2]=fRayPtrIDontOwn->n_exit2bn[2].Rotate(rotateangle, rotationaxis);

  //   for (int ilayer=0; ilayer < fSettingsPtrIDontOwn->NLAYERS; ilayer++) { // loop over layers on the payload
  //     // ifold loops over phi
  //     for (int ifold=0;ifold<this->NRX_PHI[ilayer];ifold++) {
  // 	Vector rotationaxis2=fRayPtrIDontOwn->n_exit2bn_eachboresight[2][ilayer][ifold].Cross(n_pol_eachboresight[ilayer][ifold]);
  // 	double rotateangle2=constants::PI/2.-fRayPtrIDontOwn->n_exit2bn_eachboresight[2][ilayer][ifold].Dot(n_pol_eachboresight[ilayer][ifold]);
  // 	fRayPtrIDontOwn->n_exit2bn_eachboresight[2][ilayer][ifold].Rotate(rotateangle2, rotationaxis2);
  //     } // end loop over phi
  //   } // end loop over layers
  // }//end if MAKEVERTICAL

  globalTrigger.volts_rx_rfcm_trigger.assign(16,  vector <vector <double> >(3,  vector <double>(0)));
  this->rms_rfcm_e_single_event = 0;

  Vector n_eplane = constants::const_z;
  Vector n_hplane = -constants::const_y;
  Vector n_normal = constants::const_x;

  
  // variable declarations for functions GetEcompHcompEvector and GetEcompHcompkvector - oindree
  double e_component=0; // E comp along polarization
  double h_component=0; // H comp along polarization
  double n_component=0; // normal comp along polarization
  double e_component_kvector=0; // component of e-field along the rx e-plane
  double h_component_kvector=0; // component of the e-field along the rx h-plane
  double n_component_kvector=0; // component of the e-field along the normal
  double bwslice_vnoise_thislayer[4];// for filling tree6b,  noise for each bandwidth on each layer
  
  double rx0_signal_eachband[2][5];
  double rx0_threshold_eachband[2][5];
  double rx0_noise_eachband[2][5];
  int rx0_passes_eachband[2][5];

  Vector ant_max_normal0; //Vector normal to the face of the antenna with the maximum signal for a single neutrino,  top layer
  Vector ant_max_normal1; //Vector normal to the face of the antenna with the maximum signal for a single neutrino,  middle layer
  Vector ant_max_normal2; //Vector normal to the face of the antenna with the maximum signal for a single neutrino,  bottom layer

  Vector ant_normal; //Vector normal to the face of the antenna
  int nchannels_triggered = 0; // total number of channels triggered
  int nchannels_perrx_triggered[48]; // total number of channels triggered

  // these are declared here so that they can be stuck into trees
  int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //counting how many pass trigger requirement

  int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX]; //counting how many pass trigger requirement

  double thresholdsAnt[48][2][5] = {{{0}}};
  // double thresholdsAntPass[48][2][5] = {{{0}}};
  
  
  // if (!fSettingsPtrIDontOwn->BORESIGHTS) {
  //   this->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal,  fRayPtrIDontOwn->n_exit2bn[2], e_component_kvector,  h_component_kvector,  n_component_kvector);
  //   const icemc::Vector& n_pol = fWaveformsRX.back().polarization;
  //   this->GetEcompHcompEvector(fSettingsPtrIDontOwn,  n_eplane,  n_hplane,  n_pol,  e_component,  h_component,  n_component);
  // }

  int count_rx = 0;
  VoltsRX voltsRX;
  for (int ilayer=0; ilayer < fSettingsPtrIDontOwn->NLAYERS; ilayer++) { // loop over layers on the payload
    for (int ifold=0;ifold<this->NRX_PHI[ilayer];ifold++) { // ifold loops over phi
          
      ChanTrigger chantrig1;
      // chantrig1.InitializeEachBand(anita1);
      chantrig1.InitializeEachBand(this);

      // this->GetAntennaOrientation(fSettingsPtrIDontOwn,  anita1,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);
      this->GetAntennaOrientation(fSettingsPtrIDontOwn,  this,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);
 
      // if (fSettingsPtrIDontOwn->BORESIGHTS){ // i.e. if BORESIGHTS is true
      // 	this->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal,  fRayPtrIDontOwn->n_exit2bn_eachboresight[2][ilayer][ifold],  e_component_kvector,  h_component_kvector,  n_component_kvector);
      // 	this->GetEcompHcompEvector(fSettingsPtrIDontOwn,  n_eplane,  n_hplane,  n_pol_eachboresight[ilayer][ifold], e_component,  h_component,  n_component);
      // 	Log().fslac_hitangles << ilayer << "\t" << ifold << "\t" << hitangle_e << "\t" << hitangle_h << "\t" << e_component_kvector << "\t" << h_component_kvector << "\t" << fresnel1_eachboresight[ilayer][ifold] << " " << mag1_eachboresight[ilayer][ifold] << "\n";
      // }
      
      this->GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector, hitangle_e, hitangle_h);
      // store hitangles for plotting
      hitangle_h_all[count_rx]=hitangle_h;
      hitangle_e_all[count_rx]=hitangle_e;
      // for debugging
      // if (ro.h6.GetEntries()<fSettingsPtrIDontOwn->HIST_MAX_ENTRIES && !fSettingsPtrIDontOwn->ONLYFINAL && fSettingsPtrIDontOwn->HIST==1){
      // 	ro.h6.Fill(hitangle_h, fRayPtrIDontOwn->n_exit2bn[2].Dot(this->n_bn));
      // }
      
      int antNum = this->GetRxTriggerNumbering(ilayer, ifold);

      // chantrig1.ApplyAntennaGain(fSettingsPtrIDontOwn, anita1, bn1, panel1, antNum, n_eplane, n_hplane, n_normal);
      // chantrig1.ApplyAntennaGain(fSettingsPtrIDontOwn, this, this, panel1, antNum, n_eplane, n_hplane, n_normal);
      // chantrig1.ApplyAntennaGain(fSettingsPtrIDontOwn, this, this, panel1, antNum, n_eplane, n_hplane, n_normal);
      chantrig1.ApplyAntennaGain(fSettingsPtrIDontOwn, this, this, fScreenPtrIDontOwn, antNum, n_eplane, n_hplane, n_normal);

      // chantrig1.TriggerPath(fSettingsPtrIDontOwn, anita1, antNum, bn1);
      // chantrig1.DigitizerPath(fSettingsPtrIDontOwn, anita1, antNum, bn1);
      chantrig1.TriggerPath(fSettingsPtrIDontOwn, this, antNum, this);
      chantrig1.DigitizerPath(fSettingsPtrIDontOwn, this, antNum, this);

      // chantrig1.TimeShiftAndSignalFluct(fSettingsPtrIDontOwn, anita1, ilayer, ifold, voltsRX.rfcm_lab_e_all,  voltsRX.rfcm_lab_h_all);
      chantrig1.TimeShiftAndSignalFluct(fSettingsPtrIDontOwn, this, ilayer, ifold, voltsRX.rfcm_lab_e_all,  voltsRX.rfcm_lab_h_all);	  

      // chantrig1.saveTriggerWaveforms(anita1, justSignal_trig[0][antNum], justSignal_trig[1][antNum], justNoise_trig[0][antNum], justNoise_trig[1][antNum]);
      // chantrig1.saveDigitizerWaveforms(anita1, justSignal_dig[0][antNum], justSignal_dig[1][antNum], justNoise_dig[0][antNum], justNoise_dig[1][antNum]);
      chantrig1.saveTriggerWaveforms(this, justSignal_trig[0][antNum], justSignal_trig[1][antNum], justNoise_trig[0][antNum], justNoise_trig[1][antNum]);
      chantrig1.saveDigitizerWaveforms(this, justSignal_dig[0][antNum], justSignal_dig[1][antNum], justNoise_dig[0][antNum], justNoise_dig[1][antNum]);
	  
      Tools::Zero(sumsignal, 5);

      if (this->WHICHPATH==4 && ilayer==this->GetLayer(this->rx_minarrivaltime) && ifold==this->GetIfold(this->rx_minarrivaltime)) {
	for (int ibw=0;ibw<5;ibw++) {
	  std::cout << "Just after Taper,  sumsignal is " << sumsignal_aftertaper[ibw] << "\n";
	  std::cout << "Just after antennagain,  sumsignal is " << sumsignal[ibw] << "\n";
	}
      }

      // for energy reconstruction studies
      if (count_rx==this->rx_minarrivaltime) {
	undogaintoheight_e/=(double)Anita::NFREQ;
	undogaintoheight_h/=(double)Anita::NFREQ;
	for (int k=0;k<4;k++) {
	  undogaintoheight_e_array[k]/=(double)nbins_array[k];
	  undogaintoheight_h_array[k]/=(double)nbins_array[k];
	}
      }
      if (fSettingsPtrIDontOwn->SCALEDOWNLCPRX1){
	globalTrigger.volts[0][ilayer][0]=globalTrigger.volts[0][ilayer][0]/sqrt(2.);
      }

      if (fSettingsPtrIDontOwn->RCPRX2ZERO){
	globalTrigger.volts[1][ilayer][1]=0.;
      }
	  
      if (fSettingsPtrIDontOwn->LCPRX2ZERO){
	globalTrigger.volts[0][ilayer][1]=0.;
      }
	  
      if (fSettingsPtrIDontOwn->SIGNAL_FLUCT) {
	if (fSettingsPtrIDontOwn->WHICH==0) {
	  globalTrigger.volts[ilayer][ifold][0]+=gRandom->Gaus(0., this->VNOISE_ANITALITE[ifold]);
	  globalTrigger.volts[ilayer][ifold][1]+=gRandom->Gaus(0., this->VNOISE_ANITALITE[ifold]);
	} //else
      } //if adding noise
      if (count_rx==this->rx_minarrivaltime) {
	rec_efield=sqrt(pow(globalTrigger.volts_original[0][ilayer][ifold]/(undogaintoheight_e*0.5), 2)+pow(globalTrigger.volts_original[1][ilayer][ifold]/(undogaintoheight_h*0.5), 2));
	for (int ibw=0;ibw<4;ibw++) {
	  rec_efield_array[ibw]=sqrt(pow(chantrig1.bwslice_volts_pole[ibw]/(undogaintoheight_e_array[ibw]*0.5), 2)+pow(chantrig1.bwslice_volts_polh[ibw]/(undogaintoheight_h_array[ibw]*0.5), 2));
	  bwslice_vnoise_thislayer[ibw]=this->bwslice_vnoise[ilayer][ibw];// this is just for filling into a tree
	} // end loop over bandwidth slices
      } // end if this is the closest antenna

      //+++++//+++++//+++++//+++++//+++++//+++++//+++++

      // chantrig1.WhichBandsPass(fSettingsPtrIDontOwn, anita1, globaltrig1, bn1, ilayer, ifold,  viewangle-askFreqGen.GetChangle(), emfrac, hadfrac, thresholdsAnt[antNum]);
      // chantrig1.WhichBandsPass(fSettingsPtrIDontOwn, anita1, globaltrig1, bn1, ilayer, ifold,  viewangle-askFreqGen.GetChangle(), showerProps.emFrac, showerProps.hadFrac, thresholdsAnt[antNum]);
      // chantrig1.WhichBandsPass(fSettingsPtrIDontOwn, this, globaltrig1, this, ilayer, ifold,  viewangle-askFreqGen.GetChangle(), showerProps.emFrac, showerProps.hadFrac, thresholdsAnt[antNum]);
      // chantrig1.WhichBandsPass(fSettingsPtrIDontOwn, this, &globalTrigger, this, ilayer, ifold,  viewangle-askFreqGen.GetChangle(), showerProps.emFrac, showerProps.hadFrac, thresholdsAnt[antNum]);
      chantrig1.WhichBandsPass(fSettingsPtrIDontOwn, this, &globalTrigger, this, ilayer, ifold, thresholdsAnt[antNum]);

	  
      if (Anita::GetAntennaNumber(ilayer, ifold)==this->rx_minarrivaltime) {
	for (int iband=0;iband<5;iband++) {
	  for (int ipol=0;ipol<2;ipol++) {
	    rx0_signal_eachband[ipol][iband]=chantrig1.signal_eachband[ipol][iband];
	    rx0_threshold_eachband[ipol][iband]=chantrig1.threshold_eachband[ipol][iband];
	    rx0_noise_eachband[ipol][iband]=chantrig1.noise_eachband[ipol][iband];
	    rx0_passes_eachband[ipol][iband]=chantrig1.passes_eachband[ipol][iband];
	  }
	}
      }

      //For verification plots: find antenna with max signal - added by Stephen
      if (ilayer == 0 && globalTrigger.volts[0][ilayer][ifold] > max_antenna_volts0) {
	max_antenna0 = count_rx;
	max_antenna_volts0 = globalTrigger.volts[0][ilayer][ifold];
	max_antenna_volts0_em=globalTrigger.volts_em[0][ilayer][ifold];
	ant_max_normal0 = ant_normal;
	e_comp_max1 = e_component;
	h_comp_max1 = h_component;
      }
      else if (ilayer == 0 && globalTrigger.volts[0][ilayer][ifold] == max_antenna_volts0 && globalTrigger.volts[0][ilayer][ifold] != 0){
	std::cout<<"Equal voltage on two antennas!  Event : "<<inu<<std::endl;
      }
      else if (ilayer == 1 && globalTrigger.volts[0][ilayer][ifold] > max_antenna_volts1) {
	max_antenna1 = count_rx;
	max_antenna_volts1 = globalTrigger.volts[0][ilayer][ifold];
	ant_max_normal1 = ant_normal;
	e_comp_max2 = e_component;
	h_comp_max2 = h_component;
      }
      else if (ilayer == 1 && globalTrigger.volts[0][ilayer][ifold] == max_antenna_volts1 && globalTrigger.volts[0][ilayer][ifold] != 0){
	std::cout<<"Equal voltage on two antennas!  Event : "<<inu<<std::endl;
      }
      else if (ilayer == 2 && globalTrigger.volts[0][ilayer][ifold] > max_antenna_volts2) {
	max_antenna2 = count_rx;
	max_antenna_volts2 = globalTrigger.volts[0][ilayer][ifold];
	ant_max_normal2 = ant_normal;
	e_comp_max3 = e_component;
	h_comp_max3 = h_component;
      }
      else if (ilayer == 2 && globalTrigger.volts[0][ilayer][ifold] == max_antenna_volts2 && globalTrigger.volts[0][ilayer][ifold] != 0){
	std::cout<<"Equal voltage on two antennas!  Event : "<<inu<<std::endl;
      }
      voltagearray[count_rx] = globalTrigger.volts[0][ilayer][ifold];
      //End verification plot block

      count_rx++; // counting antennas that we loop through,  for indexing

      // if (fSettingsPtrIDontOwn->TRIGTYPE==0 && ifold==1 && count_pass>=fSettingsPtrIDontOwn->NFOLD) { //added djg --line below fills "direct" voltage output file
      // 	Log().al_voltages_direct<<"0 0 0"<<"   "<<"    "<<globalTrigger.volts_original[1][0][0]<<"    "<<(globalTrigger.volts_original[0][0][0]/sqrt(2.))<<"     "<<globalTrigger.volts_original[1][0][1]<<"     "<<globalTrigger.volts_original[0][0][1]<<"      "<<this->VNOISE[0]<<"     "<<this->VNOISE[0]<<"     "<<this->VNOISE[0]<<"     "<<this->VNOISE[0]<<"  "<<fNeutrinoPath->weight<<std::endl;
      // }
    } //loop through the phi-fold antennas
  }  //loop through the layers of antennas


  this->rms_rfcm_e_single_event = sqrt(this->rms_rfcm_e_single_event / (this->HALFNFOUR * fSettingsPtrIDontOwn->NANTENNAS));

  // if(!fSettingsPtrIDontOwn->ROUGHNESS){
    //    if (fSettingsPtrIDontOwn->DISCONES==1)  {
    //   // loop through discones
    //   for (int idiscone=0;NDISCONES;idiscone++) {
    // 	ChanTrigger chantrig1;
    // 	volts_discone=0.;
    // 	polarfactor_discone=n_pol.Dot(this->n_bn); // beam pattern
    // 	for (int k=0;k<Anita::NFREQ;k++) {
    // 	  if (this->freq[k]>=FREQ_LOW_DISCONES && this->freq[k]<=FREQ_HIGH_DISCONES) {
    // 	    thislambda = constants::CLIGHT/askFreqGen.N_AIR/this->freq[k];
    // 	    heff_discone = thislambda*sqrt(2*constants::Zr*gain_dipole/constants::Z0/4/constants::PI*askFreqGen.N_AIR);   // effective height of dipole,  using formula from Ped's note

    // 	    volts_discone+=fScreenPtrIDontOwn->GetVmmhz_freq(k)*0.5*heff_discone*((fSettingsPtrIDontOwn->BW/1E6)/(double)Anita::NFREQ)*polarfactor_discone;
    // 	  }
    // 	}// end for k loop

    // 	vnoise_discone=this->VNOISE[0]*sqrt(BW_DISCONES/fSettingsPtrIDontOwn->BW_SEAVEYS);

    // 	if (fSettingsPtrIDontOwn->SIGNAL_FLUCT) {
    // 	  volts_discone+=gRandom->Gaus(0., vnoise_discone); // here I'm using the noise seen by an antenna pointed with a 10 degree cant.  Should be different for a discone but we'll change it later.
    // 	}

    // 	if (fabs(volts_discone)/vnoise_discone>this->maxthreshold){
    // 	  discones_passing++;
    // 	}
    //   } // end looping through discones
    // }
    //end if settings discones==1
  // }
  for (int irx=0;irx<fSettingsPtrIDontOwn->NANTENNAS;irx++) {
    nchannels_perrx_triggered[irx]=globalTrigger.nchannels_perrx_triggered[irx];
  }

  nchannels_triggered = Tools::iSum(globalTrigger.nchannels_perrx_triggered, fSettingsPtrIDontOwn->NANTENNAS); // find total number of antennas that were triggered.
  voltsRX.ave = GetAverageVoltageFromAntennasHit(fSettingsPtrIDontOwn, globalTrigger.nchannels_perrx_triggered, voltagearray, voltsRX.sum);





  




  


  // globalTrigger.PassesTrigger(fSettingsPtrIDontOwn, anita1, discones_passing, 2, l3trig, l2trig, l1trig, fSettingsPtrIDontOwn->antennaclump, loctrig, loctrig_nadironly, inu,
  // 				 thispasses);
  globalTrigger.PassesTrigger(fSettingsPtrIDontOwn, this, discones_passing, 2, l3trig, l2trig, l1trig, fSettingsPtrIDontOwn->antennaclump, loctrig, loctrig_nadironly, inu, thispasses);

  // for (int i=0;i<2;i++) {
  //   for (int j=0;j<16;j++) {
  //     for (int k=0;k<this->HALFNFOUR;k++) {
  // 	fCountPtrIDontOwn->nl1triggers[i][whichray]+=this->l1trig_anita3and4_inanita[i][j][k];
  //     }
  //   }
  // }  


  ///////////////////////////////////////
  //       Require that it passes      //
  //            global trigger         //
  ///////////////////////////////////////
  // for Anita-lite,  Anita Hill, just L1 requirement on 2 antennas. This option is currently disabled
  // Save events that generate an RF trigger or that are part of the min bias sample
  // Minimum bias sample: save all events that we could see at the payload
  // Independentely from the fact that they generated an RF trigger

  bool eventPassesTrigger = false;  
  if ( (thispasses[0]==1 && this->pol_allowed[0]==1)
       || (thispasses[1]==1 && this->pol_allowed[1]==1)
       || (fSettingsPtrIDontOwn->TRIGTYPE==0 && count_pass>=fSettingsPtrIDontOwn->NFOLD)
       || (fSettingsPtrIDontOwn->MINBIAS==1)){
    eventPassesTrigger = true;
  }

  return eventPassesTrigger;
}





double icemc::ANITA::GetAverageVoltageFromAntennasHit(const Settings *settings1, int *nchannels_perrx_triggered, double *voltagearray, double& volts_rx_sum) const {
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