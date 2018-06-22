#include "ANITA.h"
#include "Settings.h"
#include "GlobalTrigger.h"
#include "ChanTrigger.h"
#include "Tools.h"
#include "counting.hh"
#include "Settings.h"
#include "IcemcLog.h"
#include "Constants.h"
#include "screen.hh"
#include "RayTracer.h"
#include "VoltsRX.h"
#include "Seavey.h"
#include "RootOutput.h"
#include <memory>

#include "TFile.h" ///@todo remove after done debugging


icemc::ANITA::ANITA(const Settings* settings, RayTracer* sillyRay, Screen* sillyPanel, const RootOutput* ro)
  : fSettingsPtrIDontOwn(settings), fRayPtrIDontOwn(sillyRay), fScreenPtrIDontOwn(sillyPanel),
    fAnitaOutput(this, settings, ro->getOutputDir(), ro->getRun())
{
  const int n = 48; ///@todo don't hardcode this 
  for(int rx=0; rx < n; rx++){
    fSeaveys.emplace_back(icemc::Seavey(fSettingsPtrIDontOwn));
  }
}


icemc::ANITA::~ANITA(){
  
}


icemc::Position icemc::ANITA::getCenterOfDetector(UInt_t unixTime){
  (void) unixTime;

  // UInt_t theUnixTime = unixTime ? 1 : 0;
  for(int rx=0; rx < static_cast<int>(fSeaveys.size()); rx++){
    int layer, fold;
    getLayerFoldFromRX(rx, layer, fold);
  }
  
  return r_bn;
}


icemc::Vector icemc::ANITA::getPositionRX(Int_t rx) const {  

  return Position();
}


void icemc::ANITA::getLayerFoldFromRX(int rx, int& ilayer, int& ifold) const {
  int antNum = rx;
  
  // This is NOT how to do things...
  ilayer = -1;
  ifold = -1;
  for (int ilayerTemp=0 ;ilayerTemp < fSettingsPtrIDontOwn->NLAYERS; ilayerTemp++) { // loop over layers on the payload
    for (int ifoldTemp=0;ifoldTemp<this->NRX_PHI[ilayerTemp];ifoldTemp++) { // ifold loops over phi
      Int_t antNum2 = this->GetRxTriggerNumbering(ilayerTemp, ifoldTemp);
      if(antNum==antNum2){
	ilayer = ilayerTemp;
	ifold = ifoldTemp;
	break;
      }
    }
    if(ilayer > -1){
      break;
    }
  }
}



void icemc::ANITA::addSignalToRX(const icemc::PropagatingSignal& signal, int rx, int inu){

  int ifold, ilayer;
  getLayerFoldFromRX(rx, ilayer, ifold);
  
  static bool firstTime = true;
  if(inu == 397 && firstTime){
    for(int i=0; i < fSeaveys.size(); i++){
      if(i==34){
	fSeaveys.at(i).setDebug(true);
      }
    }
    if(rx==fSeaveys.size()-1){
      firstTime = false;
    }
  }
  else{
    for(int i=0; i < fSeaveys.size(); i++){
      fSeaveys.at(i).setDebug(false);
    }
  }
  
  if(rx >= 0 && rx < fSeaveys.size()){
    // @todo It makes much more sense to do this when the balloon position is updated!
    this->GetAntennaOrientation(fSettingsPtrIDontOwn,  this,  ilayer,  ifold,
				fSeaveys.at(rx).fEPlane, fSeaveys.at(rx).fHPlane, fSeaveys.at(rx).fNormal);

    // actually we add the signal to the new Seavey class
    fSeaveys.at(rx).addSignal(signal);
  }
}





// bool icemc::ANITA::applyTrigger(const std::vector<TGraph>& pureSignalVoltageTimeGraphs, const TVector& poyntingVector, const TVector& polarizationVector){
bool icemc::ANITA::applyTrigger(int inu){
  
  //////////////////////////////////////
  //       EVALUATE GLOBAL TRIGGER    //
  //          FOR VPOL AND HPOL       //
  //////////////////////////////////////

  int thispasses[Anita::NPOL]={0,0};

  auto globalTrigger = std::unique_ptr<GlobalTrigger>(new GlobalTrigger(fSettingsPtrIDontOwn, dynamic_cast<Anita*>(this)));

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

  // for comparing with peter
  double sumsignal[5]={0.};
  double sumsignal_aftertaper[5]={0.};
  

  constexpr int nAnt = 48;
  std::vector<double> justNoise_trig[NPOL][nAnt];
  std::vector<double> justSignal_trig[NPOL][nAnt];
  std::vector<double> justNoise_dig[NPOL][nAnt];
  std::vector<double> justSignal_dig[NPOL][nAnt];

  for(int pol=0;  pol < NPOL; pol++){
    for(int ant=0; ant < nAnt; ant++){
      justNoise_trig[pol][ant].resize(Anita::HALFNFOUR);
      justSignal_trig[pol][ant].resize(Anita::HALFNFOUR);
      justNoise_dig[pol][ant].resize(Anita::HALFNFOUR);
      justSignal_dig[pol][ant].resize(Anita::HALFNFOUR);
    }
  }
  // double justNoise_trig[2][48][Anita::HALFNFOUR] = {{{0}}};
  // double justSignal_trig[2][48][Anita::HALFNFOUR] = {{{0}}};
  // double justNoise_dig[2][48][Anita::HALFNFOUR] = {{{0}}};
  // double justSignal_dig[2][48][Anita::HALFNFOUR] = {{{0}}};
  

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
    icemcLog().fslac_hitangles << this->sslacpositions[this->islacposition] << "\n";
  }

  if (fSettingsPtrIDontOwn->CENTER){
    this->CenterPayload(hitangle_e);
  }
  
  globalTrigger->volts_rx_rfcm_trigger.assign(16,  vector <vector <double> >(3,  vector <double>(0)));
  this->rms_rfcm_e_single_event = 0;

  Vector n_eplane = constants::const_z;
  Vector n_hplane = -constants::const_y;
  Vector n_normal = constants::const_x;
  double e_component=0; // E comp along polarization
  double h_component=0; // H comp along polarization
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
  int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //counting how many pass trigger requirement
  int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX]; //counting how many pass trigger requirement
  double thresholdsAnt[48][2][5] = {{{0}}};

  int count_rx = 0;
  for (int ilayer=0; ilayer < fSettingsPtrIDontOwn->NLAYERS; ilayer++) { // loop over layers on the payload
    for (int ifold=0;ifold<this->NRX_PHI[ilayer];ifold++) { // ifold loops over phi

      ChanTrigger ct;
      ct.InitializeEachBand(this);
      this->GetAntennaOrientation(fSettingsPtrIDontOwn,  this,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);
      int antNum = this->GetRxTriggerNumbering(ilayer, ifold);
      TGraph grHack(1, &inu, &inu);
      ct.ApplyAntennaGain(fSettingsPtrIDontOwn, this, this, fScreenPtrIDontOwn, antNum, n_eplane, n_hplane, n_normal, &grHack);
      ct.TriggerPath(fSettingsPtrIDontOwn, this, antNum, this);
      ct.DigitizerPath(fSettingsPtrIDontOwn, this, antNum, this);
      ct.TimeShiftAndSignalFluct(fSettingsPtrIDontOwn, this, ilayer, ifold, fVoltsRX.rfcm_lab_e_all,  fVoltsRX.rfcm_lab_h_all);	  
      ct.saveTriggerWaveforms(this, &justSignal_trig[0][antNum][0], &justSignal_trig[1][antNum][0], &justNoise_trig[0][antNum][0], &justNoise_trig[1][antNum][0]);
      ct.saveDigitizerWaveforms(this, &justSignal_dig[0][antNum][0], &justSignal_dig[1][antNum][0], &justNoise_dig[0][antNum][0], &justNoise_dig[1][antNum][0]);

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
	globalTrigger->volts[0][ilayer][0] = globalTrigger->volts[0][ilayer][0]/sqrt(2.);
      }

      if (fSettingsPtrIDontOwn->RCPRX2ZERO){
	globalTrigger->volts[1][ilayer][1]=0.;
      }

      if (fSettingsPtrIDontOwn->LCPRX2ZERO){
	globalTrigger->volts[0][ilayer][1]=0.;
      }

      if (fSettingsPtrIDontOwn->SIGNAL_FLUCT) {
	if (fSettingsPtrIDontOwn->WHICH==0) {
	  globalTrigger->volts[ilayer][ifold][0]+=gRandom->Gaus(0., this->VNOISE_ANITALITE[ifold]);
	  globalTrigger->volts[ilayer][ifold][1]+=gRandom->Gaus(0., this->VNOISE_ANITALITE[ifold]);
	} //else
      } //if adding noise
      
      if (count_rx==this->rx_minarrivaltime) {
	rec_efield=sqrt(pow(globalTrigger->volts_original[0][ilayer][ifold]/(undogaintoheight_e*0.5), 2)+pow(globalTrigger->volts_original[1][ilayer][ifold]/(undogaintoheight_h*0.5), 2));
	for (int ibw=0;ibw<4;ibw++) {
	  rec_efield_array[ibw]=sqrt(pow(ct.bwslice_volts_pole[ibw]/(undogaintoheight_e_array[ibw]*0.5), 2)+pow(ct.bwslice_volts_polh[ibw]/(undogaintoheight_h_array[ibw]*0.5), 2));
	  bwslice_vnoise_thislayer[ibw]=this->bwslice_vnoise[ilayer][ibw];// this is just for filling into a tree
	} // end loop over bandwidth slices
      } // end if this is the closest antenna

      ct.WhichBandsPass(fSettingsPtrIDontOwn, this, globalTrigger.get(), this, ilayer, ifold, thresholdsAnt[antNum]);
	  
      if (Anita::GetAntennaNumber(ilayer, ifold)==this->rx_minarrivaltime) {
	for (int iband=0;iband<5;iband++) {
	  for (int ipol=0;ipol<2;ipol++) {
	    rx0_signal_eachband[ipol][iband] = ct.signal_eachband[ipol][iband];
	    rx0_threshold_eachband[ipol][iband] = ct.threshold_eachband[ipol][iband];
	    rx0_noise_eachband[ipol][iband] = ct.noise_eachband[ipol][iband];
	    rx0_passes_eachband[ipol][iband] = ct.passes_eachband[ipol][iband];
	  }
	}
      }

      //For verification plots: find antenna with max signal - added by Stephen
      if (ilayer == 0 && globalTrigger->volts[0][ilayer][ifold] > max_antenna_volts0) {
	max_antenna0 = count_rx;
	max_antenna_volts0 = globalTrigger->volts[0][ilayer][ifold];
	max_antenna_volts0_em=globalTrigger->volts_em[0][ilayer][ifold];
	ant_max_normal0 = ant_normal;
	e_comp_max1 = e_component;
	h_comp_max1 = h_component;
      }
      else if (ilayer == 0 && globalTrigger->volts[0][ilayer][ifold] == max_antenna_volts0 && globalTrigger->volts[0][ilayer][ifold] != 0){
	std::cout<<"Equal voltage on two antennas!  Event : "<<inu<<std::endl;
      }
      else if (ilayer == 1 && globalTrigger->volts[0][ilayer][ifold] > max_antenna_volts1) {
	max_antenna1 = count_rx;
	max_antenna_volts1 = globalTrigger->volts[0][ilayer][ifold];
	ant_max_normal1 = ant_normal;
	e_comp_max2 = e_component;
	h_comp_max2 = h_component;
      }
      else if (ilayer == 1 && globalTrigger->volts[0][ilayer][ifold] == max_antenna_volts1 && globalTrigger->volts[0][ilayer][ifold] != 0){
	std::cout<<"Equal voltage on two antennas!  Event : "<<inu<<std::endl;
      }
      else if (ilayer == 2 && globalTrigger->volts[0][ilayer][ifold] > max_antenna_volts2) {
	max_antenna2 = count_rx;
	max_antenna_volts2 = globalTrigger->volts[0][ilayer][ifold];
	ant_max_normal2 = ant_normal;
	e_comp_max3 = e_component;
	h_comp_max3 = h_component;
      }
      else if (ilayer == 2 && globalTrigger->volts[0][ilayer][ifold] == max_antenna_volts2 && globalTrigger->volts[0][ilayer][ifold] != 0){
	std::cout<<"Equal voltage on two antennas!  Event : "<<inu<<std::endl;
      }
      voltagearray[count_rx] = globalTrigger->volts[0][ilayer][ifold];
      //End verification plot block

      count_rx++; // counting antennas that we loop through,  for indexing

      // if (fSettingsPtrIDontOwn->TRIGTYPE==0 && ifold==1 && count_pass>=fSettingsPtrIDontOwn->NFOLD) { //added djg --line below fills "direct" voltage output file
      // 	Log().al_voltages_direct<<"0 0 0"<<"   "<<"    "<<globalTrigger->volts_original[1][0][0]<<"    "<<(globalTrigger->volts_original[0][0][0]/sqrt(2.))<<"     "<<globalTrigger->volts_original[1][0][1]<<"     "<<globalTrigger->volts_original[0][0][1]<<"      "<<this->VNOISE[0]<<"     "<<this->VNOISE[0]<<"     "<<this->VNOISE[0]<<"     "<<this->VNOISE[0]<<"  "<<fNeutrinoPath->weight<<std::endl;
      // }
    } //loop through the phi-fold antennas
  }  //loop through the layers of antennas


  this->rms_rfcm_e_single_event = sqrt(this->rms_rfcm_e_single_event / (this->HALFNFOUR * fSettingsPtrIDontOwn->NANTENNAS));

  for (int irx=0;irx<fSettingsPtrIDontOwn->NANTENNAS;irx++) {
    nchannels_perrx_triggered[irx]=globalTrigger->nchannels_perrx_triggered[irx];
  }

  nchannels_triggered = Tools::iSum(globalTrigger->nchannels_perrx_triggered, fSettingsPtrIDontOwn->NANTENNAS); // find total number of antennas that were triggered.
  fVoltsRX.ave = GetAverageVoltageFromAntennasHit(fSettingsPtrIDontOwn, globalTrigger->nchannels_perrx_triggered, voltagearray, fVoltsRX.sum);

  globalTrigger->PassesTrigger(fSettingsPtrIDontOwn, this, discones_passing, 2, fL3trig, fL2trig, fL1trig, fSettingsPtrIDontOwn->antennaclump, loctrig, loctrig_nadironly, inu, thispasses);

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
    fEventNumber++;
  }

  if(eventPassesTrigger){
    fAnitaOutput.fillRootifiedAnitaDataTrees(*fSettingsPtrIDontOwn, fRayPtrIDontOwn, fScreenPtrIDontOwn);
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
