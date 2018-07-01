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


icemc::ANITA::ANITA(const Settings* settings, const RayTracer* sillyRay, const RootOutput* ro)
  : Balloon(settings),
    fSettings(settings), fRayPtrIDontOwn(sillyRay),
    fNumRX(settings ? settings->NANTENNAS : 48),
    fVoltsRX(settings ? settings->NANTENNAS : 0),
    fAnitaOutput(this, settings, ro->getOutputDir(), ro->getRun())
{
  initSeaveys(settings, this);
}


icemc::ANITA::~ANITA(){
  
}


icemc::Position icemc::ANITA::getCenterOfDetector(UInt_t unixTime){
  (void) unixTime;

  // UInt_t theUnixTime = unixTime ? 1 : 0;
  // for(int rx=0; rx < static_cast<int>(fSeaveys.size()); rx++){
  //   int layer, fold;
  //   getLayerFoldFromRX(rx, layer, fold);
  // }

  for(auto& s : fSeaveys){
    s.updatePosition(Balloon::position(),
		     Balloon::getHeading(),
		     Balloon::getPitch(),
		     Balloon::getRoll());
  }

  return Balloon::position();
}


icemc::Vector icemc::ANITA::getPositionRX(Int_t rx) const {
  return fSeaveys.at(rx).getPosition(Seavey::Pol::V);
}


void icemc::ANITA::getLayerFoldFromRX(int rx, int& ilayer, int& ifold) const {
  int antNum = rx;
  
  // This is NOT how to do things...
  ilayer = -1;
  ifold = -1;
  for (int ilayerTemp=0 ;ilayerTemp < fSettings->NLAYERS; ilayerTemp++) { // loop over layers on the payload
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


void icemc::ANITA::initSeaveys(const Settings *settings1, const Anita *anita1) {

  for(int rx = 0; rx < getNumRX(); rx++){
    int ilayer = -1;
    int ifold = -1;
    getLayerFoldFromRX(rx, ilayer, ifold);
    std::cout << rx << "\t" << ilayer << "\t" << ifold << std::endl;    

    Vector n_eplane;
    Vector n_hplane;
    Vector n_normal;
    
    if(settings1->WHICH==Payload::Anita1 || settings1->WHICH==Payload::Anita2 || settings1->WHICH==Payload::Anita3) {
      n_eplane = constants::const_z.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
      n_hplane = (-constants::const_y).RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
      n_normal = constants::const_x.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
    }
    else {
      n_eplane = constants::const_z.RotateY(anita1->THETA_ZENITH[ilayer] - constants::PI/2);
      n_hplane = (-constants::const_y).RotateY(anita1->THETA_ZENITH[ilayer] - constants::PI/2);
      n_normal = constants::const_x.RotateY(anita1->THETA_ZENITH[ilayer] - constants::PI/2);
    }

    double phi = 0;
    // rotate about z axis for phi
    if (settings1->CYLINDRICALSYMMETRY==1) {
      phi=(double)ifold/(double)anita1->NRX_PHI[ilayer]*2*constants::PI + anita1->PHI_OFFSET[ilayer]; // + phi_spin;
    }
    else{
      //phi=anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer] +phi_spin;
      phi=anita1->PHI_EACHLAYER[ilayer][ifold];
    }

    n_eplane = n_eplane.RotateZ(phi);
    n_hplane = n_hplane.RotateZ(phi);
    n_normal = n_normal.RotateZ(phi);
    
    Vector positionH;
    Vector positionV;    
    for(auto pol : {Seavey::Pol::V, Seavey::Pol::H}){
      int ipol = static_cast<int>(pol);
      Vector& seaveyPayloadPos = pol == Seavey::Pol::V ? positionV : positionH;

      double phi = 0;
      if (settings1->WHICH == Payload::Anita1 ||
	  settings1->WHICH == Payload::Anita2 ||
	  settings1->WHICH == Payload::Anita3 ||
	  settings1->WHICH == Payload::Anita4 ){
	///@todo Anita needs to be instantiated with these arrays filled!
	/// currently they are done in "apply settings", this needs to change
	seaveyPayloadPos = anita1->ANTENNA_POSITION_START[ipol][ilayer][ifold];
      }
      else {
	if (settings1->CYLINDRICALSYMMETRY==1){ // for timing code
	  // phi is 0 for antenna 0 (0-31) and antenna 16 (0-31)
	  // antenna 1 (1-32) and antenna 18 (1-32)
	  phi = (double) ifold / (double) anita1->NRX_PHI[ilayer] * 2 * constants::PI + anita1->PHI_OFFSET[ilayer];
	}
	else{
	  phi = anita1->PHI_EACHLAYER[ilayer][ifold] + anita1->PHI_OFFSET[ilayer];
	}
	seaveyPayloadPos = Vector(anita1->RRX[ilayer]*cos(phi) + anita1->LAYER_HPOSITION[ilayer]*cos(anita1->LAYER_PHIPOSITION[ilayer]),
			   anita1->RRX[ilayer]*sin(phi) + anita1->LAYER_HPOSITION[ilayer]*sin(anita1->LAYER_PHIPOSITION[ilayer]),
			   anita1->LAYER_VPOSITION[ilayer]);

      }
     
    } 
    fSeaveys.emplace_back(Seavey(positionV, positionH,  n_eplane,  n_hplane, n_normal, settings1));
  }
}




void icemc::ANITA::addSignalToRX(const icemc::PropagatingSignal& signal, int rx, int inu){

  int ifold, ilayer;
  getLayerFoldFromRX(rx, ilayer, ifold);
  
  static bool firstTime = true;
  if(inu == 522 && firstTime){
    for(int i=0; i < fSeaveys.size(); i++){
      // if(true || i==34){
      if(true || i==41){
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
    // this->GetAntennaOrientation(fSettings,  this,  ilayer,  ifold,
    // 				fSeaveys.at(rx).fEPlane, fSeaveys.at(rx).fHPlane, fSeaveys.at(rx).fNormal);

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

  fVoltsRX.reset();
  int thispasses[Anita::NPOL]={0,0};

  auto globalTrigger = std::unique_ptr<GlobalTrigger>(new GlobalTrigger(fSettings, dynamic_cast<Anita*>(this)));
  // make a global trigger object (but don't touch the electric fences)
  // globaltrig1 = new GlobalTrigger(fSettings, anita1);
  // globaltrig1 = new GlobalTrigger(fSettings, fDetector);

  Tools::Zero(this->arrival_times[0], Anita::NLAYERS_MAX*Anita::NPHI_MAX);
  Tools::Zero(this->arrival_times[1], Anita::NLAYERS_MAX*Anita::NPHI_MAX);
  if (!fSettings->TRIGGEREFFSCAN){
    if(fSettings->BORESIGHTS){
      this->GetArrivalTimesBoresights(fRayPtrIDontOwn->n_exit2bn_eachboresight[2]);
    }
    else{
      this->GetArrivalTimes(fRayPtrIDontOwn->n_exit2bn[2],this,fSettings);
    }
  }
  this->rx_minarrivaltime = Tools::WhichIsMin(this->arrival_times[0], fSettings->NANTENNAS);

  if(inu==522){
    for(int rx = 0; rx < fSeaveys.size(); rx++){
      int ifold, ilayer;
      getLayerFoldFromRX(rx, ilayer, ifold);
      int rx2 = GetRx(ilayer,ifold);
      std::cout << "Arrival times " << rx << "\t"
		<< arrival_times[0][rx2] - arrival_times[0][0] << std::endl;//"\t"
		// << arrival_times[1][rx] << std::endl;
    }
 }
  //For verification plots - added by Stephen
  double voltagearray[Anita::NLAYERS_MAX*Anita::NPHI_MAX] = {0}; //Records max voltages on each antenna for one neutrino
  int discones_passing = 0;  // number of discones that pass

  // double max_antenna_volts0 = 0; //Voltage on the antenna with maximum signal,  top layer
  // // double max_antenna_volts0_em = 0; //Component of voltage from em shower on the antenna with maximum signal,  top layer
  // double max_antenna_volts1 = 0; //Voltage on the antenna with maximum signal,  middle layer
  // double max_antenna_volts2 = 0; //Voltage on the antenna with maximum signal,  bottom layer

  // double e_comp_max1=0;
  // double h_comp_max1=0;
  // double e_comp_max2=0;
  // double h_comp_max2=0;
  // double e_comp_max3=0;
  // double h_comp_max3=0;

  // double hitangle_e = 0; // angle the ray hits the antenna wrt e-plane

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
  


  // start looping over antennnas.
  // ilayer loops through vertical layers

  if (fSettings->SLAC){
    icemcLog() << icemc::error << "SLAC is no longer supported!" << std::endl;
    exit(1);
    // icemcLog().fslac_hitangles << this->sslacpositions[this->islacposition] << "\n";    
  }

  double hitangle_e;
  if (fSettings->CENTER){
    this->CenterPayload(hitangle_e);
  }
  
  globalTrigger->volts_rx_rfcm_trigger.assign(16,  vector <vector <double> >(3,  vector <double>(0)));
  this->rms_rfcm_e_single_event = 0;

  double bwslice_vnoise_thislayer[4];// for filling tree6b,  noise for each bandwidth on each layer  
  double rx0_signal_eachband[2][5];
  double rx0_threshold_eachband[2][5];
  double rx0_noise_eachband[2][5];
  int rx0_passes_eachband[2][5];
  int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //counting how many pass trigger requirement
  int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX]; //counting how many pass trigger requirement
  double thresholdsAnt[48][2][5] = {{{0}}};

  int count_rx = 0;
  for (int ilayer=0; ilayer < fSettings->NLAYERS; ilayer++) { // loop over layers on the payload
    for (int ifold=0;ifold<this->NRX_PHI[ilayer];ifold++) { // ifold loops over phi

      ChanTrigger ct;
      ct.InitializeEachBand(this);

      int antNum = this->GetRxTriggerNumbering(ilayer, ifold);
      ct.readInSeavey(fSettings,  &fSeaveys.at(antNum), antNum, this, inu);

      // this->GetAntennaOrientation(fSettings,  this,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);
      // ct.ApplyAntennaGain(fSettings, this, fScreenPtrIDontOwn, antNum, n_eplane, n_hplane, n_normal, inu);

      ct.TriggerPath(fSettings, this, antNum, this);
      ct.DigitizerPath(fSettings, this, antNum, this);
      ct.TimeShiftAndSignalFluct(fSettings, this, ilayer, ifold,
				 fVoltsRX.rfcm_lab_e_all.at(count_rx).data(),
				 fVoltsRX.rfcm_lab_h_all.at(count_rx).data(), inu);
      ct.saveTriggerWaveforms(&justSignal_trig[0][antNum][0], &justSignal_trig[1][antNum][0], &justNoise_trig[0][antNum][0], &justNoise_trig[1][antNum][0]);
      ct.saveDigitizerWaveforms(&justSignal_dig[0][antNum][0], &justSignal_dig[1][antNum][0], &justNoise_dig[0][antNum][0], &justNoise_dig[1][antNum][0]);

      if(inu==522){
	int j = TMath::LocMax(Anita::HALFNFOUR, fVoltsRX.rfcm_lab_e_all.at(count_rx).data());
	std::cout << Anita::HALFNFOUR << "\t" << count_rx << "\t" << j << "\t" <<  1000*fVoltsRX.rfcm_lab_e_all[count_rx][j] << "\n";
      }

      Tools::Zero(sumsignal, 5);

      if (this->whichPath()==FlightPath::PeterEvent && ilayer==this->GetLayer(this->rx_minarrivaltime) && ifold==this->GetIfold(this->rx_minarrivaltime)) {
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
      if (fSettings->SCALEDOWNLCPRX1){
	globalTrigger->volts[0][ilayer][0] = globalTrigger->volts[0][ilayer][0]/sqrt(2.);
      }

      if (fSettings->RCPRX2ZERO){
	globalTrigger->volts[1][ilayer][1]=0.;
      }

      if (fSettings->LCPRX2ZERO){
	globalTrigger->volts[0][ilayer][1]=0.;
      }

      if (fSettings->SIGNAL_FLUCT) {
	if (fSettings->WHICH==Payload::AnitaLite) {
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

      ct.WhichBandsPass(fSettings, this, globalTrigger.get(), this, ilayer, ifold, thresholdsAnt[antNum]);
	  
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

      voltagearray[count_rx] = globalTrigger->volts[0][ilayer][ifold];
      //End verification plot block

      count_rx++; // counting antennas that we loop through,  for indexing

      // if (fSettings->TRIGTYPE==0 && ifold==1 && count_pass>=fSettings->NFOLD) { //added djg --line below fills "direct" voltage output file
      // 	Log().al_voltages_direct<<"0 0 0"<<"   "<<"    "<<globalTrigger->volts_original[1][0][0]<<"    "<<(globalTrigger->volts_original[0][0][0]/sqrt(2.))<<"     "<<globalTrigger->volts_original[1][0][1]<<"     "<<globalTrigger->volts_original[0][0][1]<<"      "<<this->VNOISE[0]<<"     "<<this->VNOISE[0]<<"     "<<this->VNOISE[0]<<"     "<<this->VNOISE[0]<<"  "<<fNeutrinoPath->weight<<std::endl;
      // }
    } //loop through the phi-fold antennas
  }  //loop through the layers of antennas


  this->rms_rfcm_e_single_event = sqrt(this->rms_rfcm_e_single_event / (this->HALFNFOUR * fSettings->NANTENNAS));

  // for (int irx=0;irx<fSettings->NANTENNAS;irx++) {
  //   nchannels_perrx_triggered[irx]=globalTrigger->nchannels_perrx_triggered[irx];
  // }

  // nchannels_triggered = Tools::iSum(globalTrigger->nchannels_perrx_triggered, fSettings->NANTENNAS); // find total number of antennas that were triggered.

  fVoltsRX.ave = GetAverageVoltageFromAntennasHit(fSettings, globalTrigger->nchannels_perrx_triggered, voltagearray, fVoltsRX.sum);

  globalTrigger->PassesTrigger(fSettings, this, discones_passing, 2, fL3trig, fL2trig, fL1trig, fSettings->antennaclump, loctrig, loctrig_nadironly, inu, thispasses);

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
       || (fSettings->TRIGTYPE==0 && count_pass>=fSettings->NFOLD)
       || (fSettings->MINBIAS==1)){
    eventPassesTrigger = true;
    fEventNumber++;
  }

  if(eventPassesTrigger){
    fAnitaOutput.fillRootifiedAnitaDataTrees();
  }

  return eventPassesTrigger;
}





double icemc::ANITA::GetAverageVoltageFromAntennasHit(const Settings *settings1, int *nchannels_perrx_triggered, const double *voltagearray, double& volts_rx_sum) const {
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
