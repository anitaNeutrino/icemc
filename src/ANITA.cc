#include "ANITA.h"
#include "Settings.h"
#include "GlobalTrigger.h"
#include "ChanTrigger.h"
#include "Tools.h"
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


icemc::ANITA::ANITA(const Settings* settings, const RootOutput* ro)
  : Balloon(settings), Anita(settings, ro->getOutputDir(), this),
    fSettings(settings),
    fNumRX(settings ? settings->NANTENNAS : 48),
    fVoltsRX(settings ? settings->NANTENNAS : 0),
    fAnitaOutput(this, settings, ro->getOutputDir(), ro->getRun())
{
  initSeaveys(settings, this);
}


icemc::ANITA::~ANITA(){
  
}


const Geoid::Position& icemc::ANITA::getPosition(double time){
  /**
   * Here I'm using TMath's QuietNaN() as a default argument.
   * That might be a bad idea, but for now...
   */

  if(TMath::IsNaN(time) && TMath::IsNaN(fLastPositionTime)){
    time = getStartTime();
  }

  if(!TMath::IsNaN(time) && time != fLastPositionTime){

    PickBalloonPosition(time, fSettings, this);

    for(auto& s : fSeaveys){
      s.updatePosition(Balloon::position(),
		       Balloon::getHeading(),
		       Balloon::getPitch(),
		       Balloon::getRoll());
    }

    fLastPositionTime = time;

  }

  return Balloon::position();
}


TVector3 icemc::ANITA::getPositionRX(Int_t rx) const {
  return fSeaveys.at(rx).getPosition(Seavey::Pol::V);
}


void icemc::ANITA::getLayerFoldFromTriggerRX(int rx, int& ilayer, int& ifold) const {
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



void icemc::ANITA::getLayerFoldFromRX(int rx, int& ilayer, int& ifold) const {
  int antNum = rx;  
  // This is NOT how to do things...
  ilayer = -1;
  ifold = -1;
  for (int ilayerTemp=0 ;ilayerTemp < fSettings->NLAYERS; ilayerTemp++) { // loop over layers on the payload
    for (int ifoldTemp=0;ifoldTemp<this->NRX_PHI[ilayerTemp];ifoldTemp++) { // ifold loops over phi
      Int_t antNum2 = this->GetRx(ilayerTemp, ifoldTemp);
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

    TVector3 n_eplane;
    TVector3 n_hplane;
    TVector3 n_normal;
    
    if(settings1->WHICH==Payload::Anita1 ||
       settings1->WHICH==Payload::Anita2 ||
       settings1->WHICH==Payload::Anita3 ||
       settings1->WHICH==Payload::Anita4) { /// @todo presumably this is also correct for ANITA-4

      n_eplane = constants::const_z;
      n_eplane.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
      n_hplane = -constants::const_y;
      n_hplane.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
      n_normal = constants::const_x;
      n_normal.RotateY(anita1->ANTENNA_DOWN[ilayer][ifold]);
    }
    else {
      n_eplane = constants::const_z;
      n_eplane.RotateY(anita1->THETA_ZENITH[ilayer] - constants::PI/2);
      n_hplane = (-constants::const_y);
      n_hplane.RotateY(anita1->THETA_ZENITH[ilayer] - constants::PI/2);
      n_normal = constants::const_x;
      n_normal.RotateY(anita1->THETA_ZENITH[ilayer] - constants::PI/2);
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

    n_eplane.RotateZ(phi);
    n_hplane.RotateZ(phi);
    n_normal.RotateZ(phi);
    
    TVector3 positionH;
    TVector3 positionV;    
    for(auto pol : {Seavey::Pol::V, Seavey::Pol::H}){
      int ipol = static_cast<int>(pol);
      TVector3& seaveyPayloadPos = pol == Seavey::Pol::V ? positionV : positionH;

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
	seaveyPayloadPos = TVector3(anita1->RRX[ilayer]*cos(phi) + anita1->LAYER_HPOSITION[ilayer]*cos(anita1->LAYER_PHIPOSITION[ilayer]),
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

  int discones_passing = 0;  // number of discones that pass
  

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

  // start looping over antennnas.
  // ilayer loops through vertical layers

  if (fSettings->SLAC){
    icemcLog() << icemc::error << "SLAC is no longer supported!" << std::endl;
    exit(1);
    // icemcLog().fslac_hitangles << this->sslacpositions[this->islacposition] << "\n";    
  }
  
  globalTrigger->volts_rx_rfcm_trigger.assign(16,  std::vector <std::vector <double> >(3,  std::vector <double>(0)));

  int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //counting how many pass trigger requirement
  int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX]; //counting how many pass trigger requirement
  double thresholdsAnt[48][2][5] = {{{0}}};

  int count_rx = 0;
  for (int ilayer=0; ilayer < fSettings->NLAYERS; ilayer++) { // loop over layers on the payload
    for (int ifold=0;ifold<this->NRX_PHI[ilayer];ifold++) { // ifold loops over phi

      ChanTrigger ct;
      ct.InitializeEachBand(this);

      // int antNum = this->GetRxTriggerNumbering(ilayer, ifold);
      int antNum = this->GetRx(ilayer, ifold);
      ct.readInSeavey(fSettings,  &fSeaveys.at(antNum), antNum, this);

      // this->GetAntennaOrientation(fSettings,  this,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);
      // ct.ApplyAntennaGain(fSettings, this, fScreenPtrIDontOwn, antNum, n_eplane, n_hplane, n_normal, inu);

      ct.TriggerPath(fSettings, this, antNum, this);
      ct.DigitizerPath(fSettings, this, antNum, this);
      ct.TimeShiftAndSignalFluct(fSettings, this, ilayer, ifold,
				 fVoltsRX.rfcm_lab_e_all.at(count_rx).data(),
				 fVoltsRX.rfcm_lab_h_all.at(count_rx).data());
      ct.saveTriggerWaveforms(&justSignal_trig[0][antNum][0], &justSignal_trig[1][antNum][0], &justNoise_trig[0][antNum][0], &justNoise_trig[1][antNum][0]);
      ct.saveDigitizerWaveforms(&justSignal_dig[0][antNum][0], &justSignal_dig[1][antNum][0], &justNoise_dig[0][antNum][0], &justNoise_dig[1][antNum][0]);

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
      
      ct.WhichBandsPass(fSettings, this, globalTrigger.get(), this, ilayer, ifold, thresholdsAnt[antNum]);
	  
      count_rx++; // counting antennas that we loop through,  for indexing
    } //loop through the phi-fold antennas
  }  //loop through the layers of antennas


  int count_pass = 0;
  globalTrigger->PassesTrigger(fSettings, this, discones_passing, 2, fL3trig, fL2trig, fL1trig, fSettings->antennaclump, loctrig, loctrig_nadironly, inu, thispasses);

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
