#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include "vector.hh"
#include "position.hh"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLine.h"
#include "TFile.h"

#include "rx.h"
#include "Constants.h"
#include "anita.hh"
#include "balloon.hh"
#include "TRandom3.h"
#include "ChanTrigger.h"
#include <cmath>
#include "Tools.h"
#include "Settings.h"
#include "screen.hh"
#include "GlobalTrigger.h"

using std::cout;

#ifdef ANITA_UTIL_EXISTS
#include "FFTtools.h"
#include "SimulatedSignal.h"
#endif




void icemc::ChanTrigger::ConvertEHtoLREfield(double e_component,double h_component,double& lcp_component,double& rcp_component) {
    
  lcp_component=sqrt((e_component*e_component+h_component*h_component)/2);
  rcp_component=lcp_component;
    
} //ConvertEHtoLREfield
void icemc::ChanTrigger::ConvertEHtoLREnergy(double e_component,double h_component,double& lcp_component,double& rcp_component) {
    
  lcp_component=(e_component+h_component)/2;
  rcp_component=lcp_component;
    
} //ConvertEHtoLREnergy

void icemc::ChanTrigger::ConvertHVtoLRTimedomain(const int nfour,double *vvolts,
					  double *hvolts,
					  double *left,double *right) {
    
  // nfour is the real and complex values from -F to F
    
  // first perform fft on each of h and v
  // find l, r polarizations
  // take fft back
    
  double hvolts_f[nfour/2];
  double vvolts_f[nfour/2];
  for (int i=0;i<nfour/2;i++) {
    hvolts_f[i]=hvolts[i];
    vvolts_f[i]=vvolts[i];
  }
    
  //  double *hvolts_f=hvolts;
  //double *vvolts_f=vvolts;
    
  Tools::realft(hvolts_f,1,nfour/2);
  Tools::realft(vvolts_f,1,nfour/2);
    
  for (int i=0;i<nfour/4;i++) {
    //right[2*i]=1/sqrt(2.)*(hvolts_f[2*i]+vvolts_f[2*i+1]); //This is what was being done, until Jacob declared it wrong
    right[2*i]=1/sqrt(2.)*(vvolts_f[2*i]-hvolts_f[2*i+1]); //The thing Jacob declared right
    left[2*i]=1/sqrt(2.)*(hvolts_f[2*i]-vvolts_f[2*i+1]);
		
    //right[2*i+1]=1/sqrt(2.)*(hvolts_f[2*i+1]-vvolts_f[2*i]); //This is what was being done, until Jacob declared it wrong
    right[2*i+1]=1/sqrt(2.)*(vvolts_f[2*i+1]+hvolts_f[2*i]); //The thing Jacob declared right
    left[2*i+1]=1/sqrt(2.)*(hvolts_f[2*i+1]+vvolts_f[2*i]);
		

    left[2*i]=left[2*i]*2./((double)nfour/2.);
    left[2*i+1]=left[2*i+1]*2./((double)nfour/2.);
		
    right[2*i]=right[2*i]*2./((double)nfour/2.);
    right[2*i+1]=right[2*i+1]*2./((double)nfour/2.);
  }
    
  Tools::realft(left,-1,nfour/2);
  Tools::realft(right,-1,nfour/2);
    
  // now take fft back
    
}


void icemc::ChanTrigger::WhichBandsPass(const Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double thresholds[2][5]){
    

  if (settings1->USETIMEDEPENDENTTHRESHOLDS==1 && settings1->WHICH==9) {
    for(int i=0;i<4;i++) thresholds[0][i] = thresholds[1][i] = anita1->powerthreshold[i];
    int iring = (ilayer<2)*0 + (ilayer==2)*1 + (ilayer==3)*2;
    int iphi = ifold;
    if (ilayer==0) iphi = ifold*2;
    else if (ilayer==1) iphi = ifold*2+1;
    // we invert VPOL and HPOL because in AnitaEventReader (0:HPOL 1:VPOL) and in icemc (0:VPOL 1:HPOL)
    // convert scalers to power thresholds using fitted function got from ANITA-1, are they too old?
    // thresholds[1][4] = ADCCountstoPowerThreshold(anita1,0,iring*16+iphi)*(-1.);
    // thresholds[0][4] = ADCCountstoPowerThreshold(anita1,1,iring*16+iphi)*(-1.);


    // Use thresholds converted from flight scalers
    thresholds[0][4] = anita1->fakeThresholds2[0][iring*16+iphi]*(-1.);
    thresholds[1][4] = anita1->fakeThresholds2[1][iring*16+iphi]*(-1.);

    //    cout << thresholds[0][4] << " " <<  thresholds[0][4] << " \n";
  } else {
    GetThresholds(settings1,anita1,ilayer,thresholds); // get the right thresholds for this layer
  }
  globaltrig1->volts[0][ilayer][ifold]=0.;
  globaltrig1->volts[1][ilayer][ifold]=0.;
    
  //  TRandom3 Rand3;
    
  // This is the old way used for ANITA 1 where we integrate in frequency
  // to get an estimate of the signal strength
  if (settings1->TRIGGERSCHEME <= 1){
    WhichBandsPassTrigger1(settings1, anita1, globaltrig1, bn1, ilayer, ifold, thresholds);
   
  } else  if (settings1->TRIGGERSCHEME >= 2){
    // this scheme is used for ANITA 2 on.
    //cout << "i'm here.\n";
    WhichBandsPassTrigger2(settings1, anita1, globaltrig1, bn1, ilayer, ifold, thresholds);

  }
} // end which bands pass


//!
/*!
 *
 *
 *
 *
 *
 */
void icemc::ChanTrigger::WhichBandsPassTrigger1(const Settings *settings1, Anita *anita1, icemc::GlobalTrigger *globaltrig1, icemc::Balloon *bn1, int ilayer, int ifold, double thresholds[2][5]){

  double volts_thischannel;
  double energy_thischannel;
  double voltagethresh_thischannel,energythresh_thischannel;

  // add noise, then find lcp, rcp components
  for (int ibw=0;ibw<anita1->NBANDS+1;ibw++) {
	
    //     cout << "ibw, bwslice_volts_pole are " << ibw << " " << bwslice_volts_pole[ibw] << "\n";
	
    if (settings1->SIGNAL_FLUCT) {// add noise fluctuations if requested
      bwslice_volts_pole[ibw] += gRandom->Gaus(0.,anita1->bwslice_vnoise[ilayer][ibw]);
      bwslice_volts_polh[ibw] += gRandom->Gaus(0.,anita1->bwslice_vnoise[ilayer][ibw]);
    }
	
    // Convert e-plane and h-plane components into lcp,rcp
    ConvertEHtoLREfield(bwslice_volts_pole[ibw],bwslice_volts_polh[ibw],bwslice_volts_pol0[ibw],bwslice_volts_pol1[ibw]);
	
    // do the same for energy
    ConvertEHtoLREnergy(bwslice_energy_pole[ibw],bwslice_energy_polh[ibw],bwslice_energy_pol0[ibw],bwslice_energy_pol1[ibw]);
	
    globaltrig1->volts[0][ilayer][ifold]+=bwslice_volts_pol0[ibw];
    globaltrig1->volts[1][ilayer][ifold]+=bwslice_volts_pol1[ibw];
  }
      
  for (int ibw=0;ibw<anita1->NBANDS+1;ibw++) { // subbands+full band
	
    // first lcp or v pol (here called pole)
    // if we're implementing masking and the channel has been masked
    if (settings1->CHMASKING && !icemc::ChanTrigger::IsItUnmasked(bn1->surfTrigBandMask,ibw,ilayer,ifold,0)) {
	  
      globaltrig1->channels_passing[ilayer][ifold][0][ibw]=0;// channel does not pass
      globaltrig1->vchannels_passing[ilayer][ifold][0][ibw]=0;// channel does not pass
    }
    // note the last element of the array is 0 because this is lcp or e pol
    else {
	  
	  
      if (settings1->LCPRCP )  {// if we're considering lcp, rcp
	volts_thischannel=bwslice_volts_pol0[ibw];
	energy_thischannel=bwslice_energy_pol0[ibw];
      }
      else {// if we're considering v and h pol
	volts_thischannel=bwslice_volts_pole[ibw];
	energy_thischannel=bwslice_energy_pole[ibw];

      }
	  
	  
	  
      //     isurf=Anita::AntennaNumbertoSurfNumber(ilayer,ifold)-1;
      //     ichan=Anita::GetSurfChannel(Anita::GetAntennaNumber(ilayer,ifold),ibw,AnitaPol::kLeft)-1; // for vertical trigger, arbitrarily use threshold scans from left polarisation
      // get the treshold in adc counts for this channel.
      //	     powerthresh_thischannel=bn1->powerthresh[isurf][ichan];
	  
      energythresh_thischannel=thresholds[0][ibw];
      //meanp_thischannel=bn1->meanp[isurf][ichan];
      //energy_thischannel+=meanp_thischannel*INTEGRATION_TIME;
	  
      voltagethresh_thischannel=anita1->bwslice_thresholds[ibw];
	  
	  
      if (settings1->ZEROSIGNAL) {
	volts_thischannel=0.;
	energy_thischannel=0.;
      }
	  
      if (settings1->TRIGGERSCHEME==0) {
	signal_eachband[0][ibw]=(double)fabs(volts_thischannel);
	threshold_eachband[0][ibw]=voltagethresh_thischannel;
	noise_eachband[0][ibw]=anita1->bwslice_vnoise[ilayer][ibw];
	    
	vsignal_eachband[0][ibw]=(double)fabs(volts_thischannel);
	vthreshold_eachband[0][ibw]=voltagethresh_thischannel;
	vnoise_eachband[0][ibw]=anita1->bwslice_vnoise[ilayer][ibw];
	    
	    
      }
      if (settings1->TRIGGERSCHEME==1) {
	signal_eachband[0][ibw]=energy_thischannel;
	threshold_eachband[0][ibw]=energythresh_thischannel;
	noise_eachband[0][ibw]=anita1->bwslice_enoise[ibw];
	    
	vsignal_eachband[0][ibw]=energy_thischannel;
	vthreshold_eachband[0][ibw]=energythresh_thischannel;
	vnoise_eachband[0][ibw]=anita1->bwslice_enoise[ibw];
      }
	  
      // compare the signal to noise and see if it passes
      //cout << "ibw, signal, noise, threshold are " << ibw << " " << signal_eachband[0][ibw] << " " << noise_eachband[0][ibw] << " " << threshold_eachband[0][ibw] << "\n";
      if (signal_eachband[0][ibw]/noise_eachband[0][ibw]>=threshold_eachband[0][ibw]) {
	//      if (fabs(volts_thischannel)/anita1->bwslice_vnoise[ilayer][ibw]>=bwslice_thresholds[ibw]) {
	if (anita1->pol_allowed[0] && anita1->bwslice_allowed[ibw]) { // are this pol and bandwidth allowed to pass
	  globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]++; //Records number of first level triggers on each antenna for a single neutrino
	  globaltrig1->channels_passing[ilayer][ifold][0][ibw]=1; // if it does then flag the element of the channels_passing array that corresponds to lcp or v pol for this band
	  globaltrig1->vchannels_passing[ilayer][ifold][0][ibw]=1; // if it does then flag the element of the channels_passing array that corresponds to lcp or v pol for this band
	  passes_eachband[0][ibw]=1;
	  vpasses_eachband[0][ibw]=1;
	      
	} // end if this polarization and bandwidth are allowed to pass
      } // does this channel pass
	  
	  
    } // end of the else
	
	
	
	
    // now rcp or h pol (here called polh)
    // if we're implementing masking and the channel has been masked
    if (!settings1->JUSTVPOL) { // if we're considering just vpol then there's not a second polarization to consider
      if (settings1->CHMASKING && !icemc::ChanTrigger::IsItUnmasked(bn1->surfTrigBandMask,ibw,ilayer,ifold,1)) {
	globaltrig1->channels_passing[ilayer][ifold][1][ibw]=0;// channel does not pass
	globaltrig1->vchannels_passing[ilayer][ifold][1][ibw]=0;// channel does not pass
	    
	    
	// note the last element of the array is 1 because this is rcp or h pol
	passes_eachband[1][ibw]=0;
      }
      else {
	    
	if (settings1->LCPRCP) {  // if we're considering lcp, rcp
	  volts_thischannel=bwslice_volts_pol1[ibw];
	  energy_thischannel=bwslice_energy_pol1[ibw];
	}
	else { // if we're considering just e and h
	  volts_thischannel=bwslice_volts_polh[ibw];
	  energy_thischannel=bwslice_energy_polh[ibw];
	}
	    
	energythresh_thischannel=thresholds[1][ibw];
	//meanp_thischannel=bn1->meanp[isurf][ichan];
	    
	//energy_thischannel+=meanp_thischannel*INTEGRATION_TIME;
	    
	voltagethresh_thischannel=anita1->bwslice_thresholds[ibw];
	    
	    
	if (settings1->ZEROSIGNAL) {
	  volts_thischannel=0.;
	  energy_thischannel=0.;
	}
	    
	    
	if (settings1->TRIGGERSCHEME==0) {
	  signal_eachband[1][ibw]=(double)fabs(volts_thischannel);
	  threshold_eachband[1][ibw]=voltagethresh_thischannel;
	  noise_eachband[1][ibw]=anita1->bwslice_vnoise[ilayer][ibw];
	}
	if (settings1->TRIGGERSCHEME==1) {
	  signal_eachband[1][ibw]=energy_thischannel;
	  threshold_eachband[1][ibw]=energythresh_thischannel;
	  noise_eachband[1][ibw]=anita1->bwslice_enoise[ibw];
	}
	    

	// compare the signal to noise and see if it passes
	if (signal_eachband[1][ibw]/noise_eachband[1][ibw]>=threshold_eachband[1][ibw]) {
	  if (anita1->pol_allowed[1] && anita1->bwslice_allowed[ibw]) {
	    globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]++; //Records number of first level triggers on each antenna for a single neutrino
	    globaltrig1->channels_passing[ilayer][ifold][1][ibw]=1; // if it does then flag the element of the channels_passing array that corresponds to rcp for this band
	    globaltrig1->vchannels_passing[ilayer][ifold][1][ibw]=1; // if it does then flag the element of the channels_passing array that corresponds to rcp for this band
	    passes_eachband[1][ibw]=1;
	  } // is this band and polarization allowed to pass
	} //does this channel pass
	    
      } // end of the else (channel isn't masked)
	  
    } // if not justvpol
  } // end loop over bands

} // end WhichBandsPassTrigger1

//!
/*!
 *
 *
 *
 *
 *
 */
void icemc::ChanTrigger::WhichBandsPassTrigger2(const Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double thresholds[2][5]){
  
  double psignal[2][5][Anita::NFOUR];
  
  double mindiodeconvl[2][5];

  double onediodeconvl[2][5];
  
  double timedomain_output[2][5][Anita::NFOUR];
  int iant=anita1->GetRxTriggerNumbering(ilayer, ifold);
  int ipol=0;

  if (settings1->NOISEFROMFLIGHTTRIGGER){
    anita1->bwslice_rmsdiode[4] = anita1->bwslice_dioderms_fullband_allchan[ipol][iant][anita1->tuffIndex];
  }
  
  // if we use the diode to perform an integral
  // this is the number of bins to the left of center where the diode function starts to be completely overlapping with the waveform in the convolution.
  int ibinshift=(anita1->NFOUR/4-(int)(anita1->maxt_diode/anita1->TIMESTEP));
  
  // now we have converted the signal to time domain waveforms for all the bands of the antenna
      
  double integrateenergy[5]={0.,0.,0.,0.,0.};
        
  for (int iband=0;iband<5;iband++) { // Only loop over allowed bands
    if (anita1->bwslice_allowed[iband]!=1) continue; 
  
    if (settings1->TRIGGERSCHEME == 2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5){
      for (int itime=0;itime<anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP);itime++) {
  	anita1->maxbin_fortotal[iband]=anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP);
  	
  	// The below line seems to shorten the length of the waveform to less than HALFNFOUR
  	int itimenoisebin=anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP)-itime;
  	
  	// this is just the straight sum of the two
  	anita1->total_vpol_inanita[iband][itime]=anita1->timedomainnoise_rfcm_banding[0][iband][itime]+anita1->signal_vpol_inanita[iband][itime];

  	integrateenergy[iband]+=anita1->timedomainnoise_rfcm_banding[0][iband][itime]*anita1->timedomainnoise_rfcm_banding[0][iband][itime]*anita1->TIMESTEP;
   	if ( settings1->SIGNAL_FLUCT && (!settings1->NOISEFROMFLIGHTTRIGGER) ) {
  	  // this reverses the noise is time, and starts with bin anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP)
	  justNoise_trigPath[0][itime] = anita1->timedomainnoise_rfcm_banding[0][iband][itimenoisebin];
	  justNoise_trigPath[1][itime] = anita1->timedomainnoise_rfcm_banding[1][iband][itimenoisebin];
  	  v_banding_rfcm_forfft[0][iband][itime]=v_banding_rfcm_forfft[0][iband][itime]+anita1->timedomainnoise_rfcm_banding[0][iband][itimenoisebin];
  	  v_banding_rfcm_forfft[1][iband][itime]=v_banding_rfcm_forfft[1][iband][itime]+anita1->timedomainnoise_rfcm_banding[1][iband][itimenoisebin];
  	  }
      }
      for (int itime=anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP);itime<anita1->NFOUR/2;itime++) {
  	anita1->total_vpol_inanita[iband][itime]=0.;
  	v_banding_rfcm_forfft[0][iband][itime]=0.;
  	v_banding_rfcm_forfft[1][iband][itime]=0.;
      }
    }
    else if ( settings1->SIGNAL_FLUCT && (!settings1->NOISEFROMFLIGHTTRIGGER) ) {
      for (unsigned int itime = 0; itime < anita1->HALFNFOUR; ++itime){
  	// this is just a straight sum
  	anita1->total_vpol_inanita[iband][itime]=anita1->timedomainnoise_rfcm_banding[0][iband][itime]+anita1->signal_vpol_inanita[iband][itime];
  	integrateenergy[iband]+=anita1->timedomainnoise_rfcm_banding[0][iband][itime]*anita1->timedomainnoise_rfcm_banding[0][iband][itime]*anita1->TIMESTEP;
  	// this one is the one actually used by the diode
  	v_banding_rfcm_forfft[0][iband][itime] += anita1->timedomainnoise_rfcm_banding[0][iband][itime];
      }
    }
  } // end loop over bands	
  
  // Translate Anita physical layer to Anita trigger layer and phi sector
  // (4 layers with 8,8,16,8 phi sector to 3 layers with 16 phi sectors each.
  // In the nadir layer, 8 of the 16 slots are empty on the trigger layer.)
  int whichlayer,whichphisector;
  globaltrig1->GetAnitaLayerPhiSector(settings1,ilayer,ifold,whichlayer,whichphisector); 
    
  for (int iband=0;iband<5;iband++) {
    if (anita1->bwslice_allowed[iband]!=1) continue; 

    if (settings1->LCPRCP 
	// &&	!(settings1->WHICH==10 && !globaltrig1->WHICHLAYERSLCPRCP[whichlayer])
	) {
      
      // Convert Horiz and Vert polarization
      // To Left and Right circular polarization
      ConvertHVtoLRTimedomain(anita1->NFOUR, v_banding_rfcm_forfft[0][iband], v_banding_rfcm_forfft[1][iband], vm_banding_rfcm_forfft[0][iband], vm_banding_rfcm_forfft[1][iband]);
      
    } else {
  
      for (int itime=0;itime<anita1->NFOUR/2;itime++) {
	    
	vm_banding_rfcm_forfft[0][iband][itime] = v_banding_rfcm_forfft[0][iband][itime];
	vm_banding_rfcm_forfft[1][iband][itime] = v_banding_rfcm_forfft[1][iband][itime];
	    
      }
    }

    for (int itime=0;itime<anita1->NFOUR/2;itime++) {
      anita1->total_diodeinput_1_inanita[iband][itime] = vm_banding_rfcm_forfft[0][iband][itime];
      anita1->total_diodeinput_2_inanita[iband][itime] = vm_banding_rfcm_forfft[1][iband][itime];
	  
    }
    //volts_fullband_trigger_path_e.push_back(;
  } // end loop over bands
      
            
  for (int itime=0;itime<Anita::NFOUR/2;itime++) {
    anita1->total_diodeinput_1_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)][itime]=anita1->total_diodeinput_1_inanita[4][itime];
    anita1->total_diodeinput_2_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)][itime]=anita1->total_diodeinput_2_inanita[4][itime];
  }

  DiodeConvolution(settings1, anita1, globaltrig1, ilayer, ifold, mindiodeconvl[0], onediodeconvl[0], psignal[0], timedomain_output[0], ibinshift, 0, thresholds);
  DiodeConvolution(settings1, anita1, globaltrig1, ilayer, ifold, mindiodeconvl[1], onediodeconvl[1], psignal[1], timedomain_output[1], ibinshift, 1, thresholds);

  // fill channels_passing
  //  } // end if the signal is big enough the be considered
  // now we've made the diode outputs
  // now step in time and find the most number of bands that pass in the
  // right time window
  // then fills channels_passing
  int npass;

  L1Trigger(anita1,timedomain_output[0],timedomain_output[1],thresholds, //inputs
	    globaltrig1->channels_passing[ilayer][ifold][0],globaltrig1->channels_passing[ilayer][ifold][1],npass); //outputs

  //if (npass==1) std::cout << "L1 trigger " << ilayer << " " << ifold << " " << npass << std::endl;

  // if it's the closest antenna,
  // save flag_e,h in anita class for writing to tsignals tree
  int startbin=TMath::MinElement(5,anita1->iminbin);
      
  if (ilayer==anita1->GetLayer(anita1->rx_minarrivaltime) && ifold==anita1->GetIfold(anita1->rx_minarrivaltime)) {
    for (int iband=0;iband<5;iband++) {
      if (anita1->bwslice_allowed[iband]!=1) continue; 
      // cout << "zeroeing here 1.\n";
      anita1->ston[iband]=0.;
      for (int i=anita1->iminbin[iband];i<anita1->imaxbin[iband];i++) {
	// 	    if (iband==0 && i==anita1->NFOUR/4) {
	// 	      cout << "output is " << anita1->inu << "\t" << timedomain_output[0][iband][i] << "\n";
	// 	    }
	// cout << "output, bwslice_rmsdiode are " << timedomain_output[0][iband][i] << "\t" << anita1->bwslice_rmsdiode[iband] << "\n";
	if (timedomain_output[0][iband][i]/anita1->bwslice_rmsdiode[iband]<anita1->ston[iband]) {
	  anita1->ston[iband]=timedomain_output[0][iband][i]/anita1->bwslice_rmsdiode[iband];
	  // if (iband==4 && anita1->ston[iband]<0.)
	  // cout << "ston is " << anita1->ston[iband] << "\n";
	}
      }


      for (int i=0;i<anita1->HALFNFOUR;i++) {
	anita1->flag_e_inanita[iband][i]=0;
	anita1->flag_h_inanita[iband][i]=0;
	anita1->timedomain_output_inanita[0][iband][i]=timedomain_output[0][iband][i];
	anita1->timedomain_output_inanita[1][iband][i]=timedomain_output[1][iband][i];

      }
      for (int i=0;i<(int)flag_e[iband].size();i++) {
	anita1->flag_e_inanita[iband][i+startbin]=flag_e[iband][i];
      }
      for (int i=0;i<(int)flag_h[iband].size();i++) {
	anita1->flag_h_inanita[iband][i+startbin]=flag_h[iband][i];
	    
      }
    }

  }
      
      
      
  for (int iband=0;iband<5;iband++) {
    if (anita1->bwslice_allowed[iband]!=1) continue; 
    // this is for the e polarization
    // if we're implementing masking and the channel has been masked
    if (settings1->CHMASKING && !icemc::ChanTrigger::IsItUnmasked(bn1->surfTrigBandMask,iband,ilayer,ifold,0)) {
      if (globaltrig1->channels_passing[ilayer][ifold][0][iband])
	globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]--;
      globaltrig1->channels_passing[ilayer][ifold][0][iband]=0;// channel does not pass
    }
	
    // if we're implementing masking and the channel has been masked
    if (settings1->CHMASKING && !icemc::ChanTrigger::IsItUnmasked(bn1->surfTrigBandMask,iband,ilayer,ifold,1)) {
      if (globaltrig1->channels_passing[ilayer][ifold][1][iband])
	globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]--;
      globaltrig1->channels_passing[ilayer][ifold][1][iband]=0;// channel does not pass
      // note the last element of the array is 0 because this is lcp or e pol
    }
  } // end loop over first four bands
      
  // this l1_passing is only used for output to tsignals file, but really need to re-evaluate npass after the masking has been imposed before deciding whether it passes an L1 trigger
  if (npass>=anita1->trigRequirements[0]) {
    anita1->l1_passing=1;
    anita1->l1_passing_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)]=1;
  } else {
    anita1->l1_passing=0;
    anita1->l1_passing_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)]=0;
  }

      
  anita1->irx=anita1->GetRx(ilayer,ifold);
      
  if (Anita::GetLayer(anita1->rx_minarrivaltime)==ilayer && Anita::GetIfold(anita1->rx_minarrivaltime)==ifold && anita1->tsignals->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1) {
    anita1->tsignals->Fill();
  }
      
      
  if (settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5){
    // If TRIGGERSCHEME == 3, allow all bands to pass here. This allows the coherent waveform sum trigger scheme to be used.
    for (unsigned int ichannel = 0; ichannel < 5; ichannel++){
      globaltrig1->channels_passing[ilayer][ifold][0][ichannel] = 1;
      globaltrig1->channels_passing[ilayer][ifold][1][ichannel] = 1;
      anita1->channels_passing[0][ichannel] = 1;
      anita1->channels_passing[1][ichannel] = 1;
    }
	
    unsigned int ilayer_trigger;
    unsigned int iphisector_trigger;
	
    if (ilayer <= 1){
      iphisector_trigger = ifold * 2 + ilayer;
      ilayer_trigger = 0;
    }else{
      iphisector_trigger = ifold;
      ilayer_trigger = ilayer - 1;
    }
    globaltrig1->volts_rx_rfcm_trigger[iphisector_trigger][ilayer_trigger].assign(anita1->total_diodeinput_1_inanita[4], anita1->total_diodeinput_1_inanita[4] + 512);
  }

  
  
 
  for (int i=0;i<2;i++) {

    for (unsigned int ibin=0;ibin<globaltrig1->arrayofhits[whichlayer][whichphisector][i][4].size();ibin++) {
      anita1->arrayofhits_inanita[whichlayer][whichphisector][i][ibin]=globaltrig1->arrayofhits[whichlayer][whichphisector][i][4][ibin];
    }
    for (unsigned int ibin=globaltrig1->arrayofhits[whichlayer][whichphisector][i][4].size();ibin<HALFNFOUR;ibin++) {
      anita1->arrayofhits_inanita[whichlayer][whichphisector][i][ibin]=0.;
    }
    //       if (iband==4 && i==0)
    // 	cout << "whichlayer, whichphisector, size of arrayofhits is " << whichlayer << "\t" << whichphisector << "\t" << anita1->arrayofhits_inanita[whichlayer][whichphisector][i][iband].size() << "\n";
  }
 
 



}// end WhichBandsPassTrigger2



void icemc::ChanTrigger::DiodeConvolution(const Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, int ilayer, int ifold, double mindiodeconvl[5], double onediodeconvl[5], double psignal[5][Anita::NFOUR],  double timedomain_output[5][Anita::NFOUR], int ibinshift, int ipol, double thresholds[2][5]){

  int tempChansPassing[5]={0,0,0,0,0};

  // Translate Anita physical layer to Anita trigger layer and phi sector
  // (4 layers with 8,8,16,8 phi sector to 3 layers with 16 phi sectors each.
  // In the nadir layer, 8 of the 16 slots are empty on the trigger layer.)
  int whichlayer,whichphisector;
  globaltrig1->GetAnitaLayerPhiSector(settings1,ilayer,ifold,whichlayer,whichphisector); 


  
  for (int iband=0;iband<5;iband++) {
    
    anita1->channels_passing[ipol][iband]=0;
    if (anita1->bwslice_allowed[iband]!=1) continue; 
    
    // myconvl
    // this performs the convolution with the diode response
    anita1->myconvlv(vm_banding_rfcm_forfft[ipol][iband],anita1->NFOUR,anita1->fdiode_real[iband],mindiodeconvl[iband],onediodeconvl[iband],psignal[iband],timedomain_output[iband]);
    // loop from the ibinshift left + some delay + 10 ns
    

    // // TEMPORARY CODE FOR PRINTING DIODE INPUT/OUTPUT
    // if (anita1->inu==1){
    //   TCanvas *c = new TCanvas("c");
    //   string name = "newnoise";
    //   TGraph *gV = new TGraph(anita1->HALFNFOUR, anita1->fTimes, vm_banding_rfcm_forfft[ipol][iband]);
    //   gV->SetTitle("Diode Input;Time [s];Amplitude [V]");
    //   gV->Draw("Al");
    //   c->Print(Form("%s_%d_%d_%d_%d_diodeInput.png",  name.c_str(), anita1->inu, ipol, ilayer, ifold));
    //   c->Print(Form("%s_%d_%d_%d_%d_diodeInput.pdf",  name.c_str(), anita1->inu, ipol, ilayer, ifold));
    //   c->Print(Form("%s_%d_%d_%d_%d_diodeInput.C",    name.c_str(), anita1->inu, ipol, ilayer, ifold));
    //   c->Print(Form("%s_%d_%d_%d_%d_diodeInput.root", name.c_str(), anita1->inu, ipol, ilayer, ifold));
    
    	
    //   TGraph *g2 = new TGraph(anita1->HALFNFOUR, anita1->fTimes, timedomain_output[iband]);
    //   g2->SetTitle("Diode Output;Time [s];Diode Output [J]");
    //   g2->Draw("Al");
    //   TLine *l = new TLine (0, thresholds[ipol][iband] * anita1->bwslice_rmsdiode[iband], 200, thresholds[ipol][iband] * anita1->bwslice_rmsdiode[iband]);
    //   l->SetLineColor(kRed);
    //   l->Draw();
    //   c->Print(Form("%s_%d_%d_%d_%d_diodeOutput.png",  name.c_str(), anita1->inu, ipol, ilayer, ifold));
    //   c->Print(Form("%s_%d_%d_%d_%d_diodeOutput.pdf",  name.c_str(), anita1->inu, ipol, ilayer, ifold));
    //   c->Print(Form("%s_%d_%d_%d_%d_diodeOutput.C",    name.c_str(), anita1->inu, ipol, ilayer, ifold));
    //   c->Print(Form("%s_%d_%d_%d_%d_diodeOutput.root", name.c_str(), anita1->inu, ipol, ilayer, ifold));
    //   delete gV;
    //   delete g2;
    //   delete c;
    // }

    
    // now shift right to account for arrival times
    // this is done inside the impulse response function normally
    // but if we don't use it, we need to apply it manually here
    if (!settings1->APPLYIMPULSERESPONSETRIGGER){
      Tools::ShiftRight(timedomain_output[iband],anita1->NFOUR,(int)(anita1->arrival_times[ipol][anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));
    }
	
    if (settings1->TRIGGERSCHEME == 2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5){
      //      if (anita1->inu==1570)
      //cout << "shifting left.\n";
      Tools::ShiftLeft(timedomain_output[iband],anita1->NFOUR,ibinshift);
    }

	
    for (int itime=0;itime<Anita::NFOUR/2;itime++) {
      anita1->timedomain_output_allantennas[ipol][anita1->GetRxTriggerNumbering(ilayer,ifold)][itime]=timedomain_output[4][itime];
      //cerr<<iband<<" "<<i<<" "<<vm_banding_rfcm_forfft[ipol] <<"  "<<timedomain_output[iband][i]<<endl;
    }
	
    
    // want the bits in arrayofhits to each represent a 2 ns interval =TRIGTIMESTEP
    // but TIMESTEP doesn't necessarily go into TRIGTIMESTEP nicely
 
    int whichtrigbin=0;
    //cout << "1 whichtrigbin is " << whichtrigbin << "\n";
    // nextbreakpoint is how many time steps in digitization that we have to go 
    // to equal the length we want for the trigger bin
    int nextbreakpoint=(int)((double)(whichtrigbin+1)*(globaltrig1->TRIGTIMESTEP/anita1->TIMESTEP));

    //    if (whichlayer==2 && whichphisector==15)
    //       cout << "1 nextbreakpoint is " << nextbreakpoint << "\n";

    //    cout << "1 nextbreakpoint is " << nextbreakpoint << "\n";
    // keep advancing in trigger timesteps until we're past the iminbin.
    while (nextbreakpoint<anita1->iminbin[iband]) {
      nextbreakpoint=(int)((double)(whichtrigbin+1)*(globaltrig1->TRIGTIMESTEP/anita1->TIMESTEP));
      whichtrigbin++;
      // 	  if (whichlayer==2 && whichphisector==15) {
      // 	    cout << "2 nextbreakpoint is " << nextbreakpoint << "\n";
      // 	    cout << "whichtrigbin is " << whichtrigbin << "\n";
      // 	  }
    }
    // keep track of whether each trigger bin has a hit
    int thisisaone=0;


    for (int ibin = anita1->iminbin[iband]; ibin < anita1->imaxbin[iband]; ibin++) {
    
      if (timedomain_output[iband][ibin] < thresholds[ipol][iband] * anita1->bwslice_rmsdiode[iband] && anita1->pol_allowed[ipol] && anita1->bwslice_allowed[iband]) { // is this polarization and bw slice allowed to pass
	//std::cout << "VPOL : " << iband << " " << timedomain_output_1[iband][ibin]  << " " <<  thresholds[0][iband] << " " <<  anita1->bwslice_rmsdiode[iband] << std::endl;
	tempChansPassing[iband] = 1;// channel passes
    
	
	thisisaone=1;
	//	  cout << "got a hit.\n";
      }
      // if (whichlayer==2 && whichphisector==15) 
      //cout << "ibin, nextbreakpoint are " << ibin << "\t" << nextbreakpoint << "\n";
      if (ibin>nextbreakpoint) {
	//if (whichlayer==2 && whichphisector==15) 
	//cout << "I'm filling. size is " << "\t" << globaltrig1->arrayofhits[whichlayer][whichphisector][0][iband].size() << "\n";

	globaltrig1->arrayofhits[whichlayer][whichphisector][ipol][iband].push_back(thisisaone);
	// 	if (thisisaone) {
	//  	  cout << "filling with 1. phi is " << whichphisector << "\n";
	//  	}
	//if (thisisaone)
	//cout << "filling with 1.\n";
	//cout << "filling arrayofhits with " << thisisaone << "\n";
	//cout << "size of arrayofhits is " << globaltrig1->arrayofhits[ilayer][ifold][0][iband].size() << "\n";
	nextbreakpoint=(int)((double)(whichtrigbin+1)*(globaltrig1->TRIGTIMESTEP/anita1->TIMESTEP));
	//cout << "3 nextbreakpoint is " << nextbreakpoint << "\n";
	whichtrigbin++;
	//cout << "3 whichtrigbin is " << whichtrigbin << "\n";
	thisisaone=0;
      }

    } // end loop over bins in the window

    anita1->channels_passing[ipol][iband]=tempChansPassing[iband];

    
    if (tempChansPassing[iband]) {
      //Records number of first level triggers on each antenna for a single neutrino
      globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]++; 
	  
    }
  } // end loop over bands

 

}


void icemc::ChanTrigger::InitializeEachBand(Anita *anita1)
{
  unwarned=1;
  for (int ipol=0;ipol<2;ipol++) {
    for (int iband=0;iband<anita1->NBANDS+1;iband++) {

      vsignal_eachband[ipol].push_back(0.);
      vthreshold_eachband[ipol].push_back(0.);
      vnoise_eachband[ipol].push_back(0.);
      vpasses_eachband[ipol].push_back(0);
			
      signal_eachband[ipol][iband]=0.;
      threshold_eachband[ipol][iband]=0.;
      noise_eachband[ipol][iband]=0.;
      passes_eachband[ipol][iband]=0;
    }
  }
    
  Tools::Zero(bwslice_volts_pol0,5);
  Tools::Zero(bwslice_volts_pol1,5);
  Tools::Zero(bwslice_energy_pol0,5);
  Tools::Zero(bwslice_energy_pol1,5);
  Tools::Zero(bwslice_volts_pol0_em,5);
  Tools::Zero(bwslice_volts_pol1_em,5);
  Tools::Zero(bwslice_energy_pole,5);
  Tools::Zero(bwslice_energy_polh,5);
  Tools::Zero(bwslice_volts_polh,5);
  Tools::Zero(bwslice_volts_pole,5);
}



void icemc::ChanTrigger::ApplyAntennaGain(const Settings *settings1, Anita *anita1, Balloon *bn1, Screen *panel1, int ant, Vector &n_eplane, Vector &n_hplane, Vector &n_normal){
  
  e_component=0;
  h_component=0;
  n_component=0;
  e_component_kvector=0;
  h_component_kvector=0;
  n_component_kvector=0;
  hitangle_e=0;
  hitangle_h=0;
  
  int numBinShift;


  double tmp_vhz[2][anita1->NFREQ];
  double tmp_volts[2][anita1->NFOUR/2];
  
  for (int iband=0;iband<5;iband++) { // loop over bands
    
    Tools::Zero(volts_rx_forfft[0][iband], anita1->NFOUR/2);
    Tools::Zero(volts_rx_forfft[1][iband], anita1->NFOUR/2);
    Tools::Zero(vhz_rx[0][iband],          anita1->NFREQ);
    Tools::Zero(vhz_rx[1][iband],          anita1->NFREQ);
    
    if (anita1->bwslice_allowed[iband]!=1) continue;
    
    anita1->iminbin[iband]=0.;
    anita1->imaxbin[iband]=anita1->NFOUR/2;
        
    for (int jpt=0; jpt<panel1->GetNvalidPoints(); jpt++){
      for (int k=0;k<Anita::NFREQ;k++) {
        if (anita1->freq[k]>=settings1->FREQ_LOW_SEAVEYS && anita1->freq[k]<=settings1->FREQ_HIGH_SEAVEYS){

          //Copy frequency amplitude to screen point
          tmp_vhz[0][k]=tmp_vhz[1][k]=panel1->GetVmmhz_freq(jpt*Anita::NFREQ + k)/sqrt(2)/(anita1->TIMESTEP*1.E6);
          // cout << tmp_vhz[0][k] << endl;
	  bn1->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal,  panel1->GetVec2bln(jpt), e_component_kvector,  h_component_kvector,  n_component_kvector);
	  bn1->GetEcompHcompEvector(settings1,  n_eplane,  n_hplane,  panel1->GetPol(jpt),  e_component,  h_component,  n_component);
	  bn1->GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector, hitangle_e, hitangle_h);

          anita1->AntennaGain(settings1, hitangle_e, hitangle_h, e_component, h_component, k, tmp_vhz[0][k], tmp_vhz[1][k]);

          if (settings1->TUFFSON==2){
            tmp_vhz[0][k]=applyButterworthFilter(anita1->freq[k], tmp_vhz[0][k], anita1->TUFFstatus);
            tmp_vhz[1][k]=applyButterworthFilter(anita1->freq[k], tmp_vhz[1][k], anita1->TUFFstatus);
          }
          
        } // end if (seavey frequencies)
        else {
          tmp_vhz[0][k]=0;
          tmp_vhz[1][k]=0;
        }
      } // end looping over frequencies.

      anita1->MakeArrayforFFT(tmp_vhz[0],tmp_volts[0], 90., true);
      anita1->MakeArrayforFFT(tmp_vhz[1],tmp_volts[1], 90., true);

      // now v_banding_rfcm_h_forfft is in the time domain
      // and now it is really in units of V
      Tools::realft(tmp_volts[0],-1,anita1->NFOUR/2);
      Tools::realft(tmp_volts[1],-1,anita1->NFOUR/2);

      // put it in normal time ording -T to T
      // instead of 0 to T, -T to 0
      Tools::NormalTimeOrdering(anita1->NFOUR/2,tmp_volts[0]);
      Tools::NormalTimeOrdering(anita1->NFOUR/2,tmp_volts[1]);

      numBinShift = int(panel1->GetDelay(jpt) / anita1->TIMESTEP);
      if(fabs(numBinShift) >= anita1->HALFNFOUR){
        //cout<<"skipping"<<"\n";
        //don't bother adding it to the total since it's shifted out of range
      }
      else{
        if( panel1->GetDelay(jpt)>0 ){
          Tools::ShiftLeft(tmp_volts[0], anita1->NFOUR/2, numBinShift );
          Tools::ShiftLeft(tmp_volts[1], anita1->NFOUR/2, numBinShift );
        }
        else if( panel1->GetDelay(jpt)<0 ){
          Tools::ShiftRight(tmp_volts[0], anita1->NFOUR/2, -1*numBinShift );
          Tools::ShiftRight(tmp_volts[1], anita1->NFOUR/2, -1*numBinShift );
        }
    
        for (int k=0;k<anita1->NFOUR/2;k++) {
          volts_rx_forfft[0][iband][k] += tmp_volts[0][k];// * panel1->GetWeight(jpt) / panel1->GetWeightNorm();
          volts_rx_forfft[1][iband][k] += tmp_volts[1][k];// * panel1->GetWeight(jpt) / panel1->GetWeightNorm();
        }
      }
      
    } // end loop over screen points

    // Now need to convert time domain to frequency domain
    for (int k=0; k<anita1->NFOUR/2;k++){
      tmp_volts[0][k]=volts_rx_forfft[0][iband][k];
      tmp_volts[1][k]=volts_rx_forfft[1][iband][k];
      // cout << anita1->fTimes[k] << " " << tmp_volts[0][k] << endl;
    }

    
    // find back the frequency domain
    Tools::realft(tmp_volts[0],1,anita1->NFOUR/2);
    Tools::realft(tmp_volts[1],1,anita1->NFOUR/2);


    // Convert FFT arrays into standard icemc frequency amplitudes array
    anita1->GetArrayFromFFT(tmp_volts[0], vhz_rx[0][iband]);
    anita1->GetArrayFromFFT(tmp_volts[1], vhz_rx[1][iband]);
    
     
  } // end loop over bands


#ifdef ANITA_UTIL_EXISTS
  if (settings1->SIGNAL_FLUCT && (settings1->NOISEFROMFLIGHTDIGITIZER || settings1->NOISEFROMFLIGHTTRIGGER) )
    getNoiseFromFlight(anita1, ant, settings1->SIGNAL_FLUCT > 0);

  if (settings1->ADDCW){
    memset(cw_digPath, 0, sizeof(cw_digPath));
    calculateCW(anita1, 250E6, 0, 0.000005);
  }
  
#endif
  
}

void icemc::ChanTrigger::TriggerPath(const Settings *settings1, Anita *anita1, int ant, Balloon *bn1){


  double integrate_energy_freq[5]={0.,0.,0.,0.,0.};
    
  for (int iband=0;iband<5;iband++) { // loop over bands
    if (anita1->bwslice_allowed[iband]!=1) continue;
    
    anita1->iminbin[iband]=0.;
    anita1->imaxbin[iband]=anita1->NFOUR/2;

    Tools::Zero(v_banding_rfcm_forfft[0][iband],anita1->NFOUR/2);
    Tools::Zero(v_banding_rfcm_forfft[1][iband],anita1->NFOUR/2);
    
    for (int ifreq=0; ifreq<anita1->NFOUR/2; ifreq++){
      v_banding_rfcm[0][iband][ifreq]=vhz_rx[0][iband][ifreq];
      v_banding_rfcm[1][iband][ifreq]=vhz_rx[1][iband][ifreq];
    }
    
    // If not applying the trigger impulse response
    // Apply banding and RFCM
    // And then convert to time domain
    if (!settings1->APPLYIMPULSERESPONSETRIGGER){
      
      anita1->Banding(iband,anita1->freq,v_banding_rfcm[0][iband],Anita::NFREQ); 
      anita1->Banding(iband,anita1->freq,v_banding_rfcm[1][iband],Anita::NFREQ);
      anita1->RFCMs(1,1,v_banding_rfcm[0][iband]);
      anita1->RFCMs(1,1,v_banding_rfcm[1][iband]);
      
      
      for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
	if (anita1->freq[ifreq]>=settings1->FREQ_LOW_SEAVEYS && anita1->freq[ifreq]<=settings1->FREQ_HIGH_SEAVEYS){
	  addToChannelSums(settings1, anita1, iband, ifreq);
	}
      } // end loop over nfreq
      
      
      anita1->MakeArrayforFFT(v_banding_rfcm[0][iband],v_banding_rfcm_forfft[0][iband], 90., true);
      anita1->MakeArrayforFFT(v_banding_rfcm[1][iband],v_banding_rfcm_forfft[1][iband], 90., true);
      
      // for some reason I'm averaging over 10 neighboring bins
      // to get rid of the zero bins
      for (int i=0;i<anita1->NFOUR/4;i++) {
        for (int ipol=0;ipol<2;ipol++){
	  
          v_banding_rfcm_forfft_temp[ipol][iband][2*i]  =0.;
          v_banding_rfcm_forfft_temp[ipol][iband][2*i+1]=0.;
          
          int tempcount = 0;
          for (int k=i;k<i+10;k++) {
            if (k<anita1->NFOUR/4) {
              v_banding_rfcm_forfft_temp[ipol][iband][2*i]  +=v_banding_rfcm_forfft[ipol][iband][2*k];
              v_banding_rfcm_forfft_temp[ipol][iband][2*i+1]+=v_banding_rfcm_forfft[ipol][iband][2*k+1];
              tempcount++;
            }
          }
          
          v_banding_rfcm_forfft[ipol][iband][2*i]  =v_banding_rfcm_forfft_temp[ipol][iband][2*i]/tempcount;
          v_banding_rfcm_forfft[ipol][iband][2*i+1]=v_banding_rfcm_forfft_temp[ipol][iband][2*i+1]/tempcount;
          
          v_banding_rfcm_forfft_temp[ipol][iband][2*i]  =v_banding_rfcm_forfft[ipol][iband][2*i];
          v_banding_rfcm_forfft_temp[ipol][iband][2*i+1]=v_banding_rfcm_forfft[ipol][iband][2*i+1];
        }
      }

      // now v_banding_rfcm_h_forfft is in the time domain
      // and now it is really in units of V
      Tools::realft(v_banding_rfcm_forfft[0][iband],-1,anita1->NFOUR/2);
      Tools::realft(v_banding_rfcm_forfft[1][iband],-1,anita1->NFOUR/2);
      
      // put it in normal time ording -T to T
      // instead of 0 to T, -T to 0
      Tools::NormalTimeOrdering(anita1->NFOUR/2,v_banding_rfcm_forfft[0][iband]);
      Tools::NormalTimeOrdering(anita1->NFOUR/2,v_banding_rfcm_forfft[1][iband]);

      for (int itime=0; itime<anita1->NFOUR/2; itime++){
	justSig_trigPath[0][itime] = v_banding_rfcm_forfft[0][iband][itime];
	justSig_trigPath[1][itime] = v_banding_rfcm_forfft[1][iband][itime];
      }

      
    } else {
      
      for (int itime=0; itime<anita1->NFOUR/2 ; itime++){
	v_banding_rfcm_forfft[0][iband][itime]=volts_rx_forfft[0][iband][itime];
	v_banding_rfcm_forfft[1][iband][itime]=volts_rx_forfft[1][iband][itime];
      }
      
#ifdef ANITA_UTIL_EXISTS
      // if applying the impulse response
      applyImpulseResponseTrigger(settings1, anita1, ant, v_banding_rfcm_forfft[0][iband], v_banding_rfcm[0][iband], 0);
      applyImpulseResponseTrigger(settings1, anita1, ant, v_banding_rfcm_forfft[1][iband], v_banding_rfcm[1][iband], 1);
#endif
    }
    
    
    if (settings1->ZEROSIGNAL) {
      Tools::Zero(v_banding_rfcm_forfft[0][iband],anita1->NFOUR/2);
      Tools::Zero(v_banding_rfcm_forfft[1][iband],anita1->NFOUR/2);
    }
  

    for (int i=0;i<anita1->NFOUR/4;i++) {
      integrate_energy_freq[iband]+=v_banding_rfcm_forfft[0][iband][2*i]*v_banding_rfcm_forfft[0][iband][2*i]+v_banding_rfcm_forfft[0][iband][2*i+1]*v_banding_rfcm_forfft[0][iband][2*i+1];
    }
  
    // write the signal events to a tree
    for (int k=0;k<anita1->NFOUR/2;k++) {
      anita1->signal_vpol_inanita[iband][k]=v_banding_rfcm_forfft[0][iband][k];
    }
    anita1->integral_vmmhz_foranita=integral_vmmhz;
  
    // Find the p2p value before adding noise
    anita1->peak_v_banding_rfcm[0][iband]=FindPeak(v_banding_rfcm_forfft[0][iband],anita1->NFOUR/2);
    anita1->peak_v_banding_rfcm[1][iband]=FindPeak(v_banding_rfcm_forfft[1][iband],anita1->NFOUR/2);

    // Find the p2p value before adding noise
    for (int iband=0;iband<5;iband++) {
      if (anita1->bwslice_allowed[iband]!=1) continue;
      for (int ipol=0; ipol<2; ipol++) {
        anita1->peak_v_banding_rfcm[ipol][iband]=FindPeak(v_banding_rfcm_forfft[ipol][iband],anita1->NFOUR/2);
      }
    }
  } // end loop over bands
  
    

}






void icemc::ChanTrigger::DigitizerPath(const Settings *settings1, Anita *anita1, int ant, Balloon *bn1)
{
  double vhz_rx_rfcm_e[Anita::NFREQ]; // V/Hz after rx, rfcm
  double vhz_rx_rfcm_h[Anita::NFREQ];

  int fNumPoints = anita1->HALFNFOUR;
  
  // Apply anita-3 measured impulse response
  if (settings1->APPLYIMPULSERESPONSEDIGITIZER){
    anita1->GetNoiseWaveforms(); // get noise waveforms
    for (int i=0;i<fNumPoints;i++){
      volts_rx_rfcm_lab[0][i] = volts_rx_forfft[0][4][i];
      volts_rx_rfcm_lab[1][i] = volts_rx_forfft[1][4][i];
    }
    
#ifdef ANITA_UTIL_EXISTS
    applyImpulseResponseDigitizer(settings1, anita1, fNumPoints, ant, anita1->fTimes, volts_rx_rfcm_lab[0], 0);
    applyImpulseResponseDigitizer(settings1, anita1, fNumPoints, ant, anita1->fTimes, volts_rx_rfcm_lab[1], 1);
#endif
    
    if (settings1->SIGNAL_FLUCT > 0 && !settings1->NOISEFROMFLIGHTDIGITIZER){
      for (int i=0;i<anita1->NFOUR/2;i++) {
     	for (int ipol=0;ipol<2;ipol++){
     	  volts_rx_rfcm_lab[ipol][i]+=anita1->timedomainnoise_lab[ipol][i]; // add noise
     	}
      }
    }
    
  } else {


    for (int ifreq=0; ifreq<anita1->NFREQ; ifreq++){
      vhz_rx_rfcm_e[ifreq]=vhz_rx[0][4][ifreq];
      vhz_rx_rfcm_h[ifreq]=vhz_rx[1][4][ifreq];
    }
    // for frequency-domain voltage-based trigger (triggerscheme==0)
    // we do not apply rfcm's
    // for other trigger types we do
      
    // apply rfcm's
    if (settings1->TRIGGERSCHEME==1 || settings1->TRIGGERSCHEME==2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5) {
      anita1->RFCMs(1,1,vhz_rx_rfcm_e);
      anita1->RFCMs(1,1,vhz_rx_rfcm_h);
    }
    
    double scale;
    double sumpower=0.;
    if (settings1->PULSER) { // if we are using the pulser spectrum instead of simulating neutrinos
      scale=Tools::dMax(vhz_rx_rfcm_e,Anita::NFREQ)/Tools::dMax(anita1->v_pulser,anita1->NFOUR/4);
      sumpower=0.;
      int ifour;// index for fourier transform
      for (int i=0;i<Anita::NFREQ;i++) {
        ifour=Tools::Getifreq(anita1->freq[i],anita1->freq_forfft[0],anita1->freq_forfft[anita1->NFOUR/2-1],anita1->NFOUR/4);
        vhz_rx_rfcm_e[i]=scale*anita1->v_pulser[ifour];
        vhz_rx_rfcm_h[i]=0.;
        sumpower+=vhz_rx_rfcm_e[i]*vhz_rx_rfcm_e[i];
      }
    } // end if we are just using the pulser spectrum
      
      
    for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
      anita1->avgfreq_rfcm[ifreq]+=vhz_rx_rfcm_e[ifreq];
    }
    
    // change their length from Anita::NFREQ to HALFNFOUR
    anita1->MakeArrayforFFT(vhz_rx_rfcm_e,volts_rx_rfcm[0], 90., true);
    anita1->MakeArrayforFFT(vhz_rx_rfcm_h,volts_rx_rfcm[1], 90., true);
      
          
    // now the last two are in the frequency domain
    // convert to the time domain
    // still don't have any noise
      
    Tools::realft(volts_rx_rfcm[0],-1,anita1->HALFNFOUR);
    Tools::realft(volts_rx_rfcm[1],-1,anita1->HALFNFOUR);
      
    anita1->GetNoiseWaveforms(); // get noise waveforms
      
    // find the peak right here and it might be the numerator of the horizontal axis of matt's plot
    anita1->peak_rx_rfcm_signalonly[0]=icemc::ChanTrigger::FindPeak(volts_rx_rfcm[0],anita1->HALFNFOUR); // with no noise
    anita1->peak_rx_rfcm_signalonly[1]=icemc::ChanTrigger::FindPeak(volts_rx_rfcm[1],anita1->HALFNFOUR);
      
    if (settings1->SIGNAL_FLUCT) {
      for (int i=0;i<anita1->NFOUR/2;i++) {
	for (int ipol=0;ipol<2;ipol++){
	  volts_rx_rfcm[ipol][i]+=anita1->timedomainnoise_rfcm[ipol][i]; // add noise.
	}
      } 
    }
    
    
    anita1->peak_rx_rfcm[0]=icemc::ChanTrigger::FindPeak(volts_rx_rfcm[0],anita1->HALFNFOUR); // with noise 
    anita1->peak_rx_rfcm[1]=icemc::ChanTrigger::FindPeak(volts_rx_rfcm[1],anita1->HALFNFOUR); // with noise
      

    double vhz_rx_rfcm_lab_e[Anita::NFREQ]; // V/Hz after rx, rfcm and lab
    double vhz_rx_rfcm_lab_h[Anita::NFREQ];

    for (int i=0;i<Anita::NFREQ;i++) {
      vhz_rx_rfcm_lab_e[i]=vhz_rx_rfcm_e[i];
      vhz_rx_rfcm_lab_h[i]=vhz_rx_rfcm_h[i];
    }  
    // now go back to vhz_rx_e,h and apply rfcm's and surf attn. and then we will write the time domain waveforms to tsignals for andres
    // apply surf (lab) attn.
    anita1->labAttn(vhz_rx_rfcm_lab_e);
    anita1->labAttn(vhz_rx_rfcm_lab_h);
      
    for (int i=0;i<Anita::NFREQ;i++) {
      anita1->avgfreq_rfcm_lab[i]+=vhz_rx_rfcm_lab_e[i];
    }

    // change their length from Anita::NFREQ to HALFNFOUR
    anita1->MakeArrayforFFT(vhz_rx_rfcm_lab_e,volts_rx_rfcm_lab[0], 90., true);
    anita1->MakeArrayforFFT(vhz_rx_rfcm_lab_h,volts_rx_rfcm_lab[1], 90., true);
      
    // now the last two are in the frequency domain
    // convert to the time domain
    // still don't have any noise
    Tools::realft(volts_rx_rfcm_lab[0],-1,anita1->HALFNFOUR); 
    Tools::realft(volts_rx_rfcm_lab[1],-1,anita1->HALFNFOUR);

    // put it in normal time ording -T to T
    // instead of 0 to T, -T to 0 
    Tools::NormalTimeOrdering(anita1->NFOUR/2,volts_rx_rfcm_lab[0]); // EH, why only this has NormalTimeOrdering applied? Why not before?
    Tools::NormalTimeOrdering(anita1->NFOUR/2,volts_rx_rfcm_lab[1]);

    if (settings1->SIGNAL_FLUCT > 0) { 
      for (int i=0;i<anita1->NFOUR/2;i++) {
	for (int ipol=0;ipol<2;ipol++){
	  justSig_digPath[ipol][i] = volts_rx_rfcm_lab[ipol][i];
	  volts_rx_rfcm_lab[ipol][i]+=anita1->timedomainnoise_lab[ipol][i]; // add noise
	}
      }
    }//end if signal_fluct

  } // END ELSE IMPULSE RESPONSE
}// end DigitizerPath()




void icemc::ChanTrigger::TimeShiftAndSignalFluct(const Settings *settings1, Anita *anita1, int ilayer, int ifold, double volts_rx_rfcm_lab_e_all[48][512], double volts_rx_rfcm_lab_h_all[48][512])
{   
  // int ant = anita1->GetRxTriggerNumbering(ilayer, ifold);


  // now shift right to account for arrival times
  // this is done inside the impulse response function normally
  // if we don't use it, we need to do it here
  if (!settings1->APPLYIMPULSERESPONSEDIGITIZER){
    //  for (int i=0;i<48;i++) std::cout << "Arrival times " << anita1->arrival_times[0][anita1->GetRx(ilayer,ifold)] << std::endl;
    Tools::ShiftRight(volts_rx_rfcm_lab[0],anita1->NFOUR/2, int(anita1->arrival_times[0][anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));
    Tools::ShiftRight(volts_rx_rfcm_lab[1],anita1->NFOUR/2, int(anita1->arrival_times[1][anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));
  }
  
  for (int i=0;i<anita1->NFOUR/2;i++) {
    volts_rx_rfcm_lab_e_all[anita1->GetRx(ilayer, ifold)][i] = volts_rx_rfcm_lab[0][i];
    volts_rx_rfcm_lab_h_all[anita1->GetRx(ilayer, ifold)][i] = volts_rx_rfcm_lab[1][i];      
  }

  // now vmmhz_rx_rfcm_lab_e,h_forfft are the time domain waveforms after the antenna and lab attenuation
  // now find peak voltage
  // these get written to a tree
  anita1->peak_rx_rfcm_lab[0]=icemc::ChanTrigger::FindPeak(volts_rx_rfcm_lab[0],anita1->HALFNFOUR);
  anita1->peak_rx_rfcm_lab[1]=icemc::ChanTrigger::FindPeak(volts_rx_rfcm_lab[1],anita1->HALFNFOUR);  


  // END OF DIGITIZER PATH
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  
  
} 



icemc::ChanTrigger::ChanTrigger() {
    
}


double icemc::ChanTrigger::FindPeak(double *waveform,int n) {
    
  // find peak abs(voltage) of this waveform
    
  double peak=0.; // positive peak
    
  for (int i=0;i<n;i++) {
    if (fabs(waveform[i])>peak)
      peak=fabs(waveform[i]);
  }
  return peak;
    
}


void icemc::ChanTrigger::addToChannelSums(const Settings *settings1,Anita *anita1,int ibw, int k) {
    
  double term_e= v_banding_rfcm[0][ibw][k]*    // for the integral over e-field
    (settings1->BW/(double)Anita::NFREQ/1.E6); // width of a frequency bin
  // uses e- and h-components instead of lcp,rcp, for Anita-lite
    
  double energyterm_e=v_banding_rfcm[0][ibw][k]*v_banding_rfcm[0][ibw][k]* // for the integral of energy of time T=N delta t = 1/delta f
    (settings1->BW/(double)Anita::NFREQ/1.E6)*(settings1->BW/(double)Anita::NFREQ/1.E6)/50.*
    anita1->INTEGRATIONTIME; // divide by width of a freq. bin
  // end divide by 50 ohms
    
    
  double term_h=v_banding_rfcm[1][ibw][k]* // for the integral over e-field
    (settings1->BW/(double)Anita::NFREQ/1.E6);
    
  double energyterm_h=v_banding_rfcm[1][ibw][k]*v_banding_rfcm[1][ibw][k]* // for the integral over energy
    (settings1->BW/(double)Anita::NFREQ/1.E6)*(settings1->BW/(double)Anita::NFREQ/1.E6)/50.*
    anita1->INTEGRATIONTIME; // multiply by frequency bin and divide by 50 ohms
    
    
  bwslice_volts_pole[ibw]+=term_e;
  bwslice_energy_pole[ibw]+=energyterm_e;
  bwslice_volts_polh[ibw]+=term_h;
  bwslice_energy_polh[ibw]+=energyterm_h;
    
}


double icemc::ChanTrigger::ADCCountstoPowerThreshold(Anita *anita1, int ipol, int iant) {
  // first convert threshold in adc counts to the singles rate
  // do this using threshold scans taken before the flight

  // For Anita-3 using run 11927
  // these curves were read in the Balloon constructor
  // first check if the threshold in adc counts is in the allowable range
  Float_t threshadc = (Float_t)anita1->thresholds[ipol][iant];

  // If there was a problem reading the flight threshold
  // return 5
  // 5 is the average relative threshold corresponding to
  // 500kHz scalers in a full band trigger
  if (threshadc<10){
    return 5;
  }

  // Broken channel during ANITA-3 flight
  if (ipol==1 && iant==7){
    return 5;
  }
  
  if (threshadc<anita1->minadcthresh[ipol][iant]) {
    if (unwarned) {
      cout << "Warning! ADC threshold is outside range of measured threshold scans.";
      cout << "It is below the minimum so set it to the minimum.  Will not be warned again.\n";
      unwarned=0;
    }
    threshadc=anita1->minadcthresh[ipol][iant];
  }
  if (threshadc>anita1->maxadcthresh[ipol][iant]) {
    if (unwarned) {
      cout << "Warning! ADC threshold is outside range of measured threshold scans.";
      cout << "It is higher than the maximum so set it to the maximum.  Will not be warned again.\n";
      unwarned=0;
    }
    threshadc=anita1->maxadcthresh[ipol][iant];
  }

  // Now find singles rate for this threshold
  // first sort thresholds
  int index=TMath::BinarySearch(anita1->npointThresh, anita1->threshScanThresh[ipol][iant], threshadc);
  
  thisrate=(double)anita1->threshScanScaler[ipol][iant][index]; // these scalers are in kHz

  thisrate=(double)anita1->scalers[ipol][iant];
  
  // now find threshold for this scaler.  Here, scalers have to be in MHz.
  thisrate=thisrate/1.E3; // put it in MHz

  // figure out what band we're talking about
  // int iband=Anita::SurfChanneltoBand(ichan);
  // FOR THE MOMENT JUST USING THE FULL BAND
  int iband =3;
  thispowerthresh=rateToThreshold(thisrate,iband);

  // Avoid inf and nan: 5 is the relative power threshold corresponding to 500kHz scalers
  if (thispowerthresh>999999) return 5.;
  if (thispowerthresh<0.0001) return 5.;

  // // Limit on relative power threshold to avoid thermal noise to trigger
  // if (thispowerthresh<4.5) return 4.5;

  return thispowerthresh;
    
}


double icemc::ChanTrigger::getRate() {
    
  return thisrate;
}



double icemc::ChanTrigger::rateToThreshold(double rate, int band)
{
  double constant=0.;
  double slope=0.;
  //Note:  This method is currently very specific.  It assumes ANITA 06 parameters and
  //an exponential diode model with time constant 3.75e-9 s.
  //The rate to threshold relationship is modeled as  rate = exp(constant + slope * threshold).
  //The values of the constant and the slope were found by running the diode over noise traces with
  //different thresholds, and measuring the trigger rate at each threshold.
  switch (band) {
  case 0:
    constant = 5.15629;
    slope = -0.94107;
    break;
  case 1:
    constant = 5.33636;
    slope = -1.04554;
    break;
  case 2:
    constant = 5.76685;
    slope = -1.28046;
    break;
  case 3:
    constant = 5.87953;
    slope = -1.32115;
    break;
  default:
    std::cerr<<"[AnitaHardwareTrigger::rateToThreshold] : Unknown band selected (band "<<band<<"?)!\n";
    exit(1);
    break;
  } //end switch (band)
  return (log(rate) - constant) / slope;
} //method AnitaHardwareTrigger::rateToThreshold(double rate, int band)



int icemc::ChanTrigger::IsItUnmasked(unsigned short surfTrigBandMask[9][2],int ibw,int ilayer, int ifold, int ipol) {
    
  if (ibw>3)
    return 1;
    
  int whichband= Anita::WhichBand(ibw,ipol); // from 1 to 8, which scalar channel on an antenna this band corresponds to.
    
  int surf=Anita::AntennaNumbertoSurfNumber(ilayer,ifold); // which surf, 1-9
    
  int antenna=Anita::GetAntennaNumber(ilayer,ifold); // which antenna, 1-32
    
  if (antenna%2==0) // if it's divisible by 2, increase whichband by 8 so it runs 1 to 16
    whichband+=8;
    
  int antennaonsurf=(antenna-1)%4; // number that runs 0-3 for four antennas on the surf
    
  unsigned short whichbit= 1 << (whichband-1); // which bit of surfTrigBandMask we're interested in.  Starts with 1 and multiplies it by whichband-1 powers of 2
  // the max this can be is 2^15
  // which we & with surfTrigBandMask we get whether that bit is 1 or 0
    
  if (antennaonsurf<2) {
    if ((surfTrigBandMask[surf-1][0] & whichbit)>0)
      return 0;
  }
  else {
    if ((surfTrigBandMask[surf-1][1] & whichbit)>0)
      return 0;
  }
    
  return 1;
    
    
}





void icemc::ChanTrigger::L1Trigger(Anita *anita1,double timedomain_output_1[5][Anita::NFOUR],double timedomain_output_2[5][Anita::NFOUR],double powerthreshold[2][5],
			    int *channels_passing_e_forglob,int *channels_passing_h_forglob,int &npass) {
   
  int maxsample=TMath::MaxElement(5,anita1->imaxbin);
  int minsample=TMath::MinElement(5,anita1->iminbin);
    
  for (int j=0;j<5;j++) {
    flag_e[j].clear();
    flag_h[j].clear();		
      
    for (int i=minsample;i<maxsample;i++) {
      //      std::cout << anita1->inu << " " << timedomain_output_1[j][i] << " " << powerthreshold[0][j] << " " << anita1->bwslice_rmsdiode[j] << std::endl;
      if (timedomain_output_1[j][i]<powerthreshold[0][j]*anita1->bwslice_rmsdiode[j] && anita1->bwslice_allowed[j]==1) {
	flag_e[j].push_back(1);
	  
      }
      else
	flag_e[j].push_back(0);
	
	
      if (timedomain_output_2[j][i]<powerthreshold[1][j]*anita1->bwslice_rmsdiode[j] && anita1->bwslice_allowed[j]==1)
	flag_h[j].push_back(1);
      else
	flag_h[j].push_back(0);
	
    } // end loop over samples in window where we look for single channel trigger firing
  } // end loop over bands
    
  int nstayhigh=(int)(anita1->l1window/anita1->TIMESTEP);
  // now make each flag stay high for the required amount to time
  for (int j=0;j<5;j++) {
    for (int i=minsample;i<maxsample;i++) { // loop over samples in window where we looked for a single channel trigger
	
      if (flag_e[j][i-minsample]==1) // if the flag is high (remember the second index of flag_e counts from the start of the window)
	for(int k=nstayhigh-1; k>0 && i<maxsample-1; k--) { // then for nstayhigh-1 samples after than we keep it high
	  // i<maxsample-1 makes sure that the second index of flag_e is always less than maxsample-minsample-1.
	  i++;
	  flag_e[j][i-minsample]=1;
	    
	}
	
      //      if (flag_e[j][i]==1) {
      // 	for (int k=i;k<i+nstayhigh;k++) {
      // 	  if (k<NSAMPLES)
      // 	    flag_e[j][k]=1;
      // 	} // end loop over samples where we want it to stay high
	
      //       } // end if flag is high
    }
      
    for (int i=minsample;i<maxsample;i++) {
      if (flag_h[j][i-minsample]==1)
	for(int k=nstayhigh-1; k>0 && i<maxsample-1; k--) {
	  i++;
	  flag_h[j][i-minsample]=1;
	}
	
      //      if (flag_h[j][i]==1) {
      // 	for (int k=i;k<i+nstayhigh;k++) {
      // 	  if (k<NSAMPLES)
      // 	    flag_h[j][k]=1;
      // 	} // end loop over samples where we want it to stay high
	
      //       } // end if flag is high
	
    } // end loop over samples
  } // end loop over bands
    
    // now find the sample with the highest number of bands that pass
  int maxbands=0;
  int nbands_pass=0;
  for (int i=minsample;i<maxsample;i++) {
    nbands_pass=0;
    for (int j=0;j<5;j++) {
	
      if (flag_e[j][i-minsample]==1)
	nbands_pass++;
      if (flag_h[j][i-minsample]==1)
	nbands_pass++;
    }
      
    if (nbands_pass>maxbands) {
	
      maxbands=nbands_pass;
      for (int j=0;j<5;j++) {
	channels_passing_e_forglob[j]=flag_e[j][i-minsample];
	channels_passing_h_forglob[j]=flag_h[j][i-minsample];
      }
    }
  }
  npass=maxbands;


}
  


double icemc::ChanTrigger::GetNoise(const Settings *settings1,double altitude_bn,double geoid,double theta_zenith,double bw,double temp) {
    
  int NSTEPS=1000;
  //  double VSKY=150;
  double VSKY=15;
  double VICE=240;
  double VSYSTEM=200;
    
  double sum=0;
  double theta_pos=0;
  double theta_signed=0;
  double theta_cant=theta_zenith-constants::PI/2;
  double theta_horizon=acos(geoid/(geoid+altitude_bn)); // angle between horizontal on the payload and line of sight to horizon
  //  double theta_horizon=sqrt(2*altitude_bn/geoid);  // angle between horizontal on the payload and line of sight to horizon
  //double theta_horizon=-1;
  double theta_0=50*constants::RADDEG;
  double integral=0.;
  double integral_firsthalf=0;
  double integral_secondhalf=0;
  double vnoise=0;
    
  if (settings1->WHICH != 0) {
		
    for (int i=0;i<NSTEPS;i++) {
			
      // step in theta
      theta_pos=(double)fabs(-constants::PI+(double)i/(double)NSTEPS*2*constants::PI); // this is always a positive number
      theta_signed=(double)(-constants::PI+(double)i/(double)NSTEPS*2*constants::PI); // this is allowed to be signed
			
      if (theta_signed<theta_horizon-theta_cant) {
	vnoise=VSKY+VSYSTEM;
	integral_firsthalf+=exp(-2*constants::ALOG2*pow(theta_pos/theta_0,2))*2*constants::PI/NSTEPS;
      } //if
      else {
	vnoise=VICE+VSYSTEM;
	integral_secondhalf+=exp(-2*constants::ALOG2*pow(theta_pos/theta_0,2))*2*constants::PI/NSTEPS;
      } //else
			
      sum+=vnoise*exp(-2*constants::ALOG2*pow(theta_pos/theta_0,2))*2*constants::PI/NSTEPS;
      integral+=exp(-2*constants::ALOG2*pow(theta_pos/theta_0,2))*2*constants::PI/NSTEPS;
    } //if
		
    sum=sum/integral;
		
    //    cout << "sum, KBOLTZ, bw, sqrt are " << sum << " " << KBOLTZ << " " << bw << " " << sqrt(sum*50.*KBOLTZ*bw) << "\n";
    return sqrt(sum*50.*constants::KBOLTZ*bw);
  } //if (settings1->WHICH != 0)
  else if (settings1->WHICH == 0)
    return sqrt(temp*50.*constants::KBOLTZ*bw);
    
  return 0;
}//GetNoise



void icemc::ChanTrigger::GetThresholds(const Settings *settings1, const Anita *anita1,int ilayer,double thresholds[2][5]) const {
    
  if (ilayer==3 && settings1->DISCONES==2) // if it's a nadir layer
    for (int i=0;i<5;i++) {
      thresholds[0][i]=anita1->powerthreshold_nadir[i];
      thresholds[1][i]=anita1->powerthreshold_nadir[i];
    }
  else
    for (int i=0;i<5;i++) {
      thresholds[0][i]=anita1->powerthreshold[i];
      thresholds[1][i]=anita1->powerthreshold[i];
    }
}



double icemc::ChanTrigger::applyButterworthFilter(double ff, double ampl, int notchStatus[3]){
  // Butterworth filter for two ANITA notches.
  // order0 = order1 = 1 may be closer to hardware notch
  //  but 2 gives better SW rejection in analysis spectrum.
  // f0 is fairly constant, f1 does vary with time.
  double f[3]     = {250e6, 360e6, 460e6};
  double W[3]     = { 45e6,  55e6,  50e6};
  double order[3] = {   2.,    2.,    2.};

  double denominator = 1.;
  for (int i=0; i<3; i++){
    denominator += notchStatus[i]*pow(ff*W[i]/(ff*ff-f[i]*f[i]),2.*order[i]); 
  }

  return ampl/sqrt(denominator);
}






#ifdef ANITA_UTIL_EXISTS    
void icemc::ChanTrigger::applyImpulseResponseDigitizer(const Settings *settings1, Anita *anita1, int nPoints, int ant, double *x, double y[512], bool pol){

  if (settings1->ZEROSIGNAL){
    for (int i=0;i<nPoints;i++) y[i]=0;
  }

  TGraph *graph1;
  if (settings1->TRIGGEREFFSCAN && !pol){
    graph1 = getPulserAtAMPA(anita1, ant);
  } else {
    graph1 = new TGraph(nPoints, x, y);
  }

  
  // Upsample waveform to same deltaT of the signal chain impulse response
  TGraph *graphUp = FFTtools::getInterpolatedGraph(graph1, anita1->deltaT);

  int ipol=0;
  int iring=2;
  if (pol) ipol = 1;
  if (ant<16) iring=0;
  else if (ant<32) iring=1;

  int iphi = ant - (iring*16);
  
  if (settings1->ADDCW){
    for (int i=0;i<nPoints;i++){
      y[i]+=cw_digPath[ipol][i];
    }
  }

  TGraph *surfSignal;
  
  
  //Calculate convolution
  if(!settings1->TUFFSON){
    surfSignal = FFTtools::getConvolution(graphUp, anita1->fSignalChainResponseDigitizer[ipol][iring][iphi]);
  }
  else
  {
    // keith editing 1/24/18
    surfSignal = FFTtools::getConvolution(graphUp, anita1->fSignalChainResponseDigitizerTuffs[ipol][iring][iphi][anita1->tuffIndex]);
    // end keith editing
  }// end else for ANITA-4

  int irx = ant;
  if (iring==0) {
    if (ant%2==0) irx = ant/2;
    else          irx = 8 + ant/2;
  }

  // Translate waveform according to arrival times
  TGraph *surfTrans  = FFTtools::translateGraph(surfSignal, anita1->arrival_times[ipol][irx]*1e9 ) ;
  
  //Downsample again
  TGraph *surfSignalDown = FFTtools::getInterpolatedGraph(surfTrans, 1/2.6);
  
  // add thermal noise for anita-3 flight
  if (settings1->SIGNAL_FLUCT && settings1->NOISEFROMFLIGHTDIGITIZER) { 
    for (int i=0;i<nPoints;i++){
      justSig_digPath[ipol][i] = surfSignalDown->Eval(x[i]);
      y[i] = justSig_digPath[ipol][i] + justNoise_digPath[ipol][i];
    }
  } else {
    for (int i=0;i<nPoints;i++)  justSig_digPath[ipol][i] = y[i] = surfSignalDown->Eval(x[i]);
  }
  

  // if (ant ==8 && pol==0){
  //   TCanvas *c = new TCanvas("c");
  //   graph1->Draw("Al");
  //   c->Print("DigitPath_graph1.png");
  //   graphUp->Draw("Al");
  //   c->Print("DigitPath_graphUp.png");
  //   surfSignal->Draw("Al");
  //   c->Print("DigitPath_surfSignal.png");
  //   surfSignalDown->Draw("Al");
  //   c->Print("DigitPath_surfSignalDown.png");
  // }
  
  // Cleaning up
  delete surfSignalDown;
  delete surfTrans;
  delete surfSignal;
  delete graphUp;
  delete graph1;
}

void icemc::ChanTrigger::applyImpulseResponseTrigger(const Settings *settings1, Anita *anita1, int ant, double y[512], double *vhz, bool pol){

  int nPoints = anita1->HALFNFOUR;
  double *x   = anita1->fTimes;
  
  if (settings1->ZEROSIGNAL){
    for (int i=0;i<nPoints;i++) y[i]=0;
  }
  
  TGraph *graph1;
  if (settings1->TRIGGEREFFSCAN && !pol){
    graph1 = getPulserAtAMPA(anita1, ant);
  } else {
    graph1 = new TGraph(nPoints, x, y);
  }
  
  // Upsample waveform to same deltaT of the signal chain impulse response
  TGraph *graphUp = FFTtools::getInterpolatedGraph(graph1, anita1->deltaT);

  double voltsArray[512];
  int ipol=0;
  int iring=2;
  if (pol) ipol = 1;
  if (ant<16) iring=0;
  else if (ant<32) iring=1;

  int iphi = ant - (iring*16);

  //Calculate convolution

// begin keith edits
  TGraph *surfSignal;
  if (!settings1->TUFFSON){
    surfSignal = FFTtools::getConvolution(graphUp, anita1->fSignalChainResponseTrigger[ipol][iring][iphi]);
  }
  else
  {
    // keith editing 1/24/18
    surfSignal = FFTtools::getConvolution(graphUp, anita1->fSignalChainResponseTriggerTuffs[ipol][iring][iphi][anita1->tuffIndex]); 
    // end keith editing 
  }// end else anita 4
// end keith edits

  int irx = ant;
  if (iring==0) {
    if (ant%2==0) irx = ant/2;
    else          irx = 8 + ant/2;
  }
  
  // Translate signal
  TGraph *surfTrans  = FFTtools::translateGraph(surfSignal, anita1->arrival_times[ipol][irx]*1e9 ) ;
  
  //Downsample again
  TGraph *surfSignalDown = FFTtools::getInterpolatedGraph(surfTrans, 1/2.6);

  // add thermal noise for anita-3 flight
  if (settings1->SIGNAL_FLUCT && settings1->NOISEFROMFLIGHTTRIGGER) { 
    for (int i=0;i<nPoints;i++){
      justSig_trigPath[ipol][i] = surfSignalDown->Eval(x[i]);
      y[i] = voltsArray[i] = justSig_trigPath[ipol][i] + justNoise_trigPath[ipol][i];
      //  std::cout << i << " " << justNoise_trigPath[ipol][i] << std::endl;
    }
  } else {
    for (int i=0;i<nPoints;i++)  justSig_trigPath[ipol][i] = y[i] = voltsArray[i] = surfSignalDown->Eval(x[i]);
  }
  
  // find back the frequency domain
  Tools::realft(voltsArray,1,anita1->NFOUR/2);
  
  //convert the V pol time waveform into frequency amplitudes
  anita1->GetArrayFromFFT(voltsArray, vhz);
  
  // if (anita1->inu==1 && pol==0 && (ant==16 || ant==31 || ant==32 || ant==47)){
  //  TCanvas *c = new TCanvas("c");
  //  graph1->Draw("Al");
  //  c->Print(Form("TriggerPath_ant%i_graph1.png", ant));
  //  graphUp->Draw("Al");
  //  c->Print(Form("TriggerPath_ant%i_graphUp.png", ant));
  //  surfSignal->Draw("Al");
  //  c->Print(Form("TriggerPath_ant%i_surfSignal.png", ant));
  //  surfSignalDown->Draw("Al");
  //  c->Print(Form("TriggerPath_ant%i_surfSignalDown.png", ant));
  //  TGraph *gtemp = new TGraph (nPoints, x, y);
  //  gtemp->Draw("Al");
  //  c->Print(Form("TriggerPath_ant%i_surfSignalDown_noise.png", ant));
  // // TFile *out = new TFile("Icemc_signalChainTrigger.root", "recreate");
  // // graph1->Write("gInput");
  // // graphUp->Write("gInputUp");
  // // surfSignal->Write("gImpResp");
  // // surfSignalDown->Write("gImpRespDown");
  // // gtemp->Write("gImpRespDownNoise");
  // // out->Close();
  // }
  
  // Cleaning up
  delete surfSignalDown;
  delete surfTrans;
  delete surfSignal;
  delete graphUp;
  delete graph1;
}

void icemc::ChanTrigger::saveTriggerWaveforms(Anita *anita1, double sig0[48], double sig1[48], double noise0[48], double noise1[48]){

  for (int i=0; i<anita1->NFOUR/2; i++){
    sig0[i]   = justSig_trigPath[0][i];
    noise0[i] = justNoise_trigPath[0][i];
    sig1[i]   = justSig_trigPath[1][i];
    noise1[i] = justNoise_trigPath[1][i];
  }
}

void icemc::ChanTrigger::saveDigitizerWaveforms(Anita *anita1, double sig0[48], double sig1[48], double noise0[48], double noise1[48]){

  for (int i=0; i<anita1->NFOUR/2; i++){
    sig0[i]   = justSig_digPath[0][i];
    noise0[i] = justNoise_digPath[0][i];
    sig1[i]   = justSig_digPath[1][i];
    noise1[i] = justNoise_digPath[1][i];
  }
}

void icemc::ChanTrigger::getNoiseFromFlight(Anita* anita1, int ant, bool also_digi){

  Int_t numFreqs = anita1->numFreqs;
  FFTWComplex *phasorsDig  = new FFTWComplex[numFreqs];
  FFTWComplex *phasorsTrig = new FFTWComplex[numFreqs];
  phasorsDig[0].setMagPhase(0,0);
  phasorsTrig[0].setMagPhase(0,0);
  double *freqs = anita1->freqs;
  Double_t sigma, realPart, imPart, trigNorm, digNorm;
  int iring=2;
  if (ant<16) iring=0;
  else if (ant<32) iring=1;

  int iphi = ant - (iring*16);

  for (int ipol=0; ipol<2; ipol++){

    for(int i=1;i<numFreqs;i++) {
      trigNorm       = anita1->fRatioTriggerToA3DigitizerFreqDomain[ipol][iring][iphi][anita1->tuffIndex][i];
      digNorm        = anita1->fRatioDigitizerToA3DigitizerFreqDomain[ipol][iring][iphi][anita1->tuffIndex][i];
      sigma          = anita1->RayleighFits[ipol][ant]->Eval(freqs[i])*4./TMath::Sqrt(numFreqs);
      realPart       = anita1->fRand->Gaus(0,sigma);
      imPart         = anita1->fRand->Gaus(0,sigma);
      
      phasorsDig[i]  = FFTWComplex(realPart*digNorm,  imPart*digNorm  );
      phasorsTrig[i] = FFTWComplex(realPart*trigNorm, imPart*trigNorm );
    }

    RFSignal *rfNoiseDig    = new RFSignal(numFreqs,freqs,phasorsDig,1);    
    Double_t *justNoiseDig  = rfNoiseDig->GetY();
    
    RFSignal *rfNoiseTrig   = new RFSignal(numFreqs,freqs,phasorsTrig,1);
    Double_t *justNoiseTrig = rfNoiseTrig->GetY();

    
    for (int i=0; i<anita1->HALFNFOUR; i++){
      justNoise_digPath[ipol][i]  = also_digi ? justNoiseDig[i]*anita1->THERMALNOISE_FACTOR : 0;
      justNoise_trigPath[ipol][i] = justNoiseTrig[i]*anita1->THERMALNOISE_FACTOR;
    }
    
    delete rfNoiseDig;
    delete rfNoiseTrig;

  }

  // Cleaning up
  delete[] phasorsDig;
  delete[] phasorsTrig;
  
  
}

void icemc::ChanTrigger::calculateCW(Anita *anita1, double frequency, double phase, double amplitude){

  double omega;
  
  for (int itime=0; itime<anita1->HALFNFOUR; itime++){
    omega=TMath::Pi()*2*frequency;
    cw_digPath[0][itime]+=amplitude*TMath::Sin(omega*anita1->TIMESTEP*itime + phase);
    cw_digPath[1][itime]+=amplitude*TMath::Sin(omega*anita1->TIMESTEP*itime + phase);
  }

}  

TGraph *icemc::ChanTrigger::getPulserAtAMPA(Anita *anita1, int ant){

  // phiIndex is 0 for central antenna in trigger efficiency scan
  int phiIndex = anita1->trigEffScanPhi - (ant%16);
  if (phiIndex>8) phiIndex=phiIndex-16;
  double tmp_volts[10000];
  int iring=2;
  if (ant<16) iring=0;
  else if (ant<32) iring=1;
  double norm = 0;
  double att = 0;
  // only 2 phi sectors adjecent to the central one are considered in efficiency scan
  if(TMath::Abs(phiIndex)<=2){
    att   = anita1->trigEffScanAtt[phiIndex+2]-anita1->trigEffScanAtt[2];
    norm  = (anita1->trigEffScanAtt[phiIndex+2]==0)*0 + (anita1->trigEffScanAtt[phiIndex+2]!=0)*1;
    norm *= anita1->trigEffScanRingsUsed[iring];
  }
  int n = anita1->gPulseAtAmpa->GetN();
  n=n/4;
  for (int i=0;i<n;i++){
    tmp_volts[i]=norm*anita1->gPulseAtAmpa->GetY()[i]*TMath::Power(10, att/20.);
  }
  
  return new TGraph(n,  anita1->gPulseAtAmpa->GetX(), tmp_volts);
 
}

void icemc::ChanTrigger::injectImpulseAfterAntenna(Anita *anita1, int ant){
  // phiIndex is 0 for central antenna in trigger efficiency scan
  int phiIndex = anita1->trigEffScanPhi - (ant%16);
  if (phiIndex>8) phiIndex=phiIndex-16;
  int fNumPoints=anita1->HALFNFOUR;
  double tmp_volts[2][1000];
  int iring=2;
  if (ant<16) iring=0;
  else if (ant<32) iring=1;
  // only 2 phi sectors adjecent to the central one are considered in efficiency scan
  if(TMath::Abs(phiIndex)<=2){
    double att = anita1->trigEffScanAtt[phiIndex+2]-anita1->trigEffScanAtt[2];
    double norm = (anita1->trigEffScanAtt[phiIndex+2]==0)*0 + (anita1->trigEffScanAtt[phiIndex+2]!=0)*1;
    norm*=anita1->trigEffScanRingsUsed[iring];
    for (int i=0;i<fNumPoints;i++){
      tmp_volts[0][i]=norm*anita1->trigEffScanPulseAtAmpa[i]*TMath::Power(10, att/20.);
      tmp_volts[1][i]=0;
    }
    for (int i=0;i<fNumPoints;i++){
      volts_rx_forfft[0][4][i]=tmp_volts[0][i];
      volts_rx_forfft[1][4][i]=tmp_volts[1][i];
    }
    
    // put it in normal time ording -T to T
    // instead of 0 to T, -T to 0
    Tools::NormalTimeOrdering(anita1->NFOUR/2,volts_rx_forfft[0][4]);
    Tools::NormalTimeOrdering(anita1->NFOUR/2,volts_rx_forfft[1][4]);


    // find back the frequency domain
    Tools::realft(tmp_volts[0],1,anita1->NFOUR/2);
    Tools::realft(tmp_volts[1],1,anita1->NFOUR/2);
    
    //convert the V pol time waveform into frequency amplitudes
    anita1->GetArrayFromFFT(tmp_volts[0], vhz_rx[0][4]);
    anita1->GetArrayFromFFT(tmp_volts[1], vhz_rx[1][4]);

   
    
  }else{
    for (int i=0;i<fNumPoints;i++){
      volts_rx_forfft[0][4][i]=0;
      volts_rx_forfft[1][4][i]=0;
    }
    for (int i=0;i<anita1->NFREQ;i++){
      vhz_rx[0][4][i]=vhz_rx[1][4][i]=0;
    }
  }
  
  
}

void icemc::ChanTrigger::injectImpulseAtSurf(Anita *anita1, double volts_triggerPath_e[Anita::HALFNFOUR], double volts_triggerPath_h[Anita::HALFNFOUR], int ant){
  // phiIndex is 0 for central antenna in trigger efficiency scan
  int phiIndex = anita1->trigEffScanPhi - (ant%16);
  if (phiIndex>8) phiIndex=phiIndex-16;
  int fNumPoints=anita1->HALFNFOUR;
  // only 2 phi sectors adjecent to the central one are considered in efficiency scan
  if(TMath::Abs(phiIndex)<=2){
    int randomIndex=anita1->fRand->Integer(250);
    double att = anita1->trigEffScanAtt[phiIndex+2];
    double norm = (anita1->trigEffScanAtt[phiIndex+2]==999)*0 + (anita1->trigEffScanAtt[phiIndex+2]!=999)*1;

    for (int i=0;i<fNumPoints;i++){
      volts_triggerPath_e[i]=norm*anita1->trigEffScanPulseAtSurf[randomIndex][i]*TMath::Power(10, att/20);
      volts_triggerPath_h[i]=0;
    }
  }else{
    for (int i=0;i<fNumPoints;i++){
      volts_triggerPath_e[i]=0;
      volts_triggerPath_h[i]=0;
    }
  }
}
#endif
