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

#include "rx.hpp"
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
#endif




void ChanTrigger::ConvertEHtoLREfield(double e_component,double h_component,double& lcp_component,double& rcp_component) {
    
  lcp_component=sqrt((e_component*e_component+h_component*h_component)/2);
  rcp_component=lcp_component;
    
} //ConvertEHtoLREfield
void ChanTrigger::ConvertEHtoLREnergy(double e_component,double h_component,double& lcp_component,double& rcp_component) {
    
  lcp_component=(e_component+h_component)/2;
  rcp_component=lcp_component;
    
} //ConvertEHtoLREnergy

void ChanTrigger::ConvertHVtoLRTimedomain(const int nfour,double *vvolts,
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
    right[2*i]=1/sqrt(2.)*(hvolts_f[2*i]+vvolts_f[2*i+1]);
    left[2*i]=1/sqrt(2.)*(hvolts_f[2*i]-vvolts_f[2*i+1]);
		
    right[2*i+1]=1/sqrt(2.)*(hvolts_f[2*i+1]-vvolts_f[2*i]);
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





//!
/*!
 *
 *
 *
 *
 *
 */
void ChanTrigger::WhichBandsPass(int inu,Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double dangle, double emfrac, double hadfrac){
    
  double thresholds[2][5];
  if (settings1->USETIMEDEPENDENTTHRESHOLDS==1 && settings1->WHICH==9) {
    for(int i=0;i<4;i++) thresholds[0][i] = thresholds[1][i] = anita1->powerthreshold[i];
    int iring = (ilayer<2)*0 + (ilayer==2)*1 + (ilayer==3)*2;
    int iphi = ifold;
    if (ilayer==0) iphi = ifold*2;
    else if (ilayer==1) iphi = ifold*2+1;
    // we invert VPOL and HPOL because in AnitaEventReader (0:HPOL 1:VPOL) and in icemc (0:VPOL 1:HPOL)
    // convert scalers to power thresholds using fitted function got from ANITA-1, are they too old?
    thresholds[1][4] = ADCCountstoPowerThreshold(anita1,0,iring*16+iphi)*(-1.);
    thresholds[0][4] = ADCCountstoPowerThreshold(anita1,1,iring*16+iphi)*(-1.);
    //    cout << thresholds[4] << " \n";
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
    WhichBandsPassTrigger2(inu,settings1, anita1, globaltrig1, bn1, ilayer, ifold, dangle, emfrac, hadfrac, thresholds);

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
void ChanTrigger::WhichBandsPassTrigger1(Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double thresholds[2][5]){

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
    if (settings1->CHMASKING && !ChanTrigger::IsItUnmasked(bn1->surfTrigBandMask,ibw,ilayer,ifold,0)) {
	  
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
      if (settings1->CHMASKING && !ChanTrigger::IsItUnmasked(bn1->surfTrigBandMask,ibw,ilayer,ifold,1)) {
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
void ChanTrigger::WhichBandsPassTrigger2(int inu,Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double dangle, double emfrac, double hadfrac, double thresholds[2][5]){

  double psignal_e[5][anita1->NFOUR];
  double psignal_h[5][anita1->NFOUR];

  double mindiodeconvl_e[5];
  double mindiodeconvl_h[5];
  double onediodeconvl_e[5];
  double onediodeconvl_h[5];
    
  double timedomain_output_1[5][Anita::NFOUR];
  double timedomain_output_2[5][Anita::NFOUR];

  // if we use the diode to perform an integral
  // this is the number of bins to the left of center where the diode function starts to be completely overlapping with the waveform in the convolution.
  int ibinshift=(anita1->NFOUR/4-(int)(anita1->maxt_diode/anita1->TIMESTEP));
      
      
  double integrate_energy_freq[5]={0.,0.,0.,0.,0.};
  for (int iband=0;iband<5;iband++) {
    if (anita1->bwslice_allowed[iband]!=1) continue;
    
    //cout << "arrival time is " << globaltrig1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP << "\n";
    //anita1->iminbin[j]=anita1->NFOUR/4-ibinshift+anita1->idelaybeforepeak[j]+globaltrig1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP; // we start to look for single channel triggers firing
    //anita1->iminbin[iband]=anita1->NFOUR/4-ibinshift+anita1->idelaybeforepeak[iband]; // we start to look for single channel triggers firing
    anita1->iminbin[iband]=0.; // we start to look for single channel triggers firing
    // starting with the first bin where the diode function is completely
    // overlapping plus a delay that is brought about by the diode
    // idelaybeforepeak is the bin number in the diode response function
    // where it peaks
    //
    //
    //
    //
    // E. Hong added comment in 040412 : currently the code looks at the same bin for the peak (from iminbin to imaxbin) for all antennas
    // However, as there is a code which account for the time delay between antennas (Tools::ShiftRight), I think iminbin and imaxbin also have to account for that.
    // I think L, M, H bwslice could be fine with sharing same iminbin and imaxbin as total bin we are looking (imaxbin - iminbin) is large (anita1->iwindow = 20ns), but for the Full band sharing same iminbin and imaxbin with all antennas can cause missing the peak between iminbin and imaxbin as the monitoring bin width for full band is small (anita1->iwindow = 4ns)
    // Hope someone can check the part and either fix the code or confirm that the code is fine.
    //
    //
    //anita1->imaxbin[j]=anita1->NFOUR/4-ibinshift+anita1->idelaybeforepeak[j]+anita1->iwindow[j];
    //anita1->imaxbin[j]=anita1->NFOUR/4-ibinshift+anita1->idelaybeforepeak[j]+anita1->iwindow[j]+(globaltrig1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP);
    anita1->imaxbin[iband]=anita1->NFOUR/2;
	
    //Matt's plot:  horiz. axis:  peak_signal_rfcm/rms_rfcm
    // vert. axis:  peak_signal_lab/rms_lab
    Tools::Zero(v_banding_rfcm_e_forfft[iband],anita1->NFOUR/2);
    Tools::Zero(v_banding_rfcm_h_forfft[iband],anita1->NFOUR/2);
	
    anita1->MakeArraysforFFT(v_banding_rfcm_e[iband],v_banding_rfcm_h[iband],v_banding_rfcm_e_forfft[iband],v_banding_rfcm_h_forfft[iband], 90., true);
	
    // for some reason I'm averaging over 10 neighboring bins
    // to get rid of the zero bins
    for (int i=0;i<anita1->NFOUR/4;i++) {
	  
      v_banding_rfcm_e_forfft_temp[iband][2*i]  =0.;
      v_banding_rfcm_e_forfft_temp[iband][2*i+1]=0.;
      v_banding_rfcm_h_forfft_temp[iband][2*i]  =0.;
      v_banding_rfcm_h_forfft_temp[iband][2*i+1]=0.;

      int tempcount = 0;
      for (int k=i;k<i+10;k++) {
	if (k<anita1->NFOUR/4) {
	  v_banding_rfcm_e_forfft_temp[iband][2*i]  +=v_banding_rfcm_e_forfft[iband][2*k];
	  v_banding_rfcm_e_forfft_temp[iband][2*i+1]+=v_banding_rfcm_e_forfft[iband][2*k+1];
	  v_banding_rfcm_h_forfft_temp[iband][2*i]  +=v_banding_rfcm_h_forfft[iband][2*k];
	  v_banding_rfcm_h_forfft_temp[iband][2*i+1]+=v_banding_rfcm_h_forfft[iband][2*k+1];
	  tempcount++;
	}
      }
	  
      v_banding_rfcm_e_forfft[iband][2*i]  =v_banding_rfcm_e_forfft_temp[iband][2*i]/tempcount;
      v_banding_rfcm_e_forfft[iband][2*i+1]=v_banding_rfcm_e_forfft_temp[iband][2*i+1]/tempcount;
      v_banding_rfcm_h_forfft[iband][2*i]  =v_banding_rfcm_h_forfft_temp[iband][2*i]/tempcount;
      v_banding_rfcm_h_forfft[iband][2*i+1]=v_banding_rfcm_h_forfft_temp[iband][2*i+1]/tempcount;
	  
      v_banding_rfcm_e_forfft_temp[iband][2*i]  =v_banding_rfcm_e_forfft[iband][2*i];
      v_banding_rfcm_e_forfft_temp[iband][2*i+1]=v_banding_rfcm_e_forfft[iband][2*i+1];
      v_banding_rfcm_h_forfft_temp[iband][2*i]  =v_banding_rfcm_h_forfft[iband][2*i];
      v_banding_rfcm_h_forfft_temp[iband][2*i+1]=v_banding_rfcm_h_forfft[iband][2*i+1];
    }
	
    Tools::realft(v_banding_rfcm_e_forfft[iband],-1,anita1->NFOUR/2);
    // now v_banding_rfcm_e_forfft is in the time domain
    // and now it is really in units of V
	
	
    Tools::realft(v_banding_rfcm_h_forfft[iband],-1,anita1->NFOUR/2);
    // now v_banding_rfcm_h_forfft is in the time domain
    // and now it is really in units of V

       
    // put it in normal time ording -T to T
    // instead of 0 to T, -T to 0
    Tools::NormalTimeOrdering(anita1->NFOUR/2,v_banding_rfcm_e_forfft[iband]);
    Tools::NormalTimeOrdering(anita1->NFOUR/2,v_banding_rfcm_h_forfft[iband]);
    
    if (settings1->APPLYIMPULSERESPONSETRIGGER){
      applyImpulseResponseTrigger(settings1, anita1, anita1->GetRxTriggerNumbering(ilayer, ifold), v_banding_rfcm_e_forfft[iband], 0);
      applyImpulseResponseTrigger(settings1, anita1, anita1->GetRxTriggerNumbering(ilayer, ifold), v_banding_rfcm_h_forfft[iband], 1);
    }
    

    
    if (settings1->TRIGGEREFFSCAN && (settings1->TRIGGEREFFSCAPULSE==1)){  
      injectImpulseAtSurf(anita1, v_banding_rfcm_e_forfft[iband], v_banding_rfcm_h_forfft[iband], anita1->GetRxTriggerNumbering(ilayer, ifold));
    }
    
    
    if (settings1->ZEROSIGNAL) {
	 
      Tools::Zero(v_banding_rfcm_e_forfft[iband],anita1->NFOUR/2);
      Tools::Zero(v_banding_rfcm_h_forfft[iband],anita1->NFOUR/2);
	  
    }
	
	
    for (int i=0;i<anita1->NFOUR/4;i++) {
      integrate_energy_freq[iband]+=v_banding_rfcm_e_forfft[iband][2*i]*v_banding_rfcm_e_forfft[iband][2*i]+v_banding_rfcm_e_forfft[iband][2*i+1]*v_banding_rfcm_e_forfft[iband][2*i+1];
    }
	
    // write the signal events to a tree
    for (int k=0;k<anita1->NFOUR/2;k++) {
      anita1->signal_vpol_inanita[iband][k]=v_banding_rfcm_e_forfft[iband][k];
    }
    anita1->integral_vmmhz_foranita=integral_vmmhz;
	
    // Find the p2p value before adding noise
    anita1->peak_v_banding_rfcm_e[iband]=FindPeak(v_banding_rfcm_e_forfft[iband],anita1->NFOUR/2);
    anita1->peak_v_banding_rfcm_h[iband]=FindPeak(v_banding_rfcm_h_forfft[iband],anita1->NFOUR/2);
	
  } // loop over bands
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
	anita1->total_vpol_inanita[iband][itime]=anita1->timedomainnoise_rfcm_banding_e[iband][itime]+anita1->signal_vpol_inanita[iband][itime];

	integrateenergy[iband]+=anita1->timedomainnoise_rfcm_banding_e[iband][itime]*anita1->timedomainnoise_rfcm_banding_e[iband][itime]*anita1->TIMESTEP;
	if (settings1->SIGNAL_FLUCT && (!settings1->NOISEFROMFLIGHTTRIGGER)) {
	  // this reverses the noise is time, and starts with bin anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP)
	  v_banding_rfcm_e_forfft[iband][itime]=v_banding_rfcm_e_forfft[iband][itime]+anita1->timedomainnoise_rfcm_banding_e[iband][itimenoisebin];
	  v_banding_rfcm_h_forfft[iband][itime]=v_banding_rfcm_h_forfft[iband][itime]+anita1->timedomainnoise_rfcm_banding_h[iband][itimenoisebin];
	  }
      }
      for (int itime=anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP);itime<anita1->NFOUR/2;itime++) {
	anita1->total_vpol_inanita[iband][itime]=0.;
	v_banding_rfcm_e_forfft[iband][itime]=0.;
	v_banding_rfcm_h_forfft[iband][itime]=0.;
      }
    }
    else {
      for (unsigned int itime = 0; itime < anita1->HALFNFOUR; ++itime){
	// this is just a straight sum
	anita1->total_vpol_inanita[iband][itime]=anita1->timedomainnoise_rfcm_banding_e[iband][itime]+anita1->signal_vpol_inanita[iband][itime];
	integrateenergy[iband]+=anita1->timedomainnoise_rfcm_banding_e[iband][itime]*anita1->timedomainnoise_rfcm_banding_e[iband][itime]*anita1->TIMESTEP;
	// this one is the one actually used by the diode
	v_banding_rfcm_e_forfft[iband][itime] += anita1->timedomainnoise_rfcm_banding_e[iband][itime];
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
      ConvertHVtoLRTimedomain(anita1->NFOUR, v_banding_rfcm_e_forfft[iband], v_banding_rfcm_h_forfft[iband], vm_banding_rfcm_1_forfft[iband], vm_banding_rfcm_2_forfft[iband]);
      
    } else {
  
      for (int itime=0;itime<anita1->NFOUR/2;itime++) {
	    
	vm_banding_rfcm_1_forfft[iband][itime] = v_banding_rfcm_e_forfft[iband][itime];
	vm_banding_rfcm_2_forfft[iband][itime] = v_banding_rfcm_h_forfft[iband][itime];
	    
      }
    }

    for (int itime=0;itime<anita1->NFOUR/2;itime++) {
      anita1->total_diodeinput_1_inanita[iband][itime] = vm_banding_rfcm_1_forfft[iband][itime];
      anita1->total_diodeinput_2_inanita[iband][itime] = vm_banding_rfcm_2_forfft[iband][itime];
	  
    }
    //volts_fullband_trigger_path_e.push_back(;
  } // end loop over bands
      
            
  for (int itime=0;itime<Anita::NFOUR/2;itime++) {
    anita1->total_diodeinput_1_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)][itime]=anita1->total_diodeinput_1_inanita[4][itime];
    anita1->total_diodeinput_2_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)][itime]=anita1->total_diodeinput_2_inanita[4][itime];
  }

    
      

  for (int iband=0;iband<5;iband++) {
    if (anita1->bwslice_allowed[iband]!=1) continue; 

    // myconvl
    // this performs the convolution with the diode response

    anita1->myconvlv(vm_banding_rfcm_1_forfft[iband],anita1->NFOUR,anita1->fdiode_real[iband],mindiodeconvl_e[iband],onediodeconvl_e[iband],psignal_e[iband],timedomain_output_1[iband]);
    // loop from the ibinshift left + some delay + 10 ns
	
    // now shift right to account for arrival times
    Tools::ShiftRight(timedomain_output_1[iband],anita1->NFOUR,(int)(anita1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));
	
	
    if (settings1->TRIGGERSCHEME == 2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5){
      //      if (inu==1570)
      //cout << "shifting left.\n";
      Tools::ShiftLeft(timedomain_output_1[iband],anita1->NFOUR,ibinshift);
    }

	
    for (int itime=0;itime<Anita::NFOUR/2;itime++) {
      anita1->timedomain_output_1_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)][itime]=timedomain_output_1[4][itime];
      //cerr<<iband<<" "<<i<<" "<<vm_banding_rfcm_1_forfft<<"  "<<timedomain_output_1[iband][i]<<endl;
    }
	
    anita1->channels_passing_e[iband]=0; // does not pass
	
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
    
      if (timedomain_output_1[iband][ibin] < thresholds[0][iband] * anita1->bwslice_rmsdiode[iband] && anita1->pol_allowed[0] && anita1->bwslice_allowed[iband]) { // is this polarization and bw slice allowed to pass
	//std::cout << "VPOL : " << iband << " " << timedomain_output_1[iband][ibin]  << " " <<  thresholds[0][iband] << " " <<  anita1->bwslice_rmsdiode[iband] << std::endl;
	anita1->channels_passing_e[iband] = 1;// channel passes
    
	
	thisisaone=1;
	//	  cout << "got a hit.\n";
      }
      // if (whichlayer==2 && whichphisector==15) 
      //cout << "ibin, nextbreakpoint are " << ibin << "\t" << nextbreakpoint << "\n";
      if (ibin>nextbreakpoint) {
	//if (whichlayer==2 && whichphisector==15) 
	//cout << "I'm filling. size is " << "\t" << globaltrig1->arrayofhits[whichlayer][whichphisector][0][iband].size() << "\n";

	globaltrig1->arrayofhits[whichlayer][whichphisector][0][iband].push_back(thisisaone);
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
    
    if (anita1->channels_passing_e[iband]) {
      globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]++; //Records number of first level triggers on each antenna for a single neutrino
	  
    }
  } // end loop over bands
  //  if (whichlayer==2 && whichphisector==15)
  if (globaltrig1->arrayofhits[whichlayer][whichphisector][0][4].size()==0)
    cout << "inu, whichlayer, whichphisector, size is " << whichlayer << "\t" << whichphisector << "\t" << globaltrig1->arrayofhits[whichlayer][whichphisector][0][4].size() << "\n";
   
  for (int iband=0;iband<5;iband++) {
    if (anita1->bwslice_allowed[iband]!=1) continue; 
    // Find p2p value before adding noise
    anita1->peak_v_banding_rfcm_h[iband]=FindPeak(v_banding_rfcm_h_forfft[iband],anita1->NFOUR/2);
  }
      
  for (int iband=0;iband<5;iband++) {
    if (anita1->bwslice_allowed[iband]!=1) continue; 
    anita1->myconvlv(vm_banding_rfcm_2_forfft[iband],anita1->NFOUR,anita1->fdiode_real[iband],mindiodeconvl_h[iband],onediodeconvl_h[iband],psignal_h[iband],timedomain_output_2[iband]);
    // now shift right to account for arrival times
    Tools::ShiftRight(timedomain_output_2[iband],anita1->NFOUR,(int)(anita1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));

	
    if (settings1->TRIGGERSCHEME == 2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5){
      //      if (inu==1570)
      //cout << "shifting left.\n";
      Tools::ShiftLeft(timedomain_output_2[iband],anita1->NFOUR,ibinshift); 
    }




    for (int i=0;i<Anita::NFOUR/2;i++) {
      anita1->timedomain_output_2_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)][i]=timedomain_output_2[4][i];
    }


    anita1->channels_passing_h[iband]=0;


    int whichtrigbin=0;
    //cout << "1 whichtrigbin is " << whichtrigbin << "\n";
    // nextbreakpoint is how many time steps in digitization that we have to go 
    // to equal the length we want for the trigger bin
    int nextbreakpoint=(int)((double)(whichtrigbin+1)*(globaltrig1->TRIGTIMESTEP/anita1->TIMESTEP));

    //    cout << "1 nextbreakpoint is " << nextbreakpoint << "\n";
    // keep advancing in trigger timesteps until we're past the iminbin.
    while (nextbreakpoint<anita1->iminbin[iband]) {
      nextbreakpoint=(int)((double)(whichtrigbin+1)*(globaltrig1->TRIGTIMESTEP/anita1->TIMESTEP));
      whichtrigbin++;
    }
    // keep track of whether each trigger bin has a hit
    int thisisaone=0;

    for (int ibin=anita1->iminbin[iband];ibin<anita1->imaxbin[iband];ibin++) {
      if (timedomain_output_2[iband][ibin]<thresholds[1][iband]*anita1->bwslice_rmsdiode[iband] && anita1->pol_allowed[1] && anita1->bwslice_allowed[iband]) { // if this pol and band are allowed to pass
	//	      std::cout << "HPOL : " << iband << " " << timedomain_output_2[iband][ibin]  << " " << timedomain_output_1[iband][ibin] << " " <<  thresholds[1][iband] << " " <<  anita1->bwslice_rmsdiode[iband] << std::endl;
	anita1->channels_passing_h[iband]=1;
	  
	thisisaone=1;
      } // if it's over threshold for this time step
      if (ibin>nextbreakpoint) {
	globaltrig1->arrayofhits[whichlayer][whichphisector][1][iband].push_back(thisisaone);
	// 	if (thisisaone) {
	// 	  cout << "filling with 1. phi is " << whichphisector << "\n";
	// 	}
	//cout << "filling arrayofhits with " << thisisaone << "\n";
	//cout << "size of arrayofhits is " << globaltrig1->arrayofhits[ilayer][ifold][0][iband].size() << "\n";
	nextbreakpoint=(int)((double)(whichtrigbin+1)*(globaltrig1->TRIGTIMESTEP/anita1->TIMESTEP));
	//cout << "3 nextbreakpoint is " << nextbreakpoint << "\n";
	whichtrigbin++;
	//cout << "3 whichtrigbin is " << whichtrigbin << "\n";
	thisisaone=0;
      }
    }
	
    if (anita1->channels_passing_h[iband])
      globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]++; //Records number of first level triggers on each antenna for a single neutrino


  } // end loop over bands
  
  
  // fill channels_passing
  //  } // end if the signal is big enough the be considered
  // now we've made the diode outputs
  // now step in time and find the most number of bands that pass in the
  // right time window
  // then fills channels_passing
  int npass;

  L1Trigger(anita1,timedomain_output_1,timedomain_output_2,thresholds, //inputs
	    globaltrig1->channels_passing[ilayer][ifold][0],globaltrig1->channels_passing[ilayer][ifold][1],npass); //outputs

      
  // if it's the closest antenna,
  // save flag_e,h in anita class for writing to tsignals tree
  int startbin=TMath::MinElement(5,anita1->iminbin);
      
  if (ilayer==anita1->GetLayer(anita1->rx_minarrivaltime) && ifold==anita1->GetIfold(anita1->rx_minarrivaltime)) {
    for (int iband=0;iband<5;iband++) {
      if (anita1->bwslice_allowed[iband]!=1) continue; 

      //	  cout << "zeroeing here 1.\n";
      anita1->ston[iband]=0.;
      for (int i=anita1->iminbin[iband];i<anita1->imaxbin[iband];i++) {
	// 	    if (iband==0 && i==anita1->NFOUR/4) {
	// 	      cout << "output is " << inu << "\t" << timedomain_output_1[iband][i] << "\n";
	// 	    }
	//	    cout << "output, bwslice_rmsdiode are " << timedomain_output_1[iband][i] << "\t" << anita1->bwslice_rmsdiode[iband] << "\n";
	if (timedomain_output_1[iband][i]/anita1->bwslice_rmsdiode[iband]<anita1->ston[iband]) {
	  anita1->ston[iband]=timedomain_output_1[iband][i]/anita1->bwslice_rmsdiode[iband];
	  //  if (iband==4 && anita1->ston[iband]<0.)
	  //cout << "ston is " << anita1->ston[iband] << "\n";
	}
      }


      for (int i=0;i<anita1->HALFNFOUR;i++) {
	anita1->flag_e_inanita[iband][i]=0;
	anita1->flag_h_inanita[iband][i]=0;
	anita1->timedomain_output_1_inanita[iband][i]=timedomain_output_1[iband][i];
	anita1->timedomain_output_2_inanita[iband][i]=timedomain_output_2[iband][i];

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
    if (settings1->CHMASKING && !ChanTrigger::IsItUnmasked(bn1->surfTrigBandMask,iband,ilayer,ifold,0)) {
      if (globaltrig1->channels_passing[ilayer][ifold][0][iband])
	globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]--;
      globaltrig1->channels_passing[ilayer][ifold][0][iband]=0;// channel does not pass
    }
	
    // if we're implementing masking and the channel has been masked
    if (settings1->CHMASKING && !ChanTrigger::IsItUnmasked(bn1->surfTrigBandMask,iband,ilayer,ifold,1)) {
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
      
  anita1->dangle_inanita=dangle; // viewangle - changle
  anita1->emfrac_inanita=emfrac; // viewangle - changle
  anita1->hadfrac_inanita=hadfrac; // viewangle - changle
      

      
  anita1->irx=anita1->GetRx(ilayer,ifold);
      
  //       if (anita1->tdata->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1) {
  // 	anita1->tdata->Fill();
  // 	//for (int i=0;i<512;i++) {
  // 	//cout << "i, total_diodeinput_1_inanita[iband][i] are " << i << "\t" << anita1->total_diodeinput_1_inanita[4][i] << "\n";
  // 	//cout << "inu is " << inu << "\n";
  // 	//}
  // 	}
  if (Anita::GetLayer(anita1->rx_minarrivaltime)==ilayer && Anita::GetIfold(anita1->rx_minarrivaltime)==ifold && anita1->tsignals->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1) {
    //if (Anita::GetLayer(globaltrig1->rx_minarrivaltime)==ilayer && Anita::GetIfold(globaltrig1->rx_minarrivaltime)==ifold && anita1->tsignals->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1) {
	
    // 	for (int i=0;i<anita1->NFOUR/2;i++) {
    // 	  cout << "signal_vpol_inanita is " << anita1->signal_vpol_inanita[4][i] << "\n";
    // 	}

    anita1->tsignals->Fill();
    //cout << "filling.\n";
  }
      
      
  if (settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5){
    // If TRIGGERSCHEME == 3, allow all bands to pass here. This allows the coherent waveform sum trigger scheme to be used.
    for (unsigned int ichannel = 0; ichannel < 5; ichannel++){
      globaltrig1->channels_passing[ilayer][ifold][0][ichannel] = 1;
      globaltrig1->channels_passing[ilayer][ifold][1][ichannel] = 1;
      anita1->channels_passing_e[ichannel] = 1;
      anita1->channels_passing_h[ichannel] = 1;
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


void ChanTrigger::InitializeEachBand(Anita *anita1)
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


void ChanTrigger::ConvertInputWFtoAntennaWF(Settings *settings1, Anita *anita1, Balloon *bn1, Screen *panel1, Vector n_eplane, Vector n_hplane, Vector n_normal, int ilayer, int ifold)
{
  // vmmhz is V/m/MHz at the face of the antenna
  // this gets written to a tree because it is a measure of signal strength in the frequency domain
  integral_vmmhz=0.;
  //used if roughness, to first determine total of a frequency for all screen points
  for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
    integral_vmmhz_r[ifreq] = 0.;
  }

  if(settings1->ROUGHNESS==0){
    for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
      integral_vmmhz+=panel1->GetVmmhz_freq(ifreq)*(anita1->freq[1]-anita1->freq[0])/1.E6; // integrate vmmhz
    }
  }
  else{ // for each frequency, add all screen points
    int jf;
    for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
      for (int npts=0; npts<panel1->GetNvalidPoints(); npts++){
        jf = (npts * anita1->NFREQ) + ifreq;
        integral_vmmhz_r[ifreq] += panel1->GetVmmhz_freq(jf); // integrate vmmhz
      }
    }
    for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
      integral_vmmhz += integral_vmmhz_r[ifreq] * (anita1->freq[1]-anita1->freq[0])/1.E6; // integrate vmmhz
    }
  }

  // zero the necessary variables and arrays
  e_component=0;
  h_component=0;
  n_component=0;
  e_component_kvector=0;
  h_component_kvector=0;
  n_component_kvector=0;
  hitangle_e=0;
  hitangle_h=0;

  for (int i=0; i<Anita::NFREQ; i++){
    vhz_rx_e[i]=0.;
    vhz_rx_h[i]=0.;
  }
  for (int i=0; i<Anita::HALFNFOUR; i++){
    volts_rx_e_forfft[i]=0.;
    volts_rx_h_forfft[i]=0.;
  }

  int fNumPoints = anita1->HALFNFOUR;
  int ant = anita1->GetRxTriggerNumbering(ilayer, ifold);

  double tmp_vhz_rx_e[Anita::NFREQ]; // V/Hz after antenna gains
  double tmp_vhz_rx_h[Anita::NFREQ];
  double tmp_volts_rx_e_forfft[Anita::HALFNFOUR];
  double tmp_volts_rx_h_forfft[Anita::HALFNFOUR];



  // vmmhz_rx_e,h are going to be the V/m/MHz received by the rx (after gains)
  //
  // for each screen point, assemble the waveform and FFT it with the phase information to get its time-domain waveform, which we THEN add to the running total
  for (int jpt=0; jpt < panel1->GetNvalidPoints(); jpt++){
    for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
      // Convert V/m/MHz to V/m/Hz and divide by dt to prepare for fft
      //vhz_rx_e[ifreq]=vmmhz[ifreq]/sqrt(2.)/(anita1->TIMESTEP*1.E6); // EH, 1/sqrt(2) for dividing power in half for TDA and DDA?
      //vhz_rx_h[ifreq]=vmmhz[ifreq]/sqrt(2.)/(anita1->TIMESTEP*1.E6); 
      tmp_vhz_rx_e[ifreq] = panel1->GetVmmhz_freq(jpt*Anita::NFREQ + ifreq)/sqrt(2)/(anita1->TIMESTEP*1.E6);
      tmp_vhz_rx_h[ifreq] = panel1->GetVmmhz_freq(jpt*Anita::NFREQ + ifreq)/sqrt(2)/(anita1->TIMESTEP*1.E6);

      bn1->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal, panel1->GetVec2bln(jpt), e_component_kvector,  h_component_kvector,  n_component_kvector);
      bn1->GetEcompHcompEvector(settings1,  n_eplane,  n_hplane,  panel1->GetPol(jpt),  e_component,  h_component,  n_component);
      bn1->GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector, hitangle_e, hitangle_h);

      anita1->AntennaGain(settings1, hitangle_e, hitangle_h, e_component, h_component, ifreq, tmp_vhz_rx_e[ifreq], tmp_vhz_rx_h[ifreq]);
    }

    if (settings1->TRIGGEREFFSCAN && (settings1->TRIGGEREFFSCAPULSE==0)){
      injectImpulseAmplitudeAfterAntenna(anita1, tmp_vhz_rx_e, tmp_vhz_rx_h, ant);
    }

    // change their length from Anita::NFREQ to HALFNFOUR
    anita1->MakeArraysforFFT(tmp_vhz_rx_e, tmp_vhz_rx_h, tmp_volts_rx_e_forfft, tmp_volts_rx_h_forfft, 90., false);// 90 is just a placeholder
    //need to handle phase delay explicitly here (eventually when roughness gets worked out)
    for (unsigned int ifour=0;ifour<NFOUR/4;ifour++) {
      tmp_volts_rx_e_forfft[2*ifour]*=cos((90.)*PI/180.);
      tmp_volts_rx_e_forfft[2*ifour+1]*=sin((90.)*PI/180.);
      tmp_volts_rx_h_forfft[2*ifour]*=cos((90.)*PI/180.);
      tmp_volts_rx_h_forfft[2*ifour+1]*=sin((90.)*PI/180.); 
    }


    // now the last two are in the frequency domain
    // convert to the time domain
    Tools::realft(tmp_volts_rx_e_forfft,-1,anita1->HALFNFOUR);
    Tools::realft(tmp_volts_rx_h_forfft,-1,anita1->HALFNFOUR);

    for (int ii=0; ii<Anita::HALFNFOUR; ii++){
      volts_rx_e_forfft[ii] += tmp_volts_rx_e_forfft[ii] * panel1->GetWeight(jpt) / panel1->GetWeightNorm();
      volts_rx_h_forfft[ii] += tmp_volts_rx_h_forfft[ii] * panel1->GetWeight(jpt) / panel1->GetWeightNorm();
    }
  }//end int jpt loop over screen


  // now volts_rx_e,h_forfft are the FULL time domain signals FROM ALL SCREEN POINTS on the back end of the antenna
  // this is in volts
    
  // now find peak voltage
  // these get written to a tree
  // we still don't have any noise
  anita1->peak_e_rx_signalonly=ChanTrigger::FindPeak(volts_rx_e_forfft,anita1->HALFNFOUR); // with no noise
  anita1->peak_h_rx_signalonly=ChanTrigger::FindPeak(volts_rx_h_forfft,anita1->HALFNFOUR); // with no noise
}


void ChanTrigger::PrepareTriggerPath(Settings *settings1, Anita *anita1, Screen *panel1, int ilayer, int ifold, double hitangle_e, double hitangle_h, double e_component, double h_component){

  int fNumPoints = anita1->HALFNFOUR;
  int ant = anita1->GetRxTriggerNumbering(ilayer, ifold);

  for (int iband=0;iband<5;iband++) { // loop over bands
    if (anita1->bwslice_allowed[iband]!=1) continue;
    
    anita1->iminbin[iband]=0.;
    anita1->imaxbin[iband]=anita1->NFOUR/2;
    
    for (int i=0;i<anita1->NFOUR/2;i++) {
      v_banding_rfcm_e_forfft[iband][i] = 0.;
      v_banding_rfcm_h_forfft[iband][i] = 0.;
    }
    
    for (int jpt=0; jpt<panel1->GetNvalidPoints(); jpt++){
      //get the orientation for this screen point
      for (int i=0;i<Anita::NFREQ;i++) {
	anita1->vmmhz_banding[i]=panel1->GetVmmhz_freq(jpt*Anita::NFREQ + i);
      }
                
      // impose banding on the incident signal
      if (!settings1->APPLYIMPULSERESPONSETRIGGER){
	anita1->Banding(iband,anita1->freq,anita1->vmmhz_banding,Anita::NFREQ); 
      }

      for (int i=0;i<Anita::NFREQ;i++) {
	anita1->vmmhz_banding_rfcm[i]=anita1->vmmhz_banding[i];
      }
                
      // for frequency-domain voltage-based trigger (triggerscheme==0)
      // and when using impulse response (APPLYIMPULSERESPONSETRIGGER==1)
      // we do not apply rfcm's
      // for other trigger types we do
      if ((settings1->TRIGGERSCHEME==1 || settings1->TRIGGERSCHEME==2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5) && (!settings1->APPLYIMPULSERESPONSETRIGGER))
	anita1->RFCMs(1,1,anita1->vmmhz_banding_rfcm);
      
      if (settings1->TRIGGERSCHEME >=2) { // we need to prepar the signal for the diode integration
	for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
	  anita1->vmmhz_banding_rfcm[ifreq]=anita1->vmmhz_banding_rfcm[ifreq]/sqrt(2)/(anita1->TIMESTEP*1.E6);
	}
      }
      
      for (int k=0;k<Anita::NFREQ;k++) {
	if (anita1->freq[k]>=settings1->FREQ_LOW_SEAVEYS && anita1->freq[k]<=settings1->FREQ_HIGH_SEAVEYS){
	  // need to calculate lcp and rcp components after antenna voltages are recorded.
	  v_banding_rfcm_e[iband][k]=anita1->vmmhz_banding_rfcm[k];
	  v_banding_rfcm_h[iband][k]=anita1->vmmhz_banding_rfcm[k];
	  anita1->AntennaGain(settings1, hitangle_e, hitangle_h, e_component, h_component, k, v_banding_rfcm_e[iband][k], v_banding_rfcm_h[iband][k]);
	} // end if (seavey frequencies)
      } // end looping over frequencies.
      
      if (settings1->TRIGGEREFFSCAN && (settings1->TRIGGEREFFSCAPULSE==0)){
	injectImpulseAmplitudeAfterAntenna(anita1, v_banding_rfcm_e[iband], v_banding_rfcm_h[iband], ant);
	// if not using the impulse response we need to re-apply banding and rfcms
	if (!settings1->APPLYIMPULSERESPONSETRIGGER){
	  anita1->Banding(iband, anita1->freq, v_banding_rfcm_e[iband], Anita::NFREQ);
	  anita1->RFCMs(1, 1, v_banding_rfcm_e[iband]);
	}
      }
      
      // Currently not used, but don't throw it away just yet
      //       if (settings1->APPLYIMPULSERESPONSETRIGGER){
      // 	double volts_triggerPath_e[Anita::HALFNFOUR]={0.};
      // 	double volts_triggerPath_h[Anita::HALFNFOUR]={0.};
      // 	double vhz_triggerPath_e[Anita::NFREQ] = {0.};
      // 	double vhz_triggerPath_h[Anita::NFREQ] = {0.};
	
      // 	anita1->MakeArraysforFFT(v_banding_rfcm_e[iband], v_banding_rfcm_h[iband], volts_triggerPath_e, volts_triggerPath_h, 90., true);
	
      // 	// for the ROUGHNESS case, need to apply phase factors here somehow, like in ConvertInputWFtoAntennaWF() above in the digitization path
      // 	//for (int ifour=0;ifour<Anita::NFOUR/4;ifour++) {
      // 	//  volts_triggerPath_e[2*ifour] *= cos( (90.)*PI/180.);
      // 	//  volts_triggerPath_e[2*ifour+1] *= sin( (90.)*PI/180.);
      // 	//  volts_triggerPath_h[2*ifour] *= cos( (90.)*PI/180.);
      // 	//  volts_triggerPath_h[2*ifour+1] *= sin( (90.)*PI/180.);
      // 	//}
	
      // 	Tools::realft(volts_triggerPath_e,-1,anita1->NFOUR/2);
      // 	// now v_banding_rfcm_e_forfft is in the time domain
      // 	// and now it is really in units of V
	
      // 	Tools::realft(volts_triggerPath_h,-1,anita1->NFOUR/2);
      // 	// now v_banding_rfcm_h_forfft is in the time domain
      // 	// and now it is really in units of V
	
      // 	// put it in normal time ording -T to T
      // 	// instead of 0 to T, -T to 0
      // 	Tools::NormalTimeOrdering(anita1->NFOUR/2,volts_triggerPath_e);
      // 	Tools::NormalTimeOrdering(anita1->NFOUR/2,volts_triggerPath_h);
	
      // 	if (settings1->TRIGGEREFFSCAN && (settings1->TRIGGEREFFSCAPULSE==0)){
      // 	  injectImpulseAfterAntenna(anita1, volts_triggerPath_e, volts_triggerPath_h, ant);
      // 	}
	
      // #ifdef ANITA_UTIL_EXISTS    
      // 	applyImpulseResponseTrigger(settings1, anita1, fNumPoints, ant, anita1->fTimes, volts_triggerPath_e, vhz_triggerPath_e, 0);
      // 	applyImpulseResponseTrigger(settings1, anita1, fNumPoints, ant, anita1->fTimes, volts_triggerPath_h, vhz_triggerPath_h, 1);
	
      // 	for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
      // 	  if (anita1->freq[ifreq]>=settings1->FREQ_LOW_SEAVEYS && anita1->freq[ifreq]<=settings1->FREQ_HIGH_SEAVEYS){
      // 	    v_banding_rfcm_e[iband][ifreq]=vhz_triggerPath_e[ifreq];
      // 	    v_banding_rfcm_h[iband][ifreq]=vhz_triggerPath_h[ifreq];
      // 	  }
      // 	} // end loop over nfreq
      // #endif
	
      // 	// now add the screen point's waveform to the total including the weighting 
      // 	for (int i=0;i<anita1->NFOUR/2;i++) {
      // 	  v_banding_rfcm_e_forfft[iband][i] += volts_triggerPath_e[i] * panel1->GetWeight(jpt) / panel1->GetWeightNorm();
      // 	  v_banding_rfcm_h_forfft[iband][i] += volts_triggerPath_h[i] * panel1->GetWeight(jpt) / panel1->GetWeightNorm();
      // 	}

      //       }//if (settings1->APPLYIMPULSERESPONSETRIGGER)
      
      for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
	if (anita1->freq[ifreq]>=settings1->FREQ_LOW_SEAVEYS && anita1->freq[ifreq]<=settings1->FREQ_HIGH_SEAVEYS){
	  addToChannelSums(settings1, anita1, iband, ifreq);
	}
      } // end loop over nfreq
      
      
      //
    } // end loop over screen points
    
    // write the signal events to a tree
    for (int iband=0;iband<5;iband++) {
      if (anita1->bwslice_allowed[iband]!=1) continue; 
      for (int k=0;k<anita1->NFOUR/2;k++) {
	anita1->signal_vpol_inanita[iband][k]=v_banding_rfcm_e_forfft[iband][k];
      }
    }
    anita1->integral_vmmhz_foranita=integral_vmmhz;
    
    // Find the p2p value before adding noise
    for (int iband=0;iband<5;iband++) {
      if (anita1->bwslice_allowed[iband]!=1) continue; 
      anita1->peak_v_banding_rfcm_e[iband]=FindPeak(v_banding_rfcm_e_forfft[iband],anita1->NFOUR/2);
      anita1->peak_v_banding_rfcm_h[iband]=FindPeak(v_banding_rfcm_h_forfft[iband],anita1->NFOUR/2);
    }
  } // end loop over bands
  
}

void ChanTrigger::DigitizerPath(Settings *settings1, Anita *anita1, int ilayer, int ifold)
{
  double vhz_rx_rfcm_e[Anita::NFREQ]; // V/Hz after rx, rfcm
  double vhz_rx_rfcm_h[Anita::NFREQ];

  int fNumPoints = anita1->HALFNFOUR;
  int ant = anita1->GetRxTriggerNumbering(ilayer, ifold);

  
  // if (settings1->TRIGGEREFFSCAN && (settings1->TRIGGEREFFSCAPULSE==0)){
  //   injectImpulseAfterAntenna(anita1, volts_rx_e_forfft, volts_rx_h_forfft, ant);
  // }
  
  // Apply anita-3 measured impulse response
  if (settings1->APPLYIMPULSERESPONSEDIGITIZER){

    anita1->GetNoiseWaveforms(); // get noise waveforms
    
    Tools::NormalTimeOrdering(anita1->NFOUR/2,volts_rx_e_forfft); 
    Tools::NormalTimeOrdering(anita1->NFOUR/2,volts_rx_h_forfft);
  
    for (int i=0;i<fNumPoints;i++){
      anita1->volts_rx_rfcm_lab_e[i] = volts_rx_e_forfft[i];
      anita1->volts_rx_rfcm_lab_h[i] = volts_rx_h_forfft[i];
    }
    
#ifdef ANITA_UTIL_EXISTS    
    applyImpulseResponseDigitizer(settings1, anita1, fNumPoints, ant, anita1->fTimes, anita1->volts_rx_rfcm_lab_e, 0);
    applyImpulseResponseDigitizer(settings1, anita1, fNumPoints, ant, anita1->fTimes, anita1->volts_rx_rfcm_lab_h, 1);
#endif

    if (settings1->SIGNAL_FLUCT && !settings1->NOISEFROMFLIGHTDIGITIZER){
      for (int i=0;i<anita1->NFOUR/2;i++) {
	anita1->volts_rx_rfcm_lab_e[i]+=anita1->timedomainnoise_lab_e[i]; // add noise
	anita1->volts_rx_rfcm_lab_h[i]+=anita1->timedomainnoise_lab_h[i];
      }
    }
    
    
  } else {
    
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
    if (anita1->PULSER) { // if we are using the pulser spectrum instead of simulating neutrinos
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
    anita1->MakeArraysforFFT(vhz_rx_rfcm_e,vhz_rx_rfcm_h,anita1->volts_rx_rfcm_e,anita1->volts_rx_rfcm_h, 90., true);
      
    // double volts_rx_rfcm_e_freq[anita1->HALFNFOUR];
    // for (int i=0;i<anita1->HALFNFOUR;i++) {
    //   volts_rx_rfcm_e_freq[i]=anita1->volts_rx_rfcm_e[i];
    // }
          
    // now the last two are in the frequency domain
    // convert to the time domain
    // still don't have any noise
      
    Tools::realft(anita1->volts_rx_rfcm_e,-1,anita1->HALFNFOUR);
    Tools::realft(anita1->volts_rx_rfcm_h,-1,anita1->HALFNFOUR);
      
    anita1->GetNoiseWaveforms(); // get noise waveforms
      
    // find the peak right here and it might be the numerator of the horizontal axis of matt's plot
    anita1->peak_e_rx_rfcm_signalonly=ChanTrigger::FindPeak(anita1->volts_rx_rfcm_e,anita1->HALFNFOUR); // with no noise
    anita1->peak_h_rx_rfcm_signalonly=ChanTrigger::FindPeak(anita1->volts_rx_rfcm_h,anita1->HALFNFOUR);
      
    if (settings1->SIGNAL_FLUCT) {
      for (int i=0;i<anita1->NFOUR/2;i++) {
        anita1->volts_rx_rfcm_e[i]+=anita1->timedomainnoise_rfcm_e[i]; // add noise.
        anita1->volts_rx_rfcm_h[i]+=anita1->timedomainnoise_rfcm_h[i];
      } 
    }
    
    
    anita1->peak_e_rx_rfcm=ChanTrigger::FindPeak(anita1->volts_rx_rfcm_e,anita1->HALFNFOUR); // with noise 
    anita1->peak_h_rx_rfcm=ChanTrigger::FindPeak(anita1->volts_rx_rfcm_h,anita1->HALFNFOUR); // with noise
      

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
    anita1->MakeArraysforFFT(vhz_rx_rfcm_lab_e,vhz_rx_rfcm_lab_h,anita1->volts_rx_rfcm_lab_e,anita1->volts_rx_rfcm_lab_h, 90., true);
      
    // now the last two are in the frequency domain
    // convert to the time domain
    // still don't have any noise
    Tools::realft(anita1->volts_rx_rfcm_lab_e,-1,anita1->HALFNFOUR); 
    Tools::realft(anita1->volts_rx_rfcm_lab_h,-1,anita1->HALFNFOUR);

    // put it in normal time ording -T to T
    // instead of 0 to T, -T to 0 
    Tools::NormalTimeOrdering(anita1->NFOUR/2,anita1->volts_rx_rfcm_lab_e); // EH, why only this has NormalTimeOrdering applied? Why not before?
    Tools::NormalTimeOrdering(anita1->NFOUR/2,anita1->volts_rx_rfcm_lab_h);

    if (settings1->SIGNAL_FLUCT) { 
      for (int i=0;i<anita1->NFOUR/2;i++) {
        anita1->volts_rx_rfcm_lab_e[i]+=anita1->timedomainnoise_lab_e[i]; // add noise
        anita1->volts_rx_rfcm_lab_h[i]+=anita1->timedomainnoise_lab_h[i];
      }
    }//end if signal_fluct

  } // END ELSE IMPULSE RESPONSE
}// end ImpulseResponse()


void ChanTrigger::TimeShiftAndSignalFluct(Settings *settings1, Anita *anita1, int ilayer, int ifold, double volts_rx_rfcm_lab_e_all[48][512], double volts_rx_rfcm_lab_h_all[48][512])
{   
  int ant = anita1->GetRxTriggerNumbering(ilayer, ifold);

  // now shift right to account for arrival times
  // for (int i=0;i<48;i++) std::cout << arrival_times[i] << std::endl;
  Tools::ShiftRight(anita1->volts_rx_rfcm_lab_e,anita1->NFOUR/2, int(anita1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));
  Tools::ShiftRight(anita1->volts_rx_rfcm_lab_h,anita1->NFOUR/2, int(anita1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));


  if (settings1->TRIGGEREFFSCAN && (settings1->TRIGGEREFFSCAPULSE==1)){ 
    injectImpulseAtSurf(anita1, anita1->volts_rx_rfcm_lab_e, anita1->volts_rx_rfcm_lab_h, ant);
  }


  for (int i=0;i<anita1->NFOUR/2;i++) {
    volts_rx_rfcm_lab_e_all[anita1->GetRx(ilayer, ifold)][i] = anita1->volts_rx_rfcm_lab_e[i];
    volts_rx_rfcm_lab_h_all[anita1->GetRx(ilayer, ifold)][i] = anita1->volts_rx_rfcm_lab_h[i];      
  }

  // now vmmhz_rx_rfcm_lab_e,h_forfft are the time domain waveforms after the antenna and lab attenuation
  // now find peak voltage
  // these get written to a tree
  anita1->peak_e_rx_rfcm_lab=ChanTrigger::FindPeak(anita1->volts_rx_rfcm_lab_e,anita1->HALFNFOUR);
  anita1->peak_h_rx_rfcm_lab=ChanTrigger::FindPeak(anita1->volts_rx_rfcm_lab_h,anita1->HALFNFOUR);  


  // END OF DIGITIZER PATH
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  
  
} //ChanTrigger constructor







//!
/*!
 *
 *
 *
 *
 *
 */
ChanTrigger::ChanTrigger() {
    
}


//!
/*!
 *
 *
 *
 *
 *
 */
double ChanTrigger::FindPeak(double *waveform,int n) {
    
  // find peak abs(voltage) of this waveform
    
  double peak=0.; // positive peak
    
  for (int i=0;i<n;i++) {
    if (fabs(waveform[i])>peak)
      peak=fabs(waveform[i]);
  }
  return peak;
    
}




//!	
/*!
 *
 *
 *
 */
void ChanTrigger::addToChannelSums(Settings *settings1,Anita *anita1,int ibw, int k) {
  /** Bug fix:  Changed "vm_banding_rfcm_*" to "v_banding_rfcm_*" in all four
      of term_e, term_h, energyterm_e, and energyterm_h. Stephen Hoover, 10 August 2009 **/
    
  //  cout << "v_banding, vm_banding is " << v_banding_rfcm_e[ibw][k] << " " << vm_banding_rfcm_e[ibw][k] << "\n";
    
  double term_e= v_banding_rfcm_e[ibw][k]*    // for the integral over e-field
    (settings1->BW/(double)Anita::NFREQ/1.E6); // width of a frequency bin
  // uses e- and h-components instead of lcp,rcp, for Anita-lite
  //cout << "bandwidth bin is " << (settings1->BW/(double)Anita::NFREQ/1.E6) << "\n";
    
  //  cout << "term_e is " << term_e  << "\n";
    
  double energyterm_e=v_banding_rfcm_e[ibw][k]*v_banding_rfcm_e[ibw][k]* // for the integral of energy of time T=N delta t = 1/delta f
    (settings1->BW/(double)Anita::NFREQ/1.E6)*(settings1->BW/(double)Anita::NFREQ/1.E6)/50.*
    anita1->INTEGRATIONTIME; // divide by width of a freq. bin
  // end divide by 50 ohms
    
    
  double term_h=v_banding_rfcm_h[ibw][k]* // for the integral over e-field
    (settings1->BW/(double)Anita::NFREQ/1.E6);
    
  double energyterm_h=v_banding_rfcm_h[ibw][k]*v_banding_rfcm_h[ibw][k]* // for the integral over energy
    (settings1->BW/(double)Anita::NFREQ/1.E6)*(settings1->BW/(double)Anita::NFREQ/1.E6)/50.*
    anita1->INTEGRATIONTIME; // multiply by frequency bin and divide by 50 ohms
    
    
  bwslice_volts_pole[ibw]+=term_e;
  bwslice_energy_pole[ibw]+=energyterm_e;
  bwslice_volts_polh[ibw]+=term_h;
  bwslice_energy_polh[ibw]+=energyterm_h;
    
}




//!	Takes ADC threshold as an input and outputs threshold in power
/*!
 *	"this is a place holder for now, from reading off of various plots, pointed out below"
 *	
 */
double ChanTrigger::ADCCountstoPowerThreshold(Anita *anita1, int ipol, int iant) {
  //double ChanTrigger::ADCCountstoPowerThreshold(int threshadc, int isurf,int ichan) {
  // first convert threshold in adc counts to the singles rate
  // do this using threshold scans taken before the flight
  // For Anita-3 using run 11927
  // these curves were read in the Balloon constructor
  // first check if the threshold in adc counts is in the allowable range
  int threshadc = anita1->thresholds[ipol][iant];
  if (unwarned && (threshadc<anita1->minadcthresh[ipol][iant] || threshadc>anita1->maxadcthresh[ipol][iant]))
    cout << "Warning! ADC threshold is outside range of measured threshold scans.";
  if (threshadc<anita1->minadcthresh[ipol][iant]) {
    if (unwarned) {
      cout << "It is below the minimum so set it to the minimum.  Will not be warned again.\n";
      unwarned=0;
    }
    threshadc=anita1->minadcthresh[ipol][iant];
  }
  if (threshadc>anita1->maxadcthresh[ipol][iant]) {
    if (unwarned) {
      cout << "It is higher than the maximum so set it to the maximum.  Will not be warned again.\n";
      unwarned=0;
    }
    threshadc=anita1->maxadcthresh[ipol][iant];
  }
    
    
  // Now find singles rate for this threshold
  // first sort thresholds
  // int index=TMath::BinarySearch(NPOINTS,threshold[isurf][ichan],threshadc);
  int index=TMath::BinarySearch(anita1->npointThresh, anita1->threshScanThresh[ipol][iant], threshadc);

  //cout << "rate is " << rate[isurf][ichan][index] << "\n";
  //  thisrate=(double)rate[isurf][ichan][index]; // these scalers are in kHz
  thisrate=(double)anita1->threshScanScaler[ipol][iant][index]; // these scalers are in kHz
    
  // now find threshold for this scaler.  Here, scalers have to be in MHz.
  thisrate=thisrate/1.E3; // put it in MHz
  //cout << "thisrate is " << thisrate << "\n";
  // figure out what band we're talking about
  // int iband=Anita::SurfChanneltoBand(ichan);
  // FOR THE MOMENT JUST USING THE FULL BAND
  int iband =3;
  thispowerthresh=rateToThreshold(thisrate,iband);
  //cout << "thisrate, iband, thispowerthresh are " << thisrate << " " << iband << " " << thispowerthresh << "\n";
  return thispowerthresh;
    
  //  double powerthresh=-0.3*pow(10.,logsingles)/1.E6+5.091; // this has the right slope and intercept
  //cout << "threshadc, powerthresh are " << threshadc << " " << powerthresh << "\n";
  //return powerthresh;
  // return 0;
}



//!	Returns the thisrate variable value
/*!
 *	\todo	Refactor the class so that this function is either deprecated or "thisrate"
 *			is implemented differently.
 */
double ChanTrigger::getRate() {
    
  return thisrate;
}




//!	Calculates the trigger threshold for an antenna via a fit from the "singles" rate and band identifier
/*!
 *	This code is stolen from Stephen's AnitaHardwareTrigger in his Aesop software
 *	"rate" is the desired singles rate, in MHz.
 *	"band" ranges from 0 to 3
 *	The output is the trigger threshold, in units of Diode Output / <Diode Output>
 *	The output is good for all antennas, and the curves are good fits in the
 *	region ~1 MHz to 25 MHz.  The curves look like they should continue beyond
 *	that range as well.
 *
 *	\todo	The hard-coded values should be explained and the band dependence should be more explicit,
 *			because the actual bands can be modified elsewhere in the code and this function won't be
 *			changed unless the user specifically knows to.
 */
double ChanTrigger::rateToThreshold(double rate, int band)
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





//!	Returns whether the indicated antenna and band are "masked"
/*!
 *	
 *
 *
 */
int ChanTrigger::IsItUnmasked(unsigned short surfTrigBandMask[9][2],int ibw,int ilayer, int ifold, int ipol) {
    
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




//!	The L1 trigger of the Anita trigger scheme
/*!
 *
 *	
 *	\todo	The number of parameters this function accepts is rather high and should be reduced.
 *			In many cases an argument is a member object of another argument which was already
 *			passed to this function, which leads to intent confusion. The member objects like
 *			"flag_e" and "flag_h" are manipulated here but also in other places, and it makes
 *			tracking the changes through code difficult. Also, when dealing with std::vector
 *			containers, using "clear" and "push_back" is typically much slower than accessing
 *			each element by index, since the number of allocations/deallocations is guaranteed
 *			to be zero.
 */


void ChanTrigger::L1Trigger(Anita *anita1,double timedomain_output_1[5][Anita::NFOUR],double timedomain_output_2[5][Anita::NFOUR],double powerthreshold[2][5],
			    int *channels_passing_e_forglob,int *channels_passing_h_forglob,int &npass) {
   
  int maxsample=TMath::MaxElement(5,anita1->imaxbin);
  int minsample=TMath::MinElement(5,anita1->iminbin);
    
  for (int j=0;j<5;j++) {
    flag_e[j].clear();
    flag_h[j].clear();		
      
    for (int i=minsample;i<maxsample;i++) {
	
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
  



//!	
/*!
 *
 *
 *
 */
double ChanTrigger::GetNoise(Settings *settings1,double altitude_bn,double geoid,double theta_zenith,double bw,double temp) {
    
  int NSTEPS=1000;
  //  double VSKY=150;
  double VSKY=15;
  double VICE=240;
  double VSYSTEM=200;
    
  double sum=0;
  double theta_pos=0;
  double theta_signed=0;
  double theta_cant=theta_zenith-PI/2;
  double theta_horizon=acos(geoid/(geoid+altitude_bn)); // angle between horizontal on the payload and line of sight to horizon
  //  double theta_horizon=sqrt(2*altitude_bn/geoid);  // angle between horizontal on the payload and line of sight to horizon
  //double theta_horizon=-1;
  double theta_0=50*RADDEG;
  double integral=0.;
  double integral_firsthalf=0;
  double integral_secondhalf=0;
  double vnoise=0;
    
  if (settings1->WHICH != 0) {
		
    for (int i=0;i<NSTEPS;i++) {
			
      // step in theta
      theta_pos=(double)fabs(-PI+(double)i/(double)NSTEPS*2*PI); // this is always a positive number
      theta_signed=(double)(-PI+(double)i/(double)NSTEPS*2*PI); // this is allowed to be signed
			
      if (theta_signed<theta_horizon-theta_cant) {
	vnoise=VSKY+VSYSTEM;
	integral_firsthalf+=exp(-2*ALOG2*pow(theta_pos/theta_0,2))*2*PI/NSTEPS;
      } //if
      else {
	vnoise=VICE+VSYSTEM;
	integral_secondhalf+=exp(-2*ALOG2*pow(theta_pos/theta_0,2))*2*PI/NSTEPS;
      } //else
			
      sum+=vnoise*exp(-2*ALOG2*pow(theta_pos/theta_0,2))*2*PI/NSTEPS;
      integral+=exp(-2*ALOG2*pow(theta_pos/theta_0,2))*2*PI/NSTEPS;
    } //if
		
    sum=sum/integral;
		
    //    cout << "sum, KBOLTZ, bw, sqrt are " << sum << " " << KBOLTZ << " " << bw << " " << sqrt(sum*50.*KBOLTZ*bw) << "\n";
    return sqrt(sum*50.*KBOLTZ*bw);
  } //if (settings1->WHICH != 0)
  else if (settings1->WHICH == 0)
    return sqrt(temp*50.*KBOLTZ*bw);
    
  return 0;
}//GetNoise




//!	Sets the threshold values based on which payload and where the antenna is located physically
/*!
 *	The nadir antennas had a separate threshold from the other antennas due to the way that they
 *	were "OR"'d between their two neighbors.
 *
 *	\todo	Deprecate in favor of the more robust boost::multi_array or the more specialized
 *			PayloadArray class. Both have multi-index access to the same items.
 */
void ChanTrigger::GetThresholds(Settings *settings1,Anita *anita1,int ilayer,double thresholds[2][5]) {
    
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










#ifdef ANITA_UTIL_EXISTS    
void ChanTrigger::applyImpulseResponseDigitizer(Settings *settings1, Anita *anita1, int nPoints, int ant, double *x, double y[512], bool pol){

  TGraph *graph1 = new TGraph(nPoints, x, y);
  // Upsample waveform to same deltaT of the signal chain impulse response
  TGraph *graphUp = FFTtools::getInterpolatedGraph(graph1, anita1->deltaT);

  int ipol=0;
  int iring=2;
  if (pol) ipol = 1;
  if (ant<16) iring=0;
  else if (ant<32) iring=1;
  
  //Calculate convolution
  TGraph *surfSignal = FFTtools::getConvolution(graphUp, anita1->fSignalChainResponseDigitizer[ipol][iring]);

  //Downsample again
  TGraph *surfSignalDown = FFTtools::getInterpolatedGraph(surfSignal, 1/2.6);
  
  Double_t *newy = surfSignalDown->GetY();
  if (settings1->ZEROSIGNAL){
    for (int i=0;i<nPoints;i++) newy[i]=0;
  }

  // Translate graph of 35ns
  int indexTranslation = 35*2.6;

  // add thermal noise for anita-3 flight
  if (settings1->SIGNAL_FLUCT && settings1->NOISEFROMFLIGHTDIGITIZER) { 
    double *justNoise = getNoiseFromFlight(anita1, ipol, ant);
    for (int i=0;i<nPoints;i++){
      y[i]=newy[i+indexTranslation] + justNoise[i]*anita1->THERMALNOISE_FACTOR;
      // std::cout << justNoise[i] << std::endl;
    }
  } else {
    for (int i=0;i<nPoints;i++)  y[i]=newy[i+indexTranslation];
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
  delete surfSignal;
  delete graphUp;
  delete graph1;
}


void ChanTrigger::applyImpulseResponseTrigger(Settings *settings1, Anita *anita1, int nPoints, int ant, double *x, double y[512], double *vhz, bool pol){

  TGraph *graph1 = new TGraph(nPoints, x, y);
  // Upsample waveform to same deltaT of the signal chain impulse response
  TGraph *graphUp = FFTtools::getInterpolatedGraph(graph1, anita1->deltaT);

  int ipol=0;
  int iring=2;
  if (pol) ipol = 1;
  if (ant<16) iring=0;
  else if (ant<32) iring=1;
  
  //Calculate convolution
  TGraph *surfSignal = FFTtools::getConvolution(graphUp, anita1->fSignalChainResponseTrigger[ipol][iring]);
  
  //Downsample again
  TGraph *surfSignalDown = FFTtools::getInterpolatedGraph(surfSignal, 1/2.6);
  
  Double_t *newy = surfSignalDown->GetY();
  if (settings1->ZEROSIGNAL){
    for (int i=0;i<nPoints;i++) newy[i]=0;
  }

  // Translate graph of 35ns
  int indexTranslation = 35*2.6;
    
  // add thermal noise for anita-3 flight
  if (settings1->SIGNAL_FLUCT && settings1->NOISEFROMFLIGHTTRIGGER) { 
    double *justNoise = getNoiseFromFlight(anita1, ipol, ant);
    for (int i=0;i<nPoints;i++){
      y[i]=newy[i] + justNoise[i+indexTranslation]*anita1->THERMALNOISE_FACTOR;
      // std::cout << justNoise[i] << std::endl;
    }
  } else {
    for (int i=0;i<nPoints;i++)  y[i]=newy[i+indexTranslation];
  }

  FFTWComplex *theFFT = FFTtools::doFFT(nPoints, y);
    
  for (int i=0;i<anita1->NFREQ;i++){
    vhz[i]   = theFFT[i].getAbs();    
    // if (ant ==8 && pol==0) cout << vhz[i] << endl;
  }
  // double *invVals = FFTtools::doInvFFT(nPoints,theFFT);
  // for (int i=0;i<nPoints;i++){
  //   cout << y[i] << " " << invVals[i] << endl;
  // }
  // if (ant ==8 && pol==0){
  //   TCanvas *c = new TCanvas("c");
  //   graph1->Draw("Al");
  //   c->Print("TriggerPath_graph1.png");
  //   graphUp->Draw("Al");
  //   c->Print("TriggerPath_graphUp.png");
  //   surfSignal->Draw("Al");
  //   c->Print("TriggerPath_surfSignal.png");
  //   surfSignalDown->Draw("Al");
  //   c->Print("TriggerPath_surfSignalDown.png");
  // }
  // Cleaning up
  delete []theFFT;
  delete surfSignalDown;
  delete surfSignal;
  delete graphUp;
  delete graph1;
}

void ChanTrigger::applyImpulseResponseTrigger(Settings *settings1, Anita *anita1, int ant, double y[512], bool pol){

  int nPoints=anita1->HALFNFOUR;
  double *x = anita1->fTimes;
  
  TGraph *graph1 = new TGraph(nPoints, x, y);
  // Upsample waveform to same deltaT of the signal chain impulse response
  TGraph *graphUp = FFTtools::getInterpolatedGraph(graph1, anita1->deltaT);

  int ipol=0;
  int iring=2;
  if (pol) ipol = 1;
  if (ant<16) iring=0;
  else if (ant<32) iring=1;
  
  //Calculate convolution
  TGraph *surfSignal = FFTtools::getConvolution(graphUp, anita1->fSignalChainResponseTrigger[ipol][iring]);
  
  //Downsample again
  TGraph *surfSignalDown = FFTtools::getInterpolatedGraph(surfSignal, 1/2.6);
  
  Double_t *newy = surfSignalDown->GetY();
  if (settings1->ZEROSIGNAL){
    for (int i=0;i<nPoints;i++) newy[i]=0;
  } 

  // Translate graph of 35ns
  int indexTranslation = 35*2.6;

  // add thermal noise for anita-3 flight
  if (settings1->SIGNAL_FLUCT && settings1->NOISEFROMFLIGHTTRIGGER) { 
    double *justNoise = getNoiseFromFlight(anita1, ipol, ant);
    for (int i=0;i<nPoints;i++){
      y[i]=newy[i] + justNoise[i+indexTranslation]*anita1->THERMALNOISE_FACTOR;
      // std::cout << justNoise[i] << std::endl;
    }
  } else {
    for (int i=0;i<nPoints;i++)  y[i]=newy[i+indexTranslation];
  }

  delete surfSignalDown;
  delete surfSignal;
  delete graphUp;
  delete graph1;
}

double *ChanTrigger::getNoiseFromFlight(Anita* anita1, int pol, int ant){

  Int_t numFreqs = anita1->numFreqs;
  FFTWComplex *phasors = new FFTWComplex[numFreqs];
  double *freqs = anita1->freqs;
  phasors[0].setMagPhase(0,0);
  Double_t sigma, realPart, imPart;

  for(int i=1;i<numFreqs;i++) {
    sigma      = anita1->RayleighFits[pol][ant]->Eval(freqs[i])*4/TMath::Sqrt(numFreqs);
    realPart   = anita1->fRand->Gaus(0,sigma);
    imPart     = anita1->fRand->Gaus(0,sigma);
    phasors[i] = FFTWComplex(realPart, imPart);
  }
    
  RFSignal *rfNoise = new RFSignal(numFreqs,freqs,phasors,1);
    
  Double_t *justNoise=rfNoise->GetY();

  // Cleaning up
  delete[] phasors;
  delete rfNoise;
    
  return justNoise;
  
}

void ChanTrigger::injectImpulseAfterAntenna(Anita *anita1, double volts_triggerPath_e[Anita::HALFNFOUR], double volts_triggerPath_h[Anita::HALFNFOUR], int ant){
  // phiIndex is 0 for central antenna in trigger efficiency scan
  int phiIndex = anita1->trigEffScanPhi - (ant%16);
  if (phiIndex>8) phiIndex=phiIndex-16;
  int fNumPoints=anita1->HALFNFOUR;
  // only 2 phi sectors adjecent to the central one are considered in efficiency scan
  if(TMath::Abs(phiIndex)<=2){
    double att = anita1->trigEffScanAtt[phiIndex+2]-anita1->trigEffScanAtt[2];
    double norm = (anita1->trigEffScanAtt[phiIndex+2]==999)*0 + (anita1->trigEffScanAtt[phiIndex+2]!=999)*1;
    for (int i=0;i<fNumPoints;i++){
      volts_triggerPath_e[i]=norm*anita1->trigEffScanPulseAtAmpa[i]*TMath::Power(10, att/20);
      volts_triggerPath_h[i]=0;
    }
  }else{
    for (int i=0;i<fNumPoints;i++){
      volts_triggerPath_e[i]=0;
      volts_triggerPath_h[i]=0;
    }
  }
}

void ChanTrigger::injectImpulseAmplitudeAfterAntenna(Anita *anita1, double vhz_triggerPath_e[Anita::HALFNFOUR], double vhz_triggerPath_h[Anita::HALFNFOUR], int ant){
  // phiIndex is 0 for central antenna in trigger efficiency scan
  int phiIndex = anita1->trigEffScanPhi - (ant%16);
  if (phiIndex>8) phiIndex=phiIndex-16;
  int fNumPoints=anita1->NFREQ;
  // only 2 phi sectors adjecent to the central one are considered in efficiency scan
  if(TMath::Abs(phiIndex)<=2){
    double att = anita1->trigEffScanAtt[phiIndex+2]-anita1->trigEffScanAtt[2];
    double norm = (anita1->trigEffScanAtt[phiIndex+2]==999)*0 + (anita1->trigEffScanAtt[phiIndex+2]!=999)*1;
    for (int i=0;i<fNumPoints;i++){
      vhz_triggerPath_e[i]=norm*anita1->trigEffScanAmplitudeAtAmpa[i]*TMath::Power(10, att/20);
      vhz_triggerPath_h[i]=0;
    }
  }else{
    for (int i=0;i<fNumPoints;i++){
      vhz_triggerPath_e[i]=0;
      vhz_triggerPath_h[i]=0;
    }
  }
}

void ChanTrigger::injectImpulseAtSurf(Anita *anita1, double volts_triggerPath_e[Anita::HALFNFOUR], double volts_triggerPath_h[Anita::HALFNFOUR], int ant){
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
