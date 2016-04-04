#include <vector>
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
#include "TRandom3.h"

#include "rx.hpp"
#include "Constants.h"
#include "anita.hh"
#include "balloon.hh"
#include "TRandom3.h"
#include "trigger.hh"
#include <cmath>
#include "Tools.h"
#include "Settings.h"

using std::cout;

GlobalTrigger::GlobalTrigger(Settings *settings1,Anita *anita1,UShort_t phiTrigMask_bn){
    
    Tools::Zero(triggerbits,Anita::NTRIG);
    
    phiTrigMask=phiTrigMask_bn; // set the phi mask to the input value which comes from the balloon class
    phiTrigMaskH=0; // set other variables to 0 
    l1TrigMask  =0; 
    l1TrigMaskH =0;
    
    for (int i=0;i<Anita::NLAYERS_MAX;i++) {
      for (int j=0;j<Anita::NPHI_MAX;j++) {
	for (int k=0;k<2;k++) {
	  for (int p=0;p<anita1->NBANDS+1;p++) {
	    //	for (int p=0;p<5;p++) {
	    channels_passing[i][j][k][p]=0;
	    // make vchannels_passing the proper length.
	    vchannels_passing[i][j][k].push_back(0);
	  }
	}
      }
    }
    
    
    
    for (int k=0;k<2;k++) {		
      for (int i=0;i<Anita::NLAYERS_MAX;i++) {
	for (int j=0;j<Anita::NPHI_MAX;j++) {
	  volts[k][i][j]=0.;
	  volts_em[k][i][j]=0.;
	  volts_original[k][i][j]=0.; //added djg
	}
      }
    }
    
    
    
    //Zeroing
    for (int i=0;i<settings1->NANTENNAS;i++) {
      nchannels_perrx_triggered[i] = 0;
      for (int j=0;j<8;j++) {
	nchannels_perband_triggered[i][j]=0;
      }
    } //Zero the trigger array
    
    
    
    
}

GlobalTrigger::GlobalTrigger(Settings *settings1,Anita *anita1,UShort_t phiTrigMask_bn,UShort_t phiTrigMaskH_bn,UShort_t l1TrigMask_bn,UShort_t l1TrigMaskH_bn){
    
    Tools::Zero(triggerbits,Anita::NTRIG);
    
    phiTrigMask=phiTrigMask_bn; // set the phi mask to the input value which comes from the balloon class
    phiTrigMaskH=phiTrigMaskH_bn; // set the phi mask to the input value which comes from the balloon class
    l1TrigMask=l1TrigMask_bn; // set the phi mask to the input value which comes from the balloon class
    l1TrigMaskH=l1TrigMaskH_bn; // set the phi mask to the input value which comes from the balloon class
    
    for (int i=0;i<Anita::NLAYERS_MAX;i++) {
      for (int j=0;j<Anita::NPHI_MAX;j++) {
	for (int k=0;k<2;k++) {
	  for (int p=0;p<anita1->NBANDS+1;p++) {
	    //	for (int p=0;p<5;p++) {
	    channels_passing[i][j][k][p]=0;
	    // make vchannels_passing the proper length.
	    vchannels_passing[i][j][k].push_back(0);
	  }
	}
      }
    }
    
    
    
    for (int k=0;k<2;k++) {		
      for (int i=0;i<Anita::NLAYERS_MAX;i++) {
	for (int j=0;j<Anita::NPHI_MAX;j++) {
	  volts[k][i][j]=0.;
	  volts_em[k][i][j]=0.;
	  volts_original[k][i][j]=0.; //added djg
	}
      }
    }
    
    
    
    //Zeroing
    for (int i=0;i<settings1->NANTENNAS;i++) {
      nchannels_perrx_triggered[i] = 0;
      for (int j=0;j<8;j++) {
	nchannels_perband_triggered[i][j]=0;
      }
    } //Zero the trigger array
    
    
    
    
}


void AntTrigger::ConvertEHtoLREfield(double e_component,double h_component,double& lcp_component,double& rcp_component) {
    
    lcp_component=sqrt((e_component*e_component+h_component*h_component)/2);
    rcp_component=lcp_component;
    
} //ConvertEHtoLREfield
void AntTrigger::ConvertEHtoLREnergy(double e_component,double h_component,double& lcp_component,double& rcp_component) {
    
    lcp_component=(e_component+h_component)/2;
    rcp_component=lcp_component;
    
} //ConvertEHtoLREnergy
void AntTrigger::ConvertHVtoLRTimedomain(const int nfour,double *vvolts,
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
void AntTrigger::WhichBandsPass(Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double dangle, double emfrac, double hadfrac){
    
    double thresholds[5];
    GetThresholds(settings1,anita1,ilayer,thresholds); // get the right thresholds for this layer
    
    
    double volts_thischannel;
    double energy_thischannel;
    double voltagethresh_thischannel,energythresh_thischannel;
    
    double v_banding_rfcm_e_forfft[5][anita1->HALFNFOUR]; // starts out as V/s vs. freq after banding, rfcm, after fft it is V vs. t
    double v_banding_rfcm_h_forfft[5][anita1->HALFNFOUR];
    double vm_banding_rfcm_1_forfft[5][anita1->HALFNFOUR];
    double vm_banding_rfcm_2_forfft[5][anita1->HALFNFOUR];
    
    double v_banding_rfcm_e_forfft_temp[5][anita1->HALFNFOUR];
    
    double v_banding_rfcm_h_forfft_temp[5][anita1->HALFNFOUR];
    
    
    double psignal_e[5][anita1->NFOUR];
    double psignal_h[5][anita1->NFOUR];
    
    double mindiodeconvl_e[5];
    double mindiodeconvl_h[5];
    double onediodeconvl_e[5];
    double onediodeconvl_h[5];
    
    double timedomain_output_1[5][Anita::NFOUR];
    double timedomain_output_2[5][Anita::NFOUR];
    
    globaltrig1->volts[0][ilayer][ifold]=0.;
    globaltrig1->volts[1][ilayer][ifold]=0.;
    
    //  TRandom3 Rand3;
    
    if (settings1->TRIGGERSCHEME <= 1){
		
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
	if (settings1->CHMASKING && !AntTrigger::IsItUnmasked(bn1->surfTrigBandMask,ibw,ilayer,ifold,0)) {
	  
	  globaltrig1->channels_passing[ilayer][ifold][0][ibw]=0;// channel does not pass
	  globaltrig1->vchannels_passing[ilayer][ifold][0][ibw]=0;// channel does not pass
	}
	// note the last element of the array is 0 because this is lcp or e pol
	else {
	  
	  
	  if (settings1->LCPRCP)  {// if we're considering lcp, rcp
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
	  
	  energythresh_thischannel=thresholds[ibw];
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
	  if (settings1->CHMASKING && !AntTrigger::IsItUnmasked(bn1->surfTrigBandMask,ibw,ilayer,ifold,1)) {
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
	    
	    energythresh_thischannel=thresholds[ibw];
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
    } // end if settings1->TRIGGERSCHEME!=2
    // settings1->TRIGGERSCHEME==2
    else  if (settings1->TRIGGERSCHEME >= 2){ // if we use the diode to perform an integral
      // this is the number of bins to the left of center where the diode function starts to be completely overlapping with the waveform in the convolution.
      int ibinshift=(anita1->NFOUR/4-(int)(anita1->maxt_diode/anita1->TIMESTEP));
      
      
      
      double integrate_energy_freq[5]={0.,0.,0.,0.,0.};
      for (int iband=0;iband<5;iband++) {
	//cout << "arrival time is " << globaltrig1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP << "\n";
	//anita1->iminbin[j]=anita1->NFOUR/4-ibinshift+anita1->idelaybeforepeak[j]+globaltrig1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP; // we start to look for single channel triggers firing
	//anita1->iminbin[j]=anita1->NFOUR/4-ibinshift+anita1->idelaybeforepeak[j]; // we start to look for single channel triggers firing
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
	
	anita1->MakeArraysforFFT(v_banding_rfcm_e[iband],v_banding_rfcm_h[iband],v_banding_rfcm_e_forfft[iband],v_banding_rfcm_h_forfft[iband]);
	
	// for some reason I'm averaging over 10 neighboring bins
	// to get rid of the zero bins
	for (int i=0;i<anita1->NFOUR/4;i++) {
	  
	  v_banding_rfcm_e_forfft_temp[iband][2*i]=0.;
	  v_banding_rfcm_e_forfft_temp[iband][2*i+1]=0.;
	  v_banding_rfcm_h_forfft_temp[iband][2*i]=0.;
	  v_banding_rfcm_h_forfft_temp[iband][2*i+1]=0.;
	  
	  for (int k=i;k<i+10;k++) {
	    if (k<anita1->NFOUR/4) {
	      v_banding_rfcm_e_forfft_temp[iband][2*i]+=v_banding_rfcm_e_forfft[iband][2*k];
	      v_banding_rfcm_e_forfft_temp[iband][2*i+1]+=v_banding_rfcm_e_forfft[iband][2*k+1];
	      v_banding_rfcm_h_forfft_temp[iband][2*i]+=v_banding_rfcm_h_forfft[iband][2*k];
	      v_banding_rfcm_h_forfft_temp[iband][2*i+1]+=v_banding_rfcm_h_forfft[iband][2*k+1];
	    }
	  }
	  
	  v_banding_rfcm_e_forfft[iband][2*i]=v_banding_rfcm_e_forfft_temp[iband][2*i]/10.;
	  v_banding_rfcm_e_forfft[iband][2*i+1]=v_banding_rfcm_e_forfft_temp[iband][2*i+1]/10.;
	  v_banding_rfcm_h_forfft[iband][2*i]=v_banding_rfcm_h_forfft_temp[iband][2*i]/10.;
	  v_banding_rfcm_h_forfft[iband][2*i+1]=v_banding_rfcm_h_forfft_temp[iband][2*i+1]/10.;
	  
	  v_banding_rfcm_e_forfft_temp[iband][2*i]=v_banding_rfcm_e_forfft[iband][2*i];
	  v_banding_rfcm_e_forfft_temp[iband][2*i+1]=v_banding_rfcm_e_forfft[iband][2*i+1];
	  v_banding_rfcm_h_forfft_temp[iband][2*i]=v_banding_rfcm_h_forfft[iband][2*i];
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
	
	if (settings1->ZEROSIGNAL) {
	 
	  Tools::Zero(v_banding_rfcm_e_forfft[iband],anita1->NFOUR/2);
	  Tools::Zero(v_banding_rfcm_h_forfft[iband],anita1->NFOUR/2);
	  
	}

// 	  for (int i=0;i<anita1->NFOUR/4;i++) {
// 	    cout << "v_banding is " << v_banding_rfcm_e_forfft[iband][i] << "\n";
// 	  }
	
	
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
      
      if (settings1->SIGNAL_FLUCT) {
	
	for (int iband=0;iband<5;iband++) {
	  if (settings1->TRIGGERSCHEME == 2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5){
	    for (int k=0;k<anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP);k++) {
	      anita1->maxbin_fortotal[iband]=anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP);
	      
	      // The below line seems to shorten the length of the waveform to less than HALFNFOUR
	      int knoisebin=anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP)-k;
	      
	      // this is just the straight sum of the two
	      anita1->total_vpol_inanita[iband][k]=anita1->timedomainnoise_rfcm_banding_e[iband][k]+anita1->signal_vpol_inanita[iband][k];
	      
	      integrateenergy[iband]+=anita1->timedomainnoise_rfcm_banding_e[iband][k]*anita1->timedomainnoise_rfcm_banding_e[iband][k]*anita1->TIMESTEP;
	      // this reverses the noise is time, and starts with bin anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP)
	      v_banding_rfcm_e_forfft[iband][k]=v_banding_rfcm_e_forfft[iband][k]+anita1->timedomainnoise_rfcm_banding_e[iband][knoisebin];
	    }
	    for (int k=anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP);k<anita1->NFOUR/2;k++) {
	      anita1->total_vpol_inanita[iband][k]=0.;
	      v_banding_rfcm_e_forfft[iband][k]=0.;
	    }
	  }
	  else {
	    for (unsigned k = 0; k < anita1->HALFNFOUR; ++k){
	      // this is just a straight sum
	      anita1->total_vpol_inanita[iband][k]=anita1->timedomainnoise_rfcm_banding_e[iband][k]+anita1->signal_vpol_inanita[iband][k];
	      integrateenergy[iband]+=anita1->timedomainnoise_rfcm_banding_e[iband][k]*anita1->timedomainnoise_rfcm_banding_e[iband][k]*anita1->TIMESTEP;
	      // this one is the one actually used by the diode
	      v_banding_rfcm_e_forfft[iband][k] += anita1->timedomainnoise_rfcm_banding_e[iband][k];
	    }
	  }
	} // end loop over bands
	
	for (int iband=0;iband<5;iband++) { // loop over bands
	  for (int k=0;k<anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP);k++) {
	    int knoisebin=anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP)-k;
	    //	  cout << "before adding noise, p2p is " << FindPeak(v_banding_rfcm_h_forfft[j],anita1->NFOUR/2) << "\n";
	    v_banding_rfcm_h_forfft[iband][k]=v_banding_rfcm_h_forfft[iband][k]+anita1->timedomainnoise_rfcm_banding_h[iband][knoisebin];
	  }
	  for (int k=anita1->NFOUR/2-(int)(anita1->maxt_diode/anita1->TIMESTEP);k<anita1->NFOUR/2;k++) {
	    v_banding_rfcm_h_forfft[iband][k]=0.;
	  }
	  
	} // loop over 5 bands
      } // if we require signal fluctuations
      
      for (int iband=0;iband<5;iband++) {
	if (settings1->LCPRCP) {
	  ConvertHVtoLRTimedomain(anita1->NFOUR, v_banding_rfcm_e_forfft[iband], v_banding_rfcm_h_forfft[iband], vm_banding_rfcm_1_forfft[iband], vm_banding_rfcm_2_forfft[iband]);
	} else {
	  for (int i=0;i<anita1->NFOUR/2;i++) {
	    
	    vm_banding_rfcm_1_forfft[iband][i] = v_banding_rfcm_e_forfft[iband][i];
	    vm_banding_rfcm_2_forfft[iband][i] = v_banding_rfcm_h_forfft[iband][i];
	    
	  }
	}
	for (int i=0;i<anita1->NFOUR/2;i++) {
	  anita1->total_diodeinput_1_inanita[iband][i] = vm_banding_rfcm_1_forfft[iband][i];
	  anita1->total_diodeinput_2_inanita[iband][i] = vm_banding_rfcm_2_forfft[iband][i];
	  
	}
	//volts_fullband_trigger_path_e.push_back(;
      } // end loop over bands
      
      
      
      for (int i=0;i<Anita::NFOUR/2;i++) {
	anita1->total_diodeinput_1_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)][i]=anita1->total_diodeinput_1_inanita[4][i];
	anita1->total_diodeinput_2_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)][i]=anita1->total_diodeinput_2_inanita[4][i];
      }


      
      for (int j=0;j<5;j++) {
	// myconvl
	// this performs the convolution with the diode response

// 	for (int i=0;i<NFOUR/2;i++) {
// 	  cout << "vm_banding_rfmc_1_forfft is " << vm_banding_rfcm_1_forfft[j][i] << "\n";
// 	}


	anita1->myconvlv(vm_banding_rfcm_1_forfft[j],anita1->NFOUR,anita1->fdiode_real[j],mindiodeconvl_e[j],onediodeconvl_e[j],psignal_e[j],timedomain_output_1[j]);
	// loop from the ibinshift left + some delay + 10 ns
	


	// now shift right to account for arrival times
	Tools::ShiftRight(timedomain_output_1[j],anita1->NFOUR,(int)(anita1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));
	
	
	if (settings1->TRIGGERSCHEME == 2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5){
	  Tools::ShiftLeft(timedomain_output_1[j],anita1->NFOUR,ibinshift);
	  
	}

	
	for (int i=0;i<Anita::NFOUR/2;i++) {
	  anita1->timedomain_output_1_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)][i]=timedomain_output_1[4][i];
	}
	
	anita1->channels_passing_e[j]=0; // does not pass
	
	for (int ibin = anita1->iminbin[j]; ibin < anita1->imaxbin[j]; ibin++) {
	  if (timedomain_output_1[j][ibin] < thresholds[j] * anita1->bwslice_rmsdiode[j]) {
	    if (anita1->pol_allowed[0] && anita1->bwslice_allowed[j]) { // is this polarization and bw slice allowed to pass
	      //	      std::cout << j << " " << timedomain_output_1[j][ibin]  << " " <<  thresholds[j] << " " <<  anita1->bwslice_rmsdiode[j] << std::endl;
	      anita1->channels_passing_e[j] = 1;// channel passes
	    }
	  }
	} // end loop over bins in the window
	
	if (anita1->channels_passing_e[j]) {
	  globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]++; //Records number of first level triggers on each antenna for a single neutrino
	  
	}
      } // end loop over bands
      
      for (int j=0;j<5;j++) {
	// Find p2p value before adding noise
	anita1->peak_v_banding_rfcm_h[j]=FindPeak(v_banding_rfcm_h_forfft[j],anita1->NFOUR/2);
      }
      
      for (int j=0;j<5;j++) {
	anita1->myconvlv(vm_banding_rfcm_2_forfft[j],anita1->NFOUR/2,anita1->fdiode_real[j],mindiodeconvl_h[j],onediodeconvl_h[j],psignal_h[j],timedomain_output_2[j]);


	
	// now shift right to account for arrival times


 	Tools::ShiftRight(timedomain_output_2[j],anita1->NFOUR/2,(int)(anita1->arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));
	
	
 	if (settings1->TRIGGERSCHEME == 2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5){
 	  Tools::ShiftLeft(timedomain_output_2[j],anita1->NFOUR/2,ibinshift);
 	}




 	for (int i=0;i<Anita::NFOUR/2;i++) {
 	  anita1->timedomain_output_2_allantennas[anita1->GetRxTriggerNumbering(ilayer,ifold)][i]=timedomain_output_2[4][i];
 	}


	anita1->channels_passing_h[j]=0;
	for (int ibin=anita1->iminbin[j];ibin<anita1->imaxbin[j];ibin++) {
	  
	  if (timedomain_output_2[j][ibin]<thresholds[j]*anita1->bwslice_rmsdiode[j]) {
	    if (anita1->pol_allowed[1] && anita1->bwslice_allowed[j]) { // if this pol and band are allowed to pass
	      anita1->channels_passing_h[j]=1;
	    }
	  } // if it's over threshold for this time step
	  
	}
	
	
	if (anita1->channels_passing_h[j])
	  globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]++; //Records number of first level triggers on each antenna for a single neutrino
      } // end loop over bands
      
      // fill channels_passing
      //  } // end if the signal is big enough the be considered
      // now we've made the diode outputs
      // now step in time and find the most number of bands that pass in the
      // right time window
      // then fills channels_passing
      int npass;

      L1Trigger(anita1,timedomain_output_1,timedomain_output_2,thresholds,globaltrig1->channels_passing[ilayer][ifold][0],globaltrig1->channels_passing[ilayer][ifold][1],npass);

      
      // if it's the closest antenna,
      // save flag_e,h in anita class for writing to tsignals tree
      int startbin=TMath::MinElement(5,anita1->iminbin);
      
      if (ilayer==anita1->GetLayer(anita1->rx_minarrivaltime) && ifold==anita1->GetIfold(anita1->rx_minarrivaltime)) {
	for (int j=0;j<5;j++) {

	  //	  cout << "zeroeing here 1.\n";
	  anita1->ston[j]=0.;
	  for (int i=anita1->iminbin[j];i<anita1->imaxbin[j];i++) {
// 	    if (j==0 && i==anita1->NFOUR/4) {
// 	      cout << "output is " << inu << "\t" << timedomain_output_1[j][i] << "\n";
// 	    }
//	    cout << "output, bwslice_rmsdiode are " << timedomain_output_1[j][i] << "\t" << anita1->bwslice_rmsdiode[j] << "\n";
	    if (timedomain_output_1[j][i]/anita1->bwslice_rmsdiode[j]<anita1->ston[j]) {
	      anita1->ston[j]=timedomain_output_1[j][i]/anita1->bwslice_rmsdiode[j];
	      //  if (j==4 && anita1->ston[j]<0.)
	      //cout << "ston is " << anita1->ston[j] << "\n";
	    }
	  }


	  for (int i=0;i<anita1->HALFNFOUR;i++) {
	    anita1->flag_e_inanita[j][i]=0;
	    anita1->flag_h_inanita[j][i]=0;
	    anita1->timedomain_output_1_inanita[j][i]=timedomain_output_1[j][i];
	    anita1->timedomain_output_2_inanita[j][i]=timedomain_output_2[j][i];

	  }
	  for (int i=0;i<(int)flag_e[j].size();i++) {
	    anita1->flag_e_inanita[j][i+startbin]=flag_e[j][i];
	  }
	  for (int i=0;i<(int)flag_h[j].size();i++) {
	    anita1->flag_h_inanita[j][i+startbin]=flag_h[j][i];
	    
	  }
	}

      }
      
      
      
      for (int j=0;j<4;j++) {
	// this is for the e polarization
	// if we're implementing masking and the channel has been masked
	if (settings1->CHMASKING && !AntTrigger::IsItUnmasked(bn1->surfTrigBandMask,j,ilayer,ifold,0)) {
	  if (globaltrig1->channels_passing[ilayer][ifold][0][j])
	    globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]--;
	  globaltrig1->channels_passing[ilayer][ifold][0][j]=0;// channel does not pass
	}
	
	// if we're implementing masking and the channel has been masked
	if (settings1->CHMASKING && !AntTrigger::IsItUnmasked(bn1->surfTrigBandMask,j,ilayer,ifold,1)) {
	  if (globaltrig1->channels_passing[ilayer][ifold][1][j])
	    globaltrig1->nchannels_perrx_triggered[anita1->GetRx(ilayer,ifold)]--;
	  globaltrig1->channels_passing[ilayer][ifold][1][j]=0;// channel does not pass
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
// 	//cout << "i, total_diodeinput_1_inanita[j][i] are " << i << "\t" << anita1->total_diodeinput_1_inanita[4][i] << "\n";
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
AntTrigger::AntTrigger(Settings *settings1,int ilayer,int ifold,double *vmmhz,Anita *anita1,double hitangle_e,double hitangle_h,double e_component,double h_component,double *arrival_times,double volts_rx_rfcm_lab_e_all[48][512],double volts_rx_rfcm_lab_h_all[48][512])
//inline AntTrigger::AntTrigger(int ilayer,int ifold,double *vmmhz,Anita *anita1,double hitangle_e,double hitangle_h,double e_component,double h_component,double *arrival_times,int globaltrig1->rx_minarrivaltime)
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
    
    // vmmhz is V/m/MHz at the face of the antenna
    
    // this gets written to a tree because it is a measure of signal strength in the frequency domain
    integral_vmmhz=0.;
    
    for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
		integral_vmmhz+=vmmhz[ifreq]*(anita1->freq[1]-anita1->freq[0])/1.E6; // integrate vmmhz
    }
    
    
    double vhz_rx_e[Anita::NFREQ]={0.}; // V/Hz after antenna gains
    double vhz_rx_h[Anita::NFREQ]={0.};
    // same but with binning for fft
    double volts_rx_e_forfft[Anita::HALFNFOUR]={0.};
    double volts_rx_h_forfft[Anita::HALFNFOUR]={0.};
    
    // vmmhz_rx_e,h are going to be the V/m/MHz received by the rx (after gains)
    for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
		// Convert V/m/MHz to V/m/Hz and divide by dt to prepare for fft
                vhz_rx_e[ifreq]=vmmhz[ifreq]/sqrt(2.)/(anita1->TIMESTEP*1.E6); // EH, 1/sqrt(2) for dividing power in half for TDA and DDA?
		vhz_rx_h[ifreq]=vmmhz[ifreq]/sqrt(2.)/(anita1->TIMESTEP*1.E6); 
		
		// let's find the peak voltage just after the antenna, with no banding
		anita1->AntennaGain(settings1,hitangle_e,hitangle_h,e_component,h_component,ifreq,vhz_rx_e[ifreq],vhz_rx_h[ifreq]);
		// double relativegain[4];
		// anita1->AntennaGain(settings1,hitangle_e,hitangle_h,ifreq, relativegain);
		// anita1->ApplyPol(settings1,e_component,h_component,ifreq,vhz_rx_e[ifreq],vhz_rx_h[ifreq], relativegain);

		
    }
    
    // change their length from Anita::NFREQ to HALFNFOUR
    anita1->MakeArraysforFFT(vhz_rx_e,vhz_rx_h,volts_rx_e_forfft,volts_rx_h_forfft);
    
    // now the last two are in the frequency domain
    // convert to the time domain
    Tools::realft(volts_rx_e_forfft,1,anita1->HALFNFOUR); // EH, I believe this has to the -1 for inverse FFT (1 for forward FFT which is t-domain to f-domain)
    Tools::realft(volts_rx_h_forfft,1,anita1->HALFNFOUR);
    
    // now volts_rx_e,h_forfft are the time domain signals on the back end of the antenna
    // this is in volts
    
    // now find peak voltage
    // these get written to a tree
    // we still don't have any noise
    anita1->peak_e_rx_signalonly=AntTrigger::FindPeak(volts_rx_e_forfft,anita1->HALFNFOUR); // with no noise
    anita1->peak_h_rx_signalonly=AntTrigger::FindPeak(volts_rx_h_forfft,anita1->HALFNFOUR); // with no noise
    
    
    // matt - beginning of section to be replaced
    
    
    double vhz_rx_rfcm_e[Anita::NFREQ]; // V/Hz after rx, rfcm
    double vhz_rx_rfcm_h[Anita::NFREQ];
    
    
    Tools::Zero(anita1->volts_rx_rfcm_e,anita1->HALFNFOUR);// will be volts vs. time after rx, rfcm
    Tools::Zero(anita1->volts_rx_rfcm_h,anita1->HALFNFOUR);
    
    
    for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
		vhz_rx_rfcm_e[ifreq]=vhz_rx_e[ifreq]; // start with V/Hz after rx
		vhz_rx_rfcm_h[ifreq]=vhz_rx_h[ifreq];
    }
    
    // for frequency-domain voltage-based trigger (triggerscheme==0)
    // we do not apply rfcm's
    // for other trigger types we do
    
    // apply rfcm's
    if (settings1->TRIGGERSCHEME==1 || settings1->TRIGGERSCHEME==2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5) {

		 anita1->RFCMs(1,1,vhz_rx_rfcm_e);
		 anita1->RFCMs(1,1,vhz_rx_rfcm_h);
		
      // anita1->RFCMs(ilayer,ifold,vhz_rx_rfcm_e); // Linda changed
      // anita1->RFCMs(ilayer,ifold,vhz_rx_rfcm_h);
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
    anita1->MakeArraysforFFT(vhz_rx_rfcm_e,vhz_rx_rfcm_h,anita1->volts_rx_rfcm_e,anita1->volts_rx_rfcm_h);
    
    // double volts_rx_rfcm_e_freq[anita1->HALFNFOUR];
    // for (int i=0;i<anita1->HALFNFOUR;i++) {
    //   volts_rx_rfcm_e_freq[i]=anita1->volts_rx_rfcm_e[i];
    // }
        
    // now the last two are in the frequency domain
    // convert to the time domain
    // still don't have any noise
    
    Tools::realft(anita1->volts_rx_rfcm_e,1,anita1->HALFNFOUR); // EH, again, I think this has to be -1 for invFFT
    Tools::realft(anita1->volts_rx_rfcm_h,1,anita1->HALFNFOUR);
    
    anita1->GetNoiseWaveforms(); // get noise waveforms
    
    // find the peak right here and it might be the numerator of the horizontal axis of matt's plot
    anita1->peak_e_rx_rfcm_signalonly=AntTrigger::FindPeak(anita1->volts_rx_rfcm_e,anita1->HALFNFOUR); // with no noise
    anita1->peak_h_rx_rfcm_signalonly=AntTrigger::FindPeak(anita1->volts_rx_rfcm_h,anita1->HALFNFOUR);
    
    if (settings1->SIGNAL_FLUCT) {
      for (int i=0;i<anita1->NFOUR/2;i++) {
	anita1->volts_rx_rfcm_e[i]+=anita1->timedomainnoise_rfcm_e[i]; // add noise.
	anita1->volts_rx_rfcm_h[i]+=anita1->timedomainnoise_rfcm_h[i];
      }
      
    }
    
    
    
    anita1->peak_e_rx_rfcm=AntTrigger::FindPeak(anita1->volts_rx_rfcm_e,anita1->HALFNFOUR); // with noise 
    anita1->peak_h_rx_rfcm=AntTrigger::FindPeak(anita1->volts_rx_rfcm_h,anita1->HALFNFOUR); // with noise
    
    
    
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
    anita1->MakeArraysforFFT(vhz_rx_rfcm_lab_e,vhz_rx_rfcm_lab_h,anita1->volts_rx_rfcm_lab_e,anita1->volts_rx_rfcm_lab_h);
    
    
    // now the last two are in the frequency domain
    // convert to the time domain
    // still don't have any noise
    Tools::realft(anita1->volts_rx_rfcm_lab_e,1,anita1->HALFNFOUR); // EH, again, I think this has to be -1 for invFFT
    Tools::realft(anita1->volts_rx_rfcm_lab_h,1,anita1->HALFNFOUR);
    
    // matt- end of section to be replaced
    
    
    // put it in normal time ording -T to T
    // instead of 0 to T, -T to 0
    Tools::NormalTimeOrdering(anita1->NFOUR/2,anita1->volts_rx_rfcm_lab_e); // EH, why only this has NormalTimeOrdering applied? Why not before?
    Tools::NormalTimeOrdering(anita1->NFOUR/2,anita1->volts_rx_rfcm_lab_h);
    
    
    // now shift right to account for arrival times
    // for (int i=0;i<48;i++) std::cout << arrival_times[i] << std::endl;
    Tools::ShiftRight(anita1->volts_rx_rfcm_lab_e,anita1->NFOUR/2, int(arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));
    Tools::ShiftRight(anita1->volts_rx_rfcm_lab_h,anita1->NFOUR/2, int(arrival_times[anita1->GetRx(ilayer,ifold)]/anita1->TIMESTEP));
    
    
    
    
    
    if (settings1->SIGNAL_FLUCT) {
      for (int i=0;i<anita1->NFOUR/2;i++) {
	anita1->volts_rx_rfcm_lab_e[i]+=anita1->timedomainnoise_lab_e[i]; // add noise
	anita1->volts_rx_rfcm_lab_h[i]+=anita1->timedomainnoise_lab_h[i];
      }
    }
    
    for (int i=0;i<anita1->NFOUR/2;i++) {
      volts_rx_rfcm_lab_e_all[anita1->GetRx(ilayer, ifold)][i] = anita1->volts_rx_rfcm_lab_e[i];
      volts_rx_rfcm_lab_h_all[anita1->GetRx(ilayer, ifold)][i] = anita1->volts_rx_rfcm_lab_h[i];
      
    }
    
    // now vmmhz_rx_rfcm_lab_e,h_forfft are the time domain waveforms after the antenna and lab attenuation
    // now find peak voltage
    // these get written to a tree
    // still no noise
    anita1->peak_e_rx_rfcm_lab=AntTrigger::FindPeak(anita1->volts_rx_rfcm_lab_e,anita1->HALFNFOUR);
    
    anita1->peak_h_rx_rfcm_lab=AntTrigger::FindPeak(anita1->volts_rx_rfcm_lab_h,anita1->HALFNFOUR); // with no noise // EH, again I think this is WITH noise
    
    
    for (int j=0;j<5;j++) { // loop over bands
		
		for (int i=0;i<Anita::NFREQ;i++) {
			anita1->vmmhz_banding[i]=vmmhz[i]; // now copy vmmhz to vmmhz_bak instead, which we now play with to get the time domain waveforms for each subband
			// remember vmmhz is V/m/MHz at the face of the antenna
		}
		// Don't we need to apply antenna gains here?
		
		// impose banding on the incident signal
		anita1->Banding(j,anita1->freq,anita1->vmmhz_banding,Anita::NFREQ); // impose banding whatever the trigger scheme
		
		for (int i=0;i<Anita::NFREQ;i++) {
			anita1->vmmhz_banding_rfcm[i]=anita1->vmmhz_banding[i];
		}
		
		// for frequency-domain voltage-based trigger (triggerscheme==0)
		// we do not apply rfcm's
		// for other trigger types we do
		if (settings1->TRIGGERSCHEME==1 || settings1->TRIGGERSCHEME==2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5)
			anita1->RFCMs(1,1,anita1->vmmhz_banding_rfcm);
		
		
		if (settings1->TRIGGERSCHEME >=2) { // we need to prepar the signal for the diode integration
			
			
			
			for (int i=0;i<Anita::NFREQ;i++) {
				
				
				anita1->vmmhz_banding_rfcm[i]=anita1->vmmhz_banding_rfcm[i]/sqrt(2.)/(anita1->TIMESTEP*1.E6); // vmmhz was set to account for both negative and positive frequencies
				// now it has units of volts/(meter*s) so below we copy it to vm_banding_rfcm_e,h
				
				// here let's treat it as just one side of a double sided function though
				// divide by timestep so that we can use it for discrete fourier transform.
				// as in numerical recipes, F = (Delta t) * H
			}
		}
		double integral=0.;
		for (int i=0;i<Anita::NFREQ;i++) {
			vm_banding_rfcm_e[j][i]=anita1->vmmhz_banding_rfcm[i]; // this is now Volts/(m*s) vs. frequency with banding and rfcm's applied
			vm_banding_rfcm_h[j][i]=anita1->vmmhz_banding_rfcm[i];
			integral+=vm_banding_rfcm_e[j][i]*vm_banding_rfcm_e[j][i];
		} // end loop over nfreq
		
    } // end loop over bands
    
} //AntTrigger constructor







//!
/*!
 *
 *
 *
 *
 *
 */
AntTrigger::AntTrigger() {
    
}






//!
/*!
 *
 *
 *
 *
 *
 */
double AntTrigger::FindPeak(double *waveform,int n) {
    
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
void AntTrigger::addToChannelSums(Settings *settings1,Anita *anita1,int ibw, int k) {
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
    
    
    //  for (int iband=0;iband<4;iband++) {
    //if (freq>=bwslice_min[iband] && freq<bwslice_max[iband]) {
    //bwslice_volts_pol0[iband]+=lcp_component; // adding voltage in lcp polarization
    //bwslice_energy_pol0[iband]+=lcp_component_energy; // adding energy in lcp polarization
    //bwslice_volts_pol1[iband]+=rcp_component; // rcp is latter 4 components of bwslice_volts
    //bwslice_energy_pol1[iband]+=rcp_component_energy; // add up rcp energy for each bandwidth slice
    
    
    //bwslice_volts_pol0_em[iband]+=lcp_component*vmmhz_em; // David G. - don't you want to divide this by vmmhz[k]?
    //bwslice_volts_pol1_em[iband]+=rcp_component*vmmhz_em;
    
    bwslice_volts_pole[ibw]+=term_e;
    bwslice_energy_pole[ibw]+=energyterm_e;
    bwslice_volts_polh[ibw]+=term_h;
    bwslice_energy_polh[ibw]+=energyterm_h;
    
    //}
    //}
    
    
    //for (int i=0;i<4;i++) {
    //    bwslice_vnoise[ilayer][i]=AntTrigger::GetNoise(altitude_bn,surface_under_balloon,THETA_ZENITH[ilayer],Bands::bwslice_max[i]-Bands::bwslice_min[i],0.);
    //cout << "getnoise is " <<  i << " " << AntTrigger::GetNoise(altitude_bn,surface_under_balloon,THETA_ZENITH[ilayer],Bands::bwslice_max[i]-Bands::bwslice_min[i],0.) << "\n";
    //bwslice_vnoise[ilayer][i]*=THERMALNOISE_FACTOR; // for systematic studies- scale the noise. nominally equal to 1.
    //} //for loop over bandwidth slices
    
}




//!	Takes ADC threshold as an input and outputs threshold in power
/*!
 *	"this is a place holder for now, from reading off of various plots, pointed out below"
 *	
 */
double AntTrigger::ADCCountstoPowerThreshold(int threshadc, int isurf,int ichan) {
    // first convert threshold in adc counts to the singles rate
    // do this using Ryan's threshold scans from Elog #133
    // these curves were read in the Trigger constructor
    // first check if the threshold in adc counts is in the allowable range
    if (unwarned && (threshadc<minadcthresh[isurf][ichan] || threshadc>maxadcthresh[isurf][ichan]))
		cout << "Warning! ADC threshold is outside range of measured threshold scans.";
    if (threshadc<minadcthresh[isurf][ichan]) {
		if (unwarned) {
			cout << "It is below the minimum so set it to the minimum.  Will not be warned again.\n";
			unwarned=0;
		}
		threshadc=minadcthresh[isurf][ichan];
    }
    if (threshadc>maxadcthresh[isurf][ichan]) {
		if (unwarned) {
			cout << "It is higher than the maximum so set it to the maximum.  Will not be warned again.\n";
			unwarned=0;
		}
		threshadc=maxadcthresh[isurf][ichan];
    }
    
    
    // Now find singles rate for this threshold
    // first sort thresholds
    int index=TMath::BinarySearch(NPOINTS,threshold[isurf][ichan],threshadc);
    //cout << "rate is " << rate[isurf][ichan][index] << "\n";
    thisrate=(double)rate[isurf][ichan][index]; // these scalers are in kHz
    
    
    // now find threshold for this scaler.  Here, scalers have to be in MHz.
    thisrate=thisrate/1.E3; // put it in MHz
    //cout << "thisrate is " << thisrate << "\n";
    // figure out what band we're talking about
    int iband=Anita::SurfChanneltoBand(ichan);
    thispowerthresh=rateToThreshold(thisrate,iband);
    cout << "thisrate, iband, thispowerthresh are " << thisrate << " " << iband << " " << thispowerthresh << "\n";
    return thispowerthresh;
    
    
    //  double powerthresh=-0.3*pow(10.,logsingles)/1.E6+5.091; // this has the right slope and intercept
    
    //cout << "threshadc, powerthresh are " << threshadc << " " << powerthresh << "\n";
    
    //return powerthresh;
    
    return 0;
}



//!	Returns the thisrate variable value
/*!
 *	\todo	Refactor the class so that this function is either deprecated or "thisrate"
 *			is implemented differently.
 */
double AntTrigger::getRate() {
    
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
double AntTrigger::rateToThreshold(double rate, int band)
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
int AntTrigger::IsItUnmasked(unsigned short surfTrigBandMask[9][2],int ibw,int ilayer, int ifold, int ipol) {
    
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




void AntTrigger::L1Trigger(Anita *anita1,double timedomain_output_1[5][Anita::NFOUR],double timedomain_output_2[5][Anita::NFOUR],double *powerthreshold,int *channels_passing_e_forglob,int *channels_passing_h_forglob,int &npass) {
   
  int maxsample=TMath::MaxElement(5,anita1->imaxbin);
  int minsample=TMath::MinElement(5,anita1->iminbin);
    
    for (int j=0;j<5;j++) {
      flag_e[j].clear();
      flag_h[j].clear();		
      
      for (int i=minsample;i<maxsample;i++) {
	
	if (timedomain_output_1[j][i]<powerthreshold[j]*anita1->bwslice_rmsdiode[j] && anita1->bwslice_allowed[j]==1) {
	  flag_e[j].push_back(1);
	  
	}
	else
	  flag_e[j].push_back(0);
	
	
	if (timedomain_output_2[j][i]<powerthreshold[j]*anita1->bwslice_rmsdiode[j] && anita1->bwslice_allowed[j]==1)
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
double AntTrigger::GetNoise(Settings *settings1,double altitude_bn,double geoid,double theta_zenith,double bw,double temp) {
    
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
void AntTrigger::GetThresholds(Settings *settings1,Anita *anita1,int ilayer,double *thresholds) {
    
    if (ilayer==3 && settings1->DISCONES==2) // if it's a nadir layer
		for (int i=0;i<5;i++)
			thresholds[i]=anita1->powerthreshold_nadir[i];
    else
		for (int i=0;i<5;i++)
			thresholds[i]=anita1->powerthreshold[i];
    
}




//!	Evaluate the full trigger simulation for a given payload configuration
/*!
 *	PassesTrigger handles the logistics for the triggering of all trigger systems.
 *	
 *	The triggering system begins with the calculations to determine which channels will pass the "L0" trigger.
 *	Each antenna on Anita 2 had two polarizations, and each polarization was divided up into five frequency bands.
 *	The frequency bands would be independently measured to see if the signal was over the threshold.
 *	If 3 of the 8 bands for an antenna (in some experiments, at least) exceeded the threshold then that antenna
 *	would have triggered an "L1", or antenna-level trigger.
 *
 *	From there, the different triggering systems start to diverge, each requiring some different combinatoric trigger
 *	to be satisfied before escalating to the next-higher-level trigger.
 *
 *	The trigger system using TRIGGERSCHEME==3 is a good deal different from the other triggers, because of a few thing:
 *		- All bands were assumed to have passed
 *		- The trigger scans over many possible physical propagation delays in an attempt to find the highest power
 *		- If the highest power for that antenna group happens to exceed the set threshold, then...
 *		- If two neighboring phi sectors find powers above the threshold, then a global trigger occurs
 *
 *	The antenna trigger clusters have been different configurations in the past, with early Anita3 designs using the
 *	"Coherent Sum Trigger", which uses 9 antennas per cluster. More current trigger designs use 3 antennas per
 *	coherent-summing-group, though this more recent trigger is referred to as the "Summed Power Trigger"
 *
 *
 *	\todo	This function needs to be heavily refactored into either functions which perform some smaller part of
 *			every trigger, or some functions which perform one triggering system in its entirely, or a combination
 *			of the two.
 *
 *			The longer a function gets, the less confidence can be had in it because instead of many small, well-tested
 *			modular functions, functions of this length develop into massive hacks which may appear to perform the same
 *			but no such guarantee can be made.
 *
 *			There are a lot of functions and code snippets which could be very easily replaced by standard library
 *			functions, or by embracing "RAII" ("Resource Acquisition Is Initialization") and using fewer C arrays and
 *			more of the standard library functions (particularly those in "algorithm", "functional", and "numeric"
 *			and above all the container libraries. Those libraries will help (along with the refactor) alleviate the
 *			difficult-to-follow nesting of for-loops within conditional statements.
 *			
 *			Many of the combinatoric triggers can have a logical mask made in order to very clearly define and later
 *			comprehend the requirements of the triggers, and less nesting of loops is necessary.
 *
 *			Using scoped accumulators for calculating statistics is also a good idea here, as it reduces the number of
 *			functions and the number of objects needed to produce a distribution of values to just one object, making
 *			mistakes/errors much less likely. *			
 *
 *			There is a decent amount of dead code which should be pruned, as well.
 */

int GlobalTrigger::PassesTrigger(Settings *settings1, Anita *anita1, int discones_passing, int mode, int &l3trig, int *l2trig, int *l1trig, int antennaclump, int loctrig[Anita::NLAYERS_MAX][Anita::NPHI_MAX], int *loctrig_nadironly, int inu, bool ishpol) {
  double this_threshold= -4.34495;
  PassesTrigger(settings1,anita1,discones_passing,mode,l3trig,l2trig,l1trig,antennaclump,loctrig,loctrig_nadironly,inu,this_threshold, ishpol);   

}
int GlobalTrigger::PassesTrigger(Settings *settings1, Anita *anita1, int discones_passing, int mode, int &l3trig, int *l2trig, int *l1trig, int antennaclump, int loctrig[Anita::NLAYERS_MAX][Anita::NPHI_MAX], int *loctrig_nadironly, int inu,double this_threshold, bool ishpol) {

  //bool ishpol should only be used for anita3, by default do vpol
  
  int ltsum=0;
  int channsum=0;
  int ihit=0;
  int thispasses=0;
  
  
  if (settings1->TRIGGERSCHEME < 3) {
    // this is an array with 1=pass and 0=fail for each channel
    // the first two layers on the payload are "compacted"
    // into one trigger layer of 16 antennas.
    // layer number, antenna, polarization, bandwidth slice
    int channels_compacted_passing[Anita::NLAYERS_MAX][Anita::NPHI_MAX][2][5] = {{{{0}}}};
    
    l3trig=0;
    Tools::Zero(l1trig,Anita::NTRIGGERLAYERS_MAX);
    Tools::Zero(l2trig,Anita::NTRIGGERLAYERS_MAX); // for all three layers
    Tools::Zero(loctrig_nadironly,Anita::NPHI_MAX);
    // which clumps of 3 antennas
    //within the nadir layer pass
    // 0th and 1st element= antenna at phi=0,
    // 2nd and 3rd element=next antenna, etc.
    
    
    for (int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
      
      // for (int j=0;j<Anita::NPHI_MAX;j++) {
	Tools::Zero(loctrig[ilayer],Anita::NPHI_MAX);
	// which clumps of 3 antennas
	//within each trigger layer pass L2 trigger
	// layer, phi location
      // } //for
    } // end of initializing to zero
    
    int whichlayer=0;
    int whichphisector=0;
    for (int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
      for (int iphi=0;iphi<anita1->NRX_PHI[ilayer];iphi++) {
	for (int ipolar=0;ipolar<2;ipolar++) {

	  if (settings1->WHICH==9 && ishpol && ipolar==0) continue; // Anita3 : only do the polarisation required
	  else if (settings1->WHICH==9 && !ishpol && ipolar==1) continue; // Anita3 : only do the polarisation required
	  
	  for (int iband=0;iband<5;iband++) {
	    
	    GetAnitaLayerPhiSector(settings1,ilayer,iphi,whichlayer,whichphisector); // Translate Anita physical layer to Anita trigger layer and phi sector (4 layers with 8,8,16,8 phi sector to 3 layers with 16 phi sectors each.  In the nadir layer, 8 of the 16 slots are empty on the trigger layer.)
	    
	    
	    
	    // combining top two layers on the payload into one trigger layer
	    // this means the 0th antenna in the first trigger layer is
	    // physically higher than the 1st antenna in the first trigger layer
	    // which means that the nadirs are aligned with the antennas with indices 1,3,5 etc.
	    // we will still use indices 0-7 for them though
	    //	  channels_compacted_passing[0][2*iphi+ilayer][ipolar][iband]+=channels_passing[ilayer][iphi][ipolar][iband];
	    channels_compacted_passing[whichlayer][whichphisector][ipolar][iband]+=channels_passing[ilayer][iphi][ipolar][iband];
	  } //for
	} //for
      } //for
    } //for
    
    int antsum=0; // counter for how many channels on an antenna pass
    /* antenna triggers */
    
    int NBAND=5; // number of bandwidth slices
    int Nreq_l2[Anita::NLAYERS_MAX];
    // number of antennas within each clump
    //that need to pass for that clump to pass L2
    // 1st layer, 2nd layer,
    //nadir layer (to be "OR"ed with upper layers),
    // and nadir layer (for nadir-only) trigger
    //int Ntrig_chann[4]={6,6,6,8}; // NOT USED:  number of channels
    //per clump that need to pass.
    for (int i=0;i<Anita::NLAYERS_MAX;i++) {
      Nreq_l2[i]=anita1->trigRequirements[1];
    }
    
    
    int ant[Anita::NLAYERS_MAX][Anita::NPHI_MAX]; // which antennas pass at L1
    // trigger layer, phi position
    
    for (int i=0;i<Anita::NLAYERS_MAX;i++) {
      for (int j=0;j<Anita::NPHI_MAX;j++) {
	ant[i][j]=0;
      } //for
    } //for
    
    int antpass=0; // count how many antennas pass
    int iphitrig=0;  // which trigger phi sector
    int ihittrig=0;
    // local level 1 trigger at the antenna
    for(int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++) { // loop over layers
      // this is nlayers-1 because NLAYERS counts top 16 antennas as 2 layers
      
      for(int iphi=0;iphi<anita1->PHITRIG[iloc];iphi++) { // loop over phi position
	iphitrig=GetPhiSector(settings1,iloc,iphi); // get trigger phi sector
	// counts from 0
	antsum = 0; // start antenna sum

	for(int ipolar=0;ipolar<2;ipolar++) {

	  if (settings1->WHICH==9 && ishpol && ipolar==0) continue; // Anita3 : only do the polarisation required
	  else if (settings1->WHICH==9 && !ishpol && ipolar==1) continue; // Anita3 : only do the polarisation required

	  for(int iband=0;iband<NBAND;iband++) {
	    if(channels_compacted_passing[iloc][iphitrig][ipolar][iband] == 1
	       && anita1->bwslice_allowed[iband]==1 && anita1->pol_allowed[ipolar]==1) { // only increment if it's one of the allowed bands.
	      
	      if (settings1->PHIMASKING && settings1->WHICH==9) // only applying channel masking like this if it's 
		if ((ipolar==0 && (1<<iphitrig & l1TrigMask)) || (ipolar==1 && (1<<iphitrig & l1TrigMaskH)) )
		  continue; // was this channel masked?
	      antsum = antsum +1; // sum channels that pass for this antenna

	    }
	  } // loop over bands
	} // end loop over polarizations
	
	
	for (int ipolar=0;ipolar<2;ipolar++) {
	  if (settings1->WHICH==9 && ishpol && ipolar==0) continue; // Anita3 : only do the polarisation required
	  else if (settings1->WHICH==9 && !ishpol && ipolar==1) continue; // Anita3 : only do the polarisation required
	  for (int iband=0;iband<NBAND;iband++) { // notice sum over 5 bands now
	    if(channels_compacted_passing[iloc][iphitrig][ipolar][iband] == 0
	       && anita1->bwslice_required[iband]==1 && anita1->pol_allowed[ipolar]==1) { // if this band was required to pass and it didn't,
	      antsum = 0; // fatal for this antenna
	    }
	  }
	}

	// if the required bands didn't pass then set antsum=0 so the antenna doesn't pass	  
	if(antsum >= anita1->trigRequirements[0])  { // do more than 3 channels out of 8 pass?
	  ant[iloc][iphitrig] = 1; // then this antenna passes L1.
	  antpass++;
	} //if      

      } //for (phi position)
    } //for (layers)
    
    if (settings1->DISCONES==2) {
      // fill the slots in  between the actual antennas with the "or"
      // of the neighboring antennas
      FillInNadir(anita1,ant[2]);
    }
    
    for(int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++){ //This is NLAYERS-1 because NLAYERS counts top 16 antennas as 2 layers
      for (int iphi=0;iphi<anita1->PHITRIG[iloc];iphi++){
	iphitrig=GetPhiSector(settings1,iloc,iphi);
	if (ant[iloc][iphitrig]==1) {
	  l1trig[iloc] += (1<<iphitrig); // this keeps track of which antennas pass the l1 trigger
	}
      }//for (phi position)
    } //for (layers) // Hmm this is not used
    
    if (mode == 1 && antpass >= 2) {// TRIGTYPE=2 -> just require ANTtrig channels pass on 2 antennas
      //return 1;
      thispasses = 1;
    } else if (mode == 2) { // TRIGTYPE=1 ->
      
      // int iloc=0; // top array
      //int inad=0;  // phi position for nadir antennas,
      //only goes to 8, wheras iphi goes to 16
      int x;//for number of antenna in clump
      int y;//for number of antenna in clump
      x = (antennaclump-1)/2;//setting boundries for clump dependent on clump size(antennaclump)
      y = (antennaclump-1)/2 +1;
      
      for (int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++) {
	for(int iphi=0;iphi<16;iphi++) {
	  iphitrig=iphi;
	  ltsum = 0; // zero counter for number of antennas that pass per clump
	  channsum=0; // zero counter for number of channels that pass per clump
	  
	  // loop over clump -- with boundary check -- ugly implement
	  // even uglier - we want to require 2/3 antennas per clump
	  // but the middle antenna has to be one of the two
	  // so we multiply by ant[iloc][iphi] so that when the middle one
	  // passes we are multiplying by 1 and when it doesn't we are
	  // multiplying by zero.
	  for(int inbr=iphi-x;inbr<iphi+y;inbr++) {
	    ihit = (inbr + 16) % 16; // make any negative indices positive
	    ihittrig=ihit;
	    
	    if (anita1->REQUIRE_CENTRE) {
	      ltsum+=ant[iloc][ihittrig]*ant[iloc][iphitrig]; // increment if this antenna passes
	    }
	    // we multipy by ant[iloc][iphi] because this ensures that
	    // one of the triggered antennas is the center one in the clump
	    else
	      ltsum+=ant[iloc][ihittrig];
	    
	    for (int ipolar=0;ipolar<2;ipolar++) {
	      if (settings1->WHICH==9 && ishpol && ipolar==0) continue; // Anita3 : only do the polarisation required
	      else if (settings1->WHICH==9 && !ishpol && ipolar==1) continue; // Anita3 : only do the polarisation required
	      for (int iband=0;iband<4;iband++) {
		channsum+=channels_compacted_passing[iloc][ihittrig][ipolar][iband]; // increment if this channel passes
	      } //for
	    } //for
	  } //for
	  
	  if(ltsum >= Nreq_l2[iloc]){ // if enough antennas in this clump pass
	    loctrig[iloc][iphitrig] = 1;
	    l2trig[iloc] += (1<<iphi);
	  }//if
	} //for (phi position)
      } // end loop over trigger layers
      
      // fill the slots in  between the actual antennas with the "or"
      // of the neighboring antennas
      if (settings1->DISCONES==2)
	FillInNadir(settings1,anita1,l2trig[2]);
      
      if (settings1->PHIMASKING && settings1->WHICH!=9) { // not for Anita3
	// nadir antennas are aligned with the second physical layer of antennas
	for (int iphi=0;iphi<anita1->NTRIGPHISECTORS;iphi++) { // loop over phi sectors
	  if ((1<<iphi) & phiTrigMask) { // if this phi sector is masked
	    for (int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++) {
	      loctrig[iloc][iphi] = 0;
	      if(l2trig[iloc] & (1<<iphi)) l2trig[iloc] -= (1<<iphi);
	    }
	  } // if this phi sector is masked
	}
      } else if (settings1->PHIMASKING && settings1->WHICH==9) { // only for Anita3
	// nadir antennas are aligned with the second physical layer of antennas
	for (int iphi=0;iphi<anita1->NTRIGPHISECTORS;iphi++) { // loop over phi sectors
	  if ((((1<<iphi) & phiTrigMask) && !ishpol) || (((1<<iphi) & phiTrigMaskH) && ishpol))  { // if this phi sector is masked
	    for (int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++) {
	      loctrig[iloc][iphi] = 0;
	      if(l2trig[iloc] & (1<<iphi)) l2trig[iloc] -= (1<<iphi);
	    }
	  }// if this phi sector is masked
	}
      }
      
      if (anita1->trigRequirements[2]==0) {
	for (int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++) {
	  for (int iphi=0;iphi<anita1->NTRIGPHISECTORS;iphi++) {
	    if (l2trig[iloc] & (1<<iphi)) { // if we are not making a L3 trigger requirement and there was a l2 trigger then call it a pass
	      thispasses=1;
	      //cout << "This one passes.  inu is " << inu << " " << "iloc is " << iloc << "\n";
	    }
	  }
	}
      }
      
      thispasses=L3Trigger(settings1,anita1,loctrig,loctrig_nadironly,discones_passing,l3trig);
      
    } //else if (mode 2)
  } else if (settings1->TRIGGERSCHEME == 3) {
    // This is the case for (settings1->TRIGGERSCHEME == 3), that coherent waveform sum is being used to trigger, based on Andres' idea. This is specific to ANITA III's antenna geometry.
    // Written by Paul Schellin 2011
    
    // Method:
    //		All bands have passed so far, we now perform the following for the first-hit phi sector and its neighbors:
    //		    - Convert the relevent waveforms into a useful format
    //		    - Calculate all arrival time delays relevent to the 3 phi sectors (9 antennas)
    //		    f- For each of the direction hypotheses, do:
    //			- Shift the waveforms according to their delay offsets
    //			- Sum the 9 waveforms, resulting in one coherently-summed waveform
    //			- Square each element of the summed waveform, resulting in the power of the coherently-summed waveform
    //			- For each of the 16 sample windows (with length of 32 time-samples) do:
    //			    - Sum the powers of those samples
    //			    - Compare the sum to the threshold, if greater, trigger on this event (return 1)
    // Notes:
    //		There will be several things hardcoded into the following method, feel free to change these to be variables within the Settings class, but for the sake of sanity, PLEASE no more global variables!
    //		This will be made to implement all types of payloads shortly, what exists below is only a temporary specialization for ANITA III.
    double timesteps[anita1->HALFNFOUR];

    for (unsigned int i = 0; i < anita1->HALFNFOUR; i++){
      timesteps[i] = i;
    }

    for (int center_phi_sector_offset = -1; center_phi_sector_offset <= 1; center_phi_sector_offset++){
      int center_phi_sector_index = first_phi_sector_hit + center_phi_sector_offset;
      if (center_phi_sector_index > 15){center_phi_sector_index = 0;}
      if (center_phi_sector_index < 0){center_phi_sector_index = 15;}
      
      // Method:
      //		All bands have passed so far, we now perform the following for the first-hit phi sector and its neighbors:
      //		    - Convert the relevant waveforms into a useful format
      //		    - Calculate all arrival time delays relevant to the 3 phi sectors (9 antennas)
      //		    - For each of the direction hypotheses, do:
      //			- Shift the waveforms according to their delay offsets
      //			- Sum the 9 waveforms, resulting in one coherently-summed waveform
      //			- Square each element of the summed waveform, resulting in the power of the coherently-summed waveform
      //			- For each of the 16 sample windows (with length of 32 time-samples) do:
      //			    - Sum the powers of those samples
      //			    - Compare the sum to the threshold, if greater, trigger on this event (return 1)
      // Notes:
      //		There will be several things hard-coded into the following method, feel free to change these to be variables within the Settings class, but for the sake of sanity, PLEASE no more global variables!
      //		This will be made to implement all types of payloads shortly, what exists below is only a temporary specialization for ANITA III.
      
      bool coherent_trigger_passes = false;
      
      unsigned N_STEP_PHI = 16;	// The number of different phi hypotheses for a given center phi sector
      unsigned N_STEP_THETA = 81;	// The number of different theta hypotheses for a given center phi sector
      unsigned N_PHI_SECTORS = 16;	// The number of phi sectors
      //unsigned N_LAYERS_PHYSICAL = 4;   // The number of physical layers on the payload
      unsigned N_LAYERS_TRIGGER = 3;	// The number of layers as seen by the trigger
      //unsigned N_SUMMED_SECTORS = 3;    // The number of phi sectors clumped together by the coherent sum trigger
      
      double highest_power = 0;
      unsigned hi_pow_center = 0;
      unsigned hi_pow_phi_index = 0;
      unsigned hi_pow_theta_index = 0;
     
      
      // Here the 48 antennas are filled according to the waveforms passed to PassesTrigger(...).
      unsigned fill_index = 0;
      for (unsigned fill_index_phi_sector = 0; fill_index_phi_sector < 16; ++fill_index_phi_sector) {
	for (unsigned fill_index_layer = 0; fill_index_layer < 3; ++fill_index_layer) {
	  anita1->cwst_RXs[fill_index].phi_sector = fill_index_phi_sector;
	  anita1->cwst_RXs[fill_index].layer = fill_index_layer;
	  
	  unsigned physical_layer_index = ((fill_index_phi_sector%2) ? (fill_index_layer + 1) : (0));   // Maps the trigger layers {0, 1, 2} to the physical layers {0, 1, 2, 3}.
	  if (fill_index_phi_sector%2 == 0) {
	    physical_layer_index = (fill_index_layer == 0) ? 0 : (fill_index + 1);
	  }
	  // If phi sector is even then trigger layer zero maps to physical layer zero. If phi odd, maps to layer 1.
	  unsigned physical_phi_index (fill_index_phi_sector);
	  if (fill_index_layer == 0) {
	    physical_phi_index = unsigned(fill_index_phi_sector / 2.);   // Map {0, ..., 15} to {0, ..., 7}.
	  }
	  anita1->cwst_RXs[fill_index].x = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][0];
	  anita1->cwst_RXs[fill_index].y = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][1];
	  anita1->cwst_RXs[fill_index].z = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][2];
	    
	  anita1->cwst_RXs[fill_index].waveform->assign(volts_rx_rfcm_trigger[fill_index_phi_sector][fill_index_layer].begin(),volts_rx_rfcm_trigger[fill_index_phi_sector][fill_index_layer].end());
	  
	  for (unsigned fill_index_timestep = 0; fill_index_timestep < anita1->HALFNFOUR; ++fill_index_timestep) {
	    anita1->cwst_RXs[fill_index].digitized->at(fill_index_timestep) = three_bit_round(anita1->cwst_RXs[fill_index].waveform->at(fill_index_timestep)/ anita1->rms_rfcm_e_single_event, (anita1->summed_power_trigger_digitizer_zero_random->Rndm() >= 0.5), false);
	  }
	  ++fill_index;
	}
      }
      
      /*
	double hypothesis_offset[3][3];
	if (inu == 2564) {
	for (double deg_phi = 0.; deg_phi < 360.; deg_phi += 1.) {
	anita1->calculate_single_offset(12, deg_phi, -14.082, hypothesis_offset);
	
	vector <double> summed_wfm(anita1->HALFNFOUR, 0.);
	unsigned center_phi_sector_index = 12;
	
	for (int fill_index_phi_sector_offset = -1; fill_index_phi_sector_offset <= 1; ++fill_index_phi_sector_offset) {
	unsigned fill_index_phi_sector = (center_phi_sector_index + fill_index_phi_sector_offset + N_PHI_SECTORS)%N_PHI_SECTORS;
	unsigned fill_index = 0;
	for (unsigned fill_index_layer = 0; fill_index_layer < N_LAYERS_TRIGGER; ++fill_index_layer) {
	unsigned rx_index = fill_index_phi_sector * 3 + fill_index_layer;
	anita1->cwst_aligned_wfms[fill_index].phi_sector = fill_index_phi_sector;
	anita1->cwst_aligned_wfms[fill_index].layer = fill_index_layer;
	
	unsigned physical_layer_index = ((fill_index_phi_sector%2) ? (fill_index_layer + 1) : (0));   // Maps the trigger layers {0, 1, 2} to the physical layers {0, 1, 2, 3}.
	// If phi sector is even then trigger layer zero maps to physical layer zero. If phi odd, maps to layer 1.
	unsigned physical_phi_index;
	if (fill_index_layer == 0) {
	physical_phi_index = int(fill_index_phi_sector / 2.);   // Map {0, ..., 15} to {0, ..., 7}.
	}
	anita1->cwst_aligned_wfms[fill_index].x = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][0];
	anita1->cwst_aligned_wfms[fill_index].y = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][1];
	anita1->cwst_aligned_wfms[fill_index].z = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][2];
	
	unsigned time_offset = hypothesis_offset[fill_index_phi_sector_offset + 1][fill_index_layer];
	
	for (unsigned fill_index_timestep = 0; fill_index_timestep < anita1->HALFNFOUR - time_offset; ++fill_index_timestep) {
	anita1->cwst_aligned_wfms[fill_index].digitized->at(fill_index_timestep) = anita1->cwst_RXs[rx_index].digitized->at(fill_index_timestep + time_offset);
	summed_wfm.at(fill_index_timestep) += anita1->cwst_aligned_wfms[fill_index].digitized->at(fill_index_timestep);
	}
	
	for (unsigned fill_index_timestep = anita1->HALFNFOUR - time_offset; fill_index_timestep < anita1->HALFNFOUR; ++fill_index_timestep) {
	// Fill the ends of the waveforms with zeros. This should not negatively affect the simulation at all.
	// This is an attempt to fix the power = 648 issue, where the bins were all +0.5 and summing to 9*0.5, summed over 32 bins. Bah.
	anita1->cwst_aligned_wfms[fill_index].digitized->at(fill_index_timestep) = 0.;
	}
	++fill_index;
	
	} // for fill layer indices
	} // for fill phi sector indices
	
	vector <double> power_of_summed_wfm;
	square_waveform_elements(summed_wfm, power_of_summed_wfm);
	
	//for (unsigned window_index = 0; window_index < power_of_summed_wfm.size() - 32; window_index += 16) {
	//double power = summed_power_window(power_of_summed_wfm, window_index, 32);
	//if (highest_power < power){
	//highest_power = power;
	//}
	//if (true){
	anita1->fill_coherent_waveform_sum_tree(inu, center_phi_sector_index, settings1, anita1->rms_rfcm_e_single_event, 0., 0, 0 + 32, -14.082, deg_phi, 33., 33.,summed_wfm,	power_of_summed_wfm, 0.);
	//return 1;	// Return that this waveform passes.
	//}
	//} // for window indices
	}
	return 1;
	}
      */
      
      for (unsigned center_phi_sector_index = 0; center_phi_sector_index < N_PHI_SECTORS; ++center_phi_sector_index) {
	// Loop over the hypotheses to align waveforms.
	for (unsigned index_phi = 0; index_phi < N_STEP_PHI; ++index_phi) {
	  for (unsigned index_theta = 0; index_theta < N_STEP_THETA; ++index_theta) {
	    unsigned fill_index = 0;
	    vector <double> summed_wfm(anita1->HALFNFOUR, 0.);
	    
	    for (int fill_index_phi_sector_offset = -1; fill_index_phi_sector_offset <= 1; ++fill_index_phi_sector_offset) {
	      unsigned fill_index_phi_sector = (center_phi_sector_index + fill_index_phi_sector_offset + N_PHI_SECTORS)%N_PHI_SECTORS;
	      
	      for (unsigned fill_index_layer = 0; fill_index_layer < N_LAYERS_TRIGGER; ++fill_index_layer) {
		unsigned rx_index = fill_index_phi_sector * 3 + fill_index_layer;
		anita1->cwst_aligned_wfms[fill_index].phi_sector = fill_index_phi_sector;
		anita1->cwst_aligned_wfms[fill_index].layer = fill_index_layer;
		
		unsigned physical_layer_index = ((fill_index_phi_sector%2) ? (fill_index_layer + 1) : (0));   // Maps the trigger layers {0, 1, 2} to the physical layers {0, 1, 2, 3}.
		// If phi sector is even then trigger layer zero maps to physical layer zero. If phi odd, maps to layer 1.
		unsigned physical_phi_index (fill_index_phi_sector);
		if (fill_index_layer == 0) {
		  physical_phi_index = unsigned(fill_index_phi_sector / 2.);   // Map {0, ..., 15} to {0, ..., 7}.
		}
		anita1->cwst_aligned_wfms[fill_index].x = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][0];
		anita1->cwst_aligned_wfms[fill_index].y = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][1];
		anita1->cwst_aligned_wfms[fill_index].z = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][2];
			  
		unsigned time_offset = anita1->hypothesis_offsets[center_phi_sector_index][index_phi][index_theta][fill_index_phi_sector_offset + 1][fill_index_layer];
		
		for (unsigned fill_index_timestep = 0; fill_index_timestep < anita1->HALFNFOUR - time_offset; ++fill_index_timestep) {
		  anita1->cwst_aligned_wfms[fill_index].digitized->at(fill_index_timestep) = anita1->cwst_RXs[rx_index].digitized->at(fill_index_timestep + time_offset);
		  summed_wfm.at(fill_index_timestep) += anita1->cwst_aligned_wfms[fill_index].digitized->at(fill_index_timestep);
		}
		
		for (unsigned fill_index_timestep = anita1->HALFNFOUR - time_offset; fill_index_timestep < anita1->HALFNFOUR; ++fill_index_timestep) {
		  // Fill the ends of the waveforms with zeros. This should not negatively affect the simulation at all.
		  // This is an attempt to fix the power = 648 issue, where the bins were all +0.5 and summing to 9*0.5, summed over 32 bins.
		  anita1->cwst_aligned_wfms[fill_index].digitized->at(fill_index_timestep) = 0.;
		}
		++fill_index;
		
	      } // for fill layer indices
	    } // for fill phi sector indices
	    
	    vector <double> power_of_summed_wfm;
	    square_waveform_elements(summed_wfm, power_of_summed_wfm);
	    
	    for (unsigned window_index = 0; window_index < power_of_summed_wfm.size() - 32; window_index += 16) {
	      double power = summed_power_window(power_of_summed_wfm, window_index, 32);
	      if (highest_power < power){
		highest_power = power;
		hi_pow_center = center_phi_sector_index;
		hi_pow_phi_index = index_phi;
		hi_pow_theta_index = index_theta;
	      }
	      
	      if (power >= settings1->COHERENT_THRESHOLD){
		coherent_trigger_passes = true;
		thispasses += 2;
	      }
	      
	      if (power >= settings1->COHERENT_THRESHOLD){
		anita1->fill_coherent_waveform_sum_tree(inu, center_phi_sector_index, settings1, anita1->rms_rfcm_e_single_event, 0., window_index, window_index + 32, index_theta - 45., index_phi, 33., 33.,summed_wfm, power_of_summed_wfm, power);
		//thispasses = true;
		//return 1;	// Return that this waveform passes.
	      }
	      
	      if (inu == 1808 && center_phi_sector_index == 6) {
		anita1->fill_coherent_waveform_sum_tree(inu, center_phi_sector_index, settings1, anita1->rms_rfcm_e_single_event, 0., window_index, window_index + 32, index_theta - 45., index_phi, 33., 33.,summed_wfm, power_of_summed_wfm, power);
		//return 1;
	      }
	      
	      if (coherent_trigger_passes){
		return thispasses;
	      }
	      
	    } // for window indices
	  } // for hypothesis theta indices
	} // for hypothesis phi indices
      } // for hypothesis center phi sector indices
      
      /*
	cout << "\nEvent number\t" << inu;
	cout << "\nHighest power was\t" << highest_power;
	cout << "\nCenter phi sector\t" << hi_pow_center;
	cout << "\nHighest phi index\t" << hi_pow_phi_index;
	cout << "\ntHighest theta index\t" << hi_pow_theta_index;
      */
    }
  } else if (settings1->TRIGGERSCHEME == 4) {
    //	TRIGGERSCHEME == 4 is the Summed Power Trigger.
    //	For every window, all of the phi sectors find the maximum coherently summed power.
    //	These max powers are summed with the neighboring phi sectors (for three phi sectors in total) and compared to a threshold.
    
    double		SummedPowerThreshold	= settings1->COHERENT_THRESHOLD;	// The threshold against which all of the events will be compared.
    double		SummedPowerThetaMin		= -45.;		// DEGREE!	The minimum theta (elevation) angle to try in the hypotheses.
    double		SummedPowerThetaMax		= +25.;		// DEGREE!	The maximum theta (elevation) angle to try in the hypotheses.
    double		SummedPowerThetaStep	= 1.;		// DEGREE!	The angle to advance from minimum theta to maximum theta.
    
    // unsigned	N_STEP_THETA			= 81;		// The number of different theta hypotheses for a given center phi sector
    // unsigned	N_PHI_SECTORS			= 16;		// The number of phi sectors
    // unsigned	N_LAYERS_TRIGGER		= 3;		// The number of layers as seen by the trigger
    
    // double		highest_power			= 0.;
    // unsigned	hi_pow_center			= 0;
    // unsigned	hi_pow_phi_index		= 0;
    // unsigned	hi_pow_theta_index		= 0;
    
    // Here the 48 antennas are filled according to the waveforms passed to PassesTrigger(...).
    unsigned	fill_index				= 0;
    
    for (unsigned fill_index_phi_sector = 0; fill_index_phi_sector < 16; ++fill_index_phi_sector) {
      for (unsigned fill_index_layer = 0; fill_index_layer < 3; ++fill_index_layer) {
	anita1->cwst_RXs[fill_index].phi_sector = fill_index_phi_sector;
	anita1->cwst_RXs[fill_index].layer = fill_index_layer;
	
	
	unsigned physical_phi_index = fill_index_phi_sector;
	unsigned physical_layer_index;
	
	
	// If the phi sector is an odd one, start with physical layer 1. If even, start with layer 0
	if (fill_index_phi_sector%2) {
	  physical_layer_index = fill_index_layer + 1;
	} else
	  if (fill_index_layer == 0) {
	    physical_layer_index = 0;
	  } else {
	    physical_layer_index = fill_index_layer + 1;
	  }
	
	if (physical_layer_index == 0) {
	  physical_phi_index = unsigned(fill_index_phi_sector / 2);
	}
	
	if (physical_layer_index == 1) {
	  physical_phi_index = unsigned(fill_index_phi_sector / 2);
	}
	
	//	Set antenna positions
	anita1->cwst_RXs[fill_index].x = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][0];
	anita1->cwst_RXs[fill_index].y = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][1];
	anita1->cwst_RXs[fill_index].z = anita1->ANTENNA_POSITION_START[physical_layer_index][physical_phi_index][2];
		  
	//	Fill the waveforms
	anita1->cwst_RXs[fill_index].waveform->assign(volts_rx_rfcm_trigger[fill_index_phi_sector][fill_index_layer].begin(),volts_rx_rfcm_trigger[fill_index_phi_sector][fill_index_layer].end());
	
	for (unsigned fill_index_timestep = 0; fill_index_timestep < anita1->HALFNFOUR; ++fill_index_timestep) {
	  anita1->cwst_RXs[fill_index].digitized->at(fill_index_timestep) = three_bit_round(anita1->cwst_RXs[fill_index].waveform->at(fill_index_timestep)/ anita1->rms_rfcm_e_single_event, false);
	}
	++fill_index;
      }
    }
    
    /*
      For whatever reason ANTENNA_POSITION_START is a c-array, with large dimensions (5 by 400).
      Most of these are initialized as a Vector with (0.,0.,1) for (x,y,z). These are non-existant antennas.
      
      Make sure that the antennas positions being copied are actually from a real antenna, and pay attention to the physical
      antenna layers 0 & 1, as they map to the 0th trigger layer.
    */
    
    //	Calculate the theta hypothesis angles
    vector <double> theta_hypotheses_deg;
    double theta_deg = SummedPowerThetaMin;
    
    do {
      theta_hypotheses_deg.push_back (theta_deg);
      theta_deg += SummedPowerThetaStep;
    } while (theta_deg <= SummedPowerThetaMax);
    
    //	Now we begin the trigger.
    
    //	First steps over windows:
    for (unsigned index_window = 0; index_window < 31; ++index_window) {
      double max_power_each_sector[16] = {0.};
      double summed_power_each_sector[16] = {0.};
      
      for (unsigned index_phi_sector = 0; index_phi_sector < 16; ++index_phi_sector) {
	double max_power = 0.;
	
	//	Use the current and opposite phi sectors to find an even horizontal line.
	TVector3 theta_zero (anita1->cwst_RXs[index_phi_sector * 3 + 1].x - anita1->cwst_RXs[(index_phi_sector + 8)%16 * 3 + 1].x, anita1->cwst_RXs[index_phi_sector * 3 + 1].y - anita1->cwst_RXs[(index_phi_sector + 8)%16 * 3 + 1].y, anita1->cwst_RXs[index_phi_sector * 3 + 1].z - anita1->cwst_RXs[(index_phi_sector + 8)%16 * 3 + 1].z);
	//	Normalizes that line.
	theta_zero = theta_zero.Unit();
	
	//	Loop over theta hypotheses
	for (unsigned index_theta_hypothesis = 0; index_theta_hypothesis < theta_hypotheses_deg.size(); ++index_theta_hypothesis) {
	  double power = 0.;
	  
	  //	If theta_zero.Phi() is NaN, then print out the problem.
	  /*
	    if (theta_zero.Phi() != theta_zero.Phi()) {
	    std::cout << theta_zero << std::endl;
	    std::cout << anita1->cwst_RXs[index_phi_sector * 3].x << "\t" << anita1->cwst_RXs[index_phi_sector * 3].y << "\t" << anita1->cwst_RXs[index_phi_sector * 3].z << std::endl;
	    std::cout << anita1->cwst_RXs[index_phi_sector * 3 + 1].x << "\t" << anita1->cwst_RXs[index_phi_sector * 3 + 1].y << "\t" << anita1->cwst_RXs[index_phi_sector * 3 + 1].z << std::endl;
	    std::cout << anita1->cwst_RXs[index_phi_sector * 3 + 2].x << "\t" << anita1->cwst_RXs[index_phi_sector * 3 + 2].y << "\t" << anita1->cwst_RXs[index_phi_sector * 3 + 3].z << std::endl;
	    }
	  */
	  
	  //	Calculate the vector for the hypothesis' direction.
	  double hypoth_theta_rad = (theta_hypotheses_deg[index_theta_hypothesis] + 45.) * RADDEG;
	  double hypoth_phi_rad = theta_zero.Phi() * RADDEG;
	  TVector3 hypoth_vector (sin (hypoth_theta_rad) * cos (hypoth_phi_rad), sin (hypoth_theta_rad) * sin (hypoth_phi_rad), cos (hypoth_theta_rad));
	  
	  //	There will be three antenna time offsets per phi sector:
	  double dist[3] = {0.};
	  double offset[3] = {0.};
	  unsigned offset_steps[3] = {0};
	  TVector3 middle_rx (anita1->cwst_RXs[index_phi_sector * 3 + 1].x, anita1->cwst_RXs[index_phi_sector * 3 + 1].y, anita1->cwst_RXs[index_phi_sector * 3 + 1].z);
	  
	  
	  //	Calculate the relative path length distances
	  dist[0] = hypoth_vector * (TVector3(anita1->cwst_RXs[index_phi_sector * 3].x, anita1->cwst_RXs[index_phi_sector * 3].y, anita1->cwst_RXs[index_phi_sector * 3].z) - middle_rx);
	  dist[1] = 0.;	// Just let this one be defined as "zero", the other ones are relative.
	  dist[2] = hypoth_vector * (TVector3(anita1->cwst_RXs[index_phi_sector * 3 + 2].x, anita1->cwst_RXs[index_phi_sector * 3 + 2].y, anita1->cwst_RXs[index_phi_sector * 3 + 2].z) - middle_rx);
	  
	  
	  //	Find time offsets from the relative path lengths
	  for (unsigned index_antenna = 0; index_antenna < 3; ++index_antenna) {
	    offset[index_antenna] = dist[index_antenna] / CLIGHT;
	  }
	  
	  //	Make the lowest time offset be the minoffset
	  double minoffset = offset[0];
	  if (offset[1] <= minoffset) {minoffset = offset[1];}
	  if (offset[2] <= minoffset) {minoffset = offset[2];}
	  
	  
	  //	Subtract the minoffset from the others, "normalizing" them so that the lowest value is zero
	  for (unsigned index_antenna = 0; index_antenna < 3; ++index_antenna) {
	    offset[index_antenna] -= minoffset;
	    offset_steps[index_antenna] = int(Tools::round(offset[index_antenna] / anita1->TIMESTEP));
	  }
	  
	  //	Check to make sure that all of the values will be in bounds
	  if (offset_steps[0] + index_window * 16 + 32 > 511) {continue;}
	  if (offset_steps[1] + index_window * 16 + 32 > 511) {continue;}
	  if (offset_steps[2] + index_window * 16 + 32 > 511) {continue;}
	  
	  //	Calculate the powers for every time in the window and then add it to the window power
	  for (unsigned time_index = index_window * 16; time_index < index_window * 16 + 32; ++time_index) {
	    double thispower = (anita1->cwst_RXs[index_phi_sector * 3].digitized->at(time_index + offset_steps[0]) + anita1->cwst_RXs[index_phi_sector * 3 + 1].digitized->at(time_index + offset_steps[1]) + anita1->cwst_RXs[index_phi_sector * 3 + 2].digitized->at(time_index + offset_steps[2]));
	    thispower *= thispower;
	    power += thispower;
	  }
	  
	  //	Make sure that the power is sane for the windows we care about
	  if (power == 0 && index_window < 21) {
	    std::cout << "Trigger error: a trigger window which should not have had zero power had zero power. Window was\t" << index_window << "\n";
	  }
	  
	  //	If this is the highest-power hypothesis for this phi sector, save it to max_power
	  if (power > max_power) {max_power = power;}
	}
	
	max_power_each_sector[index_phi_sector] = max_power;
      }
      
      //	Calculate the power for the neighboring phi sectors
      summed_power_each_sector[0] = max_power_each_sector[15] + max_power_each_sector[0] + max_power_each_sector[1];
      for (unsigned index_phi_sector = 1; index_phi_sector < 15; ++index_phi_sector) {
	summed_power_each_sector[index_phi_sector] = max_power_each_sector[index_phi_sector - 1] + max_power_each_sector[index_phi_sector] + max_power_each_sector[index_phi_sector + 1];
      }
      summed_power_each_sector[15] = max_power_each_sector[14] + max_power_each_sector[15] + max_power_each_sector[1];
      
      
      bool temporary_passing_variable = false;
      for (unsigned index_phi_sector = 0; index_phi_sector < 16; ++index_phi_sector) {
	//std::cout << "Max Power, sector\t" << index_phi_sector << "\t:\t" << summed_power_each_sector[index_phi_sector] << std::endl;
	if (summed_power_each_sector[index_phi_sector] >= SummedPowerThreshold) {
	  std::cout << "Passed with SummedPower =\t" << summed_power_each_sector[index_phi_sector] << "\ton sector\t" << index_phi_sector << "\n";
	  //return 1;
	  temporary_passing_variable = true;
	}
      }
      if (temporary_passing_variable) {return 1;}
    }
  }
  else if (settings1->TRIGGERSCHEME == 5) {
    
    
    
    
    double threshold=this_threshold;
    
    // need to find how many in a set of 6 pass
    
    int maxsample=TMath::MaxElement(5,anita1->imaxbin);
    int minsample=TMath::MinElement(5,anita1->iminbin);
    
    int nstayhigh=(int)(anita1->l1window/anita1->TIMESTEP);
    // now make each flag stay high for the required amount of time
    

    minsample=(int)(anita1->maxt_diode/anita1->TIMESTEP)+(anita1->NFOUR/4-(int)(anita1->maxt_diode/anita1->TIMESTEP))+(int)(anita1->arrival_times[anita1->rx_minarrivaltime]/anita1->TIMESTEP);
    maxsample=anita1->NFOUR/2-(int)(anita1->arrival_times[anita1->rx_minarrivaltime]/anita1->TIMESTEP);

    for (int i=0;i<5;i++) {
      anita1->iminbin[i]=minsample;
      anita1->imaxbin[i]=maxsample;
    }

        for (unsigned center_phi_sector_index = 0; center_phi_sector_index < 16; ++center_phi_sector_index) {
    // for (unsigned center_phi_sector_index = 0; center_phi_sector_index < 16; center_phi_sector_index+=2) {
      

      for (unsigned int index_hyp=0;index_hyp<anita1->vdifferent_offsets.size();index_hyp++) {


      for (unsigned i_layer = 0; i_layer < anita1->N_SUMMED_LAYERS; ++i_layer) {
	for (unsigned i_sector = 0; i_sector < anita1->N_SUMMED_PHI_SECTORS; ++i_sector) {

	  int rx=anita1->GetRxTriggerNumbering(i_layer,(center_phi_sector_index+i_sector)%16);
	  
	  
	  double timedomain_output_1_corrected[Anita::NFOUR/2]; // these are corrected for their delays so that they should line up in time
	  double timedomain_output_2_corrected[Anita::NFOUR/2];
	  // if we're doing an anita 3 trigger, adjust for delays in the diode outputs so that we can do a time coincidence trigger	
	  
	  for (int i=0;i<Anita::NFOUR/2;i++) {
	    timedomain_output_1_corrected[i]=anita1->timedomain_output_1_allantennas[rx][i];
	    timedomain_output_2_corrected[i]=anita1->timedomain_output_2_allantennas[rx][i];
	  }
      
      //      Tools::ShiftLeft(timedomain_output_1_corrected,anita1->NFOUR/2,anita1->arrival_times[rx]);
      //Tools::ShiftLeft(timedomain_output_2_corrected,anita1->NFOUR/2,anita1->arrival_times[rx]);

	  Tools::ShiftLeft(timedomain_output_1_corrected,anita1->NFOUR/2,anita1->vdifferent_offsets[index_hyp][anita1->N_SUMMED_PHI_SECTORS*i_layer+i_sector]);
	  Tools::ShiftLeft(timedomain_output_2_corrected,anita1->NFOUR/2,anita1->vdifferent_offsets[index_hyp][anita1->N_SUMMED_PHI_SECTORS*i_layer+i_sector]);
      

	  //for (int k=0;k<5;k++) {
	  //anita1->ston[k]=anita1->peak_v_banding_rfcm_e[k]/anita1->bwslice_vrms[k];
	  //} // end loop over bands
	  
	  if (rx==anita1->rx_minarrivaltime) {

	for (int i=0;i<Anita::NFOUR/2;i++) {
	  anita1->timedomain_output_1_corrected_forplotting[0][i]=timedomain_output_1_corrected[i];
	  anita1->timedomain_output_2_corrected_forplotting[0][i]=timedomain_output_2_corrected[i];
	}
	  }

      flag_e_L1[rx].clear();
      flag_h_L1[rx].clear();		
      

      for (int i=minsample;i<maxsample;i++) {	      
	      //      for (int i=minsample;i<maxsample;i++) {
	
	//cout << "timedomain_output_1_corrected[i], threshold, bwslice_rmsdiode[4] are " << timedomain_output_1_corrected[i] << "\t" << threshold*anita1->bwslice_rmsdiode[4] << "\n";
	if (timedomain_output_1_corrected[i]<threshold*anita1->bwslice_rmsdiode[4]) {
	  flag_e_L1[rx].push_back(1);
	  //cout << "passes.\n";
	}
	else {
	  flag_e_L1[rx].push_back(0);
	  //cout << "doesn't pass.\n";
	}
	
	if (timedomain_output_2_corrected[i]<threshold*anita1->bwslice_rmsdiode[4])
	  flag_h_L1[rx].push_back(1);
	else
	  flag_h_L1[rx].push_back(0);
	
      } // end loop over samples in window where we look for single channel trigger firing

      for (int i=minsample;i<maxsample;i++) {	      
		    //	for (int i=minsample;i<maxsample;i++) { // loop over samples in window where we looked for a single channel trigger
	  
	  if (flag_e_L1[rx][i-minsample]==1) // if the flag is high (remember the second index of flag_e counts from the start of the window)
	    for(int k=nstayhigh-1; k>0 && i<maxsample-1; k--) { // then for nstayhigh-1 samples after than we keep it high
	      // i<maxsample-1 makes sure that the second index of flag_e is always less than maxsample-minsample-1.
	      i++;
	      flag_e_L1[rx][i-minsample]=1;
	      
	    }
	  
	  //      if (flag_e[j][i]==1) {
	  // 	for (int k=i;k<i+nstayhigh;k++) {
	  // 	  if (k<NSAMPLES)
	  // 	    flag_e[j][k]=1;
	  // 	} // end loop over samples where we want it to stay high
	  
	  //       } // end if flag is high
	}
	
	for (int i=minsample;i<maxsample;i++) {
	  if (flag_h_L1[rx][i-minsample]==1)
	    for(int k=nstayhigh-1; k>0 && i<maxsample-1; k--) {
	      i++;
	      flag_h_L1[rx][i-minsample]=1;
	    }
	  
	  //      if (flag_h[j][i]==1) {
	  // 	for (int k=i;k<i+nstayhigh;k++) {
	  // 	  if (k<NSAMPLES)
	  // 	    flag_h[j][k]=1;
	  // 	} // end loop over samples where we want it to stay high
	  
	  //       } // end if flag is high
	  
	} // end loop over samples
	
	} // end loop over phi sectors being considered for this L1
      } // end loop over layers being considered for this L1

	
    // now find the sample with the highest number of bands that pass
  // int maxbands=0;
  // int nbands_pass=0;
  
    
    for (int i=minsample;i<maxsample;i++) { 
      
      int nbands_pass[2]={0};
      int maxbands[2]={0};
      
      for (unsigned i_layer = 0; i_layer < anita1->N_SUMMED_LAYERS; ++i_layer) {
	for (unsigned i_sector = 0; i_sector < anita1->N_SUMMED_PHI_SECTORS; ++i_sector) {
	  
	  int rx=anita1->GetRxTriggerNumbering(i_layer,(center_phi_sector_index+i_sector)%16);
	  
	  
	  
	  if (flag_e_L1[rx][i-minsample]==1)
	    nbands_pass[0]++;
	  if (flag_h_L1[rx][i-minsample]==1)
	    nbands_pass[1]++;
	  
	  
	}
      }
      
      if (nbands_pass[0]>maxbands[0])	
	maxbands[0]=nbands_pass[0];
      if (nbands_pass[1]>maxbands[1])	
	maxbands[1]=nbands_pass[1];
      
      if (maxbands[0]>=anita1->NCH_PASS || maxbands[1]>=anita1->NCH_PASS) { 

	return 1;	
      }
      

      
    } // end loop over layers
    
      }
      
	} // end loop over center phi sectors    

    
	
	//    }
	return 0;
  }
  return thispasses;
}//PassesTrigger



//!	Virtual nadir antennas had effective L2 triggers based on neighbor antennas
/*!
 *	the nadir layer only has 8 antennas but there are 16 slots in
 *  the "trigger layer"
 *  ant[][] is 1 if the l1 trigger is fired for each antenna
 *  fill the slots in  between the actual antennas with the "or"
 *  of the neighboring antennas
 *
 *	\todo	Deprecate this function in favor of PayloadArray or boost::multi_array objects which provide the
 *			ability to have multiple indexing schemes for a given collection of objects.
 *			For this case specifically PayloadArray was designed to have iterators and accessors which
 *			"OR" the values of the neighboring antennas, and are optionally read-only.
 */
void GlobalTrigger::FillInNadir(Settings *settings1,Anita *anita1,int ant) { //overloaded function
    int antarray[Anita::NPHI_MAX]={0};
    //here the ant array is an array of binary numbers
    int whichphi;
    
    for (int iphi=0;iphi<anita1->NRX_PHI[2];iphi++) {
		
		whichphi = 1 << GetPhiSector(settings1,2,iphi); // get trigger phi sector
		if (whichphi & ant) { // if this phi sector is masked
			antarray[GetPhiSector(settings1,2,iphi)]=1;
			
		}
    }
    
    
    FillInNadir(anita1,antarray);
    
    
    ant=0;
    for (int iphi=0;iphi<anita1->NTRIGPHISECTORS;iphi++) {
		if (antarray[iphi])
			ant+=(1<<iphi);
    }
    
    
    
}



//!	Virtual nadir antennas had effective L2 triggers based on neighbor antennas
/*!
 *	the nadir layer only has 8 antennas but there are 16 slots in
 *  the "trigger layer"
 *  ant[][] is 1 if the l1 trigger is fired for each antenna
 *  fill the slots in  between the actual antennas with the "or"
 *  of the neighboring antennas
 *
 *	\todo	Deprecate this function in favor of PayloadArray or boost::multi_array objects which provide the
 *			ability to have multiple indexing schemes for a given collection of objects.
 *			For this case specifically PayloadArray was designed to have iterators and accessors which
 *			"OR" the values of the neighboring antennas, and are optionally read-only.
 */
void GlobalTrigger::FillInNadir(Anita *anita1,int *ant) {
        int ileft,iright;
    //  cout << "ntrigphisectors is " << ntrigphisectors << "\n";
    //   for (int i=0;i<ntrigphisectors;i++) {
    //   cout << "before, ant is " << ant[2][i] << "\n";
    //   }
    
    for (int i=0;i<anita1->NTRIGPHISECTORS/2;i++) {
		ileft=(2*i+anita1->NTRIGPHISECTORS-1)%anita1->NTRIGPHISECTORS;
		iright=(2*i+1)%anita1->NTRIGPHISECTORS;
		if (ant[ileft] || ant[iright])
			ant[2*i]=1;
    }
    
    //    for (int i=0;i<anita1->NTRIGPHISECTORS;i++) {
    //   cout << "after, ant is " << ant[2][i] << "\n";
    //   }
    
    
}



//!	Provides a mapping between the 4 layers and 16 phi sectors physically to the 3 layers and 16 sectors logically
/*!
 *	Because the top two layers of antennas had only 8 antennas each and were rotated perfectly out of phase of one
 *	another, such a conversion was needed to keep the trigger sane.
 *
 *	\todo	Deprecate this function in favor of PayloadArray or boost::multi_array objects which provide the
 *			ability to have multiple indexing schemes for a given collection of objects.
 */
int GlobalTrigger::GetPhiSector(Settings *settings1,int i,int j) { // given trigger layer and index, get phi sector.
    // for the upper two layers, the phi sector is just j
    // for the nadir layer, the phi sector is 2*j+1
    // warning this counts from 0
    if (i<2)
		return j;
    if (i==2) {
		if (settings1->DISCONES==2)
			return 2*j+1;
		else
			return j;
    }
    else
    {
		cout << "Input non-existent layer!\n";
		return -1;
		
    }
    return -1;
}


//!	Provides a mapping between the 4 layers and 16 phi sectors physically to the 3 layers and 16 sectors logically
/*!
 *	Because every payload had a different number of antennas, this provides a relation from
 *	where the antennas were located physically and how they were considered to by positioned
 *	logically by the triggering system.
 *
 *	\todo	Deprecate this function in favor of PayloadArray or boost::multi_array objects which provide the
 *			ability to have multiple indexing schemes for a given collection of objects.
 */
void GlobalTrigger::GetAnitaLayerPhiSector(Settings *settings1,int i,int j,int &whichlayer,int &whichphisector) {
    
    
    if (settings1->WHICH==6 || settings1->WHICH==2 || settings1->WHICH==8) {// If Anita 1 or Anita 2
		if (i==0) {
			whichlayer=0;
			whichphisector=2*j+1;
		}
		else if (i==1) {
			whichlayer=0;
			whichphisector=2*j;
		}
		else if (i==2) {
			whichlayer=1;
			whichphisector=j;
		}
		else if (i==3) {
			whichlayer=2;
			whichphisector=2*j+1;
		}
    } // end anita 1 or anita 2
    else if (settings1->WHICH==9) { // anita 3
		if (i==0) {
			whichlayer=0;
			whichphisector=2*j+1;
		}
		else if (i==1) {
			whichlayer=0;
			whichphisector=2*j;
		}
		else if (i==2) {
			whichlayer=1;
			whichphisector=j;
		}
		else if (i==3) {
			whichlayer=2;
			whichphisector=j;
		}
    }  // end anita 3
    else {
		whichlayer=i;
		whichphisector=j;
		
    }
    
}


//!	Level 3 Trigger Function
/*!
 *	Given a trigger scheme (represented via the Settings::WHICH value), the results of the
 *	previous trigger (the L2 trigger) are compared against one of the various combinatoric
 *	triggers.
 *
 *	For Anita 2, the trigger required two layers out of the three L2 triggers in a phi sector,
 *	but for sectors which did not have a nadir antennas, these "phantom nadirs" were a boolean
 *	"OR" of the two nadir antennas closest to that phi sector. Finally, two neighboring phi
 *	sectors must have accomplished this in order to produce a triggered L3.
 *
 *	\todo	The organization of this function needs to be changed so that a single variables
 *			turns on (or off) a single component of the trigger/simulation. The organization
 *			as it stands makes it difficult to know what features a particular trigger should
 *			be exhibiting.
 *			
 *			Along this topic, the values for settings should be either strings, values mapped
 *			to strings (such as boost::variable_map), or enumerations with a proper descriptive
 *			declared name. That way when settings1->WHICH would have a value of 9, it would be
 *			written settings1->WHICH==ANITA3.
 *
 *			Lastly, the boolean functions which are chained together to represent the full state
 *			machine of the trigger electronics are an impossible-to-find bug waiting to happen.
 *			"Short circuit evaluation" and less-than-obvious operator fixity rules make the
 *			given representation one of the most difficult to manage.
 *
 *			The use of standard library functions "all_of", "any_of", "none_of", "count", and
 *			"count_if" on an intelligent container like boost::multi_array or PayloadArray would
 *			make the expression of the trigger requirements much more clear and less prone to
 *			error.
 *			
 *			Finally, all modern compilers will automatically optimize power-of-two multiplication
 *			of integers into a bit shift, so nothing is lost by switching from explicit shifts to
 *			multiplication, and readability is gained.
 */
int GlobalTrigger::L3Trigger(Settings *settings1,Anita *anita1,int loctrig[Anita::NLAYERS_MAX][Anita::NPHI_MAX],int loctrig_nadironly[Anita::NPHI_MAX],int discones_passing,int &l3trig) {
    int thispasses=0;
    int whichphipass[Anita::NPHI_MAX]={0};
    
    for (int i=0;i<Anita::NTRIG;i++)
		triggerbits[i]=0;
    
    for (int i=0;i<anita1->PHITRIG[0];i++) {
		
		
		if (settings1->WHICH==9) { // anita 3
			if (loctrig[0][i]>0 && loctrig[1][i]>0 && loctrig[2][i]>0) {
				thispasses=1;
				whichphipass[i]=1;
				triggerbits[3]+=(1<<i);
			}
		}
		else { // anita 1 or anita 2
			
			for (int ilayer=0;ilayer<anita1->NTRIGGERLAYERS-1;ilayer++) {
				
				if (loctrig[ilayer][i]>0 && loctrig[ilayer+1][i]>0) { // upper two layers
					
					//l3trig += (1<<i);
					thispasses=1;
					whichphipass[i]=1;
					
					if(loctrig[0][i]>0 && loctrig[1][i]>0)
						triggerbits[0] += (1<<i);
					
					
					
				}
			}
			// anita 1, anita 1 kurt, anita 2 kurt, anita 3
			if (settings1->WHICH==2 || settings1->WHICH==6 || settings1->WHICH==8 || settings1->WHICH==9) { // Anita 1, 2 or 3
				
				
				if (((settings1->DISCONES==2) && (loctrig[1][i]>0 && loctrig[2][i]>0)) || // top and nadir
					((settings1->DISCONES==2) && (loctrig[0][i]>0 && loctrig[2][i]>0)) || // second layer and nadir
					((settings1->DISCONES==2) && loctrig_nadironly[i]>0) || // just nadir
					(settings1->DISCONES==1 && (loctrig[1][i]>0 && discones_passing>=settings1->NDISCONES_PASS)) || // top and discones
					(settings1->DISCONES==1 && (loctrig[0][i]>0 && discones_passing>=settings1->NDISCONES_PASS))   // second layer and discones
					) {
					
					thispasses=1;
					whichphipass[i]=1;
					// this is just for nadir studies->
					//how many of each trigger condition
					
					if(loctrig[1][i]>0 && loctrig[2][i]>0) {
						triggerbits[1] += (1<<i);
						
					}
					if(loctrig[0][i]>0 && loctrig[2][i]>0) {
						triggerbits[2]+=(1<<i);
						
					}
					// 	if (loctrig_nadironly[i]>0)
					// 	  triggerbits[3]+=weight2;
					// 	if ((loctrig_nadironly[i]>0 ||
					// 	     (loctrig[1][i]>0 && loctrig[2][i]>0) ||
					// 	     (loctrig[0][i]>0 && loctrig[2][i]>0)) &&
					// 	    !(loctrig[0][i]>0 && loctrig[1][i]>0))
					// 	  triggerbits[4]+=weight2;
					
					
					thispasses=1;
					
				}
				
			} //if it's the anita 1 payload
			
			
		} // if it's not anita 3
		
		
		
		if (whichphipass[i])
		  l3trig+=(1<<i);
		
    } // end looping over phi
    
    return thispasses;
}


//!	Calculate the difference in each arrival time due to propagation delay of the signal
/*!
 *	Method void GetArrivalTimes:
 *  input:
    direction of rf propagation when it arrives at balloon (earth frame),
    balloon location (earth frame)
    rotation of balloon around axis of symmetry
    
    output:
    array of doubles giving time at which rf reaches antenna feedpoint, where the first antenna is hit at t=0.
    
    
    This method takes as input the rf direction, the location of the balloon, the phi spin of the balloon, and
    an array of vectors initialized to point to the location of the antennas, assuming the balloon is upright at the
    south pole, with no spin.
    
    The method transforms the antenna location vectors to point from the current center of the balloon (actually the center
    of the top ring of antennas) to the antennas, assuming that the balloon is vertical over the surface of the ice,
    and taking into account a possible phi spin.
    
    The method projects the antenna direction vector onto the rf direction to find the distance between the center
    of the balloon and each antenna, along the direction of rf propagation.
    
    Finally, arrival times are shifted so that the first antenna triggered is at time 0.
 */
// void GlobalTrigger::GetArrivalTimes(int inu, Anita *anita1, const Vector& rf_direction) {
    
    
//     for (int antenna_index = 0; antenna_index < (anita1->number_all_antennas); antenna_index++) { //loop over layers on the payload
// 		arrival_times[antenna_index] = (anita1->antenna_positions[antenna_index] * rf_direction) / CLIGHT;
//     } // for: loop over antenna layers
    
//     double first_trigger_time = Tools::dMin(arrival_times,(anita1->number_all_antennas));
    
//     for (int i=0;i<(anita1->number_all_antennas);i++){
// 		arrival_times[i] -= first_trigger_time;
// 		if (arrival_times[i] == 0){
// 			first_phi_sector_hit = floor(i / 16);
// 		}
//     }
    
// } // GetArrivalTimes



//!	Shifts waveforms in time to simulate relative signal delays caused by propagation
/*!
 *	This function accepts nested vectors of waveforms along with nested vectors of delays
 *	and then iterates through nested for-loops whose bounds iterate with the surrounding scope.
 *
 *	The output is assembled and then truncated to prevent uneven array lengths post-shifting.
 *
 *
 *	\todo	Instead of accepting nested vectors, boost::multi_array/_views or
 *			PayloadArray objects should be passed in and out.
 *			
 *			Something to keep in mind is that RVO (Return Value Optimization) can
 *			make it actually faster to return a large object by value if it can
 *			be used to immediately construct another object.
 *			
 *			Also, instead of worrying about resizing the arrays it'd probably
 *			be better to store the time index of the last signal bin, and then
 *			simply pass that bin as the end iterator offset for the truncated array.
 *			
 *			The use of vector<>::push_back may cause performance issues but more
 *			importantly the use of clear() and element-by-element assignment causes
 *			more concern than they're worth, so they should be replaced with
 *			std::copy() for the temporaries which need it and just use functions from
 *			the standard library's #include <algorithm>.
 *
 *			Realistically though, the delays should just be used to advance the iterators.
 */
void GlobalTrigger::delay_align_antenna_waveforms(const vector< vector < vector <double> > >& waveforms, const vector < vector <unsigned int> >& delays, vector < vector <double> >& output){
    // The waveforms vector is accessed by waveforms[phi_sector_id][antenna_id][sample_id]
    // the delays are accessed by delays[phi_sector_id][antenna_id]
    
    // the output is in the same format as the input, but will only include samples which have all elements
    output.clear();
    
    unsigned int shortest_waveform = waveforms[0][0].size();
    vector <double> temp_antenna_waveform;
    for (unsigned int phi_sector_index = 0; phi_sector_index < 3; phi_sector_index++){
		for (unsigned int antenna_index = 0; antenna_index < waveforms[phi_sector_index].size(); antenna_index++){
			temp_antenna_waveform.clear();
			for (unsigned int sample_index = delays[phi_sector_index][antenna_index]; sample_index < waveforms[phi_sector_index][antenna_index].size(); sample_index++){
				temp_antenna_waveform.push_back(waveforms[phi_sector_index][antenna_index][sample_index]);
			}
			if (temp_antenna_waveform.size() < shortest_waveform){
				shortest_waveform = temp_antenna_waveform.size();
			}
			output.push_back(temp_antenna_waveform);
		}
    }
    // Must now make all have the same length
    for (unsigned int antenna_index = 0; antenna_index < output.size(); antenna_index++){
		output[antenna_index].resize(shortest_waveform);
    }
}


//!	Given a number of waveforms (usually 3) which are delay-aligned, it sums them.
/*!
 *	\todo	This should not take nested vectors as a parameter but rather
 *			a boost::multi_array or boost::multi_array_view, or an instance of
 *			the PayloadArray class, which is a multi_array container modified to
 *			have circular indices.
 *			Also, the bounds of the for-loop should not be separately declared
 *			as it just adds another point of failure.
 */
void GlobalTrigger::sum_aligned_waveforms(const vector < vector <double> >& waveforms, vector <double>& output){
    output.clear();
    unsigned waveform_size = waveforms[0].size();
    unsigned waveforms_size = waveforms.size();
    for (unsigned element_index = 0; element_index < waveform_size; ++element_index) {
		double element = 0.;
		for (unsigned antenna_index = 0; antenna_index < waveforms_size; ++antenna_index) {
			element += waveforms[antenna_index][element_index];
		}
		output.push_back(element);
    }
}




//!	Performs an element-wise squaring of the values of the input array
/*!
 *	\todo	This function should be deprecated because the same functionality
 *			is provided by the standard library:
 *
 *		std::transform(wfm.begin(), wfm.end(), result.begin(), [](auto& x){x *= x});
 *		
 *			or by using the boost library:
 *		
 *		result = boost::transform(wfm, [](auto& x){x *= x});
 *
 */
void GlobalTrigger::square_waveform_elements(const vector <double>& waveform, vector <double>& output){
    output.clear();
    
    for (unsigned int element_index = 0; element_index < waveform.size(); element_index++){
		output.push_back(waveform[element_index] * waveform[element_index]);
    }
}


//!	Sum a window from the specified starting index
/*!
 *	\todo	This function should be replaced as it is unnecessary -- the same
 *			thing can be done more legibly with std::accumulate(...) :
 *
 *		double sum = std::accumulate(wfm.begin(), wfm.end(), 0.); 
 *
 */
double GlobalTrigger::summed_power_window(const vector <double>& waveform, unsigned int start_index, unsigned int length){
    double result = 0;
    for (unsigned int index = start_index; index < start_index + length; index++){
		result += waveform[index];
    }
    return result;
}


//!	Three bit rounding function
/*!
 *	This function takes an input and digitizes it into one of 8 (2^3 == 8) intervals,
 *	with unit spacing from -3.5 to +3.5.
 *
 *	This function takes a number and will round it to the closest integer multiple
 *	plus 0.5, using randomness to get distribute exact integers.
 */
double GlobalTrigger::three_bit_round(double input, bool round_zero_up, bool allow_zero) {
    // 
    double result = 0.;
    if (input < -3.5){
		result = -3.5;
    }else if (input > 3.5){
		result = 3.5;
    }else if (floor(input) == input){
		if (!allow_zero) {
			if (round_zero_up) {
				result = input + 0.5;
			} else {
				result = input + -0.5;
			}
		} else {
			result = input;
		}
    }else if (floor(input) + 0.5 > input){
		result = ceil(input) - 0.5;
    }else if (floor(input) + 0.5 < input){
		result = floor(input) + 0.5;
    }
    
    return result;
}

//! Converts a waveform by sampling it into 3 bits, after it is normalized to the RMS
/*	
 *	This function takes a waveform array as a parameter and returns a vector <double> of 3 bit waveforms
 */
void GlobalTrigger::convert_wfm_to_3_bit(const vector <double>& wfm, double rms, vector <double>& output){
    output.clear();
    for (unsigned int index = 0; index < wfm.size(); index++){
		output.push_back(three_bit_round(wfm[index]/rms));
    }
}
