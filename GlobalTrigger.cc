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
#include "trigger.hh"
#include <cmath>
#include "Tools.h"
#include "Settings.h"
#include "screen.hh"
#include "GlobalTrigger.h"

using std::cout;

#include "GlobalTrigger.h"

using std::cout;


GlobalTrigger::GlobalTrigger(Settings *settings1,Anita *anita1){
    
  // 2=top,1=middle,0=bottom
  WHICHLAYERSLCPRCP[0]=0;
  WHICHLAYERSLCPRCP[1]=1;
  WHICHLAYERSLCPRCP[2]=0;


  // time step between sampling tunnel diode output for the trigger
  TRIGTIMESTEP=2.E-9; 
      
   
  L3_COINCIDENCE=22.5e-9;

  // L1 coincidence window, in seconds  
  L1_COINCIDENCE_ANITA3[0]=16.E-9; // B->M or T
  L1_COINCIDENCE_ANITA3[1]=22.E-9; // M->B or T
  L1_COINCIDENCE_ANITA3[2]=5.E-9; // T->B or M 

  // in this scenario B->M is the same as M->B for example
  // this needs to be generalized- using this same thing for all scenarios which isn't right
  LASTTIMETOTESTL1_ANITA3=((double)anita1->NFOUR/2)*anita1->TIMESTEP-Tools::dMax(L1_COINCIDENCE_ANITA3,3); // can't test L1 after this point because the l1_coincidence windows go past the end of the waveform.




  // layers indexed in reverse order!
  // B->M
  // B->T
  L1_COINCIDENCE_LR_SCA[0]=16.E-9;
  L1_COINCIDENCE_LR_SCA[1]=16.E-9;

  LASTTIMETOTESTL1_ANITA4LR_SCA=Tools::dMax(L1_COINCIDENCE_LR_SCA,2);
  LASTTIMETOTESTL1_ANITA4LR_SCA=((double)anita1->NFOUR/2)*anita1->TIMESTEP-LASTTIMETOTESTL1_ANITA4LR_SCA; // can't test L1 after this point because the l1_coincidence windows go past the end of the waveform.

  L1_COINCIDENCE_MOREGENERAL[0][0]=16.E-9;
  L1_COINCIDENCE_MOREGENERAL[0][1]=4.E-9;

  L1_COINCIDENCE_MOREGENERAL[1][0]=16.E-9;
  L1_COINCIDENCE_MOREGENERAL[1][1]=4.E-9;
   
  L1_COINCIDENCE_MOREGENERAL[2][0]=4.E-9; // L1 coincidence window, in seconds  ;
  L1_COINCIDENCE_MOREGENERAL[2][1]=4.E-9; // L1 coincidence window, in seconds  ;

  LASTTIMETOTESTL1_ANITA4=0.;
  for (int i=0;i<3;i++) {
    for (int j=0;j<2;j++) {
      if (L1_COINCIDENCE_MOREGENERAL[i][j]>LASTTIMETOTESTL1_ANITA4)
	LASTTIMETOTESTL1_ANITA4=L1_COINCIDENCE_MOREGENERAL[i][j];
    }
  }
  
  LASTTIMETOTESTL1_ANITA4=((double)anita1->NFOUR/2)*anita1->TIMESTEP-LASTTIMETOTESTL1_ANITA4;
  // can't test L1 after this point because the l1_coincidence windows go past the end of the waveform.

  nstepback=(int)(2.E-9/TRIGTIMESTEP);
    
  // coincidence between lcp,rcp within an antenna
  L1_COINCIDENCE_ANITA4LR_SCB=4.e-9;
  // B->T, B->M, M->T
  L2_COINCIDENCE_ANITA4LR_SCB[0]=12.e-9;
  L2_COINCIDENCE_ANITA4LR_SCB[1]=4.e-9;
  L2_COINCIDENCE_ANITA4LR_SCB[2]=8.e-9;
  LASTTIMETOTESTL1_ANITA4LR_SCB=((double)anita1->NFOUR/2)*anita1->TIMESTEP-L1_COINCIDENCE_ANITA4LR_SCB; // can't test L1 after this point because the l1_coincidence windows go past the end of the waveform.;
  LASTTIMETOTESTL2_ANITA4LR_SCB=((double)anita1->NFOUR/2)*anita1->TIMESTEP-Tools::dMax(L2_COINCIDENCE_ANITA4LR_SCB,3); // can't test L1 after this point because the l1_coincidence windows go past the end of the waveform.;
  L3_COINCIDENCE_ANITA4LR_SCB=12.e-9;
				   
  DELAYS[0]=0.; // delay top
  DELAYS[1]=4.e-9; // delay middle    
  DELAYS[2]=4.e-9; // delay bottom   

  Tools::Zero(triggerbits,Anita::NTRIG);
    
  phiTrigMask[0]=anita1->phiTrigMask; // set the phi mask to the input value which comes from the balloon class
  phiTrigMask[1]=anita1->phiTrigMaskH; // set the phi mask to the input value which comes from the balloon class
  l1TrigMask[0]=anita1->l1TrigMask; // set the phi mask to the input value which comes from the balloon class
  l1TrigMask[1]=anita1->l1TrigMaskH; // set the phi mask to the input value which comes from the balloon class    
  


  
  for (int i=0;i<Anita::NLAYERS_MAX;i++) {
    for (int j=0;j<Anita::NPHI_MAX;j++) {
      for (int k=0;k<2;k++) {
	for (int p=0;p<anita1->NBANDS+1;p++) {
	  channels_passing[i][j][k][p]=0;
	  // make vchannels_passing the proper length.
	  vchannels_passing[i][j][k].push_back(0);
	}


      }
    }
  }

  for (int i=0;i<3;i++) {
    for (int j=0;j<16;j++) {
      for (int k=0;k<2;k++) {
	for (int p=0;p<5;p++) {
	  arrayofhits[i][j][k][p].clear();
	}
      }
    }
  }    
    
  for (int k=0;k<2;k++) {		
    for (int i=0;i<Anita::NLAYERS_MAX;i++) {
      for (int j=0;j<Anita::NPHI_MAX;j++) {
	volts[k][i][j]=0.;
	volts_em[k][i][j]=0.;
	volts_original[k][i][j]=0.;
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

void GlobalTrigger::PassesTrigger(Settings *settings1,Anita *anita1,int discones_passing,int mode,int *l3trig,int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int antennaclump,int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX],int inu,
				  int *thispasses) {

  double this_threshold= anita1->powerthreshold[4]; //-4.34495;
  return PassesTrigger(settings1,anita1,discones_passing,mode,l3trig,l2trig,l1trig,antennaclump,loctrig,loctrig_nadironly,inu,this_threshold, thispasses);   
}




void GlobalTrigger::PassesTrigger(Settings *settings1,Anita *anita1,int discones_passing,int mode,int *l3trig,int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int antennaclump,int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX],int inu,double this_threshold,
				  int *thispasses) {


  //bool ishpol should only be used for anita3, by default do only vpol

  if (settings1->TRIGGERSCHEME < 3) {

    // Basic trigger refers to  frequency domain voltage (0) frequency domain energy (1) timedomain diode integration (2) 

    PassesTriggerBasic(settings1, anita1, discones_passing, mode, l3trig, l2trig, l1trig, antennaclump, loctrig, loctrig_nadironly, thispasses);

  }

  else if (settings1->TRIGGERSCHEME == 3) {

    PassesTriggerCoherentSum(settings1, anita1, inu, thispasses);
     
  }
  else if (settings1->TRIGGERSCHEME == 4) {

    PassesTriggerSummedPower(settings1, anita1);

  }
  else if (settings1->TRIGGERSCHEME == 5) {
    
    // Don't know the name of this one so I'll call it 5
    PassesTriggerScheme5(anita1, this_threshold, thispasses);
   
  }


}//PassesTrigger



void  GlobalTrigger::PassesTriggerBasic(Settings *settings1,Anita *anita1,int discones_passing,int mode,int *l3trig,int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int antennaclump,int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX], int *thispasses){

  int ltsum=0;
  int channsum=0;
  int ihit=0;
  //  int thispasses[2]={0,0};
  int required_bands_failed[2]={0,0}; // keep track of whether bands that were required to pass did not, for each polarization

    
  // this is an array with 1=pass and 0=fail for each channel
  // the first two layers on the payload are "compacted"
  // into one trigger layer of 16 antennas.
  // layer number, antenna, polarization, bandwidth slice
  int channels_compacted_passing[Anita::NLAYERS_MAX][Anita::NPHI_MAX][2][5] = {{{{0}}}};
    
  for (int k=0;k<Anita::NPOL;k++) {
    l3trig[k]=0;
    Tools::Zero(l1trig[k],Anita::NTRIGGERLAYERS_MAX);
    Tools::Zero(l2trig[k],Anita::NTRIGGERLAYERS_MAX); // for all three layers
    Tools::Zero(loctrig_nadironly[k],Anita::NPHI_MAX);
  }
  // which clumps of 3 antennas
  //within the nadir layer pass
  // 0th and 1st element= antenna at phi=0,
  // 2nd and 3rd element=next antenna, etc.
    
  for (int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
    for (int k=0;k<Anita::NPOL;k++) {
      // for (int j=0;j<Anita::NPHI_MAX;j++) {
      Tools::Zero(loctrig[k][ilayer],Anita::NPHI_MAX);
    }
    // which clumps of 3 antennas
    //within each trigger layer pass L2 trigger
    // layer, phi location
    // } //for
  } // end of initializing to zero

  int NBAND=5; // number of bandwidth slices
  int whichlayer=0;
  int whichphisector=0;
  for (int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
    for (int iphi=0;iphi<anita1->NRX_PHI[ilayer];iphi++) {
      for (int ipolar=0;ipolar<2;ipolar++) {

	if (anita1->pol_allowed[ipolar]==0) continue;// if this polarization is not allowed to contribute then don't do it
	//	  if (ishpol && ipolar==0) continue; // Anita3 : only do the polarisation required
	//else if ((!ishpol) && ipolar==1) continue; // Anita3 : only do the polarisation required
	  
	for (int iband=0;iband<NBAND;iband++) {
	    
	  GetAnitaLayerPhiSector(settings1,ilayer,iphi,whichlayer,whichphisector); // Translate Anita physical layer to Anita trigger layer and phi sector (4 layers with 8,8,16,8 phi sector to 3 layers with 16 phi sectors each.  In the nadir layer, 8 of the 16 slots are empty on the trigger layer.)
	    
	    
	    
	  // combining top two layers on the payload into one trigger layer
	  // this means the 0th antenna in the first trigger layer is
	  // physically higher than the 1st antenna in the first trigger layer
	  // which means that the nadirs are aligned with the antennas with indices 1,3,5 etc.
	  // we will still use indices 0-7 for them though
	  channels_compacted_passing[whichlayer][whichphisector][ipolar][iband]+=channels_passing[ilayer][iphi][ipolar][iband];
	} //for
      } //for
    } //for
  } //for

  int antsum[2]={0,0}; // counter for how many channels of each polarization on an antenna pass
  /* antenna triggers */
    
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
    
    
  int ant[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX] = {{{0}}}; // which antennas and which polarizations pass at L1
  // trigger layer, phi position
    
  int antpass[2]={0,0}; // count how many antennas pass
  int iphitrig=0;  // which trigger phi sector
  int ihittrig=0;
  // local level 1 trigger at the antenna
  for(int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++) { // loop over layers
    // this is nlayers-1 because NLAYERS counts top 16 antennas as 2 layers
      
    for(int iphi=0;iphi<anita1->PHITRIG[iloc];iphi++) { // loop over phi position
      iphitrig=GetPhiSector(settings1,iloc,iphi); // get trigger phi sector
      // counts from 0
	
      for(int ipolar=0;ipolar<2;ipolar++) {

	antsum[ipolar] = 0; // start sum for this antenna for each polarization
	
	for(int iband=0;iband<NBAND;iband++) {
	  if(channels_compacted_passing[iloc][iphitrig][ipolar][iband] == 1
	     && anita1->bwslice_allowed[iband]==1 && anita1->pol_allowed[ipolar]==1) { // only increment if it's one of the allowed bands and allowed polarizations.
	      
	    if (settings1->PHIMASKING && settings1->WHICH==9){ // only applying channel masking like this if it's Anita-3
	      //		if ((ipolar==0 && (1<<iphitrig & l1TrigMask[0])) || (ipolar==1 && (1<<iphitrig & l1TrigMask[1])) ){
	      if  (1<<iphitrig & l1TrigMask[ipolar])  {
		continue; // was this channel masked?
	      }
	    }
	    antsum[ipolar] = antsum[ipolar] +1; // sum channels that pass for this antenna, polarization

	  }
	} // loop over bands
      } // end loop over polarizations
	
	
      for (int ipolar=0;ipolar<2;ipolar++) {
	required_bands_failed[ipolar]=0;
	for (int iband=0;iband<NBAND;iband++) { // notice sum over 5 bands now
	  if(channels_compacted_passing[iloc][iphitrig][ipolar][iband] == 0
	     && anita1->bwslice_required[iband]==1 ) { // if this band was required to pass and it didn't,
	    required_bands_failed[ipolar] = 1; // fatal for this antenna, polarization
	  }
	    
	} // end loop over bands
      } // end loop over polarizations

	// if the required bands didn't pass then set antsum=0 so the antenna doesn't pass
	

      for (int ipolar=0;ipolar<2;ipolar++) {	  
	ant[ipolar][iloc][iphitrig] = 1; // start with every antenna passing L1.	  
	if( (anita1->pol_required[ipolar]==1 && 
	     (antsum[ipolar] < anita1->trigRequirements[0]	    
	      || required_bands_failed[ipolar]==1 ) )
	    || (anita1->pol_required[ipolar]==0)
	    )  { // if this polarization is required and it doesn't pass then make the antenna fail
	  ant[ipolar][iloc][iphitrig]=0;
	} //if
	else{
	  antpass[ipolar]+=1;// increment if it passed
	}
      }

    } //for (phi position)
  } //for (layers)
    
  if (settings1->DISCONES==2) {
    // fill the slots in  between the actual antennas with the "or"
    // of the neighboring antennas
    FillInNadir(anita1,ant[0][2]);
    FillInNadir(anita1,ant[1][2]);
  }
    
  for(int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++){ //This is NLAYERS-1 because NLAYERS counts top 16 antennas as 2 layers
    for (int iphi=0;iphi<anita1->PHITRIG[iloc];iphi++){
      iphitrig=GetPhiSector(settings1,iloc,iphi);
      for (int ipolar=0;ipolar<2;ipolar++) {	  
	if (ant[ipolar][iloc][iphitrig]==1) 
	  l1trig[ipolar][iloc] += (1<<iphitrig); // this keeps track of which antennas pass the l1 trigger
      }
    }//for (phi position)
  } //for (layers) // Hmm this is not used
    
  if (mode == 1) {// TRIGTYPE=2 -> just require ANTtrig channels pass on 2 antennas
    for (int ipolar=0;ipolar<Anita::NPOL;ipolar++) {
      if (antpass[ipolar] >= 2) 
	thispasses[ipolar] = 1;
    }
   
  } else if (mode == 2) { // TRIGTYPE=1 ->

    // ANITA-3
    if (settings1->WHICH==9) {
      // for each of these, need to set 
      // loctrig, l2trig

      std::array<std::array<std::vector<int>,16>,2> vl1trig;
	
      L1Anita3_AllPhiSectors(anita1,vl1trig);
      // just have to modify this function so it's looping over all relevant trigtimesteps
	
      for (int ipol=0;ipol<2;ipol++) {
	for (int iphi=0;iphi<16;iphi++) {

	  for (unsigned int ibin=0;ibin<vl1trig[ipol][iphi].size();ibin++) {
	    anita1->l1trig_anita3and4_inanita[ipol][iphi][ibin]=vl1trig[ipol][iphi][ibin];
	    // if (anita1->l1trig_anita3_inanita[ipol][iphi][ibin])
	    // 		cout << "This l1 is 1.\n";
	  }
	  for (int ibin=vl1trig[ipol][iphi].size();ibin<anita1->HALFNFOUR;ibin++) {
	    anita1->l1trig_anita3and4_inanita[ipol][iphi][ibin]=0;
	  }
	}
      }
      //	cout << "vl1trig size is " << vl1trig[0][0].size() << "\n";

      std::array<std::array<std::vector<int>,16>,2> vl2trig;
	
      L2Anita3and4(anita1,vl1trig,
		   vl2trig);

      if (settings1->PHIMASKING){ // only applying channel masking like this if it's Anita-3
	for (int ipol=0;ipol<2;ipol++) {
	  for (int iphi=0;iphi<16;iphi++) {
	    if  ((1<<iphi & l1TrigMask[ipol])||(1<<iphi & phiTrigMask[ipol]))  {
	      // set phi sector to 0
	      std::fill(vl2trig[ipol][iphi].begin(), vl2trig[ipol][iphi].end(), 0);
	    }
	  }
	}
      }
	
      int vl3trig[2][16];
      L3Anita3and4(anita1,vl2trig,
		   vl3trig,thispasses);

      for (int ipol=0;ipol<2;ipol++) {
	for (int iphi=0;iphi<16;iphi++) {
	  if (vl3trig[ipol][iphi]>0)    l3trig[ipol]+=(1<<iphi);
	}
      }
	
    }
    // ANITA-4
    else if (settings1->WHICH==10 && !(settings1->LCPRCP)) {

      std::array<std::array<std::vector<int>,16>,2> vl1trig;
	
      L1Anita4_AllPhiSectors(anita1,vl1trig);
      // just have to modify this function so it's looping over all relevant trigtimesteps
	
      for (int ipol=0;ipol<2;ipol++) {
	for (int iphi=0;iphi<16;iphi++) {

	  for (unsigned int ibin=0;ibin<vl1trig[ipol][iphi].size();ibin++) {
	    anita1->l1trig_anita3and4_inanita[ipol][iphi][ibin]=vl1trig[ipol][iphi][ibin];
	     
	    // 	      if (anita1->l1trig_anita3and4_inanita[ipol][iphi][ibin])
	    // 		cout << "This l1 is " << ipol << "\t" << iphi << "\t" << ibin << "\t" << anita1->l1trig_anita3and4_inanita[ipol][iphi][ibin] << "\n";
	  }
	  for (int ibin=vl1trig[ipol][iphi].size();ibin<anita1->HALFNFOUR;ibin++) {
	    anita1->l1trig_anita3and4_inanita[ipol][iphi][ibin]=0;
	  }
	}
      }
      //	cout << "vl1trig size is " << vl1trig[0][0].size() << "\n";

      std::array<std::array<std::vector<int>,16>,2> vl2trig;
	
      L2Anita3and4(anita1,vl1trig,
		   vl2trig);

      int vl3trig[2][16];
      L3Anita3and4(anita1,vl2trig,
		   vl3trig,thispasses);


    }
    else if (settings1->WHICH==10 && settings1->LCPRCP) {


      // 	std::array<std::array<std::vector<int>,16>,2> vl1trig;
      // 	//cout << "calling allphisectors.\n";
      // 	L1Anita4LR_ScA_AllPhiSectors(anita1,vl1trig);

      // 	for (int ipol=0;ipol<2;ipol++) {
      // 	  for (int iphi=0;iphi<16;iphi++) {

      // 	    for (int ibin=0;ibin<vl1trig[ipol][iphi].size();ibin++) {
      // 	      anita1->l1trig_anita3and4_inanita[ipol][iphi][ibin]=vl1trig[ipol][iphi][ibin];
	     
      // 	      //if (anita1->l1trig_anita3and4_inanita[ipol][iphi][ibin])
      // 	      //cout << "This l1 is " << ipol << "\t" << iphi << "\t" << ibin << "\t" << anita1->l1trig_anita3and4_inanita[ipol][iphi][ibin] << "\n";
      // 	    }
      // 	    for (int ibin=vl1trig[ipol][iphi].size();ibin<anita1->HALFNFOUR;ibin++) {
      // 	      anita1->l1trig_anita3and4_inanita[ipol][iphi][ibin]=0;
      // 	    }
      // 	  }
      // 	}
      // 	//	cout << "vl1trig size is " << vl1trig[0][0].size() << "\n";

      // 	std::array<std::array<std::vector<int>,16>,2> vl2trig;
	
      // 	L2Anita3and4(anita1,vl1trig,
      // 		     vl2trig);
      // 	L3Anita4LR_ScA(anita1,vl2trig,
      // 		   thispasses);

      //cout << "about to do delays.\n";
      delay_AllAntennas(anita1);

      //	cout << "did delays.\n";
      // l1 triggers for all antennas
      std::array< std::array< vector<int>,16>,3> vl1trig_anita4lr_scb;

      // keep track of whether you get a coincidence between 1, 2 or 3 antennas in a phi sector with the right windows.
      // for each phi sector, whether there is a 1, 2 or 3 coincidence
      std::array<std::array<vector<int>,3>,16> vl2_realtime_anita4_scb;
      //	cout << "did L2.\n";

      std::array<vector<int>,16> vl3trig_type0;
      std::array<vector<int>,16> vl3trig_type1;
      int thispasses_l3type0=0;
      int thispasses_l3type1=0;

	
      double time_thisbin=(double)nstepback*TRIGTIMESTEP;
      int itrigbin=nstepback;
      //cout << "about to loop.\n";
      //cout << "size of arrayofhits is " << arrayofhits[0][0][0].size() << "\n";
      // i should really do this in a different step so this function has fewer inputs.
      while (time_thisbin<LASTTIMETOTESTL1_ANITA4LR_SCB) {
	//cout << "vector size is " << vl1trig_anita4lr_scb[0][0].size() << "\n";
	int npassesl1=0;
	L1Anita4LR_ScB_AllAntennas_OneBin(itrigbin,anita1,vl1trig_anita4lr_scb,npassesl1);

	//	if (npassesl1)
	//cout << "npassesl1 is " << npassesl1 << "\n";

	itrigbin++;
	time_thisbin=(double)itrigbin*TRIGTIMESTEP;
      } // loop over time steps

	//      cout << "got out of first loop.\n";
      time_thisbin=(double)nstepback*TRIGTIMESTEP;
      itrigbin=nstepback;

      while (time_thisbin<LASTTIMETOTESTL1_ANITA4LR_SCB) {
	
	int npassesl2=0;
	int npassesl2_type0=0;
	L2Anita4LR_ScB_AllPhiSectors_OneBin(itrigbin,anita1,vl1trig_anita4lr_scb,
					    vl2_realtime_anita4_scb,npassesl2,npassesl2_type0);


	//  	if (npassesl2)
	//   	  cout << "npassesl2 is " << npassesl2 << "\n";
	//   	if (npassesl2_type0)
	//   	  cout << "npassesl2_type0 is " << npassesl2_type0 << "\n";
	
	itrigbin++;
	time_thisbin=(double)itrigbin*TRIGTIMESTEP;
      } // loop over time steps


      time_thisbin=(double)nstepback*TRIGTIMESTEP;
      itrigbin=nstepback;
      
      
      while (time_thisbin<LASTTIMETOTESTL1_ANITA4LR_SCB) {
	
	
	L3Anita4LR_ScB_OneBin(itrigbin,anita1,vl2_realtime_anita4_scb,
			      vl3trig_type0, vl3trig_type1,
			      thispasses_l3type0,thispasses_l3type1);
	
	// 	cout << "did L3.\n";
	itrigbin++;
	time_thisbin=(double)itrigbin*TRIGTIMESTEP;
      } // loop over time steps
      

      if (thispasses_l3type0)
	thispasses[0]=1;
      if (thispasses_l3type1)
	thispasses[1]=1;

      //	if (thispasses[0] || thispasses[1])
      //cout << "This passes! " << thispasses[0] << "\t" << thispasses[1] << "\n";


      for (int ilayer=0;ilayer<anita1->NTRIGGERLAYERS;ilayer++) {
	for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {
	  for (int ipolar=0;ipolar<anita1->NPOL;ipolar++) {
	    Tools::reverseTimeOrdering(anita1->HALFNFOUR,anita1->arrayofhits_inanita[ilayer][iphi][ipolar],anita1->arrayofhits_forgaryanderic[ilayer][iphi][ipolar]);
	      
	  }
	}
      }


      for (int ilayer=0;ilayer<anita1->NTRIGGERLAYERS;ilayer++) {
	for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {
	  for (unsigned int ibin=0;ibin<vl1trig_anita4lr_scb[ilayer][iphi].size();ibin++) {
	    anita1->l1trig_anita4lr_inanita[ilayer][iphi][ibin]=vl1trig_anita4lr_scb[ilayer][iphi][ibin];
	  }
	  for (int ibin=vl1trig_anita4lr_scb[ilayer][iphi].size();ibin<anita1->HALFNFOUR;ibin++) {
	    anita1->l1trig_anita4lr_inanita[ilayer][iphi][ibin]=0;
	  }

	  Tools::reverseTimeOrdering(anita1->HALFNFOUR,anita1->l1trig_anita4lr_inanita[ilayer][iphi],anita1->l1trig_anita4lr_forgaryanderic[ilayer][iphi]);

	}
      }



      for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {
	for (unsigned int ibin=0;ibin<vl2_realtime_anita4_scb[iphi][1].size();ibin++) {
	  anita1->l2trig_anita4lr_inanita[iphi][1][ibin]=vl2_realtime_anita4_scb[iphi][1][ibin];
	    
	}
	for (int ibin=vl2_realtime_anita4_scb[iphi][1].size();ibin<anita1->HALFNFOUR;ibin++) {
	  anita1->l2trig_anita4lr_inanita[iphi][1][ibin]=0;
	}
	Tools::reverseTimeOrdering(anita1->HALFNFOUR,anita1->l2trig_anita4lr_inanita[iphi][1],anita1->l2trig_anita4lr_forgaryanderic[iphi]);
	  
	for (unsigned int ibin=0;ibin<vl3trig_type0[iphi].size();ibin++) {
	  anita1->l3type0trig_anita4lr_inanita[iphi][ibin]=vl3trig_type0[iphi][ibin];
	}
	for (int ibin=vl3trig_type0[iphi].size();ibin<anita1->HALFNFOUR;ibin++) {
	  anita1->l3type0trig_anita4lr_inanita[iphi][ibin]=0;
	}
	Tools::reverseTimeOrdering(anita1->HALFNFOUR,anita1->l3type0trig_anita4lr_inanita[iphi],anita1->l3type0trig_anita4lr_forgaryanderic[iphi]);
	for (unsigned int ibin=0;ibin<vl3trig_type1[iphi].size();ibin++) {
	    
	  anita1->l3trig_anita4lr_inanita[iphi][ibin]=vl3trig_type1[iphi][ibin];
	}
	for (int ibin=vl3trig_type1[iphi].size();ibin<anita1->HALFNFOUR;ibin++) {
	  anita1->l3trig_anita4lr_inanita[iphi][ibin]=0;
	}
	Tools::reverseTimeOrdering(anita1->HALFNFOUR,anita1->l3trig_anita4lr_inanita[iphi],anita1->l3type1trig_anita4lr_forgaryanderic[iphi]);
      }
 
	
      for (int ibin=0;ibin<anita1->HALFNFOUR;ibin++) {
	anita1->time_trig[ibin]=TRIGTIMESTEP*(double)ibin;
      }














    }
    else {


      // this is the one we're using for Anita 2

      // int iloc=0; // top array
      //int inad=0;  // phi position for nadir antennas,
      //only goes to 8, wheras iphi goes to 16
      int x;//for number of antenna in clump
      int y;//for number of antenna in clump
      x = (antennaclump-1)/2;//setting boundries for clump dependent on clump size(antennaclump)
      y = (antennaclump-1)/2 +1;
      
      for (int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++) {
	for (int ipolar=0;ipolar<Anita::NPOL;ipolar++) {
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
		ltsum+=ant[ipolar][iloc][ihittrig]*ant[ipolar][iloc][iphitrig]; // increment if this antenna passes
	      }
	      // we multipy by ant[iloc][iphi] because this ensures that
	      // one of the triggered antennas is the center one in the clump
	      else
		ltsum+=ant[ipolar][iloc][ihittrig];
	    
	      for (int ipolar=0;ipolar<2;ipolar++) {
		for (int iband=0;iband<4;iband++) {
		  channsum+=channels_compacted_passing[iloc][ihittrig][ipolar][iband]; // increment if this channel passes
		} //for
	      } //for
	    } //for
	  
	    if(ltsum >= Nreq_l2[iloc]){ // if enough antennas in this clump pass
	      loctrig[ipolar][iloc][iphitrig] = 1;
	      l2trig[ipolar][iloc] += (1<<iphi);
	    }//if
	  } //for (phi position)

	} // end loop over polarizations
      } // end loop over trigger layers
      
	// fill the slots in  between the actual antennas with the "or"
	// of the neighboring antennas
      if (settings1->DISCONES==2) {
	FillInNadir(settings1,anita1,l2trig[0][2]);
	FillInNadir(settings1,anita1,l2trig[1][2]);
      }
      // only difference between these is whether they implement masking in both polarizations
      // can we simplify it
      if (settings1->PHIMASKING && settings1->WHICH!=9) { // not for Anita3
	// nadir antennas are aligned with the second physical layer of antennas
	for (int iphi=0;iphi<anita1->NTRIGPHISECTORS;iphi++) { // loop over phi sectors
	  if ((1<<iphi) & phiTrigMask[0]) { // if this phi sector is masked
	    for (int ipolar=0;ipolar<Anita::NPOL;ipolar++) {
	      for (int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++) {
		loctrig[ipolar][iloc][iphi] = 0;
		if(l2trig[ipolar][iloc] & (1<<iphi)) l2trig[ipolar][iloc] -= (1<<iphi);
	      }
	    } // end loop over polarizations
	  } // if this phi sector is masked
	}
      } else if (settings1->PHIMASKING && settings1->WHICH==9) { // only for Anita3
	// nadir antennas are aligned with the second physical layer of antennas
	for (int iphi=0;iphi<anita1->NTRIGPHISECTORS;iphi++) { // loop over phi sectors

	  for (int ipolar=0;ipolar<Anita::NPOL;ipolar++) {
	    if ( (1<<iphi) & phiTrigMask[ipolar] )   { // if this phi sector is masked
	      for (int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++) {
		loctrig[ipolar][iloc][iphi] = 0;
		if(l2trig[ipolar][iloc] & (1<<iphi)) l2trig[ipolar][iloc] -= (1<<iphi);
	      }
	    }
	  }// if this phi sector is masked
	}
      }
    }

    if (anita1->trigRequirements[2]==0) {
      for (int iloc=0;iloc<anita1->NTRIGGERLAYERS;iloc++) {
	for (int ipolar=0;ipolar<Anita::NPOL;ipolar++) {
	  for (int iphi=0;iphi<anita1->NTRIGPHISECTORS;iphi++) {
	    if (l2trig[ipolar][iloc] & (1<<iphi)) { // if we are not making a L3 trigger requirement and there was a l2 trigger then call it a pass
	      thispasses[ipolar]=1;
	      //cout << "This one passes.  inu is " << inu << " " << "iloc is " << iloc << "\n";
	    }
	  }
	} // end loop over polarizations
      }
    }
      

    L3Trigger(settings1,anita1,loctrig,loctrig_nadironly,discones_passing,l3trig,
	      thispasses);

  } //else if (mode 2)
} // end of PassesTriggerBasic




// This is the case for (settings1->TRIGGERSCHEME == 3), that coherent waveform sum is being used to trigger, based on Andres' idea. This is specific to ANITA III's antenna geometry.
    
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

void GlobalTrigger::PassesTriggerCoherentSum(Settings *settings1,Anita *anita1,int inu, int *thispasses) {
    
  for (int center_phi_sector_offset = -1; center_phi_sector_offset <= 1; center_phi_sector_offset++){
    int center_phi_sector_index = first_phi_sector_hit + center_phi_sector_offset;
    if (center_phi_sector_index > 15){center_phi_sector_index = 0;}
    if (center_phi_sector_index < 0){center_phi_sector_index = 15;}
      
    bool coherent_trigger_passes = false;
      
    unsigned N_STEP_PHI = 16;	// The number of different phi hypotheses for a given center phi sector
    unsigned N_STEP_THETA = 81;	// The number of different theta hypotheses for a given center phi sector
    unsigned N_PHI_SECTORS = 16;	// The number of phi sectors
    unsigned N_LAYERS_TRIGGER = 3;	// The number of layers as seen by the trigger
      
    double highest_power = 0;
     
      
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
	      // hi_pow_center = center_phi_sector_index;
	      // hi_pow_phi_index = index_phi;
	      // hi_pow_theta_index = index_theta;
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
	      // used to be return this_passes;
	      return;
	    }
	      
	  } // for window indices
	} // for hypothesis theta indices
      } // for hypothesis phi indices
    } // for hypothesis center phi sector indices
      

  }
} // end of PassesTriggerCoherentSum




void GlobalTrigger::PassesTriggerSummedPower(Settings *settings1,Anita *anita1){
				  
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
    
    
  //  For whatever reason ANTENNA_POSITION_START is a c-array, with large dimensions (5 by 400).
  //  Most of these are initialized as a Vector with (0.,0.,1) for (x,y,z). These are non-existant antennas.
  //  
  //  Make sure that the antennas positions being copied are actually from a real antenna, and pay attention to the physical
  //  antenna layers 0 & 1, as they map to the 0th trigger layer.
    
    
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
    if (temporary_passing_variable) {index_window=31;} //exit
  }
} // end of PassesTriggerSummedPower


void GlobalTrigger::PassesTriggerScheme5(Anita *anita1,double this_threshold, int *thispasses) {
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
      
	if (maxbands[0]>=anita1->NCH_PASS) 

	  thispasses[0]=1.;	
	if ( maxbands[1]>=anita1->NCH_PASS)  
	  thispasses[1]=1.;	

      
      } // end loop over layers
    
    }
      
  } // end loop over center phi sectors    

} // end of PassesTriggerScheme5



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
  else if (settings1->WHICH==9 || settings1->WHICH==10) { // anita 4
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
void GlobalTrigger::L3Trigger(Settings *settings1,Anita *anita1,int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX],int discones_passing,int *l3trig,
			      int *thispasses) {
  
  int whichphipass[Anita::NPOL][Anita::NPHI_MAX]={{0}};
   
  for (int i=0;i<Anita::NTRIG;i++)
    triggerbits[i]=0;
    
  // shouldn't get here for anita3 or anita4 anymore.
  //assumes the number of trigger phi sectors are the same in each layer
  for (int i=0;i<anita1->PHITRIG[0];i++) {
		
		
    //     if (settings1->WHICH==9) { // anita 3
    //       for (int ipolar=0;ipolar<2;ipolar++) {
    // 	if (loctrig[ipolar][0][i]>0 && loctrig[ipolar][1][i]>0 && loctrig[ipolar][2][i]>0) {
    // 	  thispasses[ipolar]=1;
    // 	  whichphipass[ipolar][i]=1;
    // 	  if (!(triggerbits[3] & (1<<i)))
    // 	    triggerbits[3]+=(1<<i);
    // 	}
    //       }
    //     }
    //    else { // anita 1 or anita 2
      
    for (int ilayer=0;ilayer<anita1->NTRIGGERLAYERS-1;ilayer++) {
      for (int ipolar=0;ipolar<2;ipolar++) {
	// check this - does ilayer and ilayer+1 include just 1 and 2, or also 2 and 3, and is that ok?
	if (loctrig[ipolar][ilayer][i]>0 && loctrig[ipolar][ilayer+1][i]>0) { // upper two layers
	    
	  //l3trig += (1<<i);
	  thispasses[ipolar]=1;
	  whichphipass[ipolar][i]=1;
	    
	  if(loctrig[ipolar][0][i]>0 && loctrig[ipolar][1][i]>0) {
	    if (!(triggerbits[0] & (1<<i)))
	      triggerbits[0] += (1<<i);

	  }
	    
 
	}
      }
    }
    
    // anita 1, anita 1 kurt, anita 2 kurt, anita 3
    // but this is inside a set of brackets for !anita3 
    if (settings1->WHICH==2 || settings1->WHICH==6 || settings1->WHICH==8) { // Anita 1, 2
				
      for (int ipolar=0;ipolar<2;ipolar++) {
	if (((settings1->DISCONES==2) && (loctrig[ipolar][1][i]>0 && loctrig[ipolar][2][i]>0)) || // second layer and nadir?
	    ((settings1->DISCONES==2) && (loctrig[ipolar][0][i]>0 && loctrig[ipolar][2][i]>0)) || // first layer and nadir?
	    ((settings1->DISCONES==2) && loctrig_nadironly[ipolar][i]>0) || // just nadir - does this mean having just a nadir pass makes the whole event pass?
	    (settings1->DISCONES==1 && (loctrig[ipolar][1][i]>0 && discones_passing>=settings1->NDISCONES_PASS)) || // second layer and discones
	    (settings1->DISCONES==1 && (loctrig[ipolar][0][i]>0 && discones_passing>=settings1->NDISCONES_PASS))   // top layer and discones
	    ) {
					
	  thispasses[ipolar]=1;
	  whichphipass[ipolar][i]=1;
	  // this is just for nadir studies->
	  //how many of each trigger condition
					
	  if(loctrig[ipolar][1][i]>0 && loctrig[ipolar][2][i]>0) {
	    triggerbits[1] += (1<<i);
						
	  }
	  if(loctrig[ipolar][0][i]>0 && loctrig[ipolar][2][i]>0) {
	    triggerbits[2]+=(1<<i);
						
	  }
	  // 	if (loctrig_nadironly[i]>0)
	  // 	  triggerbits[3]+=weight2;
	  // 	if ((loctrig_nadironly[i]>0 ||
	  // 	     (loctrig[1][i]>0 && loctrig[2][i]>0) ||
	  // 	     (loctrig[0][i]>0 && loctrig[2][i]>0)) &&
	  // 	    !(loctrig[0][i]>0 && loctrig[1][i]>0))
	  // 	  triggerbits[4]+=weight2;
					
					
	  thispasses[ipolar]=1;
					
	}
				
      } //if it's the anita 1 payload
    }	
			
    //    } // if it's not anita 3
		
		
    for (int ipolar=0;ipolar<2;ipolar++) {
      if (whichphipass[ipolar][i])
	l3trig[ipolar]+=(1<<i);
      // std::cout << ipolar << " " << whichphipass[ipolar][i] << " " << l3trig[ipolar] << std::endl;
    }
		
  } // end looping over phi
    

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
int GlobalTrigger::findahit(vector<int> myvector,int first,int last) {

  int yes=0;
  if (myvector.size()<=last)
    return yes;

  for (int i=first;i<=last;i++) {
    if (myvector[i]==1) {
      yes=1;
      continue;
    }
  }
  return yes;
}
int GlobalTrigger::findanl3(int *l3,int NPHISECTORS) {

  for (int j=0;j<NPHISECTORS;j++) {
    if (l3[j]==1)
      return 1;
  }
  return 0;
}
void GlobalTrigger::L2Anita3and4(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> vl1trig,
				 std::array<std::array<std::vector<int>,16>,2> &vl2trig) {

  
  for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {
    for (int ipolar=0;ipolar<2;ipolar++) {
      //cout << "iphi, ipolar, size is " <<  iphi << "\t" << ipolar << "\t" << vl1trig[ipolar][iphi].size() << "\n";
      for (unsigned int ibin=0;ibin<vl1trig[ipolar][iphi].size();ibin++) {
	//cout << "ipolar, iphi , ibin are " << ipolar << "\t" << iphi << "\t" << ibin << "\n";
	// 	if (vl1trig[ipolar][iphi][ibin])
	// 	  cout << "actually a nonzero l1trig.\n";
	vl2trig[ipolar][iphi].push_back(vl1trig[ipolar][iphi][ibin]);
      }
    }
  }

  


}
void GlobalTrigger::L3Anita3and4(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> vl2trig,
				 int vl3trig[2][16],int *thispasses) {

  int iphimod16_neighbor;

  for (int ipolar=0;ipolar<2;ipolar++) {
    for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {

      for (int iphi_neighbor=iphi-1;iphi_neighbor<=iphi+1;iphi_neighbor++) {
	if (iphi_neighbor!=iphi) {
	  iphimod16_neighbor=(iphi_neighbor+16)%16;
	  //cout << "L3_COINCIDENCE is " << L3_COINCIDENCE << "\n";
	  for (unsigned int ibin=0;ibin<vl2trig[ipolar][iphi].size();ibin++) {
	    if (ibin+(int)(L3_COINCIDENCE/TRIGTIMESTEP)<vl2trig[ipolar][iphimod16_neighbor].size()) {
	      if (vl2trig[ipolar][iphi][ibin] && 
		  findahit(vl2trig[ipolar][iphimod16_neighbor],ibin,ibin+(int)(L3_COINCIDENCE/TRIGTIMESTEP))) {
		//cout << "passes: " << ipolar << "\t" << iphi << "\t" << ibin << "\n";
		//cout << "neighbor: " << ipolar << "\t" << iphimod16_neighbor << "\n";
		vl3trig[ipolar][iphi]=1;
		thispasses[ipolar]=1;
		ibin = vl2trig[ipolar][iphi].size(); // ends this loop
		iphi_neighbor=iphi+2; // ends outer loop
	      }
	      else
		vl3trig[ipolar][iphi]=0;
 
	    } // end if so we're not going off the end of vl2trig
	    else
	      vl3trig[ipolar][iphi]=0;
	  }
	}
      }
    }
  }

}

// L1 trigger is at the antenna level again.  Just require coincidence between LCP and RCP
void GlobalTrigger::L1Anita4LR_ScB_OneBin(int IZERO,vector<int> vleft,vector<int> vright,
					  vector<int> &vl1trig) {

  if ((vleft[IZERO] && findahit(vright,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_ANITA4LR_SCB/TRIGTIMESTEP))) ||
      
      (vright[IZERO] && findahit(vleft,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_ANITA4LR_SCB/TRIGTIMESTEP)))) {
    vl1trig.push_back(1);

  }
  else
    vl1trig.push_back(0);


}
// L1 trigger is at the antenna level again.  Just require coincidence between LCP and RCP
void GlobalTrigger::L1Anita4LR_ScB_AllAntennas_OneBin(int IZERO,Anita *anita1,std::array< std::array< vector<int>,16>,3> &vl1trig_anita4lr_scb,int &npassesl1) {
 

  for (int itriglayer=0;itriglayer<anita1->NTRIGGERLAYERS;itriglayer++) {
    for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {

      //double time_thisbin=(double)nstepback*TRIGTIMESTEP;
      //int itrigbin=nstepback;
      
      // i should really do this in a different step so this function has fewer inputs.
      //while (time_thisbin<LASTTIMETOTESTL1_ANITA4LR_SCB) {
      L1Anita4LR_ScB_OneBin(IZERO,arrayofhits[itriglayer][iphi][0][4],arrayofhits[itriglayer][iphi][1][4],
			    vl1trig_anita4lr_scb[itriglayer][iphi]);
      if (vl1trig_anita4lr_scb[itriglayer][iphi][vl1trig_anita4lr_scb[itriglayer][iphi].size()-1])
	npassesl1++;

      //itrigbin++;
      //time_thisbin=(double)itrigbin*TRIGTIMESTEP;
      // } // loop over time steps
    } // loop over phi sectors
  } // loop over trigger layers
  


}


// void GlobalTrigger::L3Anita4LR_ScB(Anita *anita1,std::array<std::array<vector<int>,3>,16> vl2_realtime_anita4_scb,
// 				   std::array<vector<int>,16> &vl3trig_type0, std::array<vector<int>,16> &vl3trig_type1,
// 				   int &thispasses_l3type0,int &thispasses_l3type1) {
//   thispasses_l3type0=0;
//   thispasses_l3type1=0;


//   for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {

//     vl3trig_type0[iphi].clear();
//     vl3trig_type1[iphi].clear();

//     int iphi_next=(iphi+1+16)%16;
//     int iphi_previous=(iphi-1+16)%16;

//     if (L3or30Anita4LR_ScB_TwoPhiSectors( 
// 				     vl2_realtime_anita4_scb[iphi], // 3 neighbors, whether 1, 2 or 3 pass
// 				     vl2_realtime_anita4_scb[iphi_previous], // 3 neighbors, whether 1, 2 or 3 pass
				     
// 				     2,2,	
// 					 vl3trig_type1[iphi] ))
//       thispasses_l3type1=1;

//     if (iphi%2==0) {
//       if (L3or30Anita4LR_ScB_TwoPhiSectors( 
// 				       vl2_realtime_anita4_scb[iphi], // 3 neighbors, whether 1, 2 or 3 pass
// 				       vl2_realtime_anita4_scb[iphi_previous], // 3 neighbors, whether 1, 2 or 3 pass
				       
// 				       1,3,	
// 					   vl3trig_type0[iphi] ))
// 	thispasses_l3type0=1;
	    


//     }
//   }
  

// }
// for each phi sector, does 1, 2 or 3 pass
// void GlobalTrigger::L2Anita4LR_ScB_AllPhiSectors(Anita *anita1,std::array< std::array< vector<int>,16>,3> vl1trig_anita4lr_scb,
// 						 std::array<std::array<vector<int>,3>,16> &vl2_realtime_anita4_scb) {

//   double time_thisbin=(double)nstepback*TRIGTIMESTEP;
//   int itrigbin=nstepback;

//   while (time_thisbin<LASTTIMETOTESTL2_ANITA4LR_SCB) {

//       for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {
	
// 	//      cout << "itriglayer, iphi are " << itriglayer << "\t" << iphi << "\n";
// 	int npassesl2=0;
// 	int npassesl2type0=0;

// 	L2Anita4LR_ScB_OnePhiSector_OneBin(itrigbin,vl1trig_anita4lr_scb[2][iphi], 
// 					   vl1trig_anita4lr_scb[1][iphi], 
// 					   vl1trig_anita4lr_scb[0][iphi], 		
// 					   vl2_realtime_anita4_scb[iphi],
// 					   npassesl2,npassesl2type0);

	
//       }
      
    
//     itrigbin++;
//     time_thisbin=(double)itrigbin*TRIGTIMESTEP;
//   }


// }

// for each phi sector, does 1, 2 or 3 pass
void GlobalTrigger::L2Anita4LR_ScB_AllPhiSectors_OneBin(int IZERO,Anita *anita1,std::array< std::array< vector<int>,16>,3> vl1trig_anita4lr_scb,
							std::array<std::array<vector<int>,3>,16> &vl2_realtime_anita4_scb,int &npassesl2,int &npassesl2_type0) {

  for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {
      
    //      cout << "itriglayer, iphi are " << itriglayer << "\t" << iphi << "\n";
      
    L2Anita4LR_ScB_OnePhiSector_OneBin(IZERO,vl1trig_anita4lr_scb[2][iphi], 
				       vl1trig_anita4lr_scb[1][iphi], 
				       vl1trig_anita4lr_scb[0][iphi], 		
				       vl2_realtime_anita4_scb[iphi],npassesl2,npassesl2_type0);
      
      
  }
    
  


}

// ask if L3 type 1 (2 and 2) or L3 type 0 (3 and 1 or 1 and 3) passes
// int GlobalTrigger::L3or30Anita4LR_ScB_TwoPhiSectors(
// 						    std::array<vector<int>,3> vl2_realtime_anita4_scb, // 3 neighbors, whether 1, 2 or 3 pass
// 						    std::array<vector<int>,3> vl2_realtime_anita4_scb_other, // 3 neighbors, whether 1, 2 or 3 pass
// 						    int npass1,int npass2,	
// 						    vector<int> &vl3trig ) {

  

//   int passaone=0;



//   int thispasses=0;

  
//   if (L3or30Anita4LR_ScB_TwoPhiSectors_OneBin(itrigbin,
// 					      vl2_realtime_anita4_scb, // 3 neighbors, whether 1, 2 or 3 pass
// 					      vl2_realtime_anita4_scb_other, // 3 neighbors, whether 1, 2 or 3 pass
// 					      npass1,npass2,	
// 					      vl3trig )) {
    
    
//     thispasses=1;
//     passaone=1;
//     //cout << "got an l3 with " << npass1 << " and " << npass2 << "\n";
//   }
//   vl3trig.push_back(thispasses);
  
//   return passaone;
  
// }
// ask if L3 type 1 (2 and 2) or L3 type 0 (3 and 1 or 1 and 3) passes
// int GlobalTrigger::L3or30Anita4LR_ScB_TwoPhiSectors_OneBin( int IZERO,
// 						    std::array<vector<int>,3> vl2_realtime_anita4_scb, // 3 neighbors, whether 1, 2 or 3 pass
// 						    std::array<vector<int>,3> vl2_realtime_anita4_scb_other, // 3 neighbors, whether 1, 2 or 3 pass
// 						    int npass1,int npass2,	
// 						    vector<int> &vl3trig ) {

  

//   int passaone=0;
//   //double time_thisbin=(double)nstepback*TRIGTIMESTEP;
//   //int itrigbin=nstepback;
//   int thispasses=0;
//   //while (time_thisbin<LASTTIMETOTESTL1_ANITA4LR_SCB) {   
//     thispasses=0;


//     if ((vl2_realtime_anita4_scb[npass1-1][IZERO] &&
// 	 findahit(vl2_realtime_anita4_scb_other[npass2-1],IZERO-nstepback,IZERO-nstepback+(int)(L3_COINCIDENCE_ANITA4LR_SCB/TRIGTIMESTEP)))   ||
// 	(vl2_realtime_anita4_scb_other[npass2-1][IZERO] &&
// 	 findahit(vl2_realtime_anita4_scb[npass1-1],IZERO-nstepback,IZERO-nstepback+(int)(L3_COINCIDENCE_ANITA4LR_SCB/TRIGTIMESTEP)))) {
//       thispasses=1;
//       passaone=1;
//       //cout << "got an l3 with " << npass1 << " and " << npass2 << "\n";
//     }
//     //    cout << "adding to vl3trig.\n";
//     vl3trig.push_back(thispasses);
 
//     //itrigbin++;
//     //time_thisbin=(double)itrigbin*TRIGTIMESTEP;
//     //}
//   return passaone;
  
// }

void GlobalTrigger::L3Anita4LR_ScA(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> vl2trig,
				   int *thispasses) {


  for (int ipolar=0;ipolar<2;ipolar++) {
    for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {


      
      for (unsigned int ibin=0;ibin<vl2trig[ipolar][iphi].size();ibin++) {

	if (vl2trig[ipolar][iphi][ibin]) {

	  thispasses[ipolar]=1;
	}
      }

    
    }

  }

}
void GlobalTrigger::L1Anita3_AllPhiSectors(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> &vl1trig) {
  vector<int> vl1_realtime_vbottom;
  vector<int> vl1_realtime_vmiddle; 
  vector<int> vl1_realtime_vtop;
  vector<int> vl1_realtime_hbottom;
  vector<int> vl1_realtime_hmiddle; 
  vector<int> vl1_realtime_htop;
  double time_thisbin=(double)nstepback*TRIGTIMESTEP;
  int itrigbin=nstepback;
  //  cout << "phitrig is " << anita1->PHITRIG[0] << "\n";
  for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {
    itrigbin=nstepback;
    time_thisbin=(double)nstepback*TRIGTIMESTEP;
    while (time_thisbin<LASTTIMETOTESTL1_ANITA3) {

      vl1trig[0][iphi].push_back(L1Anita3_OnePhiSector(itrigbin,arrayofhits[2][iphi][0][4], 
						       arrayofhits[1][iphi][0][4], 
						       arrayofhits[0][iphi][0][4],
						       vl1_realtime_vbottom, vl1_realtime_vmiddle, vl1_realtime_vtop));
      vl1trig[1][iphi].push_back(L1Anita3_OnePhiSector(itrigbin,arrayofhits[2][iphi][1][4], arrayofhits[1][iphi][1][4], arrayofhits[0][iphi][1][4],
						       vl1_realtime_hbottom, vl1_realtime_hmiddle, vl1_realtime_htop));
      itrigbin++;
      //      if (vl1trig[0][iphi][vl1trig[0][iphi].size()-1]==1) {
      // 	cout << "an l1.\n";
      // 	cout << "iphi, bin is " << iphi << "\t" << vl1trig[0][iphi].size()-1 << "\n";
      //       }
      time_thisbin=(double)itrigbin*TRIGTIMESTEP;
    }
  }
  //  cout << "inside, vl1trig is " << vl1trig[0][0].size() << "\n";
}
void GlobalTrigger::L1Anita4_AllPhiSectors(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> &vl1trig) {
  vector<int> vl1_realtime_vbottom;
  vector<int> vl1_realtime_vmiddle; 
  vector<int> vl1_realtime_vtop;
  vector<int> vl1_realtime_hbottom;
  vector<int> vl1_realtime_hmiddle; 
  vector<int> vl1_realtime_htop;
  double time_thisbin=(double)nstepback*TRIGTIMESTEP;
  int itrigbin=nstepback;
  //  cout << "phitrig is " << anita1->PHITRIG[0] << "\n";
  for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {
    itrigbin=nstepback;
    time_thisbin=(double)nstepback*TRIGTIMESTEP;
    while (time_thisbin<LASTTIMETOTESTL1_ANITA4) {
      //cout << "iphi is " << iphi << "\n";
      vl1trig[0][iphi].push_back(L1Anita4_OnePhiSector(itrigbin,arrayofhits[2][iphi][0][4], arrayofhits[1][iphi][0][4], arrayofhits[0][iphi][0][4],
						       vl1_realtime_vbottom, vl1_realtime_vmiddle, vl1_realtime_vtop));
      vl1trig[1][iphi].push_back(L1Anita4_OnePhiSector(itrigbin,arrayofhits[2][iphi][1][4], arrayofhits[1][iphi][1][4], arrayofhits[0][iphi][1][4],
						       vl1_realtime_hbottom, vl1_realtime_hmiddle, vl1_realtime_htop));
   
      //      if (vl1trig[0][iphi][vl1trig[0][iphi].size()-1]==1) {
      //  	cout << "an l1.\n";
      //  	cout << "iphi, bin is " << iphi << "\t" << vl1trig[0][iphi].size()-1 << "\n";
      // 	cout << "arrayofhits are \n";
      // 	int istart=itrigbin-nstepback;
      // 	int iend=itrigbin-nstepback+(int)(L1_COINCIDENCE_MOREGENERAL[0][1]/TRIGTIMESTEP)-1;
      // 	cout << "istart, iend are " << istart << "\t" << iend << "\n";
      // 	cout << "sizes are " << arrayofhits[0][iphi][0][4].size() << "\t" << arrayofhits[1][iphi][0][4].size() << "\t" << arrayofhits[2][iphi][0][4].size() << "\n";
      // 	cout << "top:\n";
      // 	for (int i=istart;i<=iend;i++) {
      // 	  cout << arrayofhits[0][iphi][0][4][i] << "\t";
      // 	}
      // 	cout << "middle:\n";
      // 	for (int i=istart;i<=iend;i++) {
      // 	  cout << arrayofhits[1][iphi][0][4][i] << "\t";
      // 	}
      // 	cout << "bottom:\n";
      // 	for (int i=istart;i<=iend;i++) {
      // 	  cout << arrayofhits[2][iphi][0][4][i] << "\t";
      // 	}
      // 	cout << "\n";
	
      //       }
      itrigbin++;
      time_thisbin=(double)itrigbin*TRIGTIMESTEP;
    }
  }
  //  cout << "inside, vl1trig is " << vl1trig[0][0].size() << "\n";
}
void GlobalTrigger::L1Anita4LR_ScA_AllPhiSectors(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> &vl1trig) {
  vector<int> vl1_realtime_1bottom;
  vector<int> vl1_realtime_1middle; 
  vector<int> vl1_realtime_1top;
  vector<int> vl1_realtime_2bottom;
  vector<int> vl1_realtime_2middle; 
  vector<int> vl1_realtime_2top;

  double time_thisbin=(double)nstepback*TRIGTIMESTEP;
  int itrigbin=nstepback;
  int taketheor=0; // Take the "or" of the L1 with a phi sector and its neighbor to the left, and the L1 with a phi sector and its
  // neighbor to the right.
  //  cout << "phitrig is " << anita1->PHITRIG[0] << "\n";
  for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {
    itrigbin=nstepback;
    time_thisbin=(double)nstepback*TRIGTIMESTEP;
    taketheor=0;
    int iphi_rightneighbor=(iphi+1+16)%16;
    int iphi_leftneighbor=(iphi-1+16)%16;
    while (time_thisbin<LASTTIMETOTESTL1_ANITA4LR_SCA) {
      
      for (int ipolar=0;ipolar<2;ipolar++) {

	if (L1Anita4LR_ScA_TwoPhiSectors(itrigbin,ipolar,
					 arrayofhits[2][iphi][0][4],arrayofhits[2][iphi][1][4], 
					 arrayofhits[2][iphi_rightneighbor][0][4],arrayofhits[2][iphi_rightneighbor][1][4], 
					 arrayofhits[1][iphi][0][4], arrayofhits[1][iphi][1][4], 
					 arrayofhits[1][iphi_rightneighbor][0][4], arrayofhits[1][iphi_rightneighbor][1][4], 
					 arrayofhits[0][iphi][0][4], arrayofhits[0][iphi][1][4], 
					 arrayofhits[0][iphi_rightneighbor][0][4],arrayofhits[0][iphi_rightneighbor][1][4], 
					 vl1_realtime_1bottom, vl1_realtime_1middle, vl1_realtime_1top) ||
	  
	    L1Anita4LR_ScA_TwoPhiSectors(itrigbin,ipolar,
					 arrayofhits[2][iphi][0][4],arrayofhits[2][iphi][1][4], 
					 arrayofhits[2][iphi_leftneighbor][0][4],arrayofhits[2][iphi_leftneighbor][1][4], 
					 arrayofhits[1][iphi][0][4], arrayofhits[1][iphi][1][4], 
					 arrayofhits[1][iphi_leftneighbor][0][4], arrayofhits[1][iphi_leftneighbor][1][4], 
					 arrayofhits[0][iphi][0][4], arrayofhits[0][iphi][1][4], 
					 arrayofhits[0][iphi_leftneighbor][0][4],arrayofhits[0][iphi_leftneighbor][1][4], 
					 vl1_realtime_2bottom, vl1_realtime_2middle, vl1_realtime_2top))
	  vl1trig[ipolar][iphi].push_back(1);
	else
	  vl1trig[ipolar][iphi].push_back(0);

	  
	//	if (iphi==13) {
	//cout << "itrigbin is " << itrigbin << "\n";
	// cout << "first one is " << L1Anita4LR_ScA_TwoPhiSectors(iphi,itrigbin,ipolar,
	// 					   arrayofhits[2][iphi][0][4],arrayofhits[2][iphi][1][4], 
	// 					   arrayofhits[2][iphi_rightneighbor][0][4],arrayofhits[2][iphi_rightneighbor][1][4], 
	// 					   arrayofhits[1][iphi][0][4], arrayofhits[1][iphi][1][4], 
	// 					   arrayofhits[1][iphi_rightneighbor][0][4], arrayofhits[1][iphi_rightneighbor][1][4], 
	// 					   arrayofhits[0][iphi][0][4], arrayofhits[0][iphi][1][4], 
	// 					   arrayofhits[0][iphi_rightneighbor][0][4],arrayofhits[0][iphi_rightneighbor][1][4], 
	// 							      vl1_realtime_1bottom, vl1_realtime_1middle, vl1_realtime_1top) << "\n";
	// 	  cout << "second one is " << L1Anita4LR_ScA_TwoPhiSectors(iphi,itrigbin,ipolar,
	// 					   arrayofhits[2][iphi][0][4],arrayofhits[2][iphi][1][4], 
	// 					   arrayofhits[2][iphi_leftneighbor][0][4],arrayofhits[2][iphi_leftneighbor][1][4], 
	// 					   arrayofhits[1][iphi][0][4], arrayofhits[1][iphi][1][4], 
	// 					   arrayofhits[1][iphi_leftneighbor][0][4], arrayofhits[1][iphi_leftneighbor][1][4], 
	// 					   arrayofhits[0][iphi][0][4], arrayofhits[0][iphi][1][4], 
	// 					   arrayofhits[0][iphi_leftneighbor][0][4],arrayofhits[0][iphi_leftneighbor][1][4], 
	// 							       vl1_realtime_2bottom, vl1_realtime_2middle, vl1_realtime_2top) << "\n";

	//cout << "ipolar, iphi, taketheor are " << ipolar << "\t" << iphi  << "\t" << vl1trig[ipolar][iphi][itrigbin] << "\n";

	// 	}

      }

      itrigbin++;
      time_thisbin=(double)itrigbin*TRIGTIMESTEP;
    }
  }
  //  cout << "inside, vl1trig is " << vl1trig[0][0].size() << "\n";
}


int GlobalTrigger::L1Anita3_OnePhiSector(int IZERO,vector<int> &vl0_realtime_bottom, vector<int> &vl0_realtime_middle, vector<int> &vl0_realtime_top,
					 vector<int> &vl1_realtime_bottom, vector<int> &vl1_realtime_middle, vector<int> &vl1_realtime_top) {



  // ask if L1 trigger passes
  // Patrick says "check every 2 ns for 4 ns back and start 4 ns gate, or 12 ns, or 16 ns."
  if (vl0_realtime_bottom[IZERO]==1) {

    //cout << "findahit from " << IZERO-nstepback << " to " << IZERO-nstepback+(int)(L1_COINCIDENCE[0]/TRIGTIMESTEP)-1 << "\n";
    if (findahit(vl0_realtime_top,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_ANITA3[0]/TRIGTIMESTEP)) ||
	findahit(vl0_realtime_middle,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_ANITA3[0]/TRIGTIMESTEP))) {


      vl1_realtime_bottom.push_back(1);
      // cout << "got an l1 here. vl1_realtime_bottom.size() is " << vl1_realtime_bottom.size() << "\n";

    }
    else
      vl1_realtime_bottom.push_back(0);

  }
  else {
    vl1_realtime_bottom.push_back(0);
           
  }
     
  if (vl0_realtime_middle[IZERO]==1) {
    if (findahit(vl0_realtime_top,IZERO,IZERO-nstepback+(int)(L1_COINCIDENCE_ANITA3[1]/TRIGTIMESTEP)) ||
	findahit(vl0_realtime_bottom,IZERO,IZERO-nstepback+(int)(L1_COINCIDENCE_ANITA3[1]/TRIGTIMESTEP))) {
	 
      vl1_realtime_middle.push_back(1);
	 
    }
    else
      vl1_realtime_middle.push_back(0);

  }
  else {
    vl1_realtime_middle.push_back(0);
           
  }
  //     if (vl1_realtime_bottom[vl1_realtime_bottom.size()-1])
  //cout << "the bottom bit is " << vl1_realtime_bottom[vl1_realtime_bottom.size()-1] << "\n";
  if (vl0_realtime_top[IZERO]==1) {
    if (findahit(vl0_realtime_middle,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_ANITA3[2]/TRIGTIMESTEP)) ||
	findahit(vl0_realtime_bottom,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_ANITA3[2]/TRIGTIMESTEP))) {
	 
      vl1_realtime_top.push_back(1);
	 
    }
    else
      vl1_realtime_top.push_back(0);

  }
  else {
    vl1_realtime_top.push_back(0);
           
  }
     
  if (vl1_realtime_bottom[vl1_realtime_bottom.size()-1]==1 || vl1_realtime_top[vl1_realtime_top.size()-1]==1 || vl1_realtime_middle[vl1_realtime_middle.size()-1]==1) {
    //cout << "actually got an l1.\n";
    return 1;

  }
  else return 0;
}
int GlobalTrigger::L1Anita4_OnePhiSector(int IZERO,vector<int> &vl0_realtime_bottom, vector<int> &vl0_realtime_middle, vector<int> &vl0_realtime_top,
					 vector<int> &vl1_realtime_bottom, vector<int> &vl1_realtime_middle, vector<int> &vl1_realtime_top) {



  // ask if L1 trigger passes
  // Patrick says "check every 2 ns for 4 ns back and start 4 ns gate, or 12 ns, or 16 ns."
  if (vl0_realtime_bottom[IZERO]==1) {
    //cout << "findahit from " << IZERO-nstepback << " to " << IZERO-nstepback+(int)(L1_COINCIDENCE[0]/TRIGTIMESTEP)-1 << "\n";
    if (findahit(vl0_realtime_top,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_MOREGENERAL[0][0]/TRIGTIMESTEP)) ||
	findahit(vl0_realtime_middle,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_MOREGENERAL[0][1]/TRIGTIMESTEP))) {
	 
      vl1_realtime_bottom.push_back(1);
	 
    }
    else
      vl1_realtime_bottom.push_back(0);

  }
  else {
    vl1_realtime_bottom.push_back(0);
           
  }

  if (vl0_realtime_middle[IZERO]==1) {
    if (findahit(vl0_realtime_top,IZERO,IZERO-nstepback+(int)(L1_COINCIDENCE_MOREGENERAL[1][0]/TRIGTIMESTEP)) ||
	findahit(vl0_realtime_bottom,IZERO,IZERO-nstepback+(int)(L1_COINCIDENCE_MOREGENERAL[1][1]/TRIGTIMESTEP))) {
	 
      vl1_realtime_middle.push_back(1);
	 
    }
    else
      vl1_realtime_middle.push_back(0);

  }
  else {
    vl1_realtime_middle.push_back(0);
           
  }

  if (vl0_realtime_top[IZERO]==1) {
    if (findahit(vl0_realtime_middle,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_MOREGENERAL[2][0]/TRIGTIMESTEP)) ||
	findahit(vl0_realtime_bottom,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_MOREGENERAL[2][1]/TRIGTIMESTEP))) {
	 
      vl1_realtime_top.push_back(1);
	 
    }
    else
      vl1_realtime_top.push_back(0);

  }
  else {
    vl1_realtime_top.push_back(0);
           
  }

  if (vl1_realtime_bottom[vl1_realtime_bottom.size()-1]==1 || vl1_realtime_top[vl1_realtime_top.size()-1]==1 || vl1_realtime_middle[vl1_realtime_middle.size()-1]==1) 
    return 1;
  else return 0;
}
int GlobalTrigger::L1Anita4LR_ScA_TwoPhiSectors(int IZERO,int ipolar,
						vector<int> &v1l0_realtime_bottomleft, vector<int> &v2l0_realtime_bottomleft, 
						vector<int> &v1l0_realtime_bottomright, vector<int> &v2l0_realtime_bottomright, 
						vector<int> &v1l0_realtime_middleleft, vector<int> &v2l0_realtime_middleleft,
						vector<int> &v1l0_realtime_middleright, vector<int> &v2l0_realtime_middleright,
						vector<int> &v1l0_realtime_topleft, vector<int> &v2l0_realtime_topleft,
						vector<int> &v1l0_realtime_topright, vector<int> &v2l0_realtime_topright,
						vector<int> &vl1_realtime_bottom, 
						vector<int> &vl1_realtime_middle, 
						vector<int> &vl1_realtime_top) {

  
  //if (IZERO>=42 && IZERO<=44 && iphi==13)
  //cout << "inside _TwoPhiSectors. IZERO, ipolar are " <<  IZERO << "\t" << ipolar << "\n";

  vector<int> vleft;
  vector<int> vright;
  if (ipolar==0) {
    vleft=v1l0_realtime_bottomleft;
    vright=v1l0_realtime_bottomright;
  }
  else if (ipolar==1) {
    vleft=v2l0_realtime_bottomleft;
    vright=v2l0_realtime_bottomright;
  }

  if (vleft[IZERO]==1 || vright[IZERO]==1) {
    //   if (IZERO>=40 && IZERO<=44 && iphi==13)
    //cout << "got a bottom trigger.  iphi, izero, ipolar is " << iphi << "\t" << IZERO << "\t" << ipolar << "\t" << vleft[IZERO] << "\t" << vright[IZERO] << "\n";
    if (PartofL1Anita4LR_ScA_TwoPhiSectors(1,ipolar,IZERO,v1l0_realtime_middleleft, v2l0_realtime_middleleft,
					   v1l0_realtime_middleright, v2l0_realtime_middleright,
					   vl1_realtime_middle) &&
	PartofL1Anita4LR_ScA_TwoPhiSectors(0,ipolar,IZERO,v1l0_realtime_topleft, v2l0_realtime_topleft,
					   v1l0_realtime_topright, v2l0_realtime_topright,
					   vl1_realtime_top)) {
      //if (IZERO>=40 && IZERO<=44 && iphi==13)
      //cout << "other two layers triggered too.\n";
      return 1;
    }
    else {
      //       if (IZERO>=40 && IZERO<=44 && iphi==13) {
      //       cout << "middle is " << PartofL1Anita4LR_TwoPhiSectors(iphi,1,ipolar,IZERO,v1l0_realtime_middleleft, v2l0_realtime_middleleft,
      // 				       v1l0_realtime_middleright, v2l0_realtime_middleright,
      // 							     vl1_realtime_middle) << "\n";
      //       cout << "top is " << PartofL1Anita4LR_TwoPhiSectors(iphi,0,ipolar,IZERO,v1l0_realtime_topleft, v2l0_realtime_topleft,
      // 				       v1l0_realtime_topright, v2l0_realtime_topright,
      // 							  vl1_realtime_top) << "\n";
      //       cout << "other two layers did not trigger.\n";
      //       }
    }
  }
  return 0;


}
  


int GlobalTrigger::PartofL1Anita4LR_ScA_TwoPhiSectors(int ilayerreverse,int ipolar,int IZERO,
						      vector<int> &v1l0_realtime_left, vector<int> &v2l0_realtime_left, 
						      vector<int> &v1l0_realtime_right, vector<int> &v2l0_realtime_right, 
						      vector<int> &vl1_realtime) {
    
  //  if (iphi==13)
  //cout << "ilayerreverse is " << ilayerreverse << "\n";
  if (WHICHLAYERSLCPRCP[ilayerreverse]==0) {
    
    vector<int> vleft;
    vector<int> vright;

    if (ipolar==0) {
      vleft=v1l0_realtime_left;
      vright=v1l0_realtime_right;
    }
    else if (ipolar==1) {
      vleft=v2l0_realtime_left;
      vright=v2l0_realtime_right;

    }

    //    if (iphi==13 && IZERO>=40 && IZERO<=44) {
    //       cout << "start, stop are " << IZERO-nstepback << "\t" << IZERO-nstepback+(int)(L1_COINCIDENCE_LR[ilayerreverse]/TRIGTIMESTEP) << "\n";
    //       cout << "vleft is  ";
    //       for (int i=IZERO-nstepback ;i<=IZERO-nstepback+(int)(L1_COINCIDENCE_LR[ilayerreverse]/TRIGTIMESTEP);i++) {
    // 	cout << vleft[i] << " ";
    //       }
    //       cout << "\n";
    //      cout << "vright is  ";
    //       for (int i=IZERO-nstepback ;i<=IZERO-nstepback+(int)(L1_COINCIDENCE_LR[ilayerreverse]/TRIGTIMESTEP);i++) {
    // 	cout << vright[i] << " ";
    //       }
    //       cout << "\n";
    // //       cout << vleft[i] << " ";
    // //       cout << "left findahit is " << findahit(vleft,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_LR[ilayerreverse]/TRIGTIMESTEP)) << "\n";
    // //       cout << "right findahit is " << findahit(vleft,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_LR[ilayerreverse]/TRIGTIMESTEP)) << "\n";
    //     }
    if (findahit(vleft,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_LR_SCA[ilayerreverse]/TRIGTIMESTEP)) ||
	findahit(vright,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_LR_SCA[ilayerreverse]/TRIGTIMESTEP))) {
      
      vl1_realtime.push_back(1);
      //      if (IZERO>=40 && IZERO<=44)
      //cout << "same pol is triggered in " << ipolar << "\n";
    }
    else
      vl1_realtime.push_back(0);
  }
  else {
    if ( (findahit(v1l0_realtime_left,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_LR_SCA[ilayerreverse]/TRIGTIMESTEP)) &&
	  findahit(v2l0_realtime_left,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_LR_SCA[ilayerreverse]/TRIGTIMESTEP))) ||
	 (findahit(v1l0_realtime_right,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_LR_SCA[ilayerreverse]/TRIGTIMESTEP)) &&
	  findahit(v2l0_realtime_right,IZERO-nstepback,IZERO-nstepback+(int)(L1_COINCIDENCE_LR_SCA[ilayerreverse]/TRIGTIMESTEP))) ) {
      vl1_realtime.push_back(1);
      //if (IZERO>=40 && IZERO<=44)
      //cout << "lcp rcp is triggered in " << ipolar << "\n";

    }
    else
      vl1_realtime.push_back(0);

  }
  if (vl1_realtime[vl1_realtime.size()-1]) {
    //if (IZERO>=40 && IZERO<=44)
    //cout << "IZERO, ilayerreverse are " << IZERO << "\t" << ilayerreverse << " got a 1.\n";
    return 1;
  }
  //  if (IZERO>=40 && IZERO<=44)
  //cout << "returning zero.  IZERO, ilayerreverse are " << IZERO << "\t" << ilayerreverse << "\n";

  return 0;

}
void GlobalTrigger::L3Anita4LR_ScB_OneBin(int IZERO,Anita *anita1,std::array<std::array<vector<int>,3>,16> vl2_realtime_anita4_scb,
					  std::array<vector<int>,16> &vl3trig_type0, std::array<vector<int>,16> &vl3trig_type1,
					  int &thispasses_l3type0,int &thispasses_l3type1) {



  for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {


    int iphi_next=(iphi+1+16)%16;
    int iphi_previous=(iphi-1+16)%16;

    //    cout << "iphi is " << iphi << "\n";
    //    cout << "calling the first l3.\n";
    if (L3or30Anita4LR_ScB_TwoPhiSectors_OneBin(IZERO, 
						vl2_realtime_anita4_scb[iphi], // 3 neighbors, whether 1, 2 or 3 pass
						vl2_realtime_anita4_scb[iphi_previous], // 3 neighbors, whether 1, 2 or 3 pass
				     
						2,2) ||
	L3or30Anita4LR_ScB_TwoPhiSectors_OneBin(IZERO, 
						vl2_realtime_anita4_scb[iphi], // 3 neighbors, whether 1, 2 or 3 pass
						vl2_realtime_anita4_scb[iphi_next], // 3 neighbors, whether 1, 2 or 3 pass
						
						2,2)

	) {

      thispasses_l3type1=1;
      vl3trig_type1[iphi].push_back(1);

    }
    else
      vl3trig_type1[iphi].push_back(0);

    //cout << "calling the second l3.\n";
    if ((iphi%2==0 && L3or30Anita4LR_ScB_TwoPhiSectors_OneBin( IZERO,
							       vl2_realtime_anita4_scb[iphi], // 3 neighbors, whether 1, 2 or 3 pass
							       vl2_realtime_anita4_scb[iphi_next], // 3 neighbors, whether 1, 2 or 3 pass
							      
							       1,3)	  
	 ) ||
	(
	 iphi%2==1 && L3or30Anita4LR_ScB_TwoPhiSectors_OneBin( IZERO,
							       vl2_realtime_anita4_scb[iphi], // 3 neighbors, whether 1, 2 or 3 pass
							       vl2_realtime_anita4_scb[iphi_previous], // 3 neighbors, whether 1, 2 or 3 pass
							       
							       1,3)

	 )) {
      thispasses_l3type0=1;
      vl3trig_type0[iphi].push_back(1);
      
    }



    else
      vl3trig_type0[iphi].push_back(0);
  }

}


void GlobalTrigger::L2Anita4LR_ScB_OnePhiSector_OneBin(int IZERO,vector<int> vl1_bottom, 
						       vector<int> vl1_middle,
						       vector<int> vl1_top,
						       std::array<vector<int>,3> &vl2_realtime_anita4_scb,int &npassesl2,int &npassesl2_type0) {
  // keep track of whether you get a coincidence between 1, 2 or 3 antennas in a phi sector with the right windows.
  

 
  // If any of them pass l1, then the 0th element of the vpartofl2_realtime_anita4_scb array goes to one.
  //cout << "doing the one-hits.\n";
  

  if (vl1_bottom[IZERO] || vl1_middle[IZERO] || vl1_top[IZERO]) {
    vl2_realtime_anita4_scb[0].push_back(1);
    //cout << "got a one-hit.\n";
  }
  else
    vl2_realtime_anita4_scb[0].push_back(0);
  
  //cout << "doing the two-hits.\n";
  // If you get a coincidence betw. any two, then the 1st element of the vl2_realtime_anita4_scb array goes to one.
  //double time_thisbin=(double)nstepback*TRIGTIMESTEP;
  //int itrigbin=nstepback;

  //  cout << "sizes are " << vl1_bottom.size() << "\t" << vl1_middle.size() << "\t" << vl1_top.size() << "\n";
 

  //  while (time_thisbin<LASTTIMETOTESTL2_ANITA4LR_SCB) {


  if ((vl1_bottom[IZERO] &&
       (	findahit(vl1_middle,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[1]/TRIGTIMESTEP)) ||
		findahit(vl1_top,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[0]/TRIGTIMESTEP)) )) ||
	
      (vl1_middle[IZERO] &&
       (	findahit(vl1_bottom,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[1]/TRIGTIMESTEP)) ||
		findahit(vl1_top,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[2]/TRIGTIMESTEP)))) ||
	
      (vl1_top[IZERO] &&
       (	findahit(vl1_bottom,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[0]/TRIGTIMESTEP)) ||
 		findahit(vl1_middle,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[2]/TRIGTIMESTEP))))
      ) {
    vl2_realtime_anita4_scb[1].push_back(1);
    npassesl2++;
    //      cout << "start, stop bins  are " << IZERO-nstepback << "\t" << IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[1]/TRIGTIMESTEP) << "\n";
    //cout << "sizes are " << vl1_bottom.size() << "\t" << vl1_middle.size() << "\t" << vl1_top.size() << "\n";
    //cout << "got a two-hit.\n";
  }
  else
    vl2_realtime_anita4_scb[1].push_back(0);
  //itrigbin++;
  //time_thisbin=(double)itrigbin*TRIGTIMESTEP;
  //}

  //cout << "doing the three-hits.\n";
  // If you get a coincidence betw. all three, then the 2nd element of the vl2_realtime_anita4_scb array goes to one.
  //time_thisbin=(double)nstepback*TRIGTIMESTEP;
  //IZERO=nstepback;

  //  while (time_thisbin<LASTTIMETOTESTL2_ANITA4LR_SCB) {

  if ((vl1_bottom[IZERO] &&
       (	findahit(vl1_middle,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[1]/TRIGTIMESTEP)) &&
		findahit(vl1_top,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[0]/TRIGTIMESTEP)) )) ||
	
      (vl1_middle[IZERO] &&
       (	findahit(vl1_bottom,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[1]/TRIGTIMESTEP)) &&
		findahit(vl1_top,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[2]/TRIGTIMESTEP)))) ||
	
      (vl1_top[IZERO] &&
       (	findahit(vl1_bottom,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[0]/TRIGTIMESTEP)) &&
 		findahit(vl1_middle,IZERO-nstepback,IZERO-nstepback+(int)(L2_COINCIDENCE_ANITA4LR_SCB[2]/TRIGTIMESTEP))))
      ) {
    vl2_realtime_anita4_scb[2].push_back(1);
    npassesl2_type0++;
  }
  else
    vl2_realtime_anita4_scb[2].push_back(0);
  //itrigbin++;
  //time_thisbin=(double)itrigbin*TRIGTIMESTEP;
  //}  
}
// ask if L3 type 1 (2 and 2) or L3 type 0 (3 and 1 or 1 and 3) passes
int GlobalTrigger::L3or30Anita4LR_ScB_TwoPhiSectors_OneBin( int IZERO,
							    std::array<vector<int>,3> vl2_realtime_anita4_scb, // 3 neighbors, whether 1, 2 or 3 pass
							    std::array<vector<int>,3> vl2_realtime_anita4_scb_other, // 3 neighbors, whether 1, 2 or 3 pass
							    int npass1,int npass2) {

  

  int passaone=0;
  //double time_thisbin=(double)nstepback*TRIGTIMESTEP;
  //int itrigbin=nstepback;
  int thispasses=0;
  //while (time_thisbin<LASTTIMETOTESTL1_ANITA4LR_SCB) {   
  thispasses=0;


  if ((vl2_realtime_anita4_scb[npass1-1][IZERO] &&
       findahit(vl2_realtime_anita4_scb_other[npass2-1],IZERO-nstepback,IZERO-nstepback+(int)(L3_COINCIDENCE_ANITA4LR_SCB/TRIGTIMESTEP)))   ||
      (vl2_realtime_anita4_scb_other[npass2-1][IZERO] &&
       findahit(vl2_realtime_anita4_scb[npass1-1],IZERO-nstepback,IZERO-nstepback+(int)(L3_COINCIDENCE_ANITA4LR_SCB/TRIGTIMESTEP)))) {
    thispasses=1;
    passaone=1;
    //cout << "got an l3 with " << npass1 << " and " << npass2 << "\n";
  }

  return passaone;
  
}

void GlobalTrigger::delayL0(vector<int> &vl0,double delay) {
  int ndelay=(int)(delay/TRIGTIMESTEP);
  for (int i=0;i<ndelay;i++) {
    vl0.insert(vl0.begin(),0);
    vl0.erase(vl0.end()-1);
  }
  
}
void GlobalTrigger::delay_AllAntennas(Anita *anita1) {
  for (int itriglayer=0;itriglayer<anita1->NTRIGGERLAYERS;itriglayer++) {
    for (int iphi=0;iphi<anita1->PHITRIG[0];iphi++) {
      for (int ipolar=0;ipolar<anita1->NPOL;ipolar++) {
     
      
	delayL0(arrayofhits[itriglayer][iphi][ipolar][4],DELAYS[itriglayer]);
      }
    }
  }



}


