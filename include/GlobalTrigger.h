#ifndef GLOBALTRIGGER_HH_
#define GLOBALTRIGGER_HH_

////////////////////////////////////////////////////////////////////////////////////////////////
//class GlobalTrigger:
////////////////////////////////////////////////////////////////////////////////////////////////
#include "vector.hh"

namespace icemc{

  class Settings;
  class Anita;  

  //!  Global Trigger
  class GlobalTrigger {
    
  private:
    
  public:
    GlobalTrigger(Settings *settings1,Anita *anita1);
    //  GlobalTrigger(Settings *settings1,Anita *anita1,Balloon* bn1);
    void GetArrivalTimes(int inu,Anita *anita1, const Vector &rf_direction);
  
    // these are not really used now that we bin in frequency, but we keep them anyway.
  
    double TRIGTIMESTEP;
    int nstepback;  // when a trigger is fired, how many bins you step back to start a coincidence window

    // this is used for Anita 3
    //  const double L1_COINCIDENCE[3]={8.E-9,8.E-9,8.E-9}; // L1 coincidence window, in seconds  
    double L1_COINCIDENCE_ANITA3[3]; // L1 coincidence window, in seconds  
    // in this scenario B->M is the same as M->B for example
    //  const double L3_COINCIDENCE=22.5E-9; // L3 is now among neighboring phi sectors  
    double LASTTIMETOTESTL1_ANITA3; // can't test L1 after this point because the l1_coincidence windows go past the end of the waveform.
    double L3_COINCIDENCE; // L3 is now among neighboring phi sectors.  For Anita 3 and Anita 4 where there is an l3 coincidence  

    // used for an Anita 4 scenario where B->M doesn't have the same window as M->B necessarily
    double L1_COINCIDENCE_MOREGENERAL[3][2]; // L1 coincidence window, in seconds  
    // in this scenario B->M is *not* the same as M->B for example
    double LASTTIMETOTESTL1_ANITA4;
    double LASTTIMETOTESTL1_ANITA4LR_SCA;



    // In LR scenario A, a bottom antenna initiates the trigger, then require coincidences with hits in other layers
    double L1_COINCIDENCE_LR_SCA[2]; // L1 coincidence window, in seconds  
    // which layers we're considering LCP, RCP polarizations instead of V,H in scenario A
    double WHICHLAYERSLCPRCP[Anita::NTRIGGERLAYERS_MAX];


    // This if a coincidence between LCP and RCP for a single antenna for Anita 4
    double L1_COINCIDENCE_ANITA4LR_SCB; // L1 coincidence window, in seconds  
    double LASTTIMETOTESTL1_ANITA4LR_SCB;

    double L2_COINCIDENCE_ANITA4LR_SCB[3]; // L2 coincidence window, in seconds  
    double LASTTIMETOTESTL2_ANITA4LR_SCB;

    double L3_COINCIDENCE_ANITA4LR_SCB; // L3 coincidence window, in seconds  

    double DELAYS[3]; // delay top, middle, bottom

    vector<int> flag_e_L1[Anita::NPHI_MAX];
    vector<int> flag_h_L1[Anita::NPHI_MAX];
  
    UShort_t phiTrigMask[Anita::NPOL]; // which phi sector is masked for Anita 2 or Anita-3 (V-POL)
    //UShort_t phiTrigMaskH; // which phi sector is masked for Anita 3 (H-POL)
    UShort_t l1TrigMask[Anita::NPOL]; // which channel is masked for Anita-3 (V-POL)
    //UShort_t l1TrigMaskH; // which channel is masked for Anita 3 (H-POL)

    UShort_t thresholds_eachant[2][48]; // thresholds as read from the surf file: first index is pol, second is antenna number (only working for Anita3)

    double volts[2][Anita::NLAYERS_MAX][Anita::NPHI_MAX];                        // voltage (1st index=antenna,2nd index=pol., lcp=0. rcp=1)
    double volts_em[2][Anita::NLAYERS_MAX][Anita::NPHI_MAX];                        // component of voltage from em shower (1st index=antenna,2nd index=pol., lcp=0. rcp=1)
    double volts_original[2][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //added djg
    // e and h polarizations for each antenna summed across bandwidth
    
    
    int nchannels_perrx_triggered[48]; //Records number of first level triggers on each antenna for a single neutrino
    
    int nchannels_perband_triggered[48][8];  //Records individual channel triggers for each antenna.  (Either a 0 or a 1.) Channels 0-3 are lcp, 4-7 are rcp.
    
    
    
    
    
    // which channels on the payload pass.  
    // 1st component:  layers on the payload
    // 2nd component:  counting around the payload in phi
    // 3rd component:  which polarization
    // 4th component:  slices of bandwidth
    //  int channels_passing[4][16][2][5]; // keeps track of which channels pass
    int channels_passing[Anita::NLAYERS_MAX][Anita::NPHI_MAX][2][Anita::NBANDS_MAX]; // keeps track of which channels pass
    // make this an array of vectors instead so that we can have an arbitrary number of bands for satellite
    vector<int> vchannels_passing[Anita::NLAYERS_MAX][Anita::NPHI_MAX][2];
    std::array< std::array< std::array< std::array<std::vector<int>,5>, 2>, 16>, 3>  arrayofhits; 


    int triggerbits[Anita::NTRIG]; // keeps track of which trigger scenarios pass
    // for the nadir studies
    
    // this is L2 and L3 triggers
    void PassesTrigger(Settings *settings1,Anita *anita1,int discones_passing,int mode,int *l3trig,int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int antennaclump,int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX],int inu,
		       int *thispasses);
    void PassesTrigger(Settings *settings1,Anita *anita1,int discones_passing,int mode,int *l3trig,int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int antennaclump,int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX],int inu,double this_threshold,
		       int *thispasses);

  
    void PassesTriggerBasic(Settings *settings1,Anita *anita1,int discones_passing,int mode,int *l3trig,int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX],int antennaclump,int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX], int *thispasses, int inu);

    void PassesTriggerCoherentSum(Settings *settings1, Anita *anita1, int inu, int *thispasses);

    void PassesTriggerSummedPower(Settings *settings1,Anita *anita1);

    void PassesTriggerScheme5(Anita *anita1, double this_threshold, int *thispasses);
  
  
    void L3Trigger(Settings *settings1,Anita *anita1,int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX],int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX],int discones_passing,int *l3trigy,
		   int *thispasses);
    
    
    int GetPhiSector(Settings *settings1,int i,int j); // given trigger layer and index, get phi sector.
    // for the upper two layers, the phi sector is just j
    // for the nadir layer, the phi sector is 2*j+1
    void GetAnitaLayerPhiSector(Settings *settings1,int i,int j,int &whichlayer,int &whichphisector);
    void FillInNadir(Anita *anita1,int *ant);
    void FillInNadir(Settings *settings1,Anita *anita1,int ant);
    
    // The following functions are related to the coherent-sum trigger
    vector < vector < vector <double> > > volts_rx_rfcm_trigger; // This is used for easy access to the waveforms of specific phi sectors and layers, and it combines physical layers into their trigger layer.
    // Accessed by volts_rx_rfcm_trigger[phi_sector][layer][timestep]
    int first_phi_sector_hit; // This is used by the coherent waveform sum trigger scheme to make processing more efficient
    double three_bit_round(double input, bool round_zero_up = true, bool allow_zero = false);
    void convert_wfm_to_3_bit(const vector <double>& wfm, double rms, vector <double>& output);
    void delay_align_antenna_waveforms(const vector< vector < vector <double> > >& waveforms, const vector < vector <unsigned int> >& delays, vector < vector <double> >& output);
    void sum_aligned_waveforms(const vector < vector <double> >& waveforms, vector <double>& output);
    void square_waveform_elements(const vector <double>& waveform, vector <double>& output);
    double summed_power_window(const vector <double>& waveform, unsigned int start_index, unsigned int length);
    // End of functions relating to coherent-sum trigger

    int findahit(vector<int> myvector,int first,int last);
    int findanl3(int *l3,int NPHISECTORS);

    int L1Anita3_OnePhiSector(int IZERO,vector<int> &vl0_realtime_bottom, vector<int> &vl0_realtime_middle, vector<int> &vl0_realtime_top,
			      vector<int> &vl1_realtime_bottom, vector<int> &vl1_realtime_middle, vector<int> &vl1_realtime_top);

    void L1Anita3_AllPhiSectors(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> &l1trig);  


    int L1Anita4_OnePhiSector(int IZERO,vector<int> &vl0_realtime_bottom, vector<int> &vl0_realtime_middle, vector<int> &vl0_realtime_top,
			      vector<int> &vl1_realtime_bottom, vector<int> &vl1_realtime_middle, vector<int> &vl1_realtime_top);
    void L1Anita4_AllPhiSectors(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> &l1trig);

    void L2Anita3and4(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> l1trig,
		      std::array<std::array<std::vector<int>,16>,2> &l2trig);

    void L3Anita3and4(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> vl2trig,
		      int vl3trig[2][16],int *thispasses);


    int PartofL1Anita4LR_ScA_TwoPhiSectors(int ilayerreverse,int IZERO,int ipolar,
					   vector<int> &v1l0_realtime_left, vector<int> &v2l0_realtime_left, 
					   vector<int> &v1l0_realtime_right, vector<int> &v2l0_realtime_right, 
					   vector<int> &vl1_realtime);
  
    int L1Anita4LR_ScA_TwoPhiSectors(int IZERO,int ipolar,
				     vector<int> &vl0_realtime_bottomleft, vector<int> &v2l0_realtime_bottomleft, 
				     vector<int> &vl0_realtime_bottomright, vector<int> &v2l0_realtime_bottomright, 
				     vector<int> &vl0_realtime_middleleft, vector<int> &v2l0_realtime_middleleft,
				     vector<int> &vl0_realtime_middleright, vector<int> &v2l0_realtime_middleright,
				     vector<int> &vl0_realtime_topleft, vector<int> &v2l0_realtime_topleft,
				     vector<int> &vl0_realtime_topright, vector<int> &v2l0_realtime_topright,
				     vector<int> &vl1_realtime_bottom, 
				     vector<int> &vl1_realtime_middle, 
				     vector<int> &vl1_realtime_top);
  
    void L1Anita4LR_ScA_AllPhiSectors(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> &vl1trig);
    void L3Anita4LR_ScA(Anita *anita1,std::array<std::array<std::vector<int>,16>,2> vl2trig,
			int *thispasses);




    void L1Anita4LR_ScB_AllAntennas_OneBin(int IZERO,Anita *anita1,std::array< std::array< vector<int>,16>,3> &vl1trig_anita4lr_scb,int &npassesl1);
    // L1 trigger is at the antenna level again.  Just require coincidence between LCP and RCP
    void L1Anita4LR_ScB_OneBin(int IZERO,vector<int> vleft,vector<int> vright,
			       vector<int> &vl1trig);

    void L2Anita4LR_ScB_AllPhiSectors_OneBin(int IZERO,Anita *anita1,std::array< std::array< vector<int>,16>,3> vl1trig_anita4lr_scb,
					     std::array<std::array<vector<int>,3>,16> &vl2_realtime_anita4_scb, int &npassesl2, int &npassesl2_type0);

  
    // keep track of whether you get a coincidence between 1, 2 or 3 antennas in a phi sector with the right windows.

    void L2Anita4LR_ScB_OnePhiSector_OneBin(int IZERO,vector<int> vl1_bottom, 
					    vector<int> vl1_middle,
					    vector<int> vl1_top,
					    std::array<vector<int>,3> &vl2_realtime_anita4_scb,int &npassesl2,int &npassesl2_type0);

    int L3or30Anita4LR_ScB_TwoPhiSectors_OneBin(int IZERO, 
						std::array<vector<int>,3> vl2_realtime_anita4_scb, // 3 neighbors, whether 1, 2 or 3 pass
						std::array<vector<int>,3> vl2_realtime_anita4_scb_other, // 3 neighbors, whether 1, 2 or 3 pass
						int npass1,int npass2);
  
    void L3Anita4LR_ScB_OneBin(int IZERO,Anita *anita1,std::array<std::array<vector<int>,3>,16> vl2_realtime_anita4_scb,
			       std::array<vector<int>,16> &vl3trig_type0, std::array<vector<int>,16> &vl3trig_type1,
			       int &thispasses_l3type0,int &thispasses_l3type1);

  
    void delayL0(vector<int> &vl0,double delay);
    void delay_AllAntennas(Anita *anita1);




  };
}


#endif
