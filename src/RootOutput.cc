#include "RootOutput.h"
#include <iostream>

icemc::RootOutput::RootOutput() : fOutputDir(""), fRun(0){
  zeroFilePointers();
}

icemc::RootOutput::RootOutput(const char* outputDir, int run) : fOutputDir(outputDir), fRun(run) {
  zeroFilePointers();
}



icemc::RootOutput::~RootOutput(){

  if(fIceFinalFile){
    fIceFinalFile->Write();
    fIceFinalFile->Close();
    fIceFinalFile = NULL;
    fTree2 = NULL;
  }
}


void icemc::RootOutput::zeroFilePointers(){
  fIceFinalFile = NULL;
  fTree2 = NULL;
  // fTree2Output = NULL;
}


void icemc::RootOutput::make_icefinal(){

  if(fIceFinalFile){
    std::cerr << "Ice final already initialized!" << std::endl;
    return;
  }
  
  TString fileName = fOutputDir + TString::Format("/icefinal%d.root", fRun);
  fIceFinalFile = new TFile(fileName, "RECREATE", "ice");

  fTree2 = new TTree("h2000", "h2000"); // tree2 filled for each event that is beyond the horizon.
  fTree2->Branch("h2000", &fTree2Output);

  // tree2->Branch("inu", &inu, "inu/I");
  // tree2->Branch("horizcoord", &horizcoord, "horizcoord/D");
  // tree2->Branch("vertcoord", &vertcoord, "vertcoord/D");
  // tree2->Branch("scalefactor_distance", &scalefactor_distance, "scalefactor_distance/D");
  // tree2->Branch("scalefactor_attenuation", &scalefactor_attenuation, "scalefactor_attenuation/D");


  // TTree *tree3 = new TTree("h3000", "h3000"); // tree3 if signal is detectable.
  // tree3->Branch("deltheta_em_max", &deltheta_em_max, "deltheta_em_max/D");
  // tree3->Branch("deltheta_had_max", &deltheta_had_max, "deltheta_had_max/D");
  // tree3->Branch("theta_threshold_deg", &theta_threshold_deg, "theta_threshold_deg/D");
  // tree3->Branch("nsigma_em_threshold", &nsigma_em_threshold, "nsigma_em_threshold/D");
  // tree3->Branch("nsigma_had_threshold", &nsigma_had_threshold, "nsigma_had_threshold/D");
  // tree3->Branch("horizcoord", &horizcoord, "horizcoord/D");
  // tree3->Branch("vertcoord", &vertcoord, "vertcoord/D");
  // tree3->Branch("vmmhz_max", &vmmhz_max, "vmmhz_max/D");
  // tree3->Branch("vmmhz_min", &vmmhz_min, "vmmhz_min/D");
  // tree3->Branch("dviewangle_deg", &dviewangle_deg, "dviewangle_deg/D");
  // tree3->Branch("viewangle_deg", &viewangle_deg, "viewangle_deg/D");
  // tree3->Branch("changle_deg", &changle_deg, "changle_deg/D");
  // tree3->Branch("cosviewangle", &cosviewangle, "cosviewangle/D");
  // tree3->Branch("emfrac", &emfrac, "emfrac/D");
  // tree3->Branch("hadfrac", &hadfrac, "hadfrac/D");

  // TTree *tree5 = new TTree("h5000", "h5000"); // tree5 filled for each nutau.
  // tree5->Branch("vmmhz1m_max", &vmmhz1m_max, "vmmhz1m_max/D");
  // tree5->Branch("inu", &inu, "inu/I");
  // tree5->Branch("nuexitlength", &nuexitlength, "nuexitlength/D");
  // tree5->Branch("nuexitice",  &nuexitice,  "nuexitice");
  // tree5->Branch("vmmhz_max", &vmmhz_max, "vmmhz_max");
  // tree5->Branch("maxtaper", &maxtaper, "maxtaper");
  // tree5->Branch("inu", &inu, "inu/I");
  // tree5->Branch("whichray", &whichray, "whichray/I");
  // tree5->Branch("pnu", &pnu, "pnu/D");
  // tree5->Branch("costhetanu", &costhetanu, "costhetanu/D");
  // tree5->Branch("viewangle", &viewangle, "viewangle/D");
  // tree5->Branch("offaxis", &offaxis, "offaxis/D");
  // tree5->Branch("nsigma_offaxis", &nsigma_offaxis, "nsigma_offaxis/D");
  // tree5->Branch("hadfrac", &hadfrac, "hadfrac/D");
  // tree5->Branch("emfrac", &emfrac, "emfrac/D");
  // tree5->Branch("sumfrac", &sumfrac, "sumfrac/D");
  // tree5->Branch("horizcoord", &horizcoord, "horizcoord/D");
  // tree5->Branch("vertcoord", &vertcoord, "vertcoord/D");
  // tree5->Branch("weight1", &weight1, "weight1/D");
  // tree5->Branch("nearthlayers", &nearthlayers, "nearthlayers/D");
  // tree5->Branch("logchord", &logchord2, "interaction1->logchord/D");
  // tree5->Branch("diff_3tries", &diff_3tries, "diff_3tries/D");
  // tree5->Branch("fresnel2", &fresnel2, "fresnel2/D");
  // tree5->Branch("costheta_inc", &costheta_inc, "costheta_inc/D");
  // tree5->Branch("costheta_exit", &costheta_exit, "costheta_exit/D");
  // tree5->Branch("deltheta_em", &deltheta_em[0], "deltheta_em/D");
  // tree5->Branch("deltheta_had", &deltheta_had[0], "deltheta_had/F");
  // tree5->Branch("r_fromballoon", &interaction1->r_fromballoon[0], "r_fromballoon/F");
  // tree5->Branch("theta_in", &theta_in, "theta_in/D");
  // tree5->Branch("lat_in", &lat_in, "lat_in/D");


  // TTree *tree6 = new TTree("h6000", "h6000"); // tree6 filled for neutrinos that enter S of 60 deg S latitude.
  // tree6->Branch("volts_rx_0", &volts_rx_0, "volts_rx_0/D");
  // tree6->Branch("volts_rx_1", &volts_rx_1, "volts_rx_1/D");
  // tree6->Branch("horizcoord", &horizcoord, "horizcoord/D");
  // tree6->Branch("vertcoord", &vertcoord, "vertcoord/D");
  // tree6->Branch("theta_in", &theta_in, "theta_in/D");
  // tree6->Branch("chord_kgm2_bestcase", &chord_kgm2_bestcase2, "chord_kgm2_bestcase/D");
  // tree6->Branch("chord_kgm2_ice", &interaction1->chord_kgm2_ice, "chord_kgm2_ice/D");
  // tree6->Branch("costheta_nutraject", &interaction1->costheta_nutraject, "costheta_nutraject/D");
  // tree6->Branch("weight1", &weight1, "weight1/D");
  // tree6->Branch("weight_bestcase", &weight_bestcase2, "weight_bestcase/D");
  // tree6->Branch("whichray", &whichray, "whichray/I");
  // tree6->Branch("mybeta", &mybeta, "mybeta/D");
  // tree6->Branch("longitude", &longitude_this, "longitude/D");


  // TTree *tree6b = new TTree("h6001", "h6001"); // tree6b filled for the closest antenna to the interaction
  // tree6b->Branch("bwslice_vnoise", bwslice_vnoise_thislayer, "bwslice_vnoise_thislayer[4]/D");

  // TTree *tree7 = new TTree("h7000", "h7000"); // tree6 filled just after flavor is set
  // tree7->Branch("emfrac", &emfrac, "emfrac/D");
  // tree7->Branch("hadfrac", &hadfrac, "hadfrac/D");
  // tree7->Branch("current", &interaction1->currentint, "currentint/I");
  // tree7->Branch("nuflavor", &interaction1->nuflavorint, "nuflavorint/I");
  // tree7->Branch("sumfrac", &sumfrac, "sumfrac/D");
  // tree7->Branch("slopeyangle", &slopeyangle, "slopeyangle/D");

  // TTree *jaimetree=new TTree("jaimetree", "jaimetree"); // signal as it is produced at the interaction
  // jaimetree->Branch("vmmhz1m_max", &vmmhz1m_max, "vmmhz1m_max/D");
  // jaimetree->Branch("emfrac", &emfrac, "emfrac/D");
  // jaimetree->Branch("hadfrac", &hadfrac, "hadfrac/D");
  // jaimetree->Branch("deltheta_em_max", &deltheta_em_max, "deltheta_em_max/D");
  // jaimetree->Branch("deltheta_had_max", &deltheta_had_max, "deltheta_had_max/D");
  // jaimetree->Branch("sumfrac", &sumfrac, "sumfrac/D");
  // jaimetree->Branch("vmmhz1m_visible", &vmmhz1m_visible, "vmmhz1m_visible/D");

  // TTree *viewangletree=new TTree("viewangletree", "viewangletree"); // signal as it is produced at the interaction
  // viewangletree->Branch("dviewangle_deg", &dviewangle_deg, "dviewangle_deg/D");
  // viewangletree->Branch("emfrac", &emfrac, "emfrac/D");
  // viewangletree->Branch("hadfrac", &hadfrac, "hadfrac/D");
  // viewangletree->Branch("deltheta_em_max", &deltheta_em_max, "deltheta_em_max/D");
  // viewangletree->Branch("deltheta_had_max", &deltheta_had_max, "deltheta_had_max/D");
  // viewangletree->Branch("theta_threshold_deg", &theta_threshold_deg, "theta_threshold_deg/D");
  // viewangletree->Branch("dnutries", &interaction1->dnutries, "dnutries/D");
  // viewangletree->Branch("viewangle", &viewangle, "viewangle/D");
  // viewangletree->Branch("chord", &interaction1->chord, "dnutries/D");

  // TTree *neutrino_positiontree=new TTree("neutrino_positiontree", "neutrino_positiontree");
  // neutrino_positiontree->Branch("nnu", &interaction1->nnu, "nnu[3]/D");
  // neutrino_positiontree->Branch("dtryingdirection", &interaction1->dtryingdirection, "dtryingdirection/D");
  // neutrino_positiontree->Branch("bn1->dtryingposition", &bn1->dtryingposition, "bn1->dtryingposition/D");

  // //Filled just after Getchord,  where we find the neutrino's path through the Earth
  // TTree *nupathtree=new TTree("nupathtree", "nupathtree");
  // nupathtree->Branch("total_kgm2", &total_kgm2, "total_kgm2/D");
  // nupathtree->Branch("chord", &interaction1->chord, "chord/D");
  // nupathtree->Branch("crust_entered", &crust_entered, "crust_entered/I");
  // nupathtree->Branch("mantle_entered", &mantle_entered, "mantle_entered/I");
  // nupathtree->Branch("core_entered", &core_entered, "core_entered/I");
  // nupathtree->Branch("mybeta", &mybeta, "mybeta/D");
  // nupathtree->Branch("costheta_nutraject", &interaction1->costheta_nutraject, "costheta_nutraject/D");

  // TTree *finaltree = new TTree("passing_events", "passing_events"); // finaltree filled for all events that pass
  // finaltree->Branch("inu", &inu, "inu/I");
  // finaltree->Branch("vmmhz_min", &vmmhz_min, "vmmhz_min/D");
  // finaltree->Branch("vmmhz_max", &vmmhz_max, "vmmhz_max/D");
  // finaltree->Branch("thresholdsAnt", &thresholdsAnt, "thresholdsAnt[48][2][5]/D");
  // finaltree->Branch("thresholdsAntPass", &thresholdsAntPass, "thresholdsAntPass[48][2][5]/D");
  // finaltree->Branch("deadTime", &anita1->deadTime, "deadTime/D");
  // finaltree->Branch("horizcoord", &horizcoord, "horizcoord/D");
  // finaltree->Branch("vertcoord", &vertcoord, "vertcoord/D");
  // finaltree->Branch("horizcoord_bn", &bn1->horizcoord_bn, "horizcoord_bn/D");
  // finaltree->Branch("vertcoord_bn", &bn1->vertcoord_bn, "vertcoord_bn/D");
  // finaltree->Branch("r_bn", &r_bn_array, "r_bn_array[3]/D");
  // finaltree->Branch("n_bn", &n_bn_array, "n_bn_array[3]/D");
  // finaltree->Branch("longitude_bn", &longitude_this, "longitude_bn/D");
  // finaltree->Branch("heading_bn", &heading_this, "heading_bn/D");
  // finaltree->Branch("gps_offset", &gps_offset, "gps_offset/D");
  // // this one is just weight due to earth absorption
  // finaltree->Branch("weight1", &weight1, "weight1/D");
  // // this is the total weight - the one you want to use!
  // finaltree->Branch("weight", &weight, "weight/D");
  // finaltree->Branch("logweight", &logweight, "logweight/D");
  // finaltree->Branch("posnu", &posnu_array, "posnu_array[3]/D");
  // finaltree->Branch("costheta_nutraject", &costheta_nutraject2, "costheta_nutraject/D");
  // finaltree->Branch("chord_kgm2_ice",  &chord_kgm2_ice2, "chord_kgm2_ice/D");
  // finaltree->Branch("phi_nutraject", &phi_nutraject2, "phi_nutraject/D");
  // finaltree->Branch("altitude_int", &altitude_int2, "altitude_int/D");
  // finaltree->Branch("nnu", &nnu_array, "nnu_array[3]/D");
  // finaltree->Branch("n_exit2bn", &n_exit2bn_array, "n_exit2bn_array[5][3]/D");
  // finaltree->Branch("n_exit_phi", &n_exit_phi, "n_exit_phi/D");
  // finaltree->Branch("rfexit", &rfexit_array, "rfexit_array[5][3]/D");

  // finaltree->Branch("pnu", &pnu, "pnu/D");
  // finaltree->Branch("elast_y", &elast_y, "elast_y/D");
  // finaltree->Branch("emfrac", &emfrac, "emfrac/D");
  // finaltree->Branch("hadfrac", &hadfrac, "hadfrac/D");
  // finaltree->Branch("sumfrac", &sumfrac, "sumfrac/D");
  // finaltree->Branch("nuflavor", &nuflavorint2, "nuflavorint/I");//1=electron,  2=muon,  3=tau
  // finaltree->Branch("current", &currentint2, "currentint/I");//0=charged current,  1=neutral current
  // finaltree->Branch("logchord",  &logchord2,  "logchord/D");
  // finaltree->Branch("nuexitice",  &nuexitice,  "nuexitice/D");
  // finaltree->Branch("weight_bestcase",  &weight_bestcase2,  "weight_bestcase/D");
  // finaltree->Branch("chord_kgm2_bestcase",  &chord_kgm2_bestcase2,  "chord_kgm2_bestcase/D");
  // finaltree->Branch("dtryingdirection",  &dtryingdirection2,  "dtryingdirection/D");
  // finaltree->Branch("l3trig", &l3trig, "l3trig[2]/I");
  // finaltree->Branch("l2trig", &l2trig, "l2trig[2][3]/I");
  // finaltree->Branch("l1trig", &l1trig, "l1trig[2][3]/I");
  // finaltree->Branch("phiTrigMask", &anita1->phiTrigMask, "phiTrigMask/s");
  // finaltree->Branch("phiTrigMaskH", &anita1->phiTrigMaskH, "phiTrigMaskH/s");
  // finaltree->Branch("l1TrigMask", &anita1->l1TrigMask, "l1TrigMask/s");
  // finaltree->Branch("l1TrigMaskH", &anita1->l1TrigMaskH, "l1TrigMaskH/s");
  // finaltree->Branch("max_antenna0", &max_antenna0, "max_antenna0/I");
  // finaltree->Branch("max_antenna1", &max_antenna1, "max_antenna1/I");
  // finaltree->Branch("max_antenna2", &max_antenna2, "max_antenna2/I");

  // finaltree->Branch("viewangle", &viewangle, "viewangle/D");
  // finaltree->Branch("offaxis", &offaxis, "offaxis/D");
  // finaltree->Branch("rx0_signal_eachband", &rx0_signal_eachband, "rx0_signal_eachband[2][5]/D");
  // finaltree->Branch("rx0_threshold_eachband", &rx0_threshold_eachband, "rx0_threshold_eachband[2][5]/D");
  // finaltree->Branch("rx0_noise_eachband", &rx0_noise_eachband, "rx0_noise_eachband[2][5]/D");
  // finaltree->Branch("rx0_passes_eachband", &rx0_passes_eachband, "rx0_passes_eachband[2][5]/I");
  // finaltree->Branch("e_component", &e_component, "e_component/D");
  // finaltree->Branch("h_component", &h_component, "h_component/D");
  // finaltree->Branch("dist_int_bn_2d", &dist_int_bn_2d, "dist_int_bn_2d/D");
  // finaltree->Branch("d1", &d12, "d1/D");

  // finaltree->Branch("cosalpha", &cosalpha, "cosalpha/D");
  // finaltree->Branch("mytheta", &mytheta, "mytheta/D");
  // finaltree->Branch("cosbeta0", &cosbeta0, "cosbeta0/D");
  // finaltree->Branch("mybeta", &mybeta, "mybeta/D");
  // finaltree->Branch("d1", &d12, "d1/D");
  // finaltree->Branch("d2", &d22, "d2/D");

  // //Begin block added by Stephen for verification plots
  // finaltree->Branch("fresnel1", &fresnel1, "fresnel1/D");
  // finaltree->Branch("fresnel2", &fresnel2, "fresnel2/D");
  // finaltree->Branch("mag1", &mag1, "mag1/D");
  // finaltree->Branch("mag2", &mag2, "mag2/D");
  // finaltree->Branch("t_coeff_pokey", &t_coeff_pokey, "t_coeff_pokey/D");
  // finaltree->Branch("t_coeff_slappy", &t_coeff_slappy, "t_coeff_slappy/D");
  // finaltree->Branch("exponent", &settings1->EXPONENT, "EXPONENT/D");

  // finaltree->Branch("hitangle_e_all", &hitangle_e_all, "hitangle_e_all[48]/D");
  // finaltree->Branch("hitangle_h_all", &hitangle_h_all, "hitangle_h_all[48]/D");

  // finaltree->Branch("e_comp_max1", &e_comp_max1, "e_comp_max1/D");
  // finaltree->Branch("h_comp_max1", &h_comp_max1, "h_comp_max1/D");
  // finaltree->Branch("e_comp_max2", &e_comp_max2, "e_comp_max2/D");
  // finaltree->Branch("h_comp_max2", &h_comp_max2, "h_comp_max2/D");
  // finaltree->Branch("e_comp_max3", &e_comp_max3, "e_comp_max3/D");
  // finaltree->Branch("h_comp_max3", &h_comp_max3, "h_comp_max3/D");
  // finaltree->Branch("max_antenna_volts0", &max_antenna_volts0, "max_antenna_volts0/D");
  // finaltree->Branch("max_antenna_volts0_em", &max_antenna_volts0_em, "max_antenna_volts0_em/D");
  // finaltree->Branch("max_antenna_volts1", &max_antenna_volts1, "max_antenna_volts1/D");
  // finaltree->Branch("max_antenna_volts2", &max_antenna_volts2, "max_antenna_volts2/D");
  // finaltree->Branch("triggers", &nchannels_perrx_triggered, "nchannels_perrx_triggered[48]/I");
  // finaltree->Branch("nchannels_triggered", &nchannels_triggered, "nchannels_triggered/I");
  // finaltree->Branch("volts_rx_max", &volts_rx_max, "volts_rx_max/D");
  // finaltree->Branch("volts_rx_ave", &volts_rx_ave, "volts_rx_ave/D");
  // finaltree->Branch("volts_rx_sum", &volts_rx_sum, "volts_rx_sum/D");

  // finaltree->Branch("volts_rx_max_highband", &volts_rx_max_highband, "volts_rx_max_highband/D");
  // finaltree->Branch("volts_rx_max_lowband", &volts_rx_max_lowband, "volts_rx_max_lowband/D");
  // finaltree->Branch("theta_pol_measured", &theta_pol_measured, "theta_pol_measured/D");
  // finaltree->Branch("theta_rf_atbn", &theta_rf_atbn, "theta_rf_atbn/D");
  // finaltree->Branch("theta_rf_atbn_measured", &theta_rf_atbn_measured, "theta_rf_atbn_measured/D");
  // finaltree->Branch("voltage", &voltagearray, "voltagearray[48]/D");
  // finaltree->Branch("nlayers", &settings1->NLAYERS, "settings1->NLAYERS/I");

  // finaltree->Branch("vmmhz1m_max", &vmmhz1m_max, "vmmhz1m_max/D");
  // finaltree->Branch("vmmhz_lowfreq", &vmmhz_lowfreq, "vmmhz_lowfreq/D");

  // finaltree->Branch("deltheta_em_max", &deltheta_em_max, "deltheta_em_max/D");
  // finaltree->Branch("deltheta_had_max", &deltheta_had_max, "deltheta_had_max/D");
  // finaltree->Branch("r_enterice", &r_enterice_array, "r_enterice_array[3]/D");
  // finaltree->Branch("n_exit2bn_db", &n_exit2bn_db_array, "n_exit2bn_db_array[5][3]/D");

  // finaltree->Branch("rfexit_db", &rfexit_db_array, "rfexit_db_array[5][3]/D");
  // finaltree->Branch("r_in", &r_in_array, "r_in_array[3]/D");
  // finaltree->Branch("nsurf_rfexit", &nsurf_rfexit_array, "nsurf_rfexit_array[3]/D");
  // finaltree->Branch("nsurf_rfexit_db", &nsurf_rfexit_db_array, "nsurf_rfexit_db_array[3]/D");
  // finaltree->Branch("r_fromballoon", &r_fromballoon2, "r_fromballoon/D");
  // finaltree->Branch("r_fromballoon_db", &interaction1->r_fromballoon_db, "r_fromballoon_db/D");

  // finaltree->Branch("nuexitlength", &nuexitlength, "nuexitlength/D");
  // finaltree->Branch("nuentrancelength", &nuentrancelength, "nuentrancelength/D");
  // finaltree->Branch("taulength", &taulength, "taulength/D");
  // finaltree->Branch("icethickness", &icethickness, "icethickness/D");
  // finaltree->Branch("nrf_iceside", &nrf_iceside_array, "nrf_iceside_array[5][3]/D");
  // finaltree->Branch("nrf_iceside_db", &nrf_iceside_db_array, "nrf_iceside_db_array[5][3]/D");
  // finaltree->Branch("ant_normal0", &ant_max_normal0_array, "ant_max_normal0_array[3]/D");
  // finaltree->Branch("ant_normal1", &ant_max_normal1_array, "ant_max_normal1_array[3]/D");
  // finaltree->Branch("ant_normal2", &ant_max_normal2_array, "ant_max_normal2_array[3]/D");
  // finaltree->Branch("vmmhz1m_visible", &vmmhz1m_visible, "vmmhz1m_visible/D");
  // finaltree->Branch("freq_bins", &freq_bins, "freq_bins/I");
  // finaltree->Branch("vmmhz", &vmmhz, "vmmhz[freq_bins]/D");

  // finaltree->Branch("dist_int_bn_2d_chord", &dist_int_bn_2d_chord, "dist_int_bn_2d_chord/D");

  // finaltree->Branch("dviewangle_deg", &dviewangle_deg, "dviewangle_deg/D");
  // finaltree->Branch("theta_threshold_deg", &theta_threshold_deg, "theta_threshold_deg/D");
  // finaltree->Branch("total_kgm2", &total_kgm2, "total_kgm2/D");
  // finaltree->Branch("chord", &interaction1->chord, "chord/D");
  // finaltree->Branch("crust_entered", &crust_entered, "crust_entered/I");
  // finaltree->Branch("mantle_entered", &mantle_entered, "mantle_entered/I");
  // finaltree->Branch("core_entered", &core_entered, "core_entered/I");
  // finaltree->Branch("n_pol", &n_pol_array, "n_pol_array[3]/D");
  // finaltree->Branch("vmmhz_min_thatpasses", &vmmhz_min_thatpasses, "vmmhz_min_thatpasses/D");

  // finaltree->Branch("pieceofkm2sr", &pieceofkm2sr, "pieceofkm2sr/D");
  // //finaltree->Branch("volts_original", &volts_original, "volts_original[10][20][2]/D");
  // finaltree->Branch("r_exit2bn", &r_exit2bn2, "r_exit2bn/D");
  // finaltree->Branch("r_exit2bn_measured", &r_exit2bn_measured2, "r_exit2bn_measured/D");
  // finaltree->Branch("scalefactor_attenuation", &scalefactor_attenuation, "scalefactor_attenuation/D");
  // finaltree->Branch("anita1->PHI_OFFSET", &anita1->PHI_OFFSET, "anita1->PHI_OFFSET/D");
  // finaltree->Branch("igps", &bn1->igps, "igyps/I");
  // finaltree->Branch("volts_rx_rfcm_lab_e_all", &volts_rx_rfcm_lab_e_all, "volts_rx_rfcm_lab_e_all[48][512]/D");
  // finaltree->Branch("volts_rx_rfcm_lab_h_all", &volts_rx_rfcm_lab_h_all, "volts_rx_rfcm_lab_h_all[48][512]/D");
  // finaltree->Branch("ptaui", &ptaui, "ptaui/D");
  // finaltree->Branch("ptauf", &ptauf, "ptauf/D");
  // finaltree->Branch("sourceLon", &sourceLon, "sourceLon/D");
  // finaltree->Branch("sourceLat", &sourceLat, "sourceLat/D");
  // finaltree->Branch("sourceAlt", &sourceAlt, "sourceAlt/D");
  // finaltree->Branch("sourceMag", &sourceMag, "sourceMag/D");


  // TTree *mytaus_tree = new TTree("mytaus", "mytaus");
  // mytaus_tree->Branch("taus",  &TauPtr);

  // double rms_rfcm_e;
  // double rms_rfcm_h;
  // double rms_lab_e;
  // double rms_lab_h;

  // double avgfreq_rfcm[Anita::NFREQ];
  // double avgfreq_rfcm_lab[Anita::NFREQ];
  // double freq[Anita::NFREQ];

  // TTree *summarytree = new TTree("summarytree", "summarytree"); // finaltree filled for all events that pass
  // summarytree->Branch("NNU", &NNU, "NNU/I");
  // summarytree->Branch("EXPONENT", &settings1->EXPONENT, "EXPONENT/D");
  // summarytree->Branch("eventsfound_beforetrigger", &eventsfound_beforetrigger, "eventsfound_beforetrigger/D");
  // summarytree->Branch("rms_rfcm_e", &rms_rfcm_e, "rms_rfcm_e/D");
  // summarytree->Branch("rms_rfcm_h", &rms_rfcm_h, "rms_rfcm_h/D");
  // summarytree->Branch("rms_lab_e", &rms_lab_e, "rms_lab_e/D");
  // summarytree->Branch("rms_lab_h", &rms_lab_h, "rms_lab_h/D");
  // summarytree->Branch("avgfreq_rfcm", &avgfreq_rfcm, "avgfreq_rfcm[128]/D");
  // summarytree->Branch("avgfreq_rfcm_lab", &avgfreq_rfcm_lab, "avgfreq_rfcm_lab[128]/D");
  // summarytree->Branch("freq", &freq, "freq[128]/D");

  // TTree *banana_tree = new TTree("banana_tree", "banana_tree");  //To record banana plot info - Stephen
  // banana_tree->Branch("r_bn", &bn1->r_bn, "r_bn[3]/D");

  // TTree *ytree = new TTree("ytree", "ytree"); //To record y distributions
  // ytree->Branch("elast_y", &elast_y, "elast_y/D");

  // double icethck;
  // double elev;
  // double lon_ground;
  // double lat_ground;
  // double lon_ice;
  // double lat_ice;
  // double h20_depth;
  // double lon_water;
  // double lat_water;

  // TTree *icetree = new TTree("icetree", "icetree");
  // icetree->Branch("icethck", &icethck, "icethck/D");
  // icetree->Branch("lon_ice", &lon_ice, "lon_ice/D");
  // icetree->Branch("lat_ice", &lat_ice, "lat_ice/D");
  // icetree->Branch("lon_water", &lon_water, "lon_water/D");
  // icetree->Branch("lat_water", &lat_water, "lat_water/D");
  // icetree->Branch("h20_depth", &h20_depth, "h20_depth/D");

  // TTree *groundtree = new TTree("groundtree", "groundtree");
  // groundtree->Branch("elev", &elev, "elev/D");
  // groundtree->Branch("lon_ground", &lon_ground, "lon_ground/D");
  // groundtree->Branch("lat_ground", &lat_ground, "lat_ground/D");

  // //End block added by Stephen

  // TTree *tree11 = new TTree("h11000", "h11000"); // tree11
  // tree11->Branch("loctrig00", &loctrig[0][0], "loctrig0/D");
  // tree11->Branch("loctrig10", &loctrig[1][0], "loctrig0/D");
  // tree11->Branch("loctrig20", &loctrig[2][0], "loctrig0/D");
  // tree11->Branch("loctrig_nadironly0", &loctrig_nadironly[0], "loctrig_nadironly0/D");
  // tree11->Branch("loctrig01", &loctrig[0][1], "loctrig1/D");
  // tree11->Branch("loctrig11", &loctrig[1][1], "loctrig1/D");
  // tree11->Branch("loctrig21", &loctrig[2][1], "loctrig1/D");
  // tree11->Branch("loctrig_nadironly1", &loctrig_nadironly[1], "loctrig0/D");

  // TTree *tree16 = new TTree("h16000", "h16000");
  // tree16->Branch("pnu", &pnu, "pnu/D");
  // tree16->Branch("ptau", &ptau, "ptau/D");
  // tree16->Branch("taulength", &taulength, "taulength/D");
  // tree16->Branch("weight1", &weight1, "weight1/D");
  // tree16->Branch("emfrac", &emfrac, "emfrac/D");
  // tree16->Branch("hadfrac", &hadfrac, "hadfrac/D");
  // tree16->Branch("nuentrancelength", &nuentrancelength, "nuentrancelength/D");

  // int pdgcode;

  // TTree *tree18 = new TTree("h18000", "h18000");
  // tree18->Branch("emfrac",  &emfrac,  "emfrac/D");
  // tree18->Branch("hadfrac", &hadfrac, "hadfrac/D");
  // tree18->Branch("pdgcode", &pdgcode, "pdgcode/I");


  // TH1D *h1mybeta = new TH1D("betaforall", "betaforall(deg)", 180, -15, 15);
  // TH1D *h1mytheta= new TH1D("mytheta", "mytheta(deg)", 180, -90, 90);//90-incidentangle when neutrinos enter the Earth.
  // TH1F *hundogaintoheight_e=new TH1F("undogaintoheight_e", "undogaintoheight_e", 100, 0., 1.);
  // TH1F *hundogaintoheight_h=new TH1F("undogaintoheight_h", "undogaintoheight_h", 100, 0., 1.);
  // TH1F *rec_diff=new TH1F("rec_diff", "rec_diff", 100, -1., 1.);
  // TH1F *recsum_diff=new TH1F("recsum_diff", "recsum_diff", 100, -1., 1.);
  // TH1F *rec_diff0=new TH1F("rec_diff0", "rec_diff0", 100, -1., 1.);
  // TH1F *rec_diff1=new TH1F("rec_diff1", "rec_diff1", 100, -1., 1.);
  // TH1F *rec_diff2=new TH1F("rec_diff2", "rec_diff2", 100, -1., 1.);
  // TH1F *rec_diff3=new TH1F("rec_diff3", "rec_diff3", 100, -1., 1.);

  // TTree *vmmhz_tree = new TTree("vmmhz_tree", "vmmhz_tree"); //To record frequency spread at point where it is first filled
  // vmmhz_tree->Branch("freq_bins", &freq_bins, "freq_bins/I");
  // vmmhz_tree->Branch("vmmhz", &vmmhz, "vmmhz[freq_bins]/D");

  // TTree *tree1 = new TTree("h1000", "h1000"); // tree1 filled for each neutrino
  // tree1->Branch("inu", &inu, "inu/I");
  // tree1->Branch("diffexit", &diffexit, "diffexit/D");
  // tree1->Branch("diffrefr", &diffrefr, "diffrefr/D");
  // tree1->Branch("horizcoord", &horizcoord, "horizcoord/D");
  // tree1->Branch("vertcoord", &vertcoord, "vertcoord/D");
  // tree1->Branch("costhetanu", &costhetanu, "costhetanu/D");
  // tree1->Branch("vmmhz1m_max", &vmmhz1m_max, "vmmhz1m_max/D");
  // tree1->Branch("volume_thishorizon", &volume_thishorizon, "volume_thishorizon/D");
  // tree1->Branch("realtime", &realtime_this, "realtime/D");
  // tree1->Branch("longitude", &longitude_this, "longitude/D");
  // tree1->Branch("latitude", &latitude_this, "latitude/D");
  // tree1->Branch("MAXHORIZON", &bn1->MAXHORIZON, "MAXHORIZON/D");
  // tree1->Branch("igps", &bn1->igps, "igps/I");
  // tree1->Branch("passes_thisevent", &passes_thisevent, "passes_thisevent/I");
  // tree1->Branch("igps", &bn1->igps, "igps/I");
  // tree1->Branch("weight", &weight, "weight/D");
  // tree1->Branch("r_exit2bn", &interaction1->r_exit2bn, "r_exit2bn/D");
  // tree1->Branch("bn1->igps", &bn1->igps, "bn1->igps/I");

  // // set up balloontree

  // TTree *balloontree = new TTree("balloon", "balloon"); //filled for all events
  // balloontree->Branch("heading", &bn1->heading, "heading/D");
  // balloontree->Branch("pitch", &bn1->pitch, "pitch/D");
  // balloontree->Branch("roll", &bn1->roll, "roll/D");
  // balloontree->Branch("realTime_flightdata", &bn1->realTime_flightdata, "realTime_flightdata/I");
  // balloontree->Branch("latitude", &bn1->latitude, "latitude/D");
  // balloontree->Branch("longitude", &bn1->longitude, "longitude/D");
  // balloontree->Branch("altitude", &bn1->altitude, "altitude/D");
  // balloontree->Branch("horizcoord_bn", &bn1->horizcoord_bn, "horizcoord_bn/D");
  // balloontree->Branch("vertcoord_bn", &bn1->vertcoord_bn, "vertcoord_bn/D");
  

  
}
