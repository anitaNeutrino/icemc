#include "RootOutput.h"
#include <iostream>
#include "EventGenerator.h"
#include "Primaries.h"
// #include "balloon.hh"
// #include "anita.hh"
#include "ANITA.h"
#include "Settings.h"
#include "Taumodel.hh"
#include "EnvironmentVariable.h"
#include "Tools.h"
#include "ray.hh"
#include "screen.hh"
#include "position.hh"
#include "icemodel.hh"
#include "IcemcLog.h"

#include "GeneratedNeutrino.h"

#ifdef ANITA_UTIL_EXISTS
#include "UsefulAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "RawAnitaHeader.h"
#include "Adu5Pat.h"
#include "FFTtools.h"

#ifdef ANITA3_EVENTREADER
#include "TruthAnitaEvent.h"
#endif

#endif






icemc::RootOutput::RootOutput(const EventGenerator* uhen, const Settings* settings, const char* outputDir, int run)
  : fOutputDir(outputDir), fRun(run), fIceFinal(NULL),
    realEvPtr(NULL), rawHeaderPtr(NULL), Adu5PatPtr(NULL), truthEvPtr(NULL),
    fHeadFile(NULL), fGpsFile(NULL), fEventFile(NULL), fTruthFile(NULL)
{

  initIceFinal(uhen, settings);
  initRootifiedAnitaDataFiles(uhen, settings);
}



icemc::RootOutput::~RootOutput(){

  // write, close and delete all non-NULL member files.
  const int numFiles = 5;
  TFile* fs[numFiles] = {fIceFinal, fHeadFile, fGpsFile, fEventFile, fTruthFile};
  for(int i=0; i < numFiles;  i++){
    if(fs[i]){
      fs[i]->Write();
      fs[i]->Close();
      delete fs[i];
    }
  }
}



void icemc::RootOutput::initTree(TTree* t, const char* name, const char* title, TFile* f){
  t->SetName(name);
  t->SetTitle(title);
  t->SetDirectory(f);
}


void icemc::RootOutput::initHist(TH2* h, const char* name, const char* title,
				 int nx, double xMin, double xMax,
				 int ny, double yMin, double yMax){
  h->SetNameTitle(name, title);
  h->SetBins(nx, xMin, xMax, ny,  yMin, yMax);
  h->SetDirectory(fIceFinal);
}


void icemc::RootOutput::initHist(TH1* h, const char* name, const char* title,
				 int nx, double xMin, double xMax){
  h->SetNameTitle(name, title);
  h->SetBins(nx, xMin, xMax);
  h->SetDirectory(fIceFinal);
}



void icemc::RootOutput::initIceFinal(const EventGenerator* uhen2, const Settings* settings2){

  if(fIceFinal){
    Log() << icemc::warning << "IceFinal already initialized!"  << std::endl;
    return;
  }

  EventGenerator* uhen = const_cast<EventGenerator*>(uhen2);
  Settings* settings = const_cast<Settings*>(settings2);

  // first the file(s)
  TString fileName = fOutputDir + TString::Format("/icefinal%d.root", fRun);
  fIceFinal = new TFile(fileName, "RECREATE", "ice");

  TNamed* ss = settings2->makeRootSaveableSettings();
  ss->Write();
  delete ss;
  ss = NULL;

  initTree(&allTree, "allTree", "allTree", fIceFinal);
  allTree.Branch("genNu", &uhen->fGenNu);

  initTree(&passTree, "passTree", "passTree", fIceFinal);
  passTree.Branch("passNu", &uhen->fPassNu);

  // histograms
  initHist(&ref_int_coord, "ref_int_coord", "", 600, -3000, 3000, 500, -2500, 2500);
  ref_int_coord.SetMarkerSize(0.7);
  ref_int_coord.SetMarkerColor(kBlue);

  
  initHist(&dir_int_coord, "dir_int_coord", "", 600, -3000, 3000, 500, -2500, 2500);
  dir_int_coord.SetMarkerSize(0.7);
  dir_int_coord.SetMarkerStyle(30);

  
  initHist(&h1mybeta, "betaforall", "betaforall(deg)", 180, -15, 15);
  initHist(&h1mytheta, "mytheta", "mytheta(deg)", 180, -90, 90);
  initHist(&hundogaintoheight_e, "undogaintoheight_e", "undogaintoheight_e", 100, 0., 1.);
  initHist(&hundogaintoheight_h, "undogaintoheight_h", "undogaintoheight_h", 100, 0., 1.);
  initHist(&rec_diff, "rec_diff", "rec_diff", 100, -1., 1.);
  initHist(&recsum_diff, "recsum_diff", "recsum_diff", 100, -1., 1.);
  initHist(&rec_diff0, "rec_diff0", "rec_diff0", 100, -1., 1.);
  initHist(&rec_diff1, "rec_diff1", "rec_diff1", 100, -1., 1.);
  initHist(&rec_diff2, "rec_diff2", "rec_diff2", 100, -1., 1.);
  initHist(&rec_diff3, "rec_diff3", "rec_diff3", 100, -1., 1.);
  initHist(&prob_eachphi_bn, "prob_eachphi_bn", "prob_eachphi_bn", 100, 0., 6.3);
  initHist(&prob_eachilon_bn, "prob_eachilon_bn", "prob_eachilon_bn", 180, 0., 180.);
  initHist(&h6, "theta_vs_hitangle_h", "theta_vs_hitangle_h", 100, -3.14, 3.14, 100, -1.1, 1.1);
  initHist(&h10, "hitangle_e", "hitangle_e", 20, -1.6, 1.6);
  initHist(&hy, "hy", "hy", 100, 0., 1.);
  initHist(&fraction_sec_muons, "fraction_sec_muons", "fraction_sec_muons", 100, 0., 1.);
  initHist(&fraction_sec_taus, "fraction_sec_taus", "fraction_sec_taus", 100, 0., 1.);
  initHist(&n_sec_muons, "n_sec_muons", "n_sec_muons", 100, 0., 10.);
  initHist(&n_sec_taus, "n_sec_taus", "n_sec_taus", 100, 0., 10.);
  initHist(&sampleweights, "sampleweights", "sampleweights", 100, -5., 0.);
  
  
  


  // initTree(&viewangletree, "viewangletree", "viewangletree", fIceFinal); // signal as it is produced at the interaction
  // viewangletree.Branch("dviewangle_deg", &uhen->dviewangle_deg, "dviewangle_deg/D");
  // viewangletree.Branch("emfrac", &uhen->emfrac, "emfrac/D");
  // viewangletree.Branch("hadfrac", &uhen->hadfrac, "hadfrac/D");
  // viewangletree.Branch("deltheta_em_max", &uhen->deltheta_em_max, "deltheta_em_max/D");
  // viewangletree.Branch("deltheta_had_max", &uhen->deltheta_had_max, "deltheta_had_max/D");
  // viewangletree.Branch("theta_threshold_deg", &uhen->theta_threshold_deg, "theta_threshold_deg/D");
  // viewangletree.Branch("dnutries", &uhen->interaction1->dnutries, "dnutries/D");
  // viewangletree.Branch("viewangle", &uhen->viewangle, "viewangle/D");
  // viewangletree.Branch("chord", &uhen->interaction1->chord, "chord/D");

  //Filled just after Getchord,  where we find the neutrino's path through the Earth
  initTree(&nupathtree, "nupathtree", "nupathtree", fIceFinal);
  nupathtree.Branch("total_kgm2", &uhen->total_kgm2, "total_kgm2/D");
  nupathtree.Branch("chord", &uhen->interaction1->chord, "chord/D");
  nupathtree.Branch("crust_entered", &uhen->crust_entered, "crust_entered/I");
  nupathtree.Branch("mantle_entered", &uhen->mantle_entered, "mantle_entered/I");
  nupathtree.Branch("core_entered", &uhen->core_entered, "core_entered/I");
  nupathtree.Branch("mybeta", &uhen->mybeta, "mybeta/D");
  nupathtree.Branch("costheta_nutraject", &uhen->interaction1->costheta_nutraject, "costheta_nutraject/D");

  
  initTree(&finaltree, "passing_events", "passing_events", fIceFinal); // finaltree filled for all events that pass
  
  // finaltree.Branch("inu", &uhen->inu, "inu/I"); @todo TEMPORARILY COMMENT OUT DURING REFACTOR, if abandon refactor then uncomment
  finaltree.Branch("vmmhz_min", &uhen->vmmhz_min, "vmmhz_min/D");
  finaltree.Branch("vmmhz_max", &uhen->vmmhz_max, "vmmhz_max/D");
  finaltree.Branch("thresholdsAnt", &uhen->thresholdsAnt, "thresholdsAnt[48][2][5]/D");
  finaltree.Branch("thresholdsAntPass", &uhen->thresholdsAntPass, "thresholdsAntPass[48][2][5]/D");
  // finaltree.Branch("deadTime", &uhen->anita1->deadTime, "deadTime/D");
  finaltree.Branch("horizcoord", &uhen->horizcoord, "horizcoord/D");
  finaltree.Branch("vertcoord", &uhen->vertcoord, "vertcoord/D");
  // finaltree.Branch("horizcoord_bn", &uhen->bn1->horizcoord_bn, "horizcoord_bn/D");
  // finaltree.Branch("vertcoord_bn", &uhen->bn1->vertcoord_bn, "vertcoord_bn/D");
  finaltree.Branch("r_bn", &uhen->r_bn_array, "r_bn_array[3]/D");
  finaltree.Branch("n_bn", &uhen->n_bn_array, "n_bn_array[3]/D");
  // finaltree.Branch("longitude_bn", &uhen->longitude_this, "longitude_bn/D"); @todo TEMPORARILY COMMENT OUT DURING REFACTOR, if abandon refactor then uncomment
  // finaltree.Branch("heading_bn", &uhen->heading_this, "heading_bn/D"); @todo TEMPORARILY COMMENT OUT DURING REFACTOR, if abandon refactor then uncomment
  finaltree.Branch("gps_offset", &uhen->gps_offset, "gps_offset/D");
  // this one is just weight due to earth absorption
  // finaltree.Branch("weight1", &weight1, "weight1/D");
  // // this is the total weight - the one you want to use!
  // finaltree.Branch("weight", &weight, "weight/D");
  // finaltree.Branch("logweight", &logweight, "logweight/D");
  finaltree.Branch("neutrinoPath", uhen->fNeutrinoPath);
  finaltree.Branch("posnu", &uhen->posnu_array, "posnu_array[3]/D");
  finaltree.Branch("costheta_nutraject", &uhen->costheta_nutraject2, "costheta_nutraject/D");
  finaltree.Branch("chord_kgm2_ice",  &uhen->chord_kgm2_ice2, "chord_kgm2_ice/D");
  finaltree.Branch("phi_nutraject", &uhen->phi_nutraject2, "phi_nutraject/D");
  finaltree.Branch("altitude_int", &uhen->altitude_int2, "altitude_int/D");
  finaltree.Branch("nnu", &uhen->nnu_array, "nnu_array[3]/D");
  finaltree.Branch("n_exit2bn", &uhen->n_exit2bn_array, "n_exit2bn_array[5][3]/D");
  finaltree.Branch("n_exit_phi", &uhen->n_exit_phi, "n_exit_phi/D");
  finaltree.Branch("rfexit", &uhen->rfexit_array, "rfexit_array[5][3]/D");

  finaltree.Branch("pnu", &uhen->pnu, "pnu/D");
  finaltree.Branch("elast_y", &uhen->elast_y, "elast_y/D");
  // finaltree.Branch("emfrac", &uhen->emfrac, "emfrac/D");
  // finaltree.Branch("hadfrac", &uhen->hadfrac, "hadfrac/D");
  // finaltree.Branch("sumfrac", &uhen->sumfrac, "sumfrac/D");
  finaltree.Branch("nuflavor", &uhen->nuflavorint2, "nuflavorint/I");//1=electron,  2=muon,  3=tau
  finaltree.Branch("current", &uhen->currentint2, "currentint/I");//0=charged current,  1=neutral current
  finaltree.Branch("logchord",  &uhen->logchord2,  "logchord/D");
  finaltree.Branch("nuexitice",  &uhen->nuexitice,  "nuexitice/D");
  finaltree.Branch("weight_bestcase",  &uhen->weight_bestcase2,  "weight_bestcase/D");
  finaltree.Branch("chord_kgm2_bestcase",  &uhen->chord_kgm2_bestcase2,  "chord_kgm2_bestcase/D");
  finaltree.Branch("dtryingdirection",  &uhen->dtryingdirection2,  "dtryingdirection/D");
  finaltree.Branch("l3trig", &uhen->l3trig, "l3trig[2]/I");
  finaltree.Branch("l2trig", &uhen->l2trig, "l2trig[2][3]/I");
  finaltree.Branch("l1trig", &uhen->l1trig, "l1trig[2][3]/I");
  // finaltree.Branch("phiTrigMask", &uhen->anita1->phiTrigMask, "phiTrigMask/s");
  // finaltree.Branch("phiTrigMaskH", &uhen->anita1->phiTrigMaskH, "phiTrigMaskH/s");
  // finaltree.Branch("l1TrigMask", &uhen->anita1->l1TrigMask, "l1TrigMask/s");
  // finaltree.Branch("l1TrigMaskH", &uhen->anita1->l1TrigMaskH, "l1TrigMaskH/s");
  finaltree.Branch("max_antenna0", &uhen->max_antenna0, "max_antenna0/I");
  finaltree.Branch("max_antenna1", &uhen->max_antenna1, "max_antenna1/I");
  finaltree.Branch("max_antenna2", &uhen->max_antenna2, "max_antenna2/I");

  finaltree.Branch("viewangle", &uhen->viewangle, "viewangle/D");
  finaltree.Branch("offaxis", &uhen->offaxis, "offaxis/D");
  finaltree.Branch("rx0_signal_eachband", &uhen->rx0_signal_eachband, "rx0_signal_eachband[2][5]/D");
  finaltree.Branch("rx0_threshold_eachband", &uhen->rx0_threshold_eachband, "rx0_threshold_eachband[2][5]/D");
  finaltree.Branch("rx0_noise_eachband", &uhen->rx0_noise_eachband, "rx0_noise_eachband[2][5]/D");
  finaltree.Branch("rx0_passes_eachband", &uhen->rx0_passes_eachband, "rx0_passes_eachband[2][5]/I");
  finaltree.Branch("e_component", &uhen->e_component, "e_component/D");
  finaltree.Branch("h_component", &uhen->h_component, "h_component/D");
  finaltree.Branch("dist_int_bn_2d", &uhen->dist_int_bn_2d, "dist_int_bn_2d/D");
  finaltree.Branch("d1", &uhen->d12, "d1/D");

  finaltree.Branch("cosalpha", &uhen->cosalpha, "cosalpha/D");
  finaltree.Branch("mytheta", &uhen->mytheta, "mytheta/D");
  finaltree.Branch("cosbeta0", &uhen->cosbeta0, "cosbeta0/D");
  finaltree.Branch("mybeta", &uhen->mybeta, "mybeta/D");
  finaltree.Branch("d1", &uhen->d12, "d1/D");
  finaltree.Branch("d2", &uhen->d22, "d2/D");

  //Begin block added by Stephen for verification plots
  finaltree.Branch("fresnel1", &uhen->fresnel1, "fresnel1/D");
  finaltree.Branch("fresnel2", &uhen->fresnel2, "fresnel2/D");
  finaltree.Branch("mag1", &uhen->mag1, "mag1/D");
  finaltree.Branch("mag2", &uhen->mag2, "mag2/D");
  finaltree.Branch("t_coeff_pokey", &uhen->t_coeff_pokey, "t_coeff_pokey/D");
  finaltree.Branch("t_coeff_slappy", &uhen->t_coeff_slappy, "t_coeff_slappy/D");
  // finaltree.Branch("exponent", settings1.EXPONENT); //, "EXPONENT/D");
  finaltree.Branch("exponent", &settings->EXPONENT, "EXPONENT/D");

  finaltree.Branch("hitangle_e_all", &uhen->hitangle_e_all, "hitangle_e_all[48]/D");
  finaltree.Branch("hitangle_h_all", &uhen->hitangle_h_all, "hitangle_h_all[48]/D");

  finaltree.Branch("e_comp_max1", &uhen->e_comp_max1, "e_comp_max1/D");
  finaltree.Branch("h_comp_max1", &uhen->h_comp_max1, "h_comp_max1/D");
  finaltree.Branch("e_comp_max2", &uhen->e_comp_max2, "e_comp_max2/D");
  finaltree.Branch("h_comp_max2", &uhen->h_comp_max2, "h_comp_max2/D");
  finaltree.Branch("e_comp_max3", &uhen->e_comp_max3, "e_comp_max3/D");
  finaltree.Branch("h_comp_max3", &uhen->h_comp_max3, "h_comp_max3/D");
  finaltree.Branch("max_antenna_volts0", &uhen->max_antenna_volts0, "max_antenna_volts0/D");
  finaltree.Branch("max_antenna_volts0_em", &uhen->max_antenna_volts0_em, "max_antenna_volts0_em/D");
  finaltree.Branch("max_antenna_volts1", &uhen->max_antenna_volts1, "max_antenna_volts1/D");
  finaltree.Branch("max_antenna_volts2", &uhen->max_antenna_volts2, "max_antenna_volts2/D");
  finaltree.Branch("triggers", &uhen->nchannels_perrx_triggered, "nchannels_perrx_triggered[48]/I");
  finaltree.Branch("nchannels_triggered", &uhen->nchannels_triggered, "nchannels_triggered/I");
  finaltree.Branch("voltsRX", &uhen->voltsRX);
  
  // finaltree.Branch("volts_rx_max", &volts_rx_max, "volts_rx_max/D");
  // finaltree.Branch("volts_rx_ave", &volts_rx_ave, "volts_rx_ave/D");
  // finaltree.Branch("volts_rx_sum", &volts_rx_sum, "volts_rx_sum/D");

  // finaltree.Branch("volts_rx_max_highband", &volts_rx_max_highband, "volts_rx_max_highband/D");
  // finaltree.Branch("volts_rx_max_lowband", &volts_rx_max_lowband, "volts_rx_max_lowband/D");

  finaltree.Branch("theta_pol_measured", &uhen->theta_pol_measured, "theta_pol_measured/D");  
  finaltree.Branch("theta_rf_atbn", &uhen->theta_rf_atbn, "theta_rf_atbn/D");
  finaltree.Branch("theta_rf_atbn_measured", &uhen->theta_rf_atbn_measured, "theta_rf_atbn_measured/D");
  finaltree.Branch("voltage", &uhen->voltagearray, "voltagearray[48]/D");
  // finaltree.Branch("nlayers", settings1.NLAYERS, "NLAYERS/I");

  finaltree.Branch("vmmhz1m_max", &uhen->vmmhz1m_max, "vmmhz1m_max/D");
  finaltree.Branch("vmmhz_lowfreq", &uhen->vmmhz_lowfreq, "vmmhz_lowfreq/D");

  finaltree.Branch("deltheta_em_max", &uhen->deltheta_em_max, "deltheta_em_max/D");
  finaltree.Branch("deltheta_had_max", &uhen->deltheta_had_max, "deltheta_had_max/D");
  finaltree.Branch("r_enterice", &uhen->r_enterice_array, "r_enterice_array[3]/D");
  finaltree.Branch("n_exit2bn_db", &uhen->n_exit2bn_db_array, "n_exit2bn_db_array[5][3]/D");

  finaltree.Branch("rfexit_db", &uhen->rfexit_db_array, "rfexit_db_array[5][3]/D");
  finaltree.Branch("r_in", &uhen->r_in_array, "r_in_array[3]/D");
  finaltree.Branch("nsurf_rfexit", &uhen->nsurf_rfexit_array, "nsurf_rfexit_array[3]/D");
  finaltree.Branch("nsurf_rfexit_db", &uhen->nsurf_rfexit_db_array, "nsurf_rfexit_db_array[3]/D");
  finaltree.Branch("r_fromballoon", &uhen->r_fromballoon2, "r_fromballoon/D");
  finaltree.Branch("r_fromballoon_db", &uhen->interaction1->r_fromballoon_db, "r_fromballoon_db/D");

  finaltree.Branch("nuexitlength", &uhen->nuexitlength, "nuexitlength/D");
  finaltree.Branch("nuentrancelength", &uhen->nuentrancelength, "nuentrancelength/D");
  finaltree.Branch("taulength", &uhen->taulength, "taulength/D");
  finaltree.Branch("icethickness", &uhen->icethickness, "icethickness/D");
  finaltree.Branch("nrf_iceside", &uhen->nrf_iceside_array, "nrf_iceside_array[5][3]/D");
  finaltree.Branch("nrf_iceside_db", &uhen->nrf_iceside_db_array, "nrf_iceside_db_array[5][3]/D");
  finaltree.Branch("ant_normal0", &uhen->ant_max_normal0_array, "ant_max_normal0_array[3]/D");
  finaltree.Branch("ant_normal1", &uhen->ant_max_normal1_array, "ant_max_normal1_array[3]/D");
  finaltree.Branch("ant_normal2", &uhen->ant_max_normal2_array, "ant_max_normal2_array[3]/D");
  // finaltree.Branch("vmmhz1m_visible", &uhen->vmmhz1m_visible, "vmmhz1m_visible/D");
  finaltree.Branch("freq_bins", &uhen->freq_bins, "freq_bins/I");
  // finaltree.Branch("vmmhz", &uhen->vmmhz, "vmmhz[freq_bins]/D");@todo TEMPORARILY COMMENT OUT DURING REFACTOR, if abandon refactor then uncomment

  finaltree.Branch("dist_int_bn_2d_chord", &uhen->dist_int_bn_2d_chord, "dist_int_bn_2d_chord/D");

  finaltree.Branch("dviewangle_deg", &uhen->dviewangle_deg, "dviewangle_deg/D");
  finaltree.Branch("theta_threshold_deg", &uhen->theta_threshold_deg, "theta_threshold_deg/D");
  finaltree.Branch("total_kgm2", &uhen->total_kgm2, "total_kgm2/D");
  finaltree.Branch("chord", &uhen->interaction1->chord, "chord/D");
  finaltree.Branch("crust_entered", &uhen->crust_entered, "crust_entered/I");
  finaltree.Branch("mantle_entered", &uhen->mantle_entered, "mantle_entered/I");
  finaltree.Branch("core_entered", &uhen->core_entered, "core_entered/I");
  finaltree.Branch("n_pol", &uhen->n_pol_array, "n_pol_array[3]/D");
  finaltree.Branch("vmmhz_min_thatpasses", &uhen->vmmhz_min_thatpasses, "vmmhz_min_thatpasses/D");

  // finaltree.Branch("pieceofkm2sr", &pieceofkm2sr, "pieceofkm2sr/D");
  //finaltree.Branch("volts_original", &volts_original, "volts_original[10][20][2]/D");
  finaltree.Branch("r_exit2bn", &uhen->r_exit2bn2, "r_exit2bn/D");
  finaltree.Branch("r_exit2bn_measured", &uhen->r_exit2bn_measured2, "r_exit2bn_measured/D");
  // finaltree.Branch("scalefactor_attenuation", &uhen->scalefactor_attenuation, "scalefactor_attenuation/D"); @todo TEMPORARILY COMMENT OUT DURING REFACTOR, if abandon refactor then uncomment
  // finaltree.Branch("anita1->PHI_OFFSET", &uhen->anita1->PHI_OFFSET, "anita1->PHI_OFFSET/D");
  // finaltree.Branch("igps", &uhen->bn1->igps, "igyps/I");
  // finaltree.Branch("volts_rx_rfcm_lab_e_all", &uhen->volts_rx_rfcm_lab_e_all, "volts_rx_rfcm_lab_e_all[48][512]/D");
  // finaltree.Branch("volts_rx_rfcm_lab_h_all", &volts_rx_rfcm_lab_h_all, "volts_rx_rfcm_lab_h_all[48][512]/D");
  finaltree.Branch("ptaui", &uhen->ptaui, "ptaui/D");
  finaltree.Branch("ptauf", &uhen->ptauf, "ptauf/D");
  finaltree.Branch("sourceLon", &uhen->sourceLon, "sourceLon/D");
  finaltree.Branch("sourceLat", &uhen->sourceLat, "sourceLat/D");
  finaltree.Branch("sourceAlt", &uhen->sourceAlt, "sourceAlt/D");
  finaltree.Branch("sourceMag", &uhen->sourceMag, "sourceMag/D");

  initTree(&mytaus_tree, "mytaus", "mytaus", fIceFinal);
  mytaus_tree.Branch("taus", &uhen->fTauPtr); // this is how you do it!

  initTree(&summarytree, "summarytree", "summarytree", fIceFinal); // finaltree filled for all events that pass
  // summarytree.Branch("NNU", &settings1.NNU, "NNU/I");
  summarytree.Branch("NNU", &settings->NNU, "NNU/I");  
  // summarytree.Branch("EXPONENT", settings1.EXPONENT); //, "EXPONENT/D");
  summarytree.Branch("eventsfound_beforetrigger", &uhen->eventsfound_beforetrigger, "eventsfound_beforetrigger/D");
  summarytree.Branch("rms_rfcm_e", &uhen->rms_rfcm_e, "rms_rfcm_e/D");
  summarytree.Branch("rms_rfcm_h", &uhen->rms_rfcm_h, "rms_rfcm_h/D");
  summarytree.Branch("rms_lab_e", &uhen->rms_lab_e, "rms_lab_e/D");
  summarytree.Branch("rms_lab_h", &uhen->rms_lab_h, "rms_lab_h/D");
  summarytree.Branch("avgfreq_rfcm", &uhen->avgfreq_rfcm, "avgfreq_rfcm[128]/D");
  summarytree.Branch("avgfreq_rfcm_lab", &uhen->avgfreq_rfcm_lab, "avgfreq_rfcm_lab[128]/D");
  summarytree.Branch("freq", &uhen->freq, "freq[128]/D");

  initTree(&ytree, "ytree", "ytree", fIceFinal); //To record y distributions
  ytree.Branch("elast_y", &uhen->elast_y, "elast_y/D");

}






void icemc::RootOutput::initRootifiedAnitaDataFiles(const EventGenerator* uhen2, const Settings* settings1){

  EventGenerator* uhen = const_cast<EventGenerator*>(uhen2);
  
#ifdef ANITA_UTIL_EXISTS

  TString eventFileName = fOutputDir + TString::Format("/SimulatedAnitaEventFile%d.root", fRun);
  fEventFile = new TFile(eventFileName, "RECREATE");

  initTree(&eventTree, "eventTree", "eventTree", fEventFile);
  eventTree.Branch("event",             &realEvPtr           );
  // eventTree->Branch("run",               clOpts.run_no,   "run/I"   );
  eventTree.Branch("weight",            &uhen->fNeutrinoPath->weight,   "weight/D");

  TString headFileName = fOutputDir + TString::Format("/SimulatedAnitaHeadFile%d.root", fRun);
  fHeadFile = new TFile(headFileName, "RECREATE");
  
  initTree(&headTree, "headTree", "headTree", fHeadFile);
  headTree.Branch("header",  &rawHeaderPtr           );
  headTree.Branch("weight",  &uhen->fNeutrinoPath->weight,      "weight/D");

  TString gpsFileName = fOutputDir + TString::Format("/SimulatedAnitaGpsFile%d.root", fRun);
  fGpsFile = new TFile(gpsFileName, "RECREATE");


  initTree(&adu5PatTree, "adu5PatTree", "adu5PatTree", fGpsFile);
  adu5PatTree.Branch("pat",          &Adu5PatPtr                   );
  adu5PatTree.Branch("eventNumber",  &uhen->eventNumber,  "eventNumber/I");
  adu5PatTree.Branch("weight",       &uhen->fNeutrinoPath->weight,       "weight/D"     );

#ifdef ANITA3_EVENTREADER

  // Set AnitaVersion so that the right payload geometry is used
  AnitaVersion::set(settings1->ANITAVERSION);
  
  TString truthFileName = fOutputDir + TString::Format("/SimulatedAnitaTruthFile%d.root",  fRun);
  fTruthFile = new TFile(truthFileName, "RECREATE");
  TNamed* ss = settings1->makeRootSaveableSettings();
  ss->Write();
  delete ss;
  ss = NULL;

  TNamed gitversion("GitVersion", EnvironmentVariable::ICEMC_VERSION(fOutputDir).c_str());
  gitversion.Write();

  TNamed nnu("NNU", TString::Format("%d", settings1->NNU).Data());
  nnu.Write();

  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  TNamed starttime("StartTime", asctime(timeinfo));
  starttime.Write();

  Log() << icemc::info << "ICEMC GIT Repository Version: " <<  gitversion.GetTitle() << std::endl;

  initTree(&triggerSettingsTree, "triggerSettingsTree", "Trigger settings", fTruthFile);
  // Anita* anita2 = const_cast<Anita*>(uhen->anita1);
  ANITA* anita2 = const_cast<ANITA*>(uhen->fDetector);
  
  triggerSettingsTree.Branch("dioderms",  anita2->bwslice_dioderms_fullband_allchan,  "dioderms[2][48][6]/D" );
  triggerSettingsTree.Branch("diodemean", anita2->bwslice_diodemean_fullband_allchan, "diodemean[2][48][6]/D");
  triggerSettingsTree.Fill();
  
  initTree(&truthTree, "truthAnitaTree", "Truth Anita Tree", fTruthFile);
  truthTree.Branch("truth", &truthEvPtr);

#else
  Log() << icemc::warning << "Need ANITA EventReader version at least 3 to produce Truth output." << std::endl;    
#endif

#else
  Log() << icemc::warning << "Can't generate ROOTified output without satisfying ANITA_UTIL_EXISTS at compile time" << std::endl;
#endif  

}


int icemc::RootOutput::getIceMCAntfromUsefulEventAnt(const Settings *settings1,  int UsefulEventAnt){

#ifdef ANITA_UTIL_EXISTS  
  int IceMCAnt = UsefulEventAnt;
  if ((settings1->WHICH==9 || settings1->WHICH==10) && UsefulEventAnt<16) {
    IceMCAnt = (UsefulEventAnt%2==0)*UsefulEventAnt/2 + (UsefulEventAnt%2==1)*(UsefulEventAnt/2+8);
  }
  return IceMCAnt;
#else  
  return -1;
#endif
  
}


void icemc::RootOutput::fillRootifiedAnitaDataTrees(const EventGenerator* uhen, const Settings& settings1, const Ray* ray1, const Screen* panel1){  
  
#ifdef ANITA_UTIL_EXISTS
  AnitaGeomTool* geom = AnitaGeomTool::Instance();

  // const Balloon* bn1 = uhen->bn1;
  // const Anita* anita1 = uhen->anita1;
  ///@todo Make this into something proper!!!!
  const ANITA* bn1 = uhen->fDetector;
  const ANITA* anita1 = uhen->fDetector;
  
  realEvPtr     = new UsefulAnitaEvent();
  rawHeaderPtr  = new RawAnitaHeader();
  Adu5PatPtr    = new Adu5Pat();	    

  Adu5PatPtr->latitude= bn1->latitude;
  Adu5PatPtr->longitude=bn1->longitude;
  Adu5PatPtr->altitude=bn1->altitude;
  Adu5PatPtr->realTime=bn1->realTime_flightdata;
  Adu5PatPtr->heading = bn1->heading;
  Adu5PatPtr->pitch = bn1->pitch;
  Adu5PatPtr->roll = bn1->roll;
  Adu5PatPtr->run = fRun;//clOpts.run_no;

  memset(realEvPtr->fNumPoints, 0, sizeof(realEvPtr->fNumPoints) );
  memset(realEvPtr->fVolts,     0, sizeof(realEvPtr->fVolts)     );
  memset(realEvPtr->fTimes,     0, sizeof(realEvPtr->fTimes)     );

  int fNumPoints = 260;
  for (int ichan=0; ichan<108; ichan++){
    realEvPtr->fNumPoints[ichan] = fNumPoints;

    for (int j = 0; j < fNumPoints; j++) {
      // convert seconds to nanoseconds
      realEvPtr->fTimes[ichan][j] = j * anita1->TIMESTEP * 1.0E9;
    }
  }
  realEvPtr->fRFSpike = 0;// when we get as far as simulating this, we're doing well

  for (int iant = 0; iant < settings1.NANTENNAS; iant++){
    //int IceMCAnt = getIceMCAntfromUsefulEventAnt(anita1,  AnitaGeom1,  iant);
    int IceMCAnt = getIceMCAntfromUsefulEventAnt(&settings1,  iant);
    int UsefulChanIndexH = geom->getChanIndexFromAntPol(iant,  AnitaPol::kHorizontal);
    int UsefulChanIndexV = geom->getChanIndexFromAntPol(iant,  AnitaPol::kVertical);
    //              realEvPtr->fNumPoints[UsefulChanIndexV] = fNumPoints;
    //              realEvPtr->fNumPoints[UsefulChanIndexH] = fNumPoints;
    realEvPtr->chanId[UsefulChanIndexV] = UsefulChanIndexV;
    realEvPtr->chanId[UsefulChanIndexH] = UsefulChanIndexH;

    for (int j = 0; j < fNumPoints; j++) {
      // convert seconds to nanoseconds
      //                realEvPtr->fTimes[UsefulChanIndexV][j] = j * anita1->TIMESTEP * 1.0E9;
      //                realEvPtr->fTimes[UsefulChanIndexH][j] = j * anita1->TIMESTEP * 1.0E9;
      // convert volts to millivolts
      const double voltsToMilliVolts = 1000;
      realEvPtr->fVolts[UsefulChanIndexH][j] =  uhen->voltsRX.rfcm_lab_h_all[IceMCAnt][j+128]*voltsToMilliVolts;
      realEvPtr->fCapacitorNum[UsefulChanIndexH][j] = 0;
      realEvPtr->fVolts[UsefulChanIndexV][j] =  uhen->voltsRX.rfcm_lab_e_all[IceMCAnt][j+128]*voltsToMilliVolts;
      realEvPtr->fCapacitorNum[UsefulChanIndexV][j] = 0;
    }//end int j
  }// end int iant

  realEvPtr->eventNumber = uhen->eventNumber;

  rawHeaderPtr->eventNumber = uhen->eventNumber;
  rawHeaderPtr->surfSlipFlag = 0;
  rawHeaderPtr->errorFlag = 0;

  if (settings1.MINBIAS==1){
    rawHeaderPtr->trigType = 8; // soft-trigger
  }
  else{
    rawHeaderPtr->trigType = 1; // RF trigger
  }
	    
  rawHeaderPtr->run = fRun; //clOpts.run_no;
  // put the vpol only as a placeholder - these are only used in Anita-2 anyway
  rawHeaderPtr->upperL1TrigPattern = uhen->l1trig[0][0];
  rawHeaderPtr->lowerL1TrigPattern = uhen->l1trig[0][1];
  rawHeaderPtr->nadirL1TrigPattern = uhen->l1trig[0][2];

  rawHeaderPtr->upperL2TrigPattern = uhen->l2trig[0][0];
  rawHeaderPtr->lowerL2TrigPattern = uhen->l2trig[0][1];
  rawHeaderPtr->nadirL2TrigPattern = uhen->l2trig[0][2];

  if (settings1.WHICH<9){
    rawHeaderPtr->phiTrigMask  = (short) anita1->phiTrigMask;
    rawHeaderPtr->l3TrigPattern = (short) uhen->l3trig[0];
  }

  rawHeaderPtr->calibStatus = 31;
  rawHeaderPtr->realTime = bn1->realTime_flightdata;
  rawHeaderPtr->triggerTime = bn1->realTime_flightdata;
  Adu5PatPtr->latitude= bn1->latitude;
  Adu5PatPtr->longitude=bn1->longitude;
  Adu5PatPtr->altitude=bn1->altitude;
  Adu5PatPtr->realTime=bn1->realTime_flightdata;
  Adu5PatPtr->heading = bn1->heading;
  Adu5PatPtr->pitch = bn1->pitch;
  Adu5PatPtr->roll = bn1->roll;
  Adu5PatPtr->run = fRun; //clOpts.run_no;

#ifdef ANITA3_EVENTREADER
  if (settings1.WHICH==9 || settings1.WHICH==10) {
    rawHeaderPtr->setTrigPattern((short) uhen->l3trig[0], AnitaPol::kVertical);
    rawHeaderPtr->setTrigPattern((short) uhen->l3trig[1], AnitaPol::kHorizontal);
    rawHeaderPtr->setMask( (short) anita1->l1TrigMask,  (short) anita1->phiTrigMask,  AnitaPol::kVertical);
    rawHeaderPtr->setMask( (short) anita1->l1TrigMaskH, (short) anita1->phiTrigMaskH, AnitaPol::kHorizontal);
  }

  truthEvPtr                   = new TruthAnitaEvent();
  truthEvPtr->eventNumber      = uhen->eventNumber;
  truthEvPtr->realTime         = bn1->realTime_flightdata;
  truthEvPtr->run              = fRun; //clOpts.run_no;
  truthEvPtr->nuMom            = uhen->pnu;
  truthEvPtr->nu_pdg           = uhen->pdgcode;
  truthEvPtr->e_component      = uhen->e_component;
  truthEvPtr->h_component      = uhen->h_component;
  truthEvPtr->n_component      = uhen->n_component;
  truthEvPtr->e_component_k    = uhen->e_component_kvector;
  truthEvPtr->h_component_k    = uhen->h_component_kvector;
  truthEvPtr->n_component_k    = uhen->n_component_kvector;
  truthEvPtr->sourceLon        = uhen->sourceLon;
  truthEvPtr->sourceLat        = uhen->sourceLat;
  truthEvPtr->sourceAlt        = uhen->sourceAlt;
  truthEvPtr->weight           = uhen->fNeutrinoPath->weight;
  for (int i=0;i<3;i++){
    truthEvPtr->balloonPos[i]  = bn1->r_bn[i];
    truthEvPtr->balloonDir[i]  = bn1->n_bn[i];
    truthEvPtr->nuPos[i]       = uhen->interaction1->posnu[i];
    truthEvPtr->nuDir[i]       = uhen->interaction1->nnu[i];
  }
  for (int i=0;i<5;i++){
    for (int j=0;j<3;j++){
      truthEvPtr->rfExitNor[i][j] = ray1->n_exit2bn[i][j];
      truthEvPtr->rfExitPos[i][j] = ray1->rfexit[i][j];
    }
  }
  for (int i=0;i<48;i++){
    truthEvPtr->hitangle_e[i]  = uhen->hitangle_e_all[i];
    truthEvPtr->hitangle_h[i]  = uhen->hitangle_h_all[i];
  }
  if(!settings1.ROUGHNESS){
    for (int i=0;i<Anita::NFREQ;i++)
      truthEvPtr->vmmhz[i]       = panel1->GetVmmhz_freq(i);
  }

	    
  memset(truthEvPtr->SNRAtTrigger,       0, sizeof(truthEvPtr->SNRAtTrigger)       );
  memset(truthEvPtr->fSignalAtTrigger,   0, sizeof(truthEvPtr->fSignalAtTrigger)   );
  memset(truthEvPtr->fNoiseAtTrigger,    0, sizeof(truthEvPtr->fNoiseAtTrigger)    );
  memset(truthEvPtr->SNRAtDigitizer,     0, sizeof(truthEvPtr->SNRAtDigitizer)     );
  memset(truthEvPtr->thresholds,         0, sizeof(truthEvPtr->thresholds)         );
  memset(truthEvPtr->fDiodeOutput,       0, sizeof(truthEvPtr->fDiodeOutput)       );
	    
  truthEvPtr->maxSNRAtTriggerV=0;
  truthEvPtr->maxSNRAtTriggerH=0;
  truthEvPtr->maxSNRAtDigitizerV=0;
  truthEvPtr->maxSNRAtDigitizerH=0;

  for (int iant = 0; iant < settings1.NANTENNAS; iant++){
    int UsefulChanIndexH = geom->getChanIndexFromAntPol(iant,  AnitaPol::kHorizontal);
    int UsefulChanIndexV = geom->getChanIndexFromAntPol(iant,  AnitaPol::kVertical);

    truthEvPtr->SNRAtTrigger[UsefulChanIndexV] = Tools::calculateSNR(uhen->justSignal_trig[0][iant], uhen->justNoise_trig[0][iant]);
    truthEvPtr->SNRAtTrigger[UsefulChanIndexH] = Tools::calculateSNR(uhen->justSignal_trig[1][iant], uhen->justNoise_trig[1][iant]);
	      
    if (truthEvPtr->SNRAtTrigger[UsefulChanIndexV]>truthEvPtr->maxSNRAtTriggerV) truthEvPtr->maxSNRAtTriggerV=truthEvPtr->SNRAtTrigger[UsefulChanIndexV];
    if (truthEvPtr->SNRAtTrigger[UsefulChanIndexH]>truthEvPtr->maxSNRAtTriggerH) truthEvPtr->maxSNRAtTriggerH=truthEvPtr->SNRAtTrigger[UsefulChanIndexH];

    truthEvPtr->SNRAtDigitizer[UsefulChanIndexV] = Tools::calculateSNR(uhen->justSignal_dig[0][iant], uhen->justNoise_dig[0][iant]);
    truthEvPtr->SNRAtDigitizer[UsefulChanIndexH] = Tools::calculateSNR(uhen->justSignal_dig[1][iant], uhen->justNoise_dig[1][iant]);
	      
    if (truthEvPtr->SNRAtDigitizer[UsefulChanIndexV]>truthEvPtr->maxSNRAtDigitizerV) truthEvPtr->maxSNRAtDigitizerV=truthEvPtr->SNRAtDigitizer[UsefulChanIndexV];
    if (truthEvPtr->SNRAtDigitizer[UsefulChanIndexH]>truthEvPtr->maxSNRAtDigitizerH) truthEvPtr->maxSNRAtDigitizerH=truthEvPtr->SNRAtDigitizer[UsefulChanIndexH];

	      
    truthEvPtr->thresholds[UsefulChanIndexV] = uhen->thresholdsAnt[iant][0][4];
    truthEvPtr->thresholds[UsefulChanIndexH] = uhen->thresholdsAnt[iant][1][4];
    int irx = iant;
    if (iant<16){
      if (iant%2) irx = iant/2;
      else        irx = iant/2 + 1;
    }
	      
    for (int j = 0; j < fNumPoints; j++) {
      truthEvPtr->fTimes[UsefulChanIndexV][j]             = j * anita1->TIMESTEP * 1.0E9;
      truthEvPtr->fTimes[UsefulChanIndexH][j]             = j * anita1->TIMESTEP * 1.0E9;
		
      truthEvPtr->fSignalAtTrigger[UsefulChanIndexV][j]   = uhen->justSignal_trig[0][iant][j+128]*1000;
      truthEvPtr->fSignalAtTrigger[UsefulChanIndexH][j]   = uhen->justSignal_trig[1][iant][j+128]*1000;
      truthEvPtr->fNoiseAtTrigger[UsefulChanIndexV][j]    = uhen->justNoise_trig[0][iant][j+128]*1000;
      truthEvPtr->fNoiseAtTrigger[UsefulChanIndexH][j]    = uhen->justNoise_trig[1][iant][j+128]*1000;
      truthEvPtr->fSignalAtDigitizer[UsefulChanIndexV][j] = uhen->justSignal_dig[0][iant][j+128]*1000;
      truthEvPtr->fSignalAtDigitizer[UsefulChanIndexH][j] = uhen->justSignal_dig[1][iant][j+128]*1000;
      truthEvPtr->fNoiseAtDigitizer[UsefulChanIndexV][j]  = uhen->justNoise_dig[0][iant][j+128]*1000;
      truthEvPtr->fNoiseAtDigitizer[UsefulChanIndexH][j]  = uhen->justNoise_dig[1][iant][j+128]*1000;
		
      truthEvPtr->fDiodeOutput[UsefulChanIndexV][j]       = anita1->timedomain_output_allantennas[0][irx][j];
      truthEvPtr->fDiodeOutput[UsefulChanIndexH][j]       = anita1->timedomain_output_allantennas[1][irx][j];
    }//end int j
	      
  }// end int iant

  truthTree.Fill();
  delete truthEvPtr;
#endif

  headTree.Fill();
  eventTree.Fill();
  adu5PatTree.Fill();

  delete realEvPtr;
  delete rawHeaderPtr;
  delete Adu5PatPtr;
#endif


}
