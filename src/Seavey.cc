#include "Seavey.h"
#include "EnvironmentVariable.h"
#include "IcemcLog.h"
#include "balloon.hh"
#include "Constants.h"
#include "Settings.h"

#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TFile.h" // @todo for debugging

#include <algorithm>

#if defined(ANITA_UTIL_EXISTS) and defined(VECTORIZE)
#include "vectormath_trig.h"
#endif








 // number of points within bandwidth that gain is measured.
constexpr int numGainPoints = 131;
std::array<double, numGainPoints> freqHz_measured; ///< Frequency values for arrays with numGainPoints, goes from 200 MHz to 1.5 GHz with step size of 10 MHz (units in Hz)

/**
 * Ped's famous gains
 */
std::array<double, numGainPoints> gain_dBi_hh_measured; ///< gain for HPol component of signal to HPol feed (co-pol)
std::array<double, numGainPoints> gain_dBi_vv_measured; ///< gain for VPol component of signal to VPol feed (co-pol)
std::array<double, numGainPoints> gain_dBi_hv_measured; ///< gain for HPol component of signal to VPol feed (cross-pol)
std::array<double, numGainPoints> gain_dBi_vh_measured; ///< gain for VPol component of signal to HPol feed (cross-pol)


/**
 * The antenna heights (m), calculated from Ped's gains, converts an E-field (V/m) to a voltage (V)
 */
std::array<double, numGainPoints> heightHH_m; ///< converted from gain_dBi_hh_measured
std::array<double, numGainPoints> heightVV_m; ///< converted from gain_dBi_vv_measured
std::array<double, numGainPoints> heightHV_m; ///< converted from gain_dBi_hv_measured
std::array<double, numGainPoints> heightVH_m; ///< converted from gain_dBi_vh_measured

// if you come up with some geometry where some are in ice and others not or something
// refractiveIndexUsedForHeightArrays will need to change to be a per-antenna member variable
// maybe it should be anyway because this is a wee bit ugly.
double refractiveIndexUsedForHeightArrays = 0; ///< Required to find height from gains


constexpr int numAnglePoints = 7;
const std::array<double, numAnglePoints> referenceAnglesDeg {0, 5, 10, 20, 30, 45, 90}; // the off axis measurements are have these step sizes (Degrees)
const std::array<double, numAnglePoints> referenceAnglesRad {referenceAnglesDeg[0]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[1]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[2]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[3]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[4]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[5]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[6]*icemc::constants::RADDEG};
// const std::array<double,  numAnglePoints-1> invAngleBinSize {referenceAnglesDeg[1] - referenceAnglesDeg[0],
//                                                              referenceAnglesDeg[2] - referenceAnglesDeg[1],
//                                                              referenceAnglesDeg[3] - referenceAnglesDeg[2],
//                                                              referenceAnglesDeg[4] - referenceAnglesDeg[3],
//                                                              referenceAnglesDeg[5] - referenceAnglesDeg[4],
//                                                              referenceAnglesDeg[6] - referenceAnglesDeg[5]};

std::array<std::array<double, numGainPoints>, numAnglePoints> gain_v_angle_az; //[numAnglePoints][numGainPoints]
std::array<std::array<double, numGainPoints>, numAnglePoints> gain_h_angle_az; //[numAnglePoints][numGainPoints]
std::array<std::array<double, numGainPoints>, numAnglePoints> gain_v_angle_el; //[numAnglePoints][numGainPoints]
std::array<std::array<double, numGainPoints>, numAnglePoints> gain_h_angle_el; //[numAnglePoints][numGainPoints]

bool doneLoadGains = false; ///< So we only read the Seavey gains and calculate derived quantities once.







/** 
 * Utility function for loadGains
 * 
 * @param fileName the file to open
 * @param failHard exit on failure if true
 * 
 * @return the (hopefully) opened ifstream
 */
std::ifstream openCarefully(const char* fileName, bool failHard = true){
  std::ifstream dataFile;
  dataFile.open(fileName);
  if(dataFile.fail()){
    if(failHard){
      icemcLog() << icemc::error << "Quiting because I can't open " << fileName << "\n";
      exit(1);
    }
    else{
      icemcLog() << icemc::warning << "I can't open "  << fileName << "\n";
    }
  }
  return dataFile;
}





/** 
 * Read the gains into the arrays.
 * @todo This is NOT THREAD SAFE, some locking functions would be required to make it so.
 */
void loadGains(){
  if(!doneLoadGains){

    const std::string ICEMC_SRC_DIR=icemc::EnvironmentVariable::ICEMC_SRC_DIR();
    const std::string ICEMC_DATA_DIR=ICEMC_SRC_DIR+"/data/";

    std::vector<std::string> boresightFileNames {"hh_0", "vv_0",  "hv_0", "vh_0"};
    bool firstFile = true; // fill freqHz_measured first time, check against it on all others
    for(const auto& fileName : boresightFileNames){
      std::string fullPath = ICEMC_DATA_DIR + fileName;
      std::ifstream boresight_file = openCarefully(fullPath.c_str());

      auto& boresightGain = (fileName == "vv_0" ? gain_dBi_hh_measured :
			     fileName == "hh_0" ? gain_dBi_vv_measured :
			     fileName == "hv_0" ? gain_dBi_hv_measured :
			                          gain_dBi_vh_measured);
      
      for(int i = 0; i < numGainPoints; ++i) {
	double f;
	boresight_file >> f >> boresightGain.at(i);
	if(firstFile){
	  freqHz_measured.at(i) = f;
	}
	else if (f != freqHz_measured.at(i)){
	  icemcLog() << icemc::warning << "frequency = " << f << ", freqHz_measured[i] = " << freqHz_measured.at(i) << "\n";
	}
      }
      firstFile = false;
    }

    std::vector<std::string> angleFileNames {"vv_az", "hh_az", "vv_el", "hh_el"};

    for(const auto& fileName : angleFileNames){

      auto& gainVsAngle = (fileName == "vv_az" ? gain_v_angle_az :
			   fileName == "hh_az" ? gain_h_angle_az :
			   fileName == "vv_el" ? gain_v_angle_el :
			                         gain_h_angle_el );

      std::string filePath = ICEMC_DATA_DIR + fileName;
      std::ifstream angle_file = openCarefully(filePath.c_str());

      for(auto& g0 : gainVsAngle.at(0)){
	g0 = 1; // 0th bin is on boresight, so no off-axis effects... so set to 1.
      }

      for(int j = 1; j < numAnglePoints; j++){
	int k = 0;
	double f = 0;
	for(auto& g : gainVsAngle.at(j)){	  
	  angle_file >> f >> g;
	  if(f != freqHz_measured.at(k)){
	    icemcLog() << icemc::warning << "Check off-axis frequencies for " << filePath << std::endl;
	  }
	  k++;
	}
      }
      
      angle_file.close();

    }

    doneLoadGains = true;
  }
}


/** 
 * Convert a gain in dBi to antenna height, requires knowing the refractive index of the receiving material
 * 
 * @param gain_dB 
 * @param freq 
 * @param nmedium_receiver 
 * 
 * @return 
 */
inline double getAntennaHeightFromGain_dB(double gain_dB, double freq, double nmedium_receiver) {
  // from gain=4*pi*A_eff/lambda^2
  // and h_eff=2*sqrt(A_eff*Z_rx/Z_air)
  // gain_dB is in dB
  return 2*sqrt(gain_dB/4/icemc::constants::PI*icemc::constants::CLIGHT*icemc::constants::CLIGHT/(freq*freq)*icemc::constants::Zr/(icemc::constants::Z0*nmedium_receiver));      
}




/** 
 * We need the refractive index of the antenna medium to calculate the heights.
 * 
 * @param refractiveIndex is the refractive index of the medium
 * @todo This also isn't thread safe
 * 
 */
void fillHeightArrays(double refractiveIndex){
  
  loadGains();

  
  // then we recalculate everything
  if(refractiveIndex != refractiveIndexUsedForHeightArrays){
    
    for(const std::string& polString : {"vv", "hh", "hv",  "vh"}){

      auto& boresightGain = (polString == "vv" ? gain_dBi_vv_measured :
    			     polString == "hh" ? gain_dBi_hh_measured :
    			     polString == "hv" ? gain_dBi_hv_measured :
      			                         gain_dBi_vh_measured);

      auto& antennaHeight = (polString == "hh" ? heightHH_m :
    			     polString == "vv" ? heightVV_m :
    			     polString == "hv" ? heightHV_m :
        			                 heightVH_m);
      int i=0;
      for(auto g : boresightGain){
	antennaHeight.at(i) = getAntennaHeightFromGain_dB(g, freqHz_measured.at(i), refractiveIndex);
	i++;
      }
    }
  }
}













TCanvas* icemc::Seavey::plotInterpolatedGains(double freq, const int nBins){

  auto c = new TCanvas();

  double minAngle = *std::min_element(referenceAnglesDeg.begin(), referenceAnglesDeg.end());
  double maxAngle = *std::max_element(referenceAnglesDeg.begin(), referenceAnglesDeg.end());
 
  TString hName = TString::Format("hInterpGains_%d_%d", nBins, nBins);
  auto hH = new TH2D(hName, "HPol Gain;#deltaAzimuth(Degree);#deltaElevation(degree);Gain factor (no units)", nBins, minAngle, maxAngle,  nBins, minAngle, maxAngle);
  // auto hV = new TH2D("hV", "VPol Gain;#deltaAzimuth(Degree);#deltaElevation(degree);Gain factor (no units)", nBins, minAngle, maxAngle,  nBins, minAngle, maxAngle);

  for(int by=1; by <= hH->GetNbinsY(); by++){
    double dEl = hH->GetYaxis()->GetBinCenter(by);
    for(int bx=1; bx <= hH->GetNbinsX(); bx++){
      double dAz = hH->GetXaxis()->GetBinCenter(bx);
      double w = dAz + dEl;
      hH->Fill(dEl, dAz, w);
    }
  }  
  
  hH->Draw("colz");
  c->SetLogz(1);
  
  return c;
}



TCanvas* icemc::Seavey::plotGains() {
  loadGains();

  auto c = new TCanvas();

  c->Divide(1,  2);
  auto c1 = c->cd(1);

  TMultiGraph* grBoresight = new TMultiGraph();

  TGraph* gr_h =  new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_dBi_hh_measured.data());  
  TGraph* gr_v =  new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_dBi_vv_measured.data());
  TGraph* gr_hv = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_dBi_hv_measured.data());
  TGraph* gr_vh = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_dBi_vh_measured.data());

  gr_h->SetLineColor(kRed);
  gr_v->SetLineColor(kBlue);
  gr_hv->SetLineColor(kGreen);
  gr_vh->SetLineColor(kYellow);
  
  grBoresight->Add(gr_h);
  grBoresight->Add(gr_v);
  grBoresight->Add(gr_hv);
  grBoresight->Add(gr_vh);

  grBoresight->SetTitle("Seavey Boresight Gains Model in icemc;Frequency (Hz);Gain (dBi?)");
  
  auto l1 = new TLegend();
  l1->AddEntry(gr_h, "HPol Gain", "l");
  l1->AddEntry(gr_v, "VPol Gain", "l");
  l1->AddEntry(gr_hv, "H->V cross-pol Gain", "l");
  l1->AddEntry(gr_vh, "V->H cross-pol Gain", "l");  

  grBoresight->SetBit(kCanDelete);
  grBoresight->Draw("al");

  l1->SetBit(kCanDelete);
  l1->Draw();

  c1->SetLogy(1);

  auto c2 = c->cd(2);
  c2->Divide(2);
  auto c2_1 = c2->cd(1);
  
  TMultiGraph* grAz = new TMultiGraph();
  TMultiGraph* grEl = new TMultiGraph();
  auto l2_1 = new TLegend();
  l2_1->SetNColumns(2);
  auto l2_2 = new TLegend();
  l2_2->SetNColumns(2);
  for(int j=0; j < numAnglePoints; j++){
    TGraph* grH = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_h_angle_az.at(j).data());
    TGraph* grV = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_v_angle_az.at(j).data());

    TGraph* grH2 = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_h_angle_el.at(j).data());
    TGraph* grV2 = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_v_angle_el.at(j).data());
    
    grH->SetLineColor(kRed);
    grV->SetLineColor(kBlue);
    grH->SetLineStyle(j+1);
    grV->SetLineStyle(j+1);
    grH2->SetLineColor(kRed);
    grV2->SetLineColor(kBlue);
    grH2->SetLineStyle(j+1);
    grV2->SetLineStyle(j+1);
    
    grAz->Add(grV);
    grAz->Add(grH);
    grEl->Add(grV2);
    grEl->Add(grH2);

    TString legTextV = TString::Format("VPol %2.0lf^{#circ}", referenceAnglesDeg.at(j));
    TString legTextH = TString::Format("HPol %2.0lf^{#circ}", referenceAnglesDeg.at(j));
    l2_1->AddEntry(grV, legTextV, "l");    
    l2_1->AddEntry(grH, legTextH, "l");

    TString legTextV2 = TString::Format("VPol %2.0lf^{#circ}", referenceAnglesDeg.at(j));
    TString legTextH2 = TString::Format("HPol %2.0lf^{#circ}", referenceAnglesDeg.at(j));
    l2_2->AddEntry(grV2, legTextV2, "l");    
    l2_2->AddEntry(grH2, legTextH2, "l");
  }
  
  grAz->SetTitle("Off boresight gains - Azimuth;Frequency (Hz); Gain (dBi?)");
  grAz->Draw("al");
  grAz->SetBit(kCanDelete);
  l2_1->SetBit(kCanDelete);
  l2_1->Draw();
  c2_1->SetLogy(1);
  auto c2_2 = c2->cd(2);
  grEl->SetTitle("Off boresight gains - Elevation;Frequency (Hz); Gain (dBi?)");
  grEl->Draw("al");
  grEl->SetBit(kCanDelete);
  l2_2->SetBit(kCanDelete);
  l2_2->Draw();
  c2_2->SetLogy(1);
  
  return c;
}


TCanvas* icemc::Seavey::plotHeights(double refractiveIndex) {
  fillHeightArrays(refractiveIndex);

  auto c = new TCanvas();

  TMultiGraph* grBoresight = new TMultiGraph();

  TGraph* gr_h =  new TGraph(freqHz_measured.size(), freqHz_measured.data(), heightHH_m.data());  
  TGraph* gr_v =  new TGraph(freqHz_measured.size(), freqHz_measured.data(), heightVV_m.data());
  TGraph* gr_hv = new TGraph(freqHz_measured.size(), freqHz_measured.data(), heightHV_m.data());
  TGraph* gr_vh = new TGraph(freqHz_measured.size(), freqHz_measured.data(), heightVH_m.data());

  gr_h->SetLineColor(kRed);
  gr_v->SetLineColor(kBlue);
  gr_hv->SetLineColor(kGreen);
  gr_vh->SetLineColor(kYellow);
  
  grBoresight->Add(gr_h);
  grBoresight->Add(gr_v);
  grBoresight->Add(gr_hv);
  grBoresight->Add(gr_vh);

  grBoresight->SetTitle(TString::Format("Seavey Antenna Heights in icemc with n=%4.2lf;Frequency (Hz);Height (m)", refractiveIndex));
  
  auto l1 = new TLegend();
  l1->AddEntry(gr_h, "HPol Height", "l");
  l1->AddEntry(gr_v, "VPol Height", "l");
  l1->AddEntry(gr_hv, "H->V cross-pol Height", "l");
  l1->AddEntry(gr_vh, "V->H cross-pol Height", "l");  

  grBoresight->SetBit(kCanDelete);
  grBoresight->Draw("al");

  l1->SetBit(kCanDelete);
  l1->Draw();

  c->SetLogy(1);
  return c;
}


/** 
 * Linearly interpolate in one of the arrays of length numPoints
 * 
 * @param arr some type with .at() and .size() members whose values correspond to the freqHz_measured array (try std::vector or std::array)
 * @param freqHz the frequency at which to interpolate for
 */

template <typename T>
double linearInterp(const T& arr, double freqHz) {

  const double f0 = freqHz_measured.at(0);
  const double df = freqHz_measured.at(1) - f0;
  
  const int i1 = floor((freqHz - f0)/df);
  const int i2 = i1 + 1;

  const double f1 = f0 + i1*df;

  const double a1 = i1 < 0 || i1 >= arr.size() ? 0 : arr.at(i1);
  const double a2 = i2 < 0 || i2 >= arr.size() ? 0 : arr.at(i2);

  const double m = (a2 - a1)/df;
  const double a_interp = a1 + m*(freqHz - f1);

  return a_interp;
}


double icemc::Seavey::getHeight(Pol pol, double freqHz) const { 
  fillHeightArrays(fRefractiveIndex);

  switch(pol){
  case Pol::V:
    return linearInterp(heightVV_m, freqHz);
  case Pol::H:
    return linearInterp(heightHH_m, freqHz);
  default:
    icemcLog() << icemc::warning << "Requested height for unknown pol, return 0" << std::endl;
    return 0;
  }
}


double icemc::Seavey::getHeight(XPol xPol, double freqHz) const {
  fillHeightArrays(fRefractiveIndex);

  switch(xPol){
  case XPol::VtoH:
    return linearInterp(heightVH_m, freqHz);
  case XPol::HtoV:
    return linearInterp(heightHV_m, freqHz);
  default:
    icemcLog() << icemc::warning << "Requested height for unknown cross-pol, returning 0" << std::endl;
    return 0;
  }
}


/** 
 * Gets the index of the last reference angle less than angleRad
 * 
 * Whether you want Az/El doesn't matter since both are the same.
 * To be used for interpolating off-axis response.
 * If this returns j, then interpolate between j and j+1
 *  
 * @param angleRad is the off-axis angle in radians
 * 
 * @return the index of the last reference angle
 */
int getLowerAngleBin(double angleRad){

  int j1 = -1;
  angleRad = fabs(angleRad);
  for(int j=0; j < referenceAnglesRad.size(); j++){
    if(angleRad >= referenceAnglesRad.at(j)){
      j1 = j;
    }
    else {
      break;
    }
  }
  return j1;
}






bool icemc::Seavey::freqAllowedByPassBands(double freqHz) const {

  if(fPassBandsHz.size()==0){
    return true;
  }

  for(auto p : fPassBandsHz){
    // if(freqHz >= p.first && freqHz < p.second){
    if(freqHz >= p.first && freqHz <= p.second){ ///@todo make a choice here
      return true;
    }
  }
  return false;
}





icemc::Seavey::Seavey(const Settings* settings, double refractiveIndexOfMedium) : // remove include if you move this
  fRefractiveIndex(refractiveIndexOfMedium)
{
  if(settings){
    fPassBandsHz.emplace_back(std::make_pair(settings->FREQ_LOW_SEAVEYS, settings->FREQ_HIGH_SEAVEYS));
    // std::cout << "The pass band (Hz) is " << fPassBandsHz.back().first << "\t" << fPassBandsHz.back().second << std::endl;
  }
}






double icemc::Seavey::getOffAxisResponse(Pol pol,  AngleDir dir, double freqHz, double angleRad) const {
  loadGains();

  const int j1 = getLowerAngleBin(angleRad);
  const int j2 = j1 + 1;

  if(j1 < 0 || j2 >= referenceAnglesRad.size()){
    return 0;
  }
  
  double g1 = 0;
  double g2 = 0;
  switch(pol){
  case Pol::V:
    g1 = dir == AngleDir::Azimuth ? linearInterp(gain_v_angle_az.at(j1), freqHz) : linearInterp(gain_v_angle_az.at(j1), freqHz);
    g2 = dir == AngleDir::Azimuth ? linearInterp(gain_v_angle_az.at(j2), freqHz) : linearInterp(gain_v_angle_az.at(j2), freqHz);
    break;
  case Pol::H:
    g1 = dir == AngleDir::Azimuth ? linearInterp(gain_h_angle_az.at(j1), freqHz) : linearInterp(gain_h_angle_az.at(j1), freqHz);
    g2 = dir == AngleDir::Azimuth ? linearInterp(gain_h_angle_az.at(j2), freqHz) : linearInterp(gain_h_angle_az.at(j2), freqHz);
    break;
  default:
    icemcLog() << icemc::warning << "Requested off-axis gain for for unknown pol, " << (int)pol <<", returning 0\n";
    break;    
  }

  const double a1 = referenceAnglesRad.at(j1);
  const double a2 = referenceAnglesRad.at(j2);

  const double da = a2 - a1;
  const double m = (g2 - g1)/da;
  const double g = m*(fabs(angleRad) - a1) + g1;

  if(fDebug){
    std::cout << a2 << "\t" << a1 << "\t" << g2 << "\t" << g1 << "\t" << m << "\t" << fabs(angleRad) << std::endl;
  }
  
  return g;
}





void icemc::Seavey::addSignal(const icemc::PropagatingSignal& s) {  

  double e_component_kvector = 0;
  double h_component_kvector = 0;
  double n_component_kvector = 0;
  icemc::Seavey::GetEcompHcompkvector(fEPlane,  fHPlane,  fNormal, s.poynting,
				      e_component_kvector, h_component_kvector, n_component_kvector);

  double e_component = 0;
  double h_component = 0;
  double n_component = 0;
  icemc::Seavey::GetEcompHcompEvector(fEPlane, fHPlane, s.polarization, e_component, h_component, n_component);

  double hitangle_e = 0;
  double hitangle_h = 0;
  icemc::Seavey::GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector,
			      hitangle_e, hitangle_h);

  const double df_Hz = s.waveform.getDeltaF();

  // Make copies for the VPol and HPol feeds
  FTPair thisHPol = s.waveform;
  FTPair thisVPol = s.waveform;

  if(fDebug){
    static int ant = -1;
    ant++;
    const char* opt = ant == 0 ? "recreate" : "update";
    TFile* f = TFile::Open("fSeaveysDebug.root", opt);
    f->cd();

    TGraph grV = thisVPol.getTimeDomain();
    grV.SetName(TString::Format("grV_before_%d", ant));
    grV.Write();

    TGraph grH = thisHPol.getTimeDomain();
    grH.SetName(TString::Format("grH_before_%d", ant));
    grH.Write();
    
    f->Write();
    f->Close();
  }

  std::vector<std::complex<double> >& vPolFreqs = thisVPol.changeFreqDomain();
  std::vector<std::complex<double> >& hPolFreqs = thisHPol.changeFreqDomain();

  double freqHz = 0;

  /**
   * In anita.cc, they loop over 0,1 for get_gain_angle(hitangle_e)
   * then they loop over 2, 3 for get_gain_angle(hitangle_h)		 
   * the old index over looping variables goes into the gain_angle arrays
   * 0 is vv_az, 1 is hh_az, 2 is hh_el, 3 is vv_el.
   */

  /**
   * @todo So this is what I think anita.cc was doing			
   * I'm not sure it's right, but if it's wrong, it's the wrong off-axis
   * response of the cross-pol, so probably a small effect		
   */

  bool temp = fDebug;
  fDebug = false;

  for(auto& c : vPolFreqs){    

    if(freqAllowedByPassBands(freqHz)){
    
      const double heightVV = getHeight(Pol::V, freqHz); // VPol component of signal to VPol feed
      const double heightHV = getHeight(XPol::HtoV, freqHz); // HPol component of signal to VPol feed via cross-pol antenna response

      const double offAxisResponseV  = getOffAxisResponse(Pol::V, AngleDir::Azimuth, freqHz, hitangle_e);
      const double offAxisResponseHV = getOffAxisResponse(Pol::H, AngleDir::Azimuth, freqHz, hitangle_e);

      if(fDebug){
        std::cout << "Seavey     \t" << TMath::Nint(freqHz) << "\t" << std::fixed << std::setprecision(7)
      		<< heightVV << "\t" << heightHV << "\t"
      		<< offAxisResponseV << "\t" << offAxisResponseHV << "\t"
      		<< e_component << "\t" << h_component << "\t"
      		<< hitangle_e << "\t" << "\n";
      }

      // 0.5 is for voltage dividing apparently, it doesn't happen in the Seavey... but it does happen downstream... maybe
      const double totalGainFactorV = 0.5*sqrt(  heightVV*heightVV*e_component*e_component*offAxisResponseV
						 + heightHV*heightHV*h_component*h_component*offAxisResponseHV );

      c *= totalGainFactorV;
    }
    else{
      c = 0;
    }
    
    freqHz += df_Hz;
  }

  fDebug = temp;

  
  freqHz = 0; // freqHz is incremented in the loop, so reset
  for(auto& c : hPolFreqs){  

    if(freqAllowedByPassBands(freqHz)){

      // get everything going into the HPol feed... via direct and cross-pol.
      const double heightHH = getHeight(Pol::H, freqHz);
      const double heightVH = getHeight(XPol::VtoH, freqHz);

      // then you need to take acconut of how far off boresight you are... i.e. the off-axis reponse of the antennas.
      const double offAxisResponseH  = getOffAxisResponse(Pol::H, AngleDir::Elevation, freqHz, hitangle_h);
      const double offAxisResponseVH = getOffAxisResponse(Pol::V, AngleDir::Elevation, freqHz, hitangle_h);
      
      // 0.5 is for voltage dividing apparently, it doesn't happen in the Seavey... but it does happen downstream... maybe
      double totalGainFactorH = 0.5*sqrt(  heightHH*heightHH*h_component*h_component*offAxisResponseH
					   + heightVH*heightVH*e_component*e_component*offAxisResponseVH);

      if(fDebug){
        std::cout << "Seavey     \t" << TMath::Nint(freqHz) << "\t" << std::fixed << std::setprecision(7)
      		<< heightHH << "\t" << heightVH << "\t"
      		<< offAxisResponseH << "\t" << offAxisResponseVH << "\t"
      		<< e_component << "\t" << h_component << "\t"
      		<< hitangle_e << "\t" << "\n";
      }
      

      c *= totalGainFactorH;
    }
    else{
      c = 0;
    }
    
    freqHz += df_Hz;
  }


  
  /**
   * @todo In order to make SCREEN stuff work, make this additive rather than just the most recent
   * This will require doing some addition of the waveform. And FTPair doesn't do a simple +=.
   */
  fHPol = thisHPol;
  fVPol = thisVPol;

  fDebug = temp;
  if(fDebug){
    static int ant = -1;
    ant++;
    const char* opt = "update";
    TFile* f = TFile::Open("fSeaveysDebug.root", opt);
    f->cd();

    TGraph grV = fVPol.getTimeDomain();
    grV.SetName(TString::Format("grV_after_%d", ant));
    grV.Write();

    TGraph grH = fHPol.getTimeDomain();
    grH.SetName(TString::Format("grH_after_%d", ant));
    grH.Write();

    f->Write();
    f->Close();
  }

  

  // std::cout << fVPol.getTimeDomain().GetN() << std::endl;  
  
}





const icemc::FTPair& icemc::Seavey::getSignal(Pol pol){
  
  if(pol==Pol::V){
    return fVPol;
  }
  else if(pol==Pol::H){
    return fHPol;
  }
  else{
    icemcLog() << icemc::warning << "Pol for unknown pol requested, returning V.\n";
    return fVPol;
  }  
}









void icemc::Seavey::GetEcompHcompkvector(const Vector& n_eplane, const Vector& n_hplane, const Vector& n_normal,
					 const Vector n_exit2bn,
					 double& e_component_kvector, double& h_component_kvector, double& n_component_kvector) {

  // find component along e-plane for the purpose of finding hit angles, that is, in direction of k vector, direction of radio wave)
  e_component_kvector = -(n_exit2bn.Dot(n_eplane));
  // find component along h-plane for the purpose of finding hit angles, that is, in direction of k vector, direction of radio wave)
  h_component_kvector = -(n_exit2bn.Dot(n_hplane));
  // find the component normal
  n_component_kvector = -(n_exit2bn.Dot(n_normal));

} // end GetEcompHcompkvector



void icemc::Seavey::GetEcompHcompEvector(const Vector& n_eplane, const Vector& n_hplane, const Vector& n_pol,
					 double& e_component, double& h_component, double& n_component) {

  // find component along e-plane in direction of polarization, that is in direction of the E field   
  e_component = n_pol.Dot(n_eplane);
  //    std::cout << "n_component : " << n_exit2bn << " " << n_normal << " " << n_component << std::endl;
    
  // find component along h-plane in direction of polarization, that is in direction of the E field 
  h_component = n_pol.Dot(n_hplane);


  ///@todo maybe restore this at some point?
  // if (settings1->REMOVEPOLARIZATION) {
  //   //Trying to remove effects of polarization at antenna. Stephen
  //   e_component = n_pol.Dot(n_pol);
  //   h_component = 0.001;
  //   n_component = 0.001;
  // } //if
  
} // end GetEcompHcompEvector


void icemc::Seavey::GetHitAngles(double e_component_kvector, double h_component_kvector, double n_component_kvector, double& hitangle_e, double& hitangle_h) {
#if defined(ANITA_UTIL_EXISTS) and defined(VECTORIZE)
  Vec2d y(e_component_kvector, h_component_kvector); 
  Vec2d x(n_component_kvector, n_component_kvector); 
  Vec2d answer = atan2(y,x); 
  hitangle_h = answer[0]; 
  hitangle_e = answer[1]; 

#else
  hitangle_e=atan2(h_component_kvector,n_component_kvector);
  hitangle_h=atan2(e_component_kvector,n_component_kvector);
#endif

}
