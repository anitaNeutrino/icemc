#include "Seavey.h"
#include "EnvironmentVariable.h"
#include "IcemcLog.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "balloon.hh"
#include "Constants.h"

#if defined(ANITA_UTIL_EXISTS) and defined(VECTORIZE)
#include "vectormath_trig.h"
#endif











 // number of points within bandwidth that gain is measured.
constexpr int numGainPoints = 131;
std::array<double, numGainPoints> freqHz_measured; // from 200 MHz to 1.5 GHz with step size of 10 MHz (units in Hz)
std::array<double, numGainPoints> gain_h_measured;
std::array<double, numGainPoints> gain_v_measured;
std::array<double, numGainPoints> gain_hv_measured;
std::array<double, numGainPoints> gain_vh_measured;

constexpr int numAnglePoints = 7;
const std::array<double, numAnglePoints> referenceAnglesDeg {0, 5, 10, 20, 30, 45, 90}; // the off axis measurements are have these step sizes
const std::array<double, numAnglePoints> referenceAnglesRad {referenceAnglesDeg[0]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[1]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[2]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[3]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[4]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[5]*icemc::constants::RADDEG,
                                                             referenceAnglesDeg[6]*icemc::constants::RADDEG};

std::array<std::array<double, numAnglePoints>, numGainPoints> gain_v_angle_az; //[numGainPoints][numAnglePoints]
std::array<std::array<double, numAnglePoints>, numGainPoints> gain_h_angle_az; //[numGainPoints][numAnglePoints]
std::array<std::array<double, numAnglePoints>, numGainPoints> gain_v_angle_el; //[numGainPoints][numAnglePoints]
std::array<std::array<double, numAnglePoints>, numGainPoints> gain_h_angle_el; //[numGainPoints][numAnglePoints]

bool doneLoadGains = false;


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
      icemc::Log() << icemc::error << "Quiting because I can't open " << fileName << "\n";
      exit(1);
    }
    else{
      icemc::Log() << icemc::warning << "I can't open "  << fileName << "\n";
    }
  }
  return dataFile;
}



/** 
 * Read the gains into the vectors.
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

      auto& boresightGain = (fileName == "vv_0" ? gain_h_measured :
			     fileName == "hh_0" ? gain_v_measured :
			     fileName == "hv_0" ? gain_hv_measured :
			                          gain_vh_measured);

      
      for(int i = 0; i < numGainPoints; ++i) {
	double f;
	boresight_file >> f >> boresightGain.at(i);
	if(firstFile){
	  freqHz_measured.at(i) = f;
	}
	else if (f != freqHz_measured.at(i)){
	  icemc::Log() << icemc::warning << "frequency = " << f << ", freqHz_measured[i] = " << freqHz_measured.at(i) << "\n";
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
      // fileName == "hh_el" ? gain_h_angle_el : 

      std::string filePath = ICEMC_DATA_DIR + fileName;
      std::ifstream angle_file = openCarefully(filePath.c_str());

      for(int i = 0; i < numGainPoints; i++){      
	gainVsAngle.at(i).at(0) = 1; // 0th bin is on boresight, so no off-axis effects... so set to 1.
      }

      for(int j = 1; j < numAnglePoints; j++){
	for(int i = 0; i < numGainPoints; i++){
	  double f;
	  angle_file >> f >> gainVsAngle.at(i).at(j);
	  if(f != freqHz_measured.at(i)){
	    icemc::Log() << icemc::warning << "Check frequencies for " << filePath << std::endl;
	  }
	}
      }
      
      angle_file.close();
    }

    doneLoadGains = true;
  }
}





// // reads in the effect of a signal not hitting the antenna straight on
// // also reads in gainxx_measured and sets xxgaintoheight
// void icemc::Seavey::Set_gain_angle(const Settings *settings1,double nmedium_receiver) {
//   loadGains();
  
//   std::string gain_null1, gain_null2;

//   const std::string ICEMC_SRC_DIR=icemc::EnvironmentVariable::ICEMC_SRC_DIR();
//   const std::string ICEMC_DATA_DIR=ICEMC_SRC_DIR+"/data/";

    
  // double gainhv, gainhh, gainvh, gainvv;
  // double gain_step = freqHz_measured[1]-freqHz_measured[0]; // how wide are the steps in frequency;
    
  // icemc::Log() << "GAINS is " << GAINS << "\n";
  // for (int k = 0; k < NFREQ; ++k) {
  //   whichbin[k] = int((freq[k] - freqHz_measured[0]) / gain_step); // finds the gains that were measured for the frequencies closest to the frequency being considered here
  //   if((whichbin[k] >= numGainPoints || whichbin[k] < 0)) {
  //     icemc::Log() << "Set_gain_angle out of range, freq = " << freq[k] << "\twhichbin[k] = " << whichbin[k] << std::endl;
  //     // no longer exit, just set antenna gain to 0 outside band
  //     // @todo verify this works
  //     // exit(1);
  //     scalef2[k] = 0;
  //     scalef1[k] = 0;
  //   }
  //   else{

  //     //now a linear interpolation for the frequency
  //     scalef2[k] = (freq[k] - freqHz_measured[whichbin[k]]) / gain_step;
  //     // how far from the lower frequency
  //     scalef1[k] = 1. - scalef2[k]; // how far from the higher frequency
		
  //     // convert the gain at 0 degrees to the effective antenna height at 0 degrees for every frequency
  //     if(whichbin[k] == numGainPoints - 1) { // if the frequency is 1.5e9 or goes a little over due to rounding
  // 	gainhv = gainhv_measured[whichbin[k]];
  // 	gainhh = gainh_measured[whichbin[k]];
  // 	gainvh = gainvh_measured[whichbin[k]];
  // 	gainvv = gainv_measured[whichbin[k]];
  //     }
  //     else {
  // 	// These gains should be dimensionless numbers, not in dBi
  // 	gainhv = scalef1[k] * gainhv_measured[whichbin[k]] + scalef2[k] * gainhv_measured[whichbin[k] + 1];
  // 	gainhh = scalef1[k] * gainh_measured[whichbin[k]] + scalef2[k] * gainh_measured[whichbin[k] + 1];
  // 	gainvh = scalef1[k] * gainvh_measured[whichbin[k]] + scalef2[k] * gainvh_measured[whichbin[k] + 1];
  // 	gainvv = scalef1[k] * gainv_measured[whichbin[k]] + scalef2[k] * gainv_measured[whichbin[k] + 1];
  //     }
  //   }
  //   if (GAINS==0) {
			
  //     gainhv=0.;
  //     gainhh=gain[1][k];
  //     gainvh=0.;
  //     gainvv=gain[0][k];
			
  //   }
		
  //   hvGaintoHeight[k] = GaintoHeight(gainhv,freq[k],nmedium_receiver);
  //   hhGaintoHeight[k] = GaintoHeight(gainhh,freq[k],nmedium_receiver);
  //   vhGaintoHeight[k] = GaintoHeight(gainvh,freq[k],nmedium_receiver);
  //   vvGaintoHeight[k] = GaintoHeight(gainvv,freq[k],nmedium_receiver);
		
  // } // loop through frequency bins
// } // void Set_gain_angle()











TCanvas* icemc::Seavey::plotGains() {
  loadGains();

  auto c = new TCanvas();

  c->Divide(1,  2);
  auto c1 = c->cd(1);

  TMultiGraph* grBoresight = new TMultiGraph();

  TGraph* gr_h =  new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_h_measured.data());  
  TGraph* gr_v =  new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_v_measured.data());
  TGraph* gr_hv = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_hv_measured.data());
  TGraph* gr_vh = new TGraph(freqHz_measured.size(), freqHz_measured.data(), gain_vh_measured.data());

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
    TGraph* grH = new TGraph(freqHz_measured.size());
    TGraph* grV = new TGraph(freqHz_measured.size());

    TGraph* grH2 = new TGraph(freqHz_measured.size());
    TGraph* grV2 = new TGraph(freqHz_measured.size());
    
    grH->SetLineColor(kRed);
    grV->SetLineColor(kBlue);
    grH->SetLineStyle(j+1);
    grV->SetLineStyle(j+1);
    grH2->SetLineColor(kRed);
    grV2->SetLineColor(kBlue);
    grH2->SetLineStyle(j+1);
    grV2->SetLineStyle(j+1);
    
    for(int i=0; i < freqHz_measured.size(); ++i){
      grH->SetPoint(i, freqHz_measured.at(i), gain_h_angle_az.at(i).at(j));
      grV->SetPoint(i, freqHz_measured.at(i), gain_v_angle_az.at(i).at(j));
      grH2->SetPoint(i, freqHz_measured.at(i), gain_h_angle_el.at(i).at(j));
      grV2->SetPoint(i, freqHz_measured.at(i), gain_v_angle_el.at(i).at(j));
    }
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







void icemc::Seavey::applyAntennaGain(icemc::PropagatingSignal& s) const {
  
  
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
