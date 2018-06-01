#include "Seavey.h"
#include "EnvironmentVariable.h"
#include "IcemcLog.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"

std::vector<double> freq_measured;
std::vector<double> gain_h_measured;
std::vector<double> gain_v_measured;
std::vector<double> gain_hv_measured;
std::vector<double> gain_vh_measured;


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
    icemc::Log() << icemc::error << "Error! can't open "  << fileName << "\n";
    if(failHard){
      exit(1);
    }
  }
  return dataFile;
}



/** 
 * Read the gains into the vectors.
 * NOT THREAD SAFE!
 */
void loadGains(){
  if(freq_measured.size()==0){
    static const int NPOINTS_GAIN = 131; // number of points within bandwidth that gain is measured. from 200 MHz to 1.5 GHz with step size of 10 MHz
    // const double startFreqHz = 200e6;
    // const double deltaFHz = 10e6;

    const std::string ICEMC_SRC_DIR=icemc::EnvironmentVariable::ICEMC_SRC_DIR();
    const std::string ICEMC_DATA_DIR=ICEMC_SRC_DIR+"/data/";

    freq_measured.clear();
    gain_h_measured.clear();
    gain_v_measured.clear();
    gain_hv_measured.clear();
    gain_vh_measured.clear();

    freq_measured.reserve(NPOINTS_GAIN);
    gain_h_measured.reserve(NPOINTS_GAIN);
    gain_v_measured.reserve(NPOINTS_GAIN);
    gain_hv_measured.reserve(NPOINTS_GAIN);
    gain_vh_measured.reserve(NPOINTS_GAIN);

    // gains from university of hawaii measurements.
    std::ifstream gains_h_file = openCarefully((ICEMC_DATA_DIR+"/hh_0").c_str());
    for(int i = 0; i < NPOINTS_GAIN; ++i){
      double f, g;
      gains_h_file >> f >> g;
      freq_measured.emplace_back(f);
      gain_h_measured.emplace_back(g);
      // gainsfile >> freq_measured[iii] >> gain_h_measured[iii];      
    }
    gains_h_file.close();
    
    std::ifstream gains_v_file = openCarefully((ICEMC_DATA_DIR+"/vv_0").c_str()); // gains for vertical polarization
    for(int i = 0; i < NPOINTS_GAIN; ++i) {
      double f, g;
      gains_v_file >> f >> g;
      gain_v_measured.emplace_back(g);
      if(f != freq_measured[i]){
	icemc::Log() << icemc::warning << "frequency = " << f << ", freq_measured[i] = " << freq_measured[i] << "\n";
      }
    }
    gains_v_file.close();

    
    auto gains_hv_file = openCarefully((ICEMC_DATA_DIR+"/hv_0").c_str()); // gains for h-->v cross polarization
    for(int i = 0; i < NPOINTS_GAIN; ++i) {
      double f, g;
      gains_hv_file >> f >> g;//gainhv_measured[iii];
      gain_hv_measured.emplace_back(g);
      if(f != freq_measured[i]){
	icemc::Log() << icemc::warning << "frequency = " << f << ", freq_measured[i] = " << freq_measured[i] << "\n";
      }
    }
    gains_hv_file.close();
    
    auto gains_vh_file = openCarefully((ICEMC_DATA_DIR+"/vh_0").c_str()); // gains for v-->h cross polarization
    for(int i = 0; i < NPOINTS_GAIN; ++i)  {
      double f, g;
      gains_vh_file >> f >> g; //gainvh_measured[iii];
      gain_vh_measured.emplace_back(g);
      if(f != freq_measured[i]){
	icemc::Log() << icemc::warning << "freqeuncy = " << f << ", freq_measured[i] = " << freq_measured[i] << "\n";
      }
    }
    gains_vh_file.close();

  }
}


TCanvas* icemc::Seavey::plotGains() const {
  loadGains();

  auto c = new TCanvas();
  TMultiGraph* gr = new TMultiGraph();
  
  TGraph* gr_h = new TGraph(freq_measured.size(), &freq_measured[0], &gain_h_measured[0]);  
  TGraph* gr_v = new TGraph(freq_measured.size(), &freq_measured[0], &gain_v_measured[0]);
  TGraph* gr_hv = new TGraph(freq_measured.size(), &freq_measured[0], &gain_hv_measured[0]);
  TGraph* gr_vh = new TGraph(freq_measured.size(), &freq_measured[0], &gain_vh_measured[0]);

  gr_h->SetLineColor(kRed);
  gr_v->SetLineColor(kBlue);
  gr_hv->SetLineColor(kGreen);
  gr_vh->SetLineColor(kYellow);
  
  gr->Add(gr_h);
  gr->Add(gr_v);
  gr->Add(gr_hv);
  gr->Add(gr_vh);

  gr->SetTitle("Seavey Boresight Gains Model in icemc;Frequency(Hz);Gain (dBi?)");
  
  TLegend* l = new TLegend();
  l->AddEntry(gr_h, "HPol Gain", "l");
  l->AddEntry(gr_v, "VPol Gain", "l");
  l->AddEntry(gr_hv, "H->V cross-pol Gain", "l");
  l->AddEntry(gr_vh, "V->H cross-pol Gain", "l");  

  gr->SetBit(kCanDelete);
  gr->Draw("al");

  l->SetBit(kCanDelete);
  l->Draw();

  c->SetLogy(1);
  
  return c;
}

