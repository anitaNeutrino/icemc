#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
//#include "Constants.h"
#include <math.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"

#include "TText.h"

#include "TGraphAsymmErrors.h"
//#include "vector.h"
#include <vector>
#include "TMath.h"

const double PI = TMath::Pi();
const double poissonerror_plus[20] = {1};
const double poissonerror_minus[20] = {1};

TStyle* RootStyle();

void LogToLine(int N, double *Data);

using namespace std;


// copy error find function from AraSim
void findErrorOnSumWeights(double *eventsfound_binned,double &error_plus,double &error_minus);

// make plot (with info I already have)
void MakePlot();

// apply factor to array
void ApplyFactorArray( double *output, double *input, double *factor, int size );

// get average factor difference between arrays
double GetAvgFactor( double *EArray1, double *FluxArray1, int bin1, double *EArray2, double *FluxArray2, int bin2 );

// get RMS factor difference between arrays
double GetRMSFactor( double *EArray1, double *FluxArray1, int bin1, double *EArray2, double *FluxArray2, int bin2, double AvgFactor );


int main (int argc, char **argv) {





  // TStyle *color=RootStyle();
  // gStyle=color;



  // double total_weight = 0.;
  // double bin_val[10];

  // // initialize
  // for (int i=0; i<10; i++) {
  //   bin_val[i] = 0.;
  // }

  // int total_run = 0;


  // for (int run=1; run<argc; run++) {

  //   total_run++;
  //   cout<<"run : "<<run<<" total_run : "<<total_run<<endl;

  //   ifstream ifile (argv[run]);

  //   string line, label;

  //   char delim[] = " \n";
  //   char *token;
  //   char *line_char;

  //   if (ifile.is_open()) {
  //     //while ( ifile.good() ) {

  //     getline( ifile, line );

  //     label = line.substr(0,line.find_first_of("="));
  //     if (label == "Total_Weight") {  // read total weight from the file
  // 	total_weight += atof( line.substr(line.find_first_of("=")+1).c_str() );
  //     }

  //     getline( ifile, line ); // next line should have bin values
  //     line_char = &(line[0]);
  //     //memcpy(line_char,line.c_str(),line.size());

  //     token = strtok (line_char, delim);
  //     bin_val[0] += atof(token);
  //     //cout<<"bin val : "<<bin_val<<endl;

  //     for (int i=0; i<9; i++) {

  // 	token = strtok (NULL, delim);
  // 	bin_val[i+1] += atof(token);
  // 	//cout<<"bin val : "<<bin_val<<endl;
  //     }


  //     //}
  //     ifile.close();
  //   } // end if file open


  // }   // end run loop


  // int total_gt_events = 0;

  // cout<<"total_weight : "<<total_weight<<endl;
  // for (int i=0; i<10; i++) {
  //   cout<<bin_val[i]<<" ";
  //   total_gt_events += bin_val[i];
  // }
  // cout<<"\n";
  // cout<<"Total global triggered events : "<<total_gt_events<<endl;
  // double error_plus = 0;
  // double error_minus = 0;
  // findErrorOnSumWeights( bin_val, error_plus, error_minus );
  // cout<<"error plus : "<<error_plus<<" error minus : "<<error_minus<<endl;

  // //double Veff_factor_km3sr = 977.446;
  // double Veff_factor_km3sr;
  // //double IceVolume = 1.35717e+12; // get this from output log file, this one is for radius 12000m
  // double IceVolume = 3.39292e+11; // get this from output log file, this one is for radius 12000m
  // double RHOICE = 917; // ice density kg/m^3
  // double RHOH20 = 1000; // water density kg/m^3

  // double Veff_other_factors;
  // int NNU = 10000;
  // //int NNU = 50000;
  // cout<<"NNU assumed as : "<<NNU<<". CHECK!!"<<endl;

  // Veff_factor_km3sr = IceVolume * 4. * PI * RHOICE / RHOH20;
  // Veff_other_factors = 1./((double)total_run * (double) NNU);

  // //cout<<"\nVeff : "<<Veff_factor_km3sr * Veff_other_factors * total_weight<<" km3sr"<<endl;
  // cout<<"\nVeff : "<<Veff_factor_km3sr * Veff_other_factors * total_weight * 1.E-9<<" km3sr"<<endl;
  // cout<<"Veff_error_plus : "<<Veff_factor_km3sr * Veff_other_factors * error_plus * 1.E-9<<"    minus : "<<Veff_factor_km3sr * Veff_other_factors * error_minus * 1.E-9<<endl;



  MakePlot();


  return 0;
}



void findErrorOnSumWeights(double *eventsfound_binned,double &error_plus,double &error_minus) {
  // below 3 values are from AraSim counting class
  int NBINS = 10;
  double MAX_LOGWEIGHT = 0;
  double MIN_LOGWEIGHT = -3;

  for (int i=0;i<NBINS;i++) { // we are going to sum the error on the weight squared over the weight bins
    // in each bin, the error on the weight is the weight for that bin times the error on the number of events in that bin
    double thislogweight=((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT; // find log weight for this bin
    if (eventsfound_binned[i]<=20) {  // if the number of events in this bin <20, use poisson errors
      error_plus+=pow(poissonerror_plus[(int)eventsfound_binned[i]]*pow(10.,thislogweight),2);
      error_minus+=pow(poissonerror_minus[(int)eventsfound_binned[i]]*pow(10.,thislogweight),2);
    }
    else {// otherwise, use sqrt(n) errors
      error_plus+=eventsfound_binned[i]*pow(pow(10.,thislogweight),2);
      error_minus=error_plus;

    }

  }
  error_plus=sqrt(error_plus); // take the sqrt of the sum of the squares
  error_minus=sqrt(error_minus);
}


void MakePlot() {

  double x[2], xerr[2];    // x values shared

  double y_vhvh[2], yerr_vhvh[2]; // for vhvh

  double y_vhv[2], yerr_vhv[2]; // for vhv

  double y_vhvv[2], yerr_vhvv[2]; // for vhvv

  double y_vhhh[2], yerr_vhhh[2]; // for vhhh

  double y_vhh[2], yerr_vhh[2]; // for vhh


  // set shared values
  x[0] = 17.;
  x[1] = 18.;

  xerr[0] = 0.;
  xerr[1] = 0.;

  // set values are each cases
  y_vhvh[0] = 2.0231;
  y_vhvh[1] = 8.5754;
  yerr_vhvh[0] = 0.0353;
  yerr_vhvh[1] = 0.07113;

  y_vhv[0] = 1.7832;
  y_vhv[1] = 8.0907;
  yerr_vhv[0] = 0.03019;
  yerr_vhv[1] = 0.06490;

  y_vhvv[0] = 2.1033;
  y_vhvv[1] = 8.7296;
  yerr_vhvv[0] = 0.03503;
  yerr_vhvv[1] = 0.07269;

  y_vhhh[0] = 2.0313;
  y_vhhh[1] = 8.2336;
  yerr_vhhh[0] = 0.03393;
  yerr_vhhh[1] = 0.066545;

  y_vhh[0] = 1.9577;
  y_vhh[1] = 7.8734;
  yerr_vhh[0] = 0.0571;
  yerr_vhh[1] = 0.06335;

  TGraphErrors *g_vhvh = new TGraphErrors( 2, x, y_vhvh, xerr, yerr_vhvh );

  TGraphErrors *g_vhv = new TGraphErrors( 2, x, y_vhv, xerr, yerr_vhv );

  TGraphErrors *g_vhvv = new TGraphErrors( 2, x, y_vhvv, xerr, yerr_vhvv );

  TGraphErrors *g_vhhh = new TGraphErrors( 2, x, y_vhhh, xerr, yerr_vhhh );

  TGraphErrors *g_vhh = new TGraphErrors( 2, x, y_vhh, xerr, yerr_vhh );


  TCanvas *cVeff = new TCanvas("cVeff","A Simple Graph Example",200,10,1000,700);

  cVeff->cd();

  cVeff->SetLogy();
  cVeff->SetGrid();
  g_vhvh->SetTitle("Effective volume from ARA station 1 (~1M evt)");
  g_vhvh->GetHistogram()->SetXTitle("log E");
  g_vhvh->GetHistogram()->SetYTitle("km^3 sr");
  //g_vhvh->GetHistogram()->SetMaximum(100);
  //g_vhvh->GetHistogram()->SetMinimum(0.1);
  g_vhvh->GetHistogram()->SetMaximum(12.5);
  g_vhvh->GetHistogram()->SetMinimum(1.5);
  g_vhvh->Draw("al*");

  g_vhv -> SetLineColor(kRed);
  g_vhv->Draw("l*");

  g_vhvv -> SetLineColor(kGreen);
  g_vhvv->Draw("l*");

  g_vhhh -> SetLineColor(kBlue);
  g_vhhh->Draw("l*");

  g_vhh -> SetLineColor(kPink);
  g_vhh -> SetLineWidth(2);
  g_vhh -> SetLineStyle(7);
  g_vhh->Draw("l*");


  TLegend *Leg_Veff = new TLegend(1., 0.75, 0.85,0.6);
  Leg_Veff -> AddEntry(g_vhvh, "VHVH", "l");
  Leg_Veff -> AddEntry(g_vhv, "VHV", "l");
  Leg_Veff -> AddEntry(g_vhvv, "VHVV", "l");
  Leg_Veff -> AddEntry(g_vhhh, "VHHH", "l");
  Leg_Veff -> AddEntry(g_vhh, "VHH", "l");
  Leg_Veff -> Draw();

  cVeff->Print("test_Veff_ara1.pdf");



  double x_5p[5], xerr_5p[5];
  double y_vhvh_5p[5], yerr_vhvh_5p[5];
  x_5p[0] = 16;
  x_5p[1] = 17;
  x_5p[2] = 18;
  x_5p[3] = 19;
  x_5p[4] = 20;

  xerr_5p[0] = 0;
  xerr_5p[1] = 0;
  xerr_5p[2] = 0;
  xerr_5p[3] = 0;
  xerr_5p[4] = 0;


  // with signal and noise (after fix ndepth)
  y_vhvh_5p[0] = 11.0418;
  y_vhvh_5p[1] = 71.6261;
  y_vhvh_5p[2] = 266.506;
  y_vhvh_5p[3] = 504.048;
  y_vhvh_5p[4] = 729.863;
  yerr_vhvh_5p[0] = 0.29644;
  yerr_vhvh_5p[1] = 0.7438;
  yerr_vhvh_5p[2] = 1.445;
  yerr_vhvh_5p[3] = 2.001;
  yerr_vhvh_5p[4] = 2.409;


  TGraphErrors *g_vhvh_37 = new TGraphErrors( 5, x_5p, y_vhvh_5p, xerr_5p, yerr_vhvh_5p );


  // now with only signal (same threshold) (I think this is before fixing ndepth)
  double y_vhvh_5p_PS[5], yerr_vhvh_5p_PS[5];

  y_vhvh_5p_PS[0] = 5.877;
  y_vhvh_5p_PS[1] = 64.041;
  y_vhvh_5p_PS[2] = 256.477;
  y_vhvh_5p_PS[3] = 488.454;
  y_vhvh_5p_PS[4] = 724.319;
  yerr_vhvh_5p_PS[0] = 0.218;
  yerr_vhvh_5p_PS[1] = 0.704;
  yerr_vhvh_5p_PS[2] = 1.417;
  yerr_vhvh_5p_PS[3] = 1.959;
  yerr_vhvh_5p_PS[4] = 2.399;


  TGraphErrors *g_vhvh_37_PS = new TGraphErrors( 5, x_5p, y_vhvh_5p_PS, xerr_5p, yerr_vhvh_5p_PS );


  // with signal and noise (setndepth fixed, signal+noise - pure noise)
  double y_vhvh_5p_N[5], yerr_vhvh_5p_N[5];

  y_vhvh_5p_N[0] = 5.4086;
  y_vhvh_5p_N[1] = 66.34495;
  y_vhvh_5p_N[2] = 262.143;
  y_vhvh_5p_N[3] = 500.124;
  y_vhvh_5p_N[4] = 725.962;
  //yerr_vhvh_5p_N[0] = 0.299;
  //yerr_vhvh_5p_N[1] = 0.754;
  //yerr_vhvh_5p_N[2] = 1.446;
  //yerr_vhvh_5p_N[3] = 1.999;
  //yerr_vhvh_5p_N[4] = 2.419;


  //TGraphErrors *g_vhvh_37_N = new TGraphErrors( 5, x_5p, y_vhvh_5p_N, xerr_5p, yerr_vhvh_5p_N );
  TGraph *g_vhvh_37_N = new TGraph( 5, x_5p, y_vhvh_5p_N );


  // Peter's Veff plot
  double y_vhvh_5p_Peter[5];

  y_vhvh_5p_Peter[0] = 0.7;
  y_vhvh_5p_Peter[1] = 40.;
  y_vhvh_5p_Peter[2] = 250.;
  y_vhvh_5p_Peter[3] = 500.;
  y_vhvh_5p_Peter[4] = 700.;

  TGraph *g_vhvh_37_Peter = new TGraph( 5, x_5p, y_vhvh_5p_Peter );



  TCanvas *cVeff_37 = new TCanvas("cVeff_37","A Simple Graph Example",200,10,1000,700);

  cVeff_37->cd();

  cVeff_37->SetLogy();
  cVeff_37->SetGrid();
  g_vhvh_37->SetTitle("Effective volume from ARA station 37 (~1M evt)");
  g_vhvh_37->GetHistogram()->SetXTitle("log E");
  g_vhvh_37->GetHistogram()->SetYTitle("km^3 sr");
  //g_vhvh->GetHistogram()->SetMaximum(100);
  //g_vhvh->GetHistogram()->SetMinimum(0.1);
  g_vhvh_37->GetHistogram()->SetMaximum(3000);
  g_vhvh_37->GetHistogram()->SetMinimum(0.3);
  g_vhvh_37->Draw("al*");


  //g_vhvh_37_PS->SetLineColor(kRed);
  //g_vhvh_37_PS->Draw("l*");


  //g_vhvh_37_N->SetLineColor(kBlue);
  g_vhvh_37_N->SetLineColor(kRed);
  g_vhvh_37_N->Draw("l*");

  g_vhvh_37_Peter->SetLineColor(kBlue);
  g_vhvh_37_Peter->Draw("l*");


  TLegend *Leg_Veff37 = new TLegend(1., 1., 0.8,0.8);
  Leg_Veff37 -> AddEntry(g_vhvh_37, "Signal+noise", "l");
  Leg_Veff37 -> AddEntry(g_vhvh_37_N, "(Signal+noise) - Noise", "l");
  Leg_Veff37 -> AddEntry(g_vhvh_37_Peter, "Peter's", "l");
  //Leg_Veff37 -> AddEntry(g_vhvh_37_N, "Signal+noise (fixed)", "l");
  //Leg_Veff37 -> AddEntry(g_vhvh_37_PS, "Pure signal", "l");
  Leg_Veff37 -> Draw();

  cVeff_37->Print("test_Veff_ara37.pdf");






  // Veff for Testbed comparison between RARA (bug fixed?) and AraSim (old RF mode)
  //

  double Veff_RARA_x[9] = { 16.5086,     17.0005,         17.5078,         17.9921,         18.4926,         18.9853,         19.4938,         19.9946,         20.4954 };
  //double Veff_RARA_y[9] = { -5.07826,     -0.748688,         4.53638,         7.60887,         9.17314,         10.3351,         11.095,         11.2514,         11.5587 };
  double Veff_RARA_y[9] = { 0.210001,     0.63908,         2.47996,         5.34291,         7.94328,         10.6605,         12.7513,         13.4208,         14.3072 };

  TGraph *g_Veff_Testbed_RARA = new TGraph( 9, Veff_RARA_x, Veff_RARA_y );


  double Veff_AraSim_oldRF_x[9] = { 17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5, 21 };
  double Veff_AraSim_oldRF_y[9] = { 0.894408, 2.02302, 3.27728, 5.22782, 7.62459, 10.2982, 12.8196, 16.2872, 18.1278 };

  double Veff_AraSim_newRF_x[9] = { 17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5, 21 };
  double Veff_AraSim_newRF_y[9] = { 0.261292, 0.751385, 1.803, 3.51358, 5.72288, 8.54767, 11.4697, 14.3365, 17.1533 };

  // LPM effect constant
  //
  // 17: 0.8689
  // 17.5: 0.8204
  // 18: 0.7719
  // 18.5: 0.7451
  // 19: 0.7183
  // 19.5: 0.715
  // 20: 0.712
  // 20.5: 0.754
  // 21: 0.796

  double Veff_AraSim_newRF_LPMeffect[9] = { 0.8689, 0.8204, 0.7719, 0.7451, 0.7183, 0.715, 0.712, 0.754, 0.796 };
  for (int i=0; i<9; i++) {

    Veff_AraSim_newRF_y[i] = Veff_AraSim_newRF_y[i] * Veff_AraSim_newRF_LPMeffect[i];
  }



  TGraph *g_Veff_Testbed_AraSim_oldRF = new TGraph( 9, Veff_AraSim_oldRF_x, Veff_AraSim_oldRF_y );

  TGraph *g_Veff_Testbed_AraSim_newRF = new TGraph( 9, Veff_AraSim_newRF_x, Veff_AraSim_newRF_y );



  TCanvas *cVeff_Testbed_compare = new TCanvas("cVeff_Testbed_compare","",800,600);

  cVeff_Testbed_compare->cd();

  cVeff_Testbed_compare->SetLogy();
  cVeff_Testbed_compare->SetGrid();

  g_Veff_Testbed_RARA->SetTitle("Testbed Effective Volume");
  g_Veff_Testbed_RARA->GetHistogram()->SetXTitle("log E");
  g_Veff_Testbed_RARA->GetHistogram()->SetYTitle("km^{3} sr (ice volume)     ");
  g_Veff_Testbed_RARA->GetHistogram()->SetMaximum(30);
  g_Veff_Testbed_RARA->GetHistogram()->SetMinimum(1.E-1);
  g_Veff_Testbed_RARA->GetXaxis()->SetLimits(16,21);
  g_Veff_Testbed_RARA->SetLineWidth(3);

  g_Veff_Testbed_RARA->GetHistogram()->SetTitleSize  ( 0.04,"X");
  g_Veff_Testbed_RARA->GetHistogram()->SetLabelOffset( 0.006,"X");
  g_Veff_Testbed_RARA->GetHistogram()->SetLabelSize( 0.04,"X");
  g_Veff_Testbed_RARA->GetHistogram()->SetTitleOffset( 1.0,"X");
  g_Veff_Testbed_RARA->GetHistogram()->SetLabelSize( 0.04,"Y");
  g_Veff_Testbed_RARA->GetHistogram()->SetLabelOffset( 0.007,"Y");
  g_Veff_Testbed_RARA->GetHistogram()->SetTitleSize  ( 0.04,"Y");
  //g_Veff_Testbed_RARA->GetHistogram()->SetTitleOffset( 2.0,"Y");
  g_Veff_Testbed_RARA->GetHistogram()->SetTitleOffset( 1.0,"Y");

  g_Veff_Testbed_RARA->Draw("al*");

  g_Veff_Testbed_AraSim_oldRF->SetLineColor(4);
  g_Veff_Testbed_AraSim_oldRF->SetLineWidth(3);
  g_Veff_Testbed_AraSim_oldRF->Draw("l*");

  g_Veff_Testbed_AraSim_newRF->SetLineColor(2);
  g_Veff_Testbed_AraSim_newRF->SetLineWidth(3);
  g_Veff_Testbed_AraSim_newRF->Draw("l*");

  TLegend *Leg_Veff_Testbed = new TLegend(0.95, 0.15, 0.5,0.45);
  Leg_Veff_Testbed-> AddEntry(g_Veff_Testbed_RARA, "RARA", "l*");
  Leg_Veff_Testbed-> AddEntry(g_Veff_Testbed_AraSim_oldRF, "AraSim (oldRF mode w/ LPM on)", "l*");
  Leg_Veff_Testbed-> AddEntry(g_Veff_Testbed_AraSim_newRF, "AraSim (newRF mode w/ LPM account)", "l*");
  Leg_Veff_Testbed-> Draw();

  cVeff_Testbed_compare->Print("Veff_Testbed_AraSim_RARA.pdf");







  // UCL analysis efficiency plot (digitized)
  //

  // use same x bin
  double UCL_x[20] = { 2.50867,     3.51875,         4.49976,         5.5093,         6.50453,         7.52891,         8.52487,         9.5351,         10.5307,         11.4979,         12.508,         13.5181,         14.5141,         15.5099,         16.5057,         17.5303,         18.5115,         19.5072,         20.5173,         21.513 };
  double UCL_CW_y[20] = { 0.996373,     0.995104,         0.982318,         0.951096,         0.926787,         0.918605,         0.934619,         0.941415,         0.938996,         0.959623,         0.957202,         0.960542,         0.97886,         0.98681,         0.993608,         0.995794,         0.994529,         0.996718,         0.996602,         0.996487 };

  double UCL_MPROB_y[20] = { 0.95605,     0.953631,         0.975402,         0.952248,         0.925635,         0.913998,         0.927709,         0.939115,         0.930934,         0.944642,         0.938769,         0.938652,         0.922409,         0.869297,         0.787385,         0.705473,         0.534851,         0.413766,         0.391764,         0.32713 };

  double UCL_CHISQ_y[20] = { 0.245219,     0.353396,         0.802593,         0.906165,         0.92794,         0.937038,         0.943834,         0.947177,         0.950517,         0.951558,         0.949136,         0.957084,         0.952363,         0.948788,         0.944067,         0.931278,         0.961117,         0.962153,         0.962041,         0.957312 };

  double UCL_POW_y[20] = { 0.211809,     0.285425,         0.660887,         0.953398,         0.997064,         0.998098,         0.999136,         0.999019,         0.997754,         0.998791,         0.998677,         0.99856,         0.998447,         0.998331,         0.997062,         0.996948,         0.996833,         0.99787,         0.997756,         0.996489 };

  double UCL_GEOM_y[20] = { 0.403054,     0.453629,         0.57909,         0.602014,         0.636464,         0.640955,         0.637385,         0.643029,         0.652131,         0.654321,         0.676094,         0.693262,         0.688539,         0.668837,         0.645679,         0.695104,         0.713422,         0.725981,         0.715496,         0.710771 };


  double UCL_ALLCUTS_y[20] = { 0.0113479,     0.0377276,         0.278397,         0.40962,         0.436002,         0.457776,         0.470333,         0.489803,         0.517339,         0.544874,         0.555128,         0.614921,         0.566417,         0.509851,         0.441764,         0.424366,         0.349367,         0.281279,         0.258121,         0.206162 };

  double UCL_empty[20] = { 0. };
  double UCL_xbar[20] = { 0.5 };


  //TGraph *gUCL_CW = new TGraph ( 20, UCL_x, UCL_CW_y );

  TGraphAsymmErrors *g_UCL_CW = new TGraphAsymmErrors( 20, UCL_x, UCL_CW_y, UCL_xbar, UCL_xbar, UCL_empty, UCL_empty ); // only x range bar

  TGraphAsymmErrors *g_UCL_MPROB = new TGraphAsymmErrors( 20, UCL_x, UCL_MPROB_y, UCL_xbar, UCL_xbar, UCL_empty, UCL_empty ); // only x range bar

  TGraphAsymmErrors *g_UCL_CHISQ = new TGraphAsymmErrors( 20, UCL_x, UCL_CHISQ_y, UCL_xbar, UCL_xbar, UCL_empty, UCL_empty ); // only x range bar

  TGraphAsymmErrors *g_UCL_POW = new TGraphAsymmErrors( 20, UCL_x, UCL_POW_y, UCL_xbar, UCL_xbar, UCL_empty, UCL_empty ); // only x range bar

  TGraphAsymmErrors *g_UCL_GEOM = new TGraphAsymmErrors( 20, UCL_x, UCL_GEOM_y, UCL_xbar, UCL_xbar, UCL_empty, UCL_empty ); // only x range bar

  TGraphAsymmErrors *g_UCL_ALLCUTS = new TGraphAsymmErrors( 20, UCL_x, UCL_ALLCUTS_y, UCL_xbar, UCL_xbar, UCL_empty, UCL_empty ); // only x range bar


  TCanvas *cUCL_Efficiency = new TCanvas("cUCL_Efficiency","",800,600);

  cUCL_Efficiency->cd();

  g_UCL_CW->SetTitle("UCL CSW Efficiencies");
  g_UCL_CW->GetHistogram()->SetXTitle("SNR");
  g_UCL_CW->GetHistogram()->SetYTitle("Efficiency");
  g_UCL_CW->GetHistogram()->SetMaximum(1.5);
  g_UCL_CW->GetHistogram()->SetMinimum(0.);
  g_UCL_CW->GetXaxis()->SetLimits(1,23);
  g_UCL_CW->SetLineWidth(3);

  g_UCL_CW->GetHistogram()->SetTitleSize  ( 0.04,"X");
  g_UCL_CW->GetHistogram()->SetLabelOffset( 0.006,"X");
  g_UCL_CW->GetHistogram()->SetLabelSize( 0.04,"X");
  g_UCL_CW->GetHistogram()->SetTitleOffset( 1.0,"X");
  g_UCL_CW->GetHistogram()->SetLabelSize( 0.04,"Y");
  g_UCL_CW->GetHistogram()->SetLabelOffset( 0.007,"Y");
  g_UCL_CW->GetHistogram()->SetTitleSize  ( 0.04,"Y");
  //g_UCL_CW->GetHistogram()->SetTitleOffset( 2.0,"Y");
  g_UCL_CW->GetHistogram()->SetTitleOffset( 1.0,"Y");

  g_UCL_CW->Draw("ap");

  /*
    g_Veff_Testbed_AraSim_oldRF->SetLineColor(4);
    g_Veff_Testbed_AraSim_oldRF->SetLineWidth(3);
    g_Veff_Testbed_AraSim_oldRF->Draw("l*");

    g_Veff_Testbed_AraSim_newRF->SetLineColor(2);
    g_Veff_Testbed_AraSim_newRF->SetLineWidth(3);
    g_Veff_Testbed_AraSim_newRF->Draw("l*");
  */

  /*
    TLegend *Leg_Veff_Testbed = new TLegend(0.95, 0.15, 0.5,0.45);
    Leg_Veff_Testbed-> AddEntry(g_Veff_Testbed_RARA, "RARA", "l*");
    Leg_Veff_Testbed-> AddEntry(g_Veff_Testbed_AraSim_oldRF, "AraSim (oldRF mode)", "l*");
    Leg_Veff_Testbed-> AddEntry(g_Veff_Testbed_AraSim_newRF, "AraSim (newRF mode)", "l*");
    Leg_Veff_Testbed-> Draw();
  */

  cUCL_Efficiency->Print("UCL_Cut_Efficiency.pdf");













  // effective are from ARA37 stations
  //
  double vhvh_Aeff_37[5];

  vhvh_Aeff_37[0] = 1.21E-3;
  vhvh_Aeff_37[1] = 1.83E-2;
  vhvh_Aeff_37[2] = 1.66E-1;
  vhvh_Aeff_37[3] = 7.74E-1;
  vhvh_Aeff_37[4] = 2.77;

  TGraph *g_vhvh_Aeff_37 = new TGraph( 5, x_5p, vhvh_Aeff_37 );

  TCanvas *cAeff_37 = new TCanvas("cAeff_37","A Simple Graph Example",200,10,1000,700);

  cAeff_37->cd();

  cAeff_37->SetLogy();
  cAeff_37->SetGrid();

  g_vhvh_Aeff_37->SetTitle("Effective area from ARA station 37 (~1M evt)");
  g_vhvh_Aeff_37->GetHistogram()->SetXTitle("log E");
  g_vhvh_Aeff_37->GetHistogram()->SetYTitle("km^2 sr");
  //g_vhvh->GetHistogram()->SetMaximum(100);
  //g_vhvh->GetHistogram()->SetMinimum(0.1);
  g_vhvh_Aeff_37->GetHistogram()->SetMaximum(10);
  g_vhvh_Aeff_37->GetHistogram()->SetMinimum(1.E-4);
  g_vhvh_Aeff_37->Draw("al*");

  cAeff_37->Print("test_Aeff_ara37.pdf");





  // TestBed '11 from Davez -> Aeff calculated from sensitivity value
  //
  //double TestBed_Dave_Aeff_x[6] = { 16.6176, 16.8607, 17.1207, 18.2206, 19.2074, 20.9303}; // from livetime 5month assumption
  //double TestBed_Dave_Aeff_y[6] = { 2.32e-5, 0.000189, 0.0005996, 0.01487, 0.07518, 2.39};
  double TestBed_Dave_Aeff_x[6] = { 16.6176, 16.8607, 17.1207, 18.2206, 19.2074, 20.9303}; // from livetime 2 yrs assumption
  double TestBed_Dave_Aeff_y[6] = { 4.838e-6, 3.934e-5, 0.000125, 0.003098, 0.015663, 0.498 };

  TGraph *gAeff_TestBed_Dave = new TGraph (6, TestBed_Dave_Aeff_x, TestBed_Dave_Aeff_y);

  TGraph *gAeff_TestBed_Dave_2 = new TGraph (6, TestBed_Dave_Aeff_x, TestBed_Dave_Aeff_y);

  TGraph *gAeff_TestBed_Dave_3 = new TGraph (6, TestBed_Dave_Aeff_x, TestBed_Dave_Aeff_y);



  // ARA37 Aeffs
  TGraph *gAeff_ARA37[2];
  double ARA37_Aeff_Fdomain_x[9] = { 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20 };
  double ARA37_Aeff_Fdomain_y[9] = { 0.000738727, 0.00457682, 0.0211768, 0.075472, 0.202656, 0.482979, 0.957373, 1.58857, 2.93544 };
  gAeff_ARA37[0] = new TGraph (9, ARA37_Aeff_Fdomain_x, ARA37_Aeff_Fdomain_y);

  double ARA37_Aeff_Tdomain_x[6] = { 17.5, 18, 18.5, 19, 19.5, 20 };
  double ARA37_Aeff_Tdomain_y[6] = { 0.0174291, 0.0689039, 0.191027, 0.404153, 0.734779, 1.21211 };
  gAeff_ARA37[1] = new TGraph (6, ARA37_Aeff_Tdomain_x, ARA37_Aeff_Tdomain_y);


  // ARA10 Aeffs
  TGraph *gAeff_ARA10[2];
  double ARA10_Aeff_Fdomain_x[4] = { 17, 18, 19, 20 };
  double ARA10_Aeff_Fdomain_y[4] = { 0.00570696, 0.0624085, 0.319845, 1.05987 };
  gAeff_ARA10[0] = new TGraph (4, ARA10_Aeff_Fdomain_x, ARA10_Aeff_Fdomain_y);

  double ARA10_Aeff_Tdomain_x[6] = { 17.5, 18, 18.5, 19, 19.5, 20 };
  double ARA10_Aeff_Tdomain_y[6] = { 0.00466818, 0.0183012, 0.0537196, 0.122696, 0.233236, 0.399833 };
  gAeff_ARA10[1] = new TGraph (6, ARA10_Aeff_Tdomain_x, ARA10_Aeff_Tdomain_y);


  double avg_factor_37_10_Tdomain = GetAvgFactor( ARA37_Aeff_Tdomain_x, ARA37_Aeff_Tdomain_y, 6, ARA10_Aeff_Tdomain_x, ARA10_Aeff_Tdomain_y, 6 );
  double rms_factor_37_10_Tdomain = GetRMSFactor( ARA37_Aeff_Tdomain_x, ARA37_Aeff_Tdomain_y, 6, ARA10_Aeff_Tdomain_x, ARA10_Aeff_Tdomain_y, 6, avg_factor_37_10_Tdomain );

  double avg_factor_37_10_Fdomain = GetAvgFactor( ARA37_Aeff_Fdomain_x, ARA37_Aeff_Fdomain_y, 9, ARA10_Aeff_Fdomain_x, ARA10_Aeff_Fdomain_y, 4 );
  double rms_factor_37_10_Fdomain = GetRMSFactor( ARA37_Aeff_Fdomain_x, ARA37_Aeff_Fdomain_y, 9, ARA10_Aeff_Fdomain_x, ARA10_Aeff_Fdomain_y, 4, avg_factor_37_10_Fdomain );

  cout<<"\nAvgFactor from ARA37 -> 10, Tdomain : "<<avg_factor_37_10_Tdomain<<", Fdomain : "<<avg_factor_37_10_Fdomain<<endl;
  cout<<"RMSFactor from ARA37 -> 10, Tdomain : "<<rms_factor_37_10_Tdomain<<", Fdomain : "<<rms_factor_37_10_Fdomain<<endl;



  // ARA6 Aeffs
  TGraph *gAeff_ARA6[2];
  double ARA6_Aeff_Fdomain_x[4] = { 17, 18, 19, 20 };
  double ARA6_Aeff_Fdomain_y[4] = { 0.0033724, 0.0372783, 0.180956, 0.566852 };
  gAeff_ARA6[0] = new TGraph (4, ARA6_Aeff_Fdomain_x, ARA6_Aeff_Fdomain_y);

  double ARA6_Aeff_Tdomain_x[6] = { 17.5, 18, 18.5, 19, 19.5, 20 };
  double ARA6_Aeff_Tdomain_y[6] = { 0.00280024, 0.0108438, 0.033258, 0.0704812, 0.128889, 0.20954 };
  gAeff_ARA6[1] = new TGraph (6, ARA6_Aeff_Tdomain_x, ARA6_Aeff_Tdomain_y);


  double avg_factor_10_6_Tdomain = GetAvgFactor( ARA10_Aeff_Tdomain_x, ARA10_Aeff_Tdomain_y, 6, ARA6_Aeff_Tdomain_x, ARA6_Aeff_Tdomain_y, 6  );
  double rms_factor_10_6_Tdomain = GetRMSFactor( ARA10_Aeff_Tdomain_x, ARA10_Aeff_Tdomain_y, 6, ARA6_Aeff_Tdomain_x, ARA6_Aeff_Tdomain_y, 6, avg_factor_10_6_Tdomain );

  double avg_factor_10_6_Fdomain = GetAvgFactor( ARA10_Aeff_Fdomain_x, ARA10_Aeff_Fdomain_y, 4, ARA6_Aeff_Fdomain_x, ARA6_Aeff_Fdomain_y, 4 );
  double rms_factor_10_6_Fdomain = GetRMSFactor( ARA10_Aeff_Fdomain_x, ARA10_Aeff_Fdomain_y, 4, ARA6_Aeff_Fdomain_x, ARA6_Aeff_Fdomain_y, 4, avg_factor_10_6_Fdomain );

  cout<<"\nAvgFactor from ARA10 -> 6, Tdomain : "<<avg_factor_10_6_Tdomain<<", Fdomain : "<<avg_factor_10_6_Fdomain<<endl;
  cout<<"RMSFactor from ARA10 -> 6, Tdomain : "<<rms_factor_10_6_Tdomain<<", Fdomain : "<<rms_factor_10_6_Fdomain<<endl;




  // ARA23 Aeffs
  TGraph *gAeff_ARA23[2];
  double ARA23_Aeff_Fdomain_x[4] = { 17, 18, 19, 20 };
  double ARA23_Aeff_Fdomain_y[4] = { 0.00116123, 0.0135043, 0.0730713, 0.267898 };
  gAeff_ARA23[0] = new TGraph (4, ARA23_Aeff_Fdomain_x, ARA23_Aeff_Fdomain_y);

  double ARA23_Aeff_Tdomain_x[6] = { 17.5, 18, 18.5, 19, 19.5, 20 };
  double ARA23_Aeff_Tdomain_y[6] = { 0.000987022, 0.00421745, 0.0119462, 0.0297662, 0.0583245, 0.100393 };
  gAeff_ARA23[1] = new TGraph (6, ARA23_Aeff_Tdomain_x, ARA23_Aeff_Tdomain_y);


  double avg_factor_6_2_Tdomain = GetAvgFactor( ARA6_Aeff_Tdomain_x, ARA6_Aeff_Tdomain_y, 6, ARA23_Aeff_Tdomain_x, ARA23_Aeff_Tdomain_y, 6 );
  double rms_factor_6_2_Tdomain = GetRMSFactor( ARA6_Aeff_Tdomain_x, ARA6_Aeff_Tdomain_y, 6, ARA23_Aeff_Tdomain_x, ARA23_Aeff_Tdomain_y, 6, avg_factor_6_2_Tdomain );

  double avg_factor_6_2_Fdomain = GetAvgFactor( ARA6_Aeff_Fdomain_x, ARA6_Aeff_Fdomain_y, 4, ARA23_Aeff_Fdomain_x, ARA23_Aeff_Fdomain_y, 4 );
  double rms_factor_6_2_Fdomain = GetRMSFactor( ARA6_Aeff_Fdomain_x, ARA6_Aeff_Fdomain_y, 4, ARA23_Aeff_Fdomain_x, ARA23_Aeff_Fdomain_y, 4, avg_factor_6_2_Fdomain );

  cout<<"\nAvgFactor from ARA6 -> 2&3, Tdomain : "<<avg_factor_6_2_Tdomain<<", Fdomain : "<<avg_factor_6_2_Fdomain<<endl;
  cout<<"RMSFactor from ARA6 -> 2&3, Tdomain : "<<rms_factor_6_2_Tdomain<<", Fdomain : "<<rms_factor_6_2_Fdomain<<endl;





  // ARA1 200m Aeffs
  TGraph *gAeff_ARA1[2];
  double ARA1_200m_Aeff_Fdomain_x[4] = { 17, 18, 19, 20 };
  double ARA1_200m_Aeff_Fdomain_y[4] = { 0.000574066, 0.00694497, 0.0419789, 0.1501 };
  gAeff_ARA1[0] = new TGraph (4, ARA1_200m_Aeff_Fdomain_x, ARA1_200m_Aeff_Fdomain_y);

  double ARA1_200m_Aeff_Tdomain_x[4] = { 17, 18, 19, 20 };
  double ARA1_200m_Aeff_Tdomain_y[4] = { 9.38336e-05, 0.00213447, 0.0170407, 0.0618513 };
  gAeff_ARA1[1] = new TGraph (4, ARA1_200m_Aeff_Tdomain_x, ARA1_200m_Aeff_Tdomain_y);


  double avg_factor_2_1_Tdomain = GetAvgFactor( ARA23_Aeff_Tdomain_x, ARA23_Aeff_Tdomain_y, 4, ARA1_200m_Aeff_Tdomain_x, ARA1_200m_Aeff_Tdomain_y, 4 );
  double rms_factor_2_1_Tdomain = GetRMSFactor( ARA23_Aeff_Tdomain_x, ARA23_Aeff_Tdomain_y, 4, ARA1_200m_Aeff_Tdomain_x, ARA1_200m_Aeff_Tdomain_y, 4, avg_factor_2_1_Tdomain );

  double avg_factor_2_1_Fdomain = GetAvgFactor( ARA23_Aeff_Fdomain_x, ARA23_Aeff_Fdomain_y, 4, ARA1_200m_Aeff_Fdomain_x, ARA1_200m_Aeff_Fdomain_y, 4 );
  double rms_factor_2_1_Fdomain = GetRMSFactor( ARA23_Aeff_Fdomain_x, ARA23_Aeff_Fdomain_y, 4, ARA1_200m_Aeff_Fdomain_x, ARA1_200m_Aeff_Fdomain_y, 4, avg_factor_2_1_Fdomain );

  cout<<"\nAvgFactor from ARA2&3 -> 1 (200m), Tdomain : "<<avg_factor_2_1_Tdomain<<", Fdomain : "<<avg_factor_2_1_Fdomain<<endl;
  cout<<"RMSFactor from ARA2&3 -> 1 (200m), Tdomain : "<<rms_factor_2_1_Tdomain<<", Fdomain : "<<rms_factor_2_1_Fdomain<<endl;






  // ARA1 100m Aeffs
  TGraph *gAeff_ARA1_100m[2];
  double ARA1_100m_Aeff_Fdomain_x[4] = { 17, 18, 19, 20 };
  double ARA1_100m_Aeff_Fdomain_y[4] = { 0.000433444, 0.00398613, 0.0174801, 0.05533 };
  gAeff_ARA1_100m[0] = new TGraph (4, ARA1_100m_Aeff_Fdomain_x, ARA1_100m_Aeff_Fdomain_y);

  double ARA1_100m_Aeff_Tdomain_x[7] = { 17, 17.5, 18, 18.5, 19, 19.5, 20 };
  double ARA1_100m_Aeff_Tdomain_y[7] = { 8.08302e-05, 0.000380774, 0.00131882, 0.00348674, 0.00729003, 0.0128874, 0.0210407 };
  gAeff_ARA1_100m[1] = new TGraph (7, ARA1_100m_Aeff_Tdomain_x, ARA1_100m_Aeff_Tdomain_y);


  double avg_factor_200m_100m_Tdomain = GetAvgFactor( ARA1_200m_Aeff_Tdomain_x, ARA1_200m_Aeff_Tdomain_y, 4, ARA1_100m_Aeff_Tdomain_x, ARA1_100m_Aeff_Tdomain_y, 7 );
  double rms_factor_200m_100m_Tdomain = GetRMSFactor( ARA1_200m_Aeff_Tdomain_x, ARA1_200m_Aeff_Tdomain_y, 4, ARA1_100m_Aeff_Tdomain_x, ARA1_100m_Aeff_Tdomain_y, 7, avg_factor_200m_100m_Tdomain );

  double avg_factor_200m_100m_Fdomain = GetAvgFactor( ARA1_200m_Aeff_Fdomain_x, ARA1_200m_Aeff_Fdomain_y, 4, ARA1_100m_Aeff_Fdomain_x, ARA1_100m_Aeff_Fdomain_y, 4 );
  double rms_factor_200m_100m_Fdomain = GetRMSFactor( ARA1_200m_Aeff_Fdomain_x, ARA1_200m_Aeff_Fdomain_y, 4, ARA1_100m_Aeff_Fdomain_x, ARA1_100m_Aeff_Fdomain_y, 4, avg_factor_200m_100m_Fdomain );

  cout<<"\nAvgFactor from ARA1 200m -> 100m, Tdomain : "<<avg_factor_200m_100m_Tdomain<<", Fdomain : "<<avg_factor_200m_100m_Fdomain<<endl;
  cout<<"RMSFactor from ARA1 200m -> 100m, Tdomain : "<<rms_factor_200m_100m_Tdomain<<", Fdomain : "<<rms_factor_200m_100m_Fdomain<<endl;





  // ARA1 100m Aeffs with trig by only bottom 8 chs
  TGraph *gAeff_ARA1_100m_Low8Chs_Tdomain;
  double ARA1_100m_Low8Chs_Aeff_Tdomain_x[3] = { 18, 19, 20 };
  double ARA1_100m_Low8Chs_Aeff_Tdomain_y[3] = { 0.0013365, 0.00908398, 0.0346322 };
  gAeff_ARA1_100m_Low8Chs_Tdomain = new TGraph (3, ARA1_100m_Low8Chs_Aeff_Tdomain_x, ARA1_100m_Low8Chs_Aeff_Tdomain_y);




  // ARA1 100m Aeffs with trig by only bottom 8 chs, TB noise
  TGraph *gAeff_ARA1_100m_Low8Chs_TBnoise_Tdomain;
  double ARA1_100m_Low8Chs_TBnoise_Aeff_Tdomain_x[3] = { 18, 19, 20 };
  double ARA1_100m_Low8Chs_TBnoise_Aeff_Tdomain_y[3] = { 0.00132485, 0.00973719, 0.0350489 };
  gAeff_ARA1_100m_Low8Chs_TBnoise_Tdomain = new TGraph (3, ARA1_100m_Low8Chs_TBnoise_Aeff_Tdomain_x, ARA1_100m_Low8Chs_TBnoise_Aeff_Tdomain_y);




  // ARA1 100m Aeffs with trig by only bottom 8 chs, TB noise, TB thres
  TGraph *gAeff_ARA1_100m_Low8Chs_TBnoise_TBthres_Tdomain;
  double ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_x[3] = { 18, 19, 20 };
  double ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_y[3] = { 0.00195696, 0.0107611, 0.0389152 };
  gAeff_ARA1_100m_Low8Chs_TBnoise_TBthres_Tdomain = new TGraph (3, ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_x, ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_y);



  double avg_factor_TBtune_100m_Tdomain = GetAvgFactor( ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_x, ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_y, 3, ARA1_100m_Aeff_Tdomain_x, ARA1_100m_Aeff_Tdomain_y, 7 );
  double rms_factor_TBtune_100m_Tdomain = GetRMSFactor( ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_x, ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_y, 3, ARA1_100m_Aeff_Tdomain_x, ARA1_100m_Aeff_Tdomain_y, 7, avg_factor_TBtune_100m_Tdomain );

  cout<<"\nAvgFactor from ARA1 100m(TBtune) -> 100m, Tdomain : "<<avg_factor_TBtune_100m_Tdomain<<endl;
  cout<<"RMSFactor from ARA1 100m(TBtune) -> 100m, Tdomain : "<<rms_factor_TBtune_100m_Tdomain<<endl;





  // TestBed Aeffs
  TGraph *gAeff_TestBed[2];
  double TestBed_Aeff_Fdomain_x[11] = { 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21 };
  double TestBed_Aeff_Fdomain_y[11] = { 1.86736e-05, 6.38504e-05, 0.000221994, 0.000704702, 0.00177821, 0.00389573, 0.00793405, 0.0134924, 0.0263898, 0.0414937, 0.0672489 };
  gAeff_TestBed[0] = new TGraph (11, TestBed_Aeff_Fdomain_x, TestBed_Aeff_Fdomain_y);

  double TestBed_Aeff_Tdomain_x[9] = { 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21 };
  double TestBed_Aeff_Tdomain_y[9] = { 6.21577e-05, 0.000207001, 0.000653315, 0.00168405, 0.00407532, 0.00796773, 0.0140794, 0.0226091, 0.0331422 };
  gAeff_TestBed[1] = new TGraph (9, TestBed_Aeff_Tdomain_x, TestBed_Aeff_Tdomain_y);





  // TestBed Aeff with 100m deeper
  TGraph *gAeff_TestBed_100m;
  double TestBed_100m_Aeff_Tdomain_x[3] = { 18, 19, 20 };
  double TestBed_100m_Aeff_Tdomain_y[3] = { 0.001717, 0.0121, 0.0438 };
  gAeff_TestBed_100m = new TGraph (3, TestBed_100m_Aeff_Tdomain_x, TestBed_100m_Aeff_Tdomain_y);


  double avg_factor_TB_100m_org_Tdomain = GetAvgFactor( TestBed_100m_Aeff_Tdomain_x, TestBed_100m_Aeff_Tdomain_y, 3, TestBed_Aeff_Tdomain_x, TestBed_Aeff_Tdomain_y, 9 );
  double rms_factor_TB_100m_org_Tdomain = GetRMSFactor( TestBed_100m_Aeff_Tdomain_x, TestBed_100m_Aeff_Tdomain_y, 3, TestBed_Aeff_Tdomain_x, TestBed_Aeff_Tdomain_y, 9, avg_factor_TB_100m_org_Tdomain );

  cout<<"\nAvgFactor from TestBed 100m -> org, Tdomain : "<<avg_factor_TB_100m_org_Tdomain<<endl;
  cout<<"RMSFactor from TestBed 100m -> org, Tdomain : "<<rms_factor_TB_100m_org_Tdomain<<endl;




  double avg_factor_TBtune_TB100m_Tdomain = GetAvgFactor( ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_x, ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_y, 3, TestBed_100m_Aeff_Tdomain_x, TestBed_100m_Aeff_Tdomain_y, 3 );
  double rms_factor_TBtune_TB100m_Tdomain = GetRMSFactor( ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_x, ARA1_100m_Low8Chs_TBnoise_TBthres_Aeff_Tdomain_y, 3, TestBed_100m_Aeff_Tdomain_x, TestBed_100m_Aeff_Tdomain_y, 3, avg_factor_TBtune_TB100m_Tdomain );

  cout<<"\nAvgFactor from ARA1 100m(TBtune) -> TB100m, Tdomain : "<<avg_factor_TBtune_TB100m_Tdomain<<endl;
  cout<<"RMSFactor from ARA1 100m(TBtune) -> TB100m, Tdomain : "<<rms_factor_TBtune_TB100m_Tdomain<<endl;






  // Aeff plots!
  //TCanvas *cAeff_compare = new TCanvas("cAeff_compare","",200,10,800,600);
  //TCanvas *cAeff_compare = new TCanvas("cAeff_compare","",200,10,1000,1400);
  TCanvas *cAeff_compare = new TCanvas("cAeff_compare","",200,10,1000,2100);
  cAeff_compare->Divide(1,3);

  // cd 1 for Fdomain
  cAeff_compare->cd(1);
  cAeff_compare->cd(1)->SetLogy();
  cAeff_compare->cd(1)->SetGrid();

  //gAeff_TestBed_Dave->SetTitle("Effective area for ARA (km^{2} sr)");
  //gAeff_TestBed_Dave->SetTitle("Effective area for ARA (km^{2} sr) Fdomain");
  gAeff_TestBed_Dave->SetTitle("Effective area for ARA (km^{2} sr) pre phase");
  gAeff_TestBed_Dave->GetHistogram()->SetXTitle("log E");
  gAeff_TestBed_Dave->GetHistogram()->SetYTitle("km^2 sr");

  gAeff_TestBed_Dave->GetHistogram()->SetMaximum(5);
  gAeff_TestBed_Dave->GetHistogram()->SetMinimum(1.E-5);
  //
  gAeff_TestBed_Dave->SetLineStyle(8);
  gAeff_TestBed_Dave->SetLineWidth(2);
  gAeff_TestBed_Dave->Draw("al");


  int domain = 0; // 0 for Fdomain, 1 for Tdomain

  gAeff_ARA37[domain]->SetLineColor(1);
  gAeff_ARA37[domain]->SetLineWidth(2);
  gAeff_ARA37[domain]->Draw("l");

  gAeff_ARA10[domain]->SetLineColor(kGreen);
  gAeff_ARA10[domain]->SetLineWidth(2);
  gAeff_ARA10[domain]->Draw("l");

  gAeff_ARA6[domain]->SetLineColor(kBlue);
  gAeff_ARA6[domain]->SetLineWidth(2);
  gAeff_ARA6[domain]->Draw("l");

  gAeff_ARA23[domain]->SetLineColor(kRed);
  gAeff_ARA23[domain]->SetLineWidth(2);
  gAeff_ARA23[domain]->Draw("l");

  gAeff_ARA1[domain]->SetLineColor(kGreen);
  gAeff_ARA1[domain]->SetLineStyle(2);
  gAeff_ARA1[domain]->SetLineWidth(2);
  gAeff_ARA1[domain]->Draw("l");

  gAeff_ARA1_100m[domain]->SetLineColor(kBlue);
  gAeff_ARA1_100m[domain]->SetLineStyle(2);
  gAeff_ARA1_100m[domain]->SetLineWidth(2);
  gAeff_ARA1_100m[domain]->Draw("l");

  gAeff_TestBed[domain]->SetLineColor(kRed);
  gAeff_TestBed[domain]->SetLineStyle(2);
  gAeff_TestBed[domain]->SetLineWidth(2);
  gAeff_TestBed[domain]->Draw("l");


  TLegend *Leg_Aeff_1 = new TLegend(0.95, 0.15, 0.5,0.45);
  Leg_Aeff_1-> AddEntry(gAeff_TestBed_Dave, "TestBed Dave (calc)", "l");
  Leg_Aeff_1-> AddEntry(gAeff_ARA37[domain], "ARA37 AraSim", "l");
  Leg_Aeff_1-> AddEntry(gAeff_ARA10[domain], "ARA10 AraSim", "l");
  Leg_Aeff_1-> AddEntry(gAeff_ARA6[domain], "ARA6 AraSim", "l");
  Leg_Aeff_1-> AddEntry(gAeff_ARA23[domain], "ARA2&3 (200m) AraSim", "l");
  Leg_Aeff_1-> AddEntry(gAeff_ARA1[domain], "ARA1 (200m) AraSim", "l");
  Leg_Aeff_1-> AddEntry(gAeff_ARA1_100m[domain], "ARA1 (100m) AraSim", "l");
  Leg_Aeff_1-> AddEntry(gAeff_TestBed[domain], "TestBed AraSim", "l");
  //Leg_Aeff_1-> AddEntry(gAeff_TB_AS_100m, "TestBed AraSim phase 100m deep", "l");
  Leg_Aeff_1 -> Draw();





  // cd 2 for Tdomain
  cAeff_compare->cd(2);
  cAeff_compare->cd(2)->SetLogy();
  cAeff_compare->cd(2)->SetGrid();

  //gAeff_TestBed_Dave_2->SetTitle("Effective area for ARA (km^{2} sr)");
  //gAeff_TestBed_Dave_2->SetTitle("Effective area for ARA (km^{2} sr Tdomain)");
  gAeff_TestBed_Dave_2->SetTitle("Effective area for ARA (km^{2} sr) after phase)");
  gAeff_TestBed_Dave_2->GetHistogram()->SetXTitle("log E");
  gAeff_TestBed_Dave_2->GetHistogram()->SetYTitle("km^2 sr");

  gAeff_TestBed_Dave_2->GetHistogram()->SetMaximum(5);
  gAeff_TestBed_Dave_2->GetHistogram()->SetMinimum(1.E-5);
  //
  gAeff_TestBed_Dave_2->SetLineStyle(8);
  gAeff_TestBed_Dave_2->SetLineWidth(2);
  gAeff_TestBed_Dave_2->Draw("al");


  domain = 1; // 0 for Fdomain, 1 for Tdomain

  gAeff_ARA37[domain]->SetLineColor(1);
  gAeff_ARA37[domain]->SetLineWidth(2);
  gAeff_ARA37[domain]->Draw("l");

  gAeff_ARA10[domain]->SetLineColor(kGreen);
  gAeff_ARA10[domain]->SetLineWidth(2);
  gAeff_ARA10[domain]->Draw("l");

  gAeff_ARA6[domain]->SetLineColor(kBlue);
  gAeff_ARA6[domain]->SetLineWidth(2);
  gAeff_ARA6[domain]->Draw("l");

  gAeff_ARA23[domain]->SetLineColor(kRed);
  gAeff_ARA23[domain]->SetLineWidth(2);
  gAeff_ARA23[domain]->Draw("l");

  gAeff_ARA1[domain]->SetLineColor(kGreen);
  gAeff_ARA1[domain]->SetLineStyle(2);
  gAeff_ARA1[domain]->SetLineWidth(2);
  gAeff_ARA1[domain]->Draw("l");

  gAeff_ARA1_100m[domain]->SetLineColor(kBlue);
  gAeff_ARA1_100m[domain]->SetLineStyle(2);
  gAeff_ARA1_100m[domain]->SetLineWidth(2);
  gAeff_ARA1_100m[domain]->Draw("l");

  gAeff_TestBed[domain]->SetLineColor(kRed);
  gAeff_TestBed[domain]->SetLineStyle(2);
  gAeff_TestBed[domain]->SetLineWidth(2);
  gAeff_TestBed[domain]->Draw("l");


  // TestBed 100m only available for Tdomain for now
  gAeff_TestBed_100m->SetLineColor(kBlue);
  gAeff_TestBed_100m->SetLineStyle(5);
  gAeff_TestBed_100m->SetLineWidth(2);
  gAeff_TestBed_100m->Draw("l");


  TLegend *Leg_Aeff_2 = new TLegend(0.95, 0.15, 0.5,0.45);
  Leg_Aeff_2-> AddEntry(gAeff_TestBed_Dave_2, "TestBed Dave (calc)", "l");
  Leg_Aeff_2-> AddEntry(gAeff_ARA37[domain], "ARA37 AraSim", "l");
  Leg_Aeff_2-> AddEntry(gAeff_ARA10[domain], "ARA10 AraSim", "l");
  Leg_Aeff_2-> AddEntry(gAeff_ARA6[domain], "ARA6 AraSim", "l");
  Leg_Aeff_2-> AddEntry(gAeff_ARA23[domain], "ARA2&3 (200m) AraSim", "l");
  Leg_Aeff_2-> AddEntry(gAeff_ARA1[domain], "ARA1 (200m) AraSim", "l");
  Leg_Aeff_2-> AddEntry(gAeff_ARA1_100m[domain], "ARA1 (100m) AraSim", "l");
  Leg_Aeff_2-> AddEntry(gAeff_TestBed[domain], "TestBed AraSim", "l");
  Leg_Aeff_2-> AddEntry(gAeff_TestBed_100m, "TestBed AraSim 100m deep", "l");
  Leg_Aeff_2 -> Draw();






  // cd 3 for Tdomain, ARA1 100m modes
  cAeff_compare->cd(3);
  cAeff_compare->cd(3)->SetLogy();
  cAeff_compare->cd(3)->SetGrid();

  //gAeff_TestBed_Dave_3->SetTitle("Effective area for ARA (km^{2} sr)");
  //gAeff_TestBed_Dave_3->SetTitle("Effective area for ARA (km^{2} sr Tdomain)");
  gAeff_TestBed_Dave_3->SetTitle("Effective area for ARA (km^{2} sr) after phase");
  gAeff_TestBed_Dave_3->GetHistogram()->SetXTitle("log E");
  gAeff_TestBed_Dave_3->GetHistogram()->SetYTitle("km^2 sr");

  //gAeff_TestBed_Dave_3->GetHistogram()->SetMaximum(5);
  //gAeff_TestBed_Dave_3->GetHistogram()->SetMinimum(1.E-5);
  gAeff_TestBed_Dave_3->GetHistogram()->SetMaximum(1.E-1);
  gAeff_TestBed_Dave_3->GetHistogram()->SetMinimum(1.E-4);

  //
  gAeff_TestBed_Dave_3->SetLineStyle(8);
  gAeff_TestBed_Dave_3->SetLineWidth(2);
  gAeff_TestBed_Dave_3->Draw("al");


  domain = 1; // 0 for Fdomain, 1 for Tdomain


  gAeff_ARA1_100m[domain]->SetLineColor(kBlue);
  gAeff_ARA1_100m[domain]->SetLineStyle(2);
  gAeff_ARA1_100m[domain]->SetLineWidth(2);
  gAeff_ARA1_100m[domain]->Draw("l");


  gAeff_ARA1_100m_Low8Chs_Tdomain->SetLineColor(kBlack);
  //gAeff_ARA1_100m_Low8Chs_Tdomain->SetLineStyle(2);
  gAeff_ARA1_100m_Low8Chs_Tdomain->SetLineWidth(2);
  gAeff_ARA1_100m_Low8Chs_Tdomain->Draw("l");


  gAeff_ARA1_100m_Low8Chs_TBnoise_Tdomain->SetLineColor(kRed);
  //gAeff_ARA1_100m_Low8Chs_TBnoise_Tdomain->SetLineStyle(2);
  gAeff_ARA1_100m_Low8Chs_TBnoise_Tdomain->SetLineWidth(2);
  gAeff_ARA1_100m_Low8Chs_TBnoise_Tdomain->Draw("l");


  gAeff_ARA1_100m_Low8Chs_TBnoise_TBthres_Tdomain->SetLineColor(kGreen+2);
  //gAeff_ARA1_100m_Low8Chs_TBnoise_TBthres_Tdomain->SetLineStyle(2);
  gAeff_ARA1_100m_Low8Chs_TBnoise_TBthres_Tdomain->SetLineWidth(2);
  gAeff_ARA1_100m_Low8Chs_TBnoise_TBthres_Tdomain->Draw("l");



  // TestBed 100m only available for Tdomain for now
  gAeff_TestBed_100m->SetLineColor(kBlue);
  gAeff_TestBed_100m->SetLineStyle(5);
  gAeff_TestBed_100m->SetLineWidth(2);
  gAeff_TestBed_100m->Draw("l");


  TLegend *Leg_Aeff_3 = new TLegend(0.95, 0.15, 0.5,0.45);
  Leg_Aeff_3-> AddEntry(gAeff_TestBed_Dave_3, "TestBed Dave (calc)", "l");
  Leg_Aeff_3-> AddEntry(gAeff_ARA1_100m[domain], "ARA1 (100m) AraSim", "l");
  Leg_Aeff_3-> AddEntry(gAeff_ARA1_100m_Low8Chs_Tdomain, "ARA1 (100m) trig low 8chs", "l");
  Leg_Aeff_3-> AddEntry(gAeff_ARA1_100m_Low8Chs_TBnoise_Tdomain, "ARA1 (100m) trig low 8chs, TBnoise", "l");
  Leg_Aeff_3-> AddEntry(gAeff_ARA1_100m_Low8Chs_TBnoise_TBthres_Tdomain, "ARA1 (100m) trig low 8chs, TBnoise,thres", "l");
  Leg_Aeff_3-> AddEntry(gAeff_TestBed_100m, "TestBed AraSim 100m deep", "l");
  Leg_Aeff_3 -> Draw();


  cAeff_compare->Print("test_Aeff_compare.pdf");

  delete cAeff_compare;












  // just couple of Veff plots
  // ARA37 Veffs
  TGraph *gVeff_ARA37[2];
  double ARA37_Veff_Fdomain_x[9] = { 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20 };
  double ARA37_Veff_Fdomain_y[9] = { 5.22266, 20.2339, 60.3621, 142.238, 257.911, 422.559, 584.753, 686.503, 908.165 };
  gVeff_ARA37[0] = new TGraph (9, ARA37_Veff_Fdomain_x, ARA37_Veff_Fdomain_y);

  double ARA37_Veff_Tdomain_x[6] = { 17.5, 18, 18.5, 19, 19.5, 20 };
  double ARA37_Veff_Tdomain_y[6] = { 32.8477, 87.6905, 167.13, 246.852, 317.536, 375.004 };
  gVeff_ARA37[1] = new TGraph (6, ARA37_Veff_Tdomain_x, ARA37_Veff_Tdomain_y);



  // ARA1 200m Veffs
  TGraph *gVeff_ARA1[2];
  double ARA1_200m_Veff_Fdomain_x[4] = { 17, 18, 19, 20 };
  double ARA1_200m_Veff_Fdomain_y[4] = { 1.63631, 8.83851, 25.6402, 46.4379 };
  gVeff_ARA1[0] = new TGraph (4, ARA1_200m_Veff_Fdomain_x, ARA1_200m_Veff_Fdomain_y);

  double ARA1_200m_Veff_Tdomain_x[4] = { 17, 18, 19, 20 };
  double ARA1_200m_Veff_Tdomain_y[4] = { 0.267462, 2.71644, 10.4082, 19.1355 };
  gVeff_ARA1[1] = new TGraph (4, ARA1_200m_Veff_Tdomain_x, ARA1_200m_Veff_Tdomain_y);


  // TestBed Veffs
  TGraph *gVeff_TestBed[2];
  double TestBed_Veff_Fdomain_x[11] = { 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21 };
  double TestBed_Veff_Fdomain_y[11] = { 0.132019, 0.28228, 0.632769, 1.32811, 2.26304, 3.40838, 4.84602, 5.83077, 8.16448, 9.2868, 10.9906 };
  gVeff_TestBed[0] = new TGraph (11, TestBed_Veff_Fdomain_x, TestBed_Veff_Fdomain_y);

  double TestBed_Veff_Tdomain_x[9] = { 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21 };
  double TestBed_Veff_Tdomain_y[9] = { 0.177174, 0.390123, 0.831442, 1.47338, 2.48916, 3.44327, 4.35588, 5.0602, 5.4165 };
  gVeff_TestBed[1] = new TGraph (9, TestBed_Veff_Tdomain_x, TestBed_Veff_Tdomain_y);



  /*
  // Veff plots!
  TCanvas *cVeff_compare = new TCanvas("cVeff_compare","",200,10,1000,1400);
  cVeff_compare->Divide(1,2);

  // cd 1 for Fdomain
  cVeff_compare->cd(1);
  cVeff_compare->cd(1)->SetLogy();
  cVeff_compare->cd(1)->SetGrid();

  domain = 0; // 0 for Fdomain, 1 for Tdomain

  gVeff_ARA37[domain]->SetTitle("Effective volume for ARA (km^{3} sr) pre phase");
  gVeff_ARA37[domain]->GetHistogram()->SetXTitle("log E");
  gVeff_ARA37[domain]->GetHistogram()->SetYTitle("km^3 sr");

  gVeff_ARA37[domain]->GetHistogram()->SetMaximum(1.e3);
  gVeff_ARA37[domain]->GetHistogram()->SetMinimum(1.e-2);
  //
  //gVeff_ARA37[domain]->SetLineStyle(8);
  gVeff_ARA37[domain]->SetLineWidth(2);
  gVeff_ARA37[domain]->Draw("al");

  gVeff_ARA1[domain]->SetLineColor(kBlue);
  //gVeff_ARA1[domain]->SetLineStyle(2);
  gVeff_ARA1[domain]->SetLineWidth(2);
  gVeff_ARA1[domain]->Draw("l");

  gVeff_TestBed[domain]->SetLineColor(kRed);
  //gVeff_TestBed[domain]->SetLineStyle(2);
  gVeff_TestBed[domain]->SetLineWidth(2);
  gVeff_TestBed[domain]->Draw("l");


  TLegend *Leg_Veff_1 = new TLegend(0.95, 0.15, 0.5,0.45);
  Leg_Veff_1-> AddEntry(gVeff_ARA37[domain], "ARA37 AraSim", "l");
  Leg_Veff_1-> AddEntry(gVeff_ARA1[domain], "ARA1 (200m) AraSim", "l");
  Leg_Veff_1-> AddEntry(gVeff_TestBed[domain], "TestBed AraSim", "l");
  //Leg_Veff_1-> AddEntry(gVeff_TB_AS_100m, "TestBed AraSim phase 100m deep", "l");
  Leg_Veff_1 -> Draw();





  // cd 2 for Tdomain
  cVeff_compare->cd(2);
  cVeff_compare->cd(2)->SetLogy();
  cVeff_compare->cd(2)->SetGrid();


  domain = 1; // 0 for Fdomain, 1 for Tdomain

  gVeff_ARA37[domain]->SetTitle("Effective volume for ARA (km^{3} sr) after phase");
  gVeff_ARA37[domain]->GetHistogram()->SetXTitle("log E");
  gVeff_ARA37[domain]->GetHistogram()->SetYTitle("km^3 sr");

  gVeff_ARA37[domain]->GetHistogram()->SetMaximum(1.e3);
  gVeff_ARA37[domain]->GetHistogram()->SetMinimum(1.e-2);
  //
  //gVeff_ARA37[domain]->SetLineStyle(8);
  gVeff_ARA37[domain]->SetLineWidth(2);
  gVeff_ARA37[domain]->Draw("al");

  gVeff_ARA1[domain]->SetLineColor(kBlue);
  //gVeff_ARA1[domain]->SetLineStyle(2);
  gVeff_ARA1[domain]->SetLineWidth(2);
  gVeff_ARA1[domain]->Draw("l");

  gVeff_TestBed[domain]->SetLineColor(kRed);
  //gVeff_TestBed[domain]->SetLineStyle(2);
  gVeff_TestBed[domain]->SetLineWidth(2);
  gVeff_TestBed[domain]->Draw("l");


  TLegend *Leg_Veff_2 = new TLegend(0.95, 0.15, 0.5,0.45);
  Leg_Veff_2-> AddEntry(gVeff_ARA37[domain], "ARA37 AraSim", "l");
  Leg_Veff_2-> AddEntry(gVeff_ARA1[domain], "ARA1 (200m) AraSim", "l");
  Leg_Veff_2-> AddEntry(gVeff_TestBed[domain], "TestBed AraSim", "l");
  //Leg_Veff_1-> AddEntry(gVeff_TB_AS_100m, "TestBed AraSim phase 100m deep", "l");
  Leg_Veff_2 -> Draw();


  cVeff_compare->Print("test_Veff_compare.pdf");

  delete cVeff_compare;
  */



















  // plots for detected neutrino events
  // ESS baseline

  double nu_eventx[5];
  nu_eventx[0] = 16;
  nu_eventx[1] = 17;
  nu_eventx[2] = 18;
  nu_eventx[3] = 19;
  nu_eventx[4] = 20;


  double ess_base[5];
  ess_base[0] = 0.03588;
  ess_base[1] = 1.1983;
  ess_base[2] = 12.2774;
  ess_base[3] = 8.1696;
  ess_base[4] = 1.0442;

  TCanvas *c_nu_event = new TCanvas("c_nu_event","A Simple Graph Example",200,10,1000,700);
  TGraph *g_ess_base = new TGraph( 5, nu_eventx, ess_base );

  c_nu_event->cd();

  c_nu_event->SetLogy();
  g_ess_base->SetTitle("ESS baseline model events with ARA37 3yrs");
  g_ess_base->GetHistogram()->SetXTitle("log E");
  g_ess_base->GetHistogram()->SetYTitle("event no");
  g_ess_base->GetHistogram()->SetMaximum(50);
  //g_vhvh->GetHistogram()->SetMinimum(0.1);
  g_ess_base->Draw("al*");

  c_nu_event->Print("test_ESS_baseline1.pdf");


  // plots for detected neutrino events
  // WB upper, lower

  double WB_up[5], WB_low[5];
  WB_up[0] = 11.58;
  WB_up[1] = 17.52;
  WB_up[2] = 15.89;
  WB_up[3] = 7.41;
  WB_up[4] = 2.65;

  WB_low[0] = 2.37;
  WB_low[1] = 3.58;
  WB_low[2] = 3.25;
  WB_low[3] = 1.52;
  WB_low[4] = 0.54;

  TCanvas *c_nu_event2 = new TCanvas("c_nu_event2","A Simple Graph Example",200,10,1000,700);
  TGraph *g_WB_up = new TGraph( 5, nu_eventx, WB_up );
  TGraph *g_WB_low = new TGraph( 5, nu_eventx, WB_low );

  c_nu_event2->cd();

  c_nu_event2->SetLogy();
  g_WB_up->SetTitle("WB limit events with ARA37 3yrs");
  g_WB_up->GetHistogram()->SetXTitle("log E");
  g_WB_up->GetHistogram()->SetYTitle("event no");
  g_WB_up->GetHistogram()->SetMaximum(30);
  g_WB_up->GetHistogram()->SetMinimum(0.1);
  g_WB_up->Draw("al*");

  g_WB_low->SetLineColor(kRed);
  g_WB_low->Draw("l*");

  TLegend *Leg_WB = new TLegend(1., 0.75, 0.85,0.6);
  Leg_WB -> AddEntry(g_WB_up, "upper(rapid evol)", "l");
  Leg_WB -> AddEntry(g_WB_low, "lower(no evol)", "l");
  Leg_WB -> Draw();


  c_nu_event2->Print("test_WB_limit1.pdf");



  // AVE max, min models
  //
  double ave_max[5], ave_min[5];

  ave_max[0] = 0.00996;
  ave_max[1] = 0.3687;
  ave_max[2] = 4.438;
  ave_max[3] = 0.651;
  ave_max[4] = 0.016;

  ave_min[0] = 9.13E-4;
  ave_min[1] = 0.045;
  ave_min[2] = 0.013;
  ave_min[3] = 2.46E-4;
  ave_min[4] = 2.81E-8;

  TGraph *g_ave_max = new TGraph( 5, nu_eventx, ave_max );
  TGraph *g_ave_min = new TGraph( 5, nu_eventx, ave_min );


  // Kotera 2010
  //
  double kotera_fe[4];

  kotera_fe[0] = 0.0034;
  kotera_fe[1] = 0.0175;
  kotera_fe[2] = 0.0621;
  kotera_fe[3] = 0.0053;

  TGraph *g_kotera_fe = new TGraph( 4, nu_eventx, kotera_fe );


  // don't understand ess fig9?!
  /*
  // ess fig9
  //
  double ess_fig9[5];

  ess_fig9[0] = 0.185;
  ess_fig9[1] = 5.788;
  ess_fig9[2] = 38.39;
  ess_fig9[3] = 20.08;
  ess_fig9[4] = 2.03;
  */


  // event plot for all models I calculated
  TCanvas *c_nu_event_all = new TCanvas("c_nu_event_all","A Simple Graph Example",200,10,1000,700);

  c_nu_event_all->cd();

  c_nu_event_all->SetLogy();
  c_nu_event_all->SetGrid();

  g_ess_base->SetTitle("Expected events from ARA37 3yrs");
  g_ess_base->GetHistogram()->SetXTitle("log E");
  g_ess_base->GetHistogram()->SetYTitle("event no");
  g_ess_base->GetHistogram()->SetMaximum(100);
  g_ess_base->GetHistogram()->SetMinimum(0.001);
  g_ess_base->SetLineWidth(2);
  g_ess_base->Draw("al*");


  g_WB_up->SetLineColor(kRed);
  g_WB_up->SetLineWidth(2);
  g_WB_up->Draw("l*");

  g_WB_low->SetLineColor(kRed);
  g_WB_low->SetLineStyle(7);
  g_WB_low->SetLineWidth(2);
  g_WB_low->Draw("l*");


  g_ave_max->SetLineColor(kBlue);
  g_ave_max->SetLineWidth(2);
  g_ave_max->Draw("l*");

  g_ave_min->SetLineColor(kBlue);
  g_ave_min->SetLineStyle(7);
  g_ave_min->SetLineWidth(2);
  g_ave_min->Draw("l*");

  g_kotera_fe->SetLineColor(kGreen);
  g_kotera_fe->SetLineWidth(2);
  g_kotera_fe->Draw("l*");

  TLegend *Leg_models = new TLegend(1., 1., 0.8,0.8);
  Leg_models -> AddEntry(g_ess_base, "ESS baseline", "l");
  Leg_models -> AddEntry(g_WB_up, "WB upper(rapid evol)", "l");
  Leg_models -> AddEntry(g_WB_low, "WB lower(no evol)", "l");
  Leg_models -> AddEntry(g_ave_max, "AVE iron Emax10^22eV", "l");
  Leg_models -> AddEntry(g_ave_min, "AVE iron Emax10^21eV", "l");
  Leg_models -> AddEntry(g_kotera_fe, "Kotera iron only", "l");
  Leg_models -> Draw();

  c_nu_event_all->Print("test_Nu_events_all.pdf");



  // event plot for all models I calculated for ARA3
  double Veff_factor_ARA37toARA3[5];
  Veff_factor_ARA37toARA3[0] = 0.104; // from Aeff_ARA3 / Aeff_ARA37
  Veff_factor_ARA37toARA3[1] = 0.112;
  Veff_factor_ARA37toARA3[2] = 0.117;
  Veff_factor_ARA37toARA3[3] = 0.127;
  Veff_factor_ARA37toARA3[4] = 0.133;

  // calculate event N for models for ARA3
  double ess_base_ARA3[5];
  ApplyFactorArray( ess_base_ARA3, ess_base, Veff_factor_ARA37toARA3, 5 );
  TGraph *g_ess_base_ARA3 = new TGraph( 5, nu_eventx, ess_base_ARA3 );

  double WB_up_ARA3[5];
  ApplyFactorArray( WB_up_ARA3, WB_up, Veff_factor_ARA37toARA3, 5 );
  TGraph *g_WB_up_ARA3 = new TGraph( 5, nu_eventx, WB_up_ARA3 );

  double WB_low_ARA3[5];
  ApplyFactorArray( WB_low_ARA3, WB_low, Veff_factor_ARA37toARA3, 5 );
  TGraph *g_WB_low_ARA3 = new TGraph( 5, nu_eventx, WB_low_ARA3 );

  double ave_max_ARA3[5];
  ApplyFactorArray( ave_max_ARA3, ave_max, Veff_factor_ARA37toARA3, 5 );
  TGraph *g_ave_max_ARA3 = new TGraph( 5, nu_eventx, ave_max_ARA3 );

  double ave_min_ARA3[5];
  ApplyFactorArray( ave_min_ARA3, ave_min, Veff_factor_ARA37toARA3, 5 );
  TGraph *g_ave_min_ARA3 = new TGraph( 5, nu_eventx, ave_min_ARA3 );

  double kotera_fe_ARA3[4];
  ApplyFactorArray( kotera_fe_ARA3, kotera_fe, Veff_factor_ARA37toARA3, 4 );
  TGraph *g_kotera_fe_ARA3 = new TGraph( 4, nu_eventx, kotera_fe_ARA3 );

  TCanvas *c_nu_event_all_ARA3 = new TCanvas("c_nu_event_all_ARA3","A Simple Graph Example",200,10,1000,700);

  c_nu_event_all_ARA3->cd();

  c_nu_event_all_ARA3->SetLogy();
  c_nu_event_all_ARA3->SetGrid();

  g_ess_base_ARA3->SetTitle("Expected events from ARA3 3yrs");
  g_ess_base_ARA3->GetHistogram()->SetXTitle("log E");
  g_ess_base_ARA3->GetHistogram()->SetYTitle("event no");
  g_ess_base_ARA3->GetHistogram()->SetMaximum(10);
  g_ess_base_ARA3->GetHistogram()->SetMinimum(0.0001);
  g_ess_base_ARA3->SetLineWidth(2);
  g_ess_base_ARA3->Draw("al*");


  g_WB_up_ARA3->SetLineColor(kRed);
  g_WB_up_ARA3->SetLineWidth(2);
  g_WB_up_ARA3->Draw("l*");

  g_WB_low_ARA3->SetLineColor(kRed);
  g_WB_low_ARA3->SetLineStyle(7);
  g_WB_low_ARA3->SetLineWidth(2);
  g_WB_low_ARA3->Draw("l*");


  g_ave_max_ARA3->SetLineColor(kBlue);
  g_ave_max_ARA3->SetLineWidth(2);
  g_ave_max_ARA3->Draw("l*");

  g_ave_min_ARA3->SetLineColor(kBlue);
  g_ave_min_ARA3->SetLineStyle(7);
  g_ave_min_ARA3->SetLineWidth(2);
  g_ave_min_ARA3->Draw("l*");

  g_kotera_fe_ARA3->SetLineColor(kGreen);
  g_kotera_fe_ARA3->SetLineWidth(2);
  g_kotera_fe_ARA3->Draw("l*");

  TLegend *Leg_models_ARA3 = new TLegend(1., 1., 0.8,0.8);
  Leg_models_ARA3 -> AddEntry(g_ess_base_ARA3, "ESS baseline", "l");
  Leg_models_ARA3 -> AddEntry(g_WB_up_ARA3, "WB upper(rapid evol)", "l");
  Leg_models_ARA3 -> AddEntry(g_WB_low_ARA3, "WB lower(no evol)", "l");
  Leg_models_ARA3 -> AddEntry(g_ave_max_ARA3, "AVE iron Emax10^22eV", "l");
  Leg_models_ARA3 -> AddEntry(g_ave_min_ARA3, "AVE iron Emax10^21eV", "l");
  Leg_models_ARA3 -> AddEntry(g_kotera_fe_ARA3, "Kotera iron only", "l");
  Leg_models_ARA3 -> Draw();

  c_nu_event_all_ARA3->Print("test_Nu_events_all_ARA3.pdf");






  // plot flux models and ARA constraints
  //
  double ESS_base_x[42] = { 14.4286  ,
			    14.5374    ,
			    14.6689    ,
			    14.8005    ,
			    14.9773    ,
			    15.1723    ,
			    15.381 ,
			    15.5714    ,
			    15.7438    ,
			    15.9252    ,
			    16.1429    ,
			    16.3288    ,
			    16.5102    ,
			    16.6463    ,
			    16.8685    ,
			    17.0136    ,
			    17.1769    ,
			    17.3311    ,
			    17.4898    ,
			    17.6531    ,
			    17.8299    ,
			    18.0295    ,
			    18.1791    ,
			    18.3787    ,
			    18.6145    ,
			    18.8231    ,
			    19.0363    ,
			    19.2449    ,
			    19.4762    ,
			    19.6032    ,
			    19.8118    ,
			    19.966 ,
			    20.0839    ,
			    20.1746    ,
			    20.3016    ,
			    20.424 ,
			    20.542 ,
			    20.6463    ,
			    20.7551    ,
			    20.8458    ,
			    20.8957    ,
			    20.9365};

  double ESS_base_y[42] = { -17.5286 ,
			    -17.4274   ,
			    -17.3209   ,
			    -17.2144   ,
			    -17.0866   ,
			    -16.9747   ,
			    -16.8895   ,
			    -16.8415   ,
			    -16.8309   ,
			    -16.8575   ,
			    -16.8842   ,
			    -16.8735   ,
			    -16.8149   ,
			    -16.7297   ,
			    -16.6019   ,
			    -16.5166   ,
			    -16.4314   ,
			    -16.3888   ,
			    -16.3569   ,
			    -16.3569   ,
			    -16.3941   ,
			    -16.4634   ,
			    -16.5646   ,
			    -16.6977   ,
			    -16.9161   ,
			    -17.1132   ,
			    -17.3529   ,
			    -17.5872   ,
			    -17.9015   ,
			    -18.0932   ,
			    -18.4021   ,
			    -18.7057   ,
			    -18.9134   ,
			    -19.1052   ,
			    -19.3928   ,
			    -19.6538   ,
			    -19.9521   ,
			    -20.2503   ,
			    -20.5965   ,
			    -20.8522   ,
			    -21.028    ,
			    -21.1824};

  LogToLine(42, ESS_base_x);
  LogToLine(42, ESS_base_y);

  TGraph *g_ESS_base = new TGraph( 42, ESS_base_x, ESS_base_y );


  // Ahlers Emin=10^18.5 eV
  double Ahlers_x[32] = { 15.0091    ,
			  15.2721    ,
			  15.5488    ,
			  15.8209    ,
			  15.9615    ,
			  16.1338    ,
			  16.2925    ,
			  16.5102    ,
			  16.7415    ,
			  16.9637    ,
			  17.1224    ,
			  17.2857    ,
			  17.4308    ,
			  17.6122    ,
			  17.8118    ,
			  17.9932    ,
			  18.1701    ,
			  18.3016    ,
			  18.4966    ,
			  18.7098    ,
			  18.9002    ,
			  19.0771    ,
			  19.2721    ,
			  19.4218    ,
			  19.5578    ,
			  19.712 ,
			  19.907 ,
			  20.0249    ,
			  20.1655    ,
			  20.2698    ,
			  20.356 ,
			  20.4512 };
  double Ahlers_y[32] = { -16.1172   ,
			  -16.1811   ,
			  -16.2503   ,
			  -16.3142   ,
			  -16.3462   ,
			  -16.3462   ,
			  -16.3462   ,
			  -16.3569   ,
			  -16.3515   ,
			  -16.3622   ,
			  -16.3515   ,
			  -16.3249   ,
			  -16.3036   ,
			  -16.3356   ,
			  -16.3782   ,
			  -16.4953   ,
			  -16.6498   ,
			  -16.8043   ,
			  -17.0226   ,
			  -17.3103   ,
			  -17.5819   ,
			  -17.8162   ,
			  -18.1252   ,
			  -18.3435   ,
			  -18.6205   ,
			  -18.9401   ,
			  -19.3395   ,
			  -19.6858   ,
			  -20.1385   ,
			  -20.5113   ,
			  -20.8415   ,
			  -21.1664 };

  LogToLine(32, Ahlers_x);
  LogToLine(32, Ahlers_y);


  TGraph *g_Ahlers = new TGraph( 32, Ahlers_x, Ahlers_y );


  double Ave07Femix_x[43] = { 14.4286    ,
			      14.5692    ,
			      14.6961    ,
			      14.9002    ,
			      15.0136    ,
			      15.1406    ,
			      15.3492    ,
			      15.4717    ,
			      15.6213    ,
			      15.78  ,
			      15.966 ,
			      16.0975    ,
			      16.2109    ,
			      16.3424    ,
			      16.4512    ,
			      16.5646    ,
			      16.6417    ,
			      16.7778    ,
			      16.8912    ,
			      17.0227    ,
			      17.1859    ,
			      17.2812    ,
			      17.4127    ,
			      17.5125    ,
			      17.6395    ,
			      17.7347    ,
			      17.8345    ,
			      17.9705    ,
			      18.1383    ,
			      18.2834    ,
			      18.3968    ,
			      18.5374    ,
			      18.678 ,
			      18.7868    ,
			      18.941 ,
			      19.0816    ,
			      19.1587    ,
			      19.3039    ,
			      19.3946    ,
			      19.5351    ,
			      19.5986    ,
			      19.6803    ,
			      19.7528 };

  double Ave07Femix_y[43] = { -17.8003   ,
			      -17.6778   ,
			      -17.5872   ,
			      -17.4434   ,
			      -17.3848   ,
			      -17.3209   ,
			      -17.257    ,
			      -17.2304   ,
			      -17.2357   ,
			      -17.2783   ,
			      -17.3688   ,
			      -17.4168   ,
			      -17.4221   ,
			      -17.4008   ,
			      -17.3529   ,
			      -17.2943   ,
			      -17.2517   ,
			      -17.1558   ,
			      -17.0866   ,
			      -17.0013   ,
			      -16.9161   ,
			      -16.8788   ,
			      -16.8415   ,
			      -16.8256   ,
			      -16.8202   ,
			      -16.8362   ,
			      -16.8629   ,
			      -16.9268   ,
			      -17.0439   ,
			      -17.1664   ,
			      -17.3103   ,
			      -17.518    ,
			      -17.7523   ,
			      -17.9441   ,
			      -18.269    ,
			      -18.5992   ,
			      -18.8336   ,
			      -19.3395   ,
			      -19.7177   ,
			      -20.2663   ,
			      -20.5379   ,
			      -20.8629   ,
			      -21.1558 };

  LogToLine(43, Ave07Femix_x);
  LogToLine(43, Ave07Femix_y);

  TGraph *g_Ave07Femix = new TGraph( 43, Ave07Femix_x, Ave07Femix_y );

  double ESS_strong_x[41] = { 14.424 ,
			      14.5828    ,
			      14.7868    ,
			      14.9819    ,
			      15.1633    ,
			      15.3311    ,
			      15.517 ,
			      15.7029    ,
			      15.9252    ,
			      16.1156    ,
			      16.3469    ,
			      16.4739    ,
			      16.6372    ,
			      16.8594    ,
			      17.0408    ,
			      17.2449    ,
			      17.4218    ,
			      17.576 ,
			      17.7256    ,
			      17.8934    ,
			      18.0703    ,
			      18.3016    ,
			      18.551 ,
			      18.746 ,
			      18.9592    ,
			      19.1905    ,
			      19.458 ,
			      19.6757    ,
			      19.8481    ,
			      19.9977    ,
			      20.093 ,
			      20.2018    ,
			      20.2971    ,
			      20.3968    ,
			      20.5238    ,
			      20.619 ,
			      20.7098    ,
			      20.8231    ,
			      20.9048    ,
			      20.9864    ,
			      21.059 };

  double ESS_strong_y[41] = { -17.0333   ,
			      -16.9001   ,
			      -16.7297   ,
			      -16.5859   ,
			      -16.4794   ,
			      -16.4154   ,
			      -16.3515   ,
			      -16.3356   ,
			      -16.3622   ,
			      -16.3835   ,
			      -16.3675   ,
			      -16.3302   ,
			      -16.245    ,
			      -16.1065   ,
			      -16    ,
			      -15.9148   ,
			      -15.8775   ,
			      -15.8562   ,
			      -15.8668   ,
			      -15.9308   ,
			      -16    ,
			      -16.1545   ,
			      -16.3622   ,
			      -16.5379   ,
			      -16.7723   ,
			      -17.0226   ,
			      -17.3848   ,
			      -17.6937   ,
			      -17.9814   ,
			      -18.253    ,
			      -18.4288   ,
			      -18.6525   ,
			      -18.8762   ,
			      -19.0945   ,
			      -19.4301   ,
			      -19.6964   ,
			      -19.9574   ,
			      -20.2983   ,
			      -20.5965   ,
			      -20.8735   ,
			      -21.1664 };

  LogToLine(41, ESS_strong_x);
  LogToLine(41, ESS_strong_y);


  TGraph *g_ESS_strong = new TGraph( 41, ESS_strong_x, ESS_strong_y );

  // test shade area betwee ESS base and ESS strong
  TGraph *g_ESS_shade = new TGraph ( 41 + 42 ); // 42 for ESS_base and 41 for ESS_strong
  for (int i=0; i<41; i++) {
    g_ESS_shade->SetPoint(i, ESS_strong_x[i], ESS_strong_y[i]);
  }
  for (int i=0; i<42; i++) {
    g_ESS_shade->SetPoint(i+41, ESS_base_x[42-i-1], ESS_base_y[42-i-1]);
  }


  // test shade area betwee ESS base and ESS strong
  TGraph *g_ESS_shade_ver2 = new TGraph ( 41 + 42 + 2 ); // 42 for ESS_base and 41 for ESS_strong
  for (int i=0; i<41; i++) {
    g_ESS_shade_ver2->SetPoint(i, ESS_strong_x[i], ESS_strong_y[i]);
  }
  g_ESS_shade_ver2->SetPoint( 41, pow(10., 21.24), pow(10., -22.) );
  g_ESS_shade_ver2->SetPoint( 42, pow(10., 21.15), pow(10., -22.) );
  for (int i=0; i<42; i++) {
    g_ESS_shade_ver2->SetPoint(i+43, ESS_base_x[42-i-1], ESS_base_y[42-i-1]);
  }


  double ARIANNA_1296_high_x[18]= {
    1.0000E+08,
    1.4734E+08,
    2.1219E+08,
    3.0560E+08,
    4.4011E+08,
    7.1037E+08,
    1.1730E+09,
    1.8507E+09,
    3.2723E+09,
    5.2817E+09,
    8.5250E+09,
    1.5421E+10,
    2.6653E+10,
    4.3019E+10,
    6.9436E+10,
    1.0955E+11,
    1.6894E+11,
    3.1264E+11
  };

  double ARIANNA_1296_high_y[18]= {
    1.2311E-08,
    1.0309E-08,
    8.8539E-09,
    7.4139E-09,
    7.0472E-09,
    6.6987E-09,
    6.8708E-09,
    7.4139E-09,
    8.4161E-09,
    9.5538E-09,
    1.1410E-08,
    1.4335E-08,
    1.7559E-08,
    2.1508E-08,
    2.7022E-08,
    3.3100E-08,
    4.0544E-08,
    5.7825E-08
  };

  double ARIANNA_1296_low_x[17]= {
    9.7746E+07,
    1.4734E+08,
    2.2721E+08,
    3.1264E+08,
    5.5281E+08,
    9.7746E+08,
    1.8089E+09,
    3.2723E+09,
    5.9194E+09,
    1.0000E+10,
    1.6894E+10,
    2.8539E+10,
    4.5026E+10,
    7.1037E+10,
    1.2001E+11,
    1.9817E+11,
    3.0560E+11
  };

  double ARIANNA_1296_low_y[17]= {
    6.0526E-09,
    5.1983E-09,
    4.1375E-09,
    3.5535E-09,
    3.3778E-09,
    3.2108E-09,
    3.6448E-09,
    4.1375E-09,
    5.0681E-09,
    6.0526E-09,
    7.2282E-09,
    8.8539E-09,
    1.1124E-08,
    1.3285E-08,
    1.7119E-08,
    2.2060E-08,
    2.7022E-08
  };

  for (int i=0;i<18;i++) {
    ARIANNA_1296_high_y[i]=ARIANNA_1296_high_y[i]/ARIANNA_1296_high_x[i];
    ARIANNA_1296_high_x[i]=ARIANNA_1296_high_x[i]*1.e9;
  }
  for (int i=0;i<17;i++) {
    ARIANNA_1296_low_y[i]=ARIANNA_1296_low_y[i]/ARIANNA_1296_low_x[i];
    ARIANNA_1296_low_x[i]=ARIANNA_1296_low_x[i]*1.e9;
  }



  TGraph *g_ARIANNA_1296_shade = new TGraph ( 18 + 17 ); // 18 for max and 17 for low
  for (int i=0; i<18; i++) {
    g_ARIANNA_1296_shade->SetPoint(i, ARIANNA_1296_high_x[i], ARIANNA_1296_high_y[i]);
  }
  for (int i=0; i<17; i++) {
    g_ARIANNA_1296_shade->SetPoint(i+18, ARIANNA_1296_low_x[17-i-1], ARIANNA_1296_low_y[17-i-1]);
  }


  double Kotera_low_x[37] = { 14.4422  ,
			      14.5737    ,
			      14.737 ,
			      14.8685    ,
			      15.0045    ,
			      15.1633    ,
			      15.3628    ,
			      15.5714    ,
			      15.7891    ,
			      15.9569    ,
			      16.1519    ,
			      16.3061    ,
			      16.5465    ,
			      16.7914    ,
			      16.941 ,
			      17.1497    ,
			      17.3129    ,
			      17.4762    ,
			      17.6667    ,
			      17.8934    ,
			      18.0794    ,
			      18.2608    ,
			      18.4059    ,
			      18.5601    ,
			      18.6871    ,
			      18.8095    ,
			      18.9546    ,
			      19.0771    ,
			      19.1542    ,
			      19.2449    ,
			      19.3356    ,
			      19.3991    ,
			      19.4898    ,
			      19.5669    ,
			      19.6349    ,
			      19.6757    ,
			      19.712 };

  double Kotera_low_y[37] = { -15.8509   ,
			      -15.8296   ,
			      -15.8083   ,
			      -15.8296   ,
			      -15.8775   ,
			      -15.9308   ,
			      -16    ,
			      -16.1385   ,
			      -16.2557   ,
			      -16.3622   ,
			      -16.5273   ,
			      -16.6551   ,
			      -16.8522   ,
			      -16.9747   ,
			      -16.996    ,
			      -16.9854   ,
			      -16.9587   ,
			      -16.9587   ,
			      -16.996    ,
			      -17.0972   ,
			      -17.2037   ,
			      -17.3901   ,
			      -17.5712   ,
			      -17.8056   ,
			      -18.0666   ,
			      -18.3276   ,
			      -18.6631   ,
			      -18.9987   ,
			      -19.2064   ,
			      -19.4674   ,
			      -19.7603   ,
			      -20    ,
			      -20.3036   ,
			      -20.5965   ,
			      -20.8735   ,
			      -21.0546   ,
			      -21.1718 };

  /*
    for (int i=0; i<37; i++) {
    cout<<Kotera_low_x[i]<<endl;
    }

    cout<<"\n";

    for (int i=0; i<37; i++) {
    cout<<Kotera_low_y[i]<<endl;
    }
  */

  LogToLine(37, Kotera_low_x);
  LogToLine(37, Kotera_low_y);

  TGraph *g_Kotera_low = new TGraph( 37, Kotera_low_x, Kotera_low_y );


  double Kotera_mid_x[34] = { 14.4195    ,
			      14.5692    ,
			      14.7324    ,
			      14.9683    ,
			      15.2177    ,
			      15.4082    ,
			      15.6621    ,
			      15.8889    ,
			      16.1156    ,
			      16.2925    ,
			      16.483 ,
			      16.6054    ,
			      16.7868    ,
			      16.9501    ,
			      17.1315    ,
			      17.3855    ,
			      17.6803    ,
			      17.9751    ,
			      18.2608    ,
			      18.4649    ,
			      18.6417    ,
			      18.8821    ,
			      19.0952    ,
			      19.2268    ,
			      19.3583    ,
			      19.5215    ,
			      19.7438    ,
			      19.9478    ,
			      20.093 ,
			      20.2381    ,
			      20.3469    ,
			      20.4966    ,
			      20.5828    ,
			      20.6961 };

  double Kotera_mid_y[34] = { -14.4874   ,
			      -14.4607   ,
			      -14.466    ,
			      -14.5459   ,
			      -14.6897   ,
			      -14.8229   ,
			      -15.0573   ,
			      -15.2863   ,
			      -15.5579   ,
			      -15.771    ,
			      -16.0852   ,
			      -16.2184   ,
			      -16.3196   ,
			      -16.3728   ,
			      -16.4368   ,
			      -16.5326   ,
			      -16.6498   ,
			      -16.8469   ,
			      -17.0599   ,
			      -17.2517   ,
			      -17.4541   ,
			      -17.7577   ,
			      -18.0346   ,
			      -18.1944   ,
			      -18.3968   ,
			      -18.6365   ,
			      -18.9774   ,
			      -19.3023   ,
			      -19.5952   ,
			      -19.9041   ,
			      -20.1651   ,
			      -20.5433   ,
			      -20.8149   ,
			      -21.1771 };

  LogToLine(34, Kotera_mid_x);
  LogToLine(34, Kotera_mid_y);

  TGraph *g_Kotera_mid = new TGraph( 34, Kotera_mid_x, Kotera_mid_y );


  double Kotera_max_x[42] = { 14.4331    ,
			      14.5918    ,
			      14.7642    ,
			      14.932 ,
			      15.1224    ,
			      15.2857    ,
			      15.3991    ,
			      15.5488    ,
			      15.7166    ,
			      15.9025    ,
			      16.102 ,
			      16.3061    ,
			      16.5057    ,
			      16.6735    ,
			      16.8322    ,
			      16.9864    ,
			      17.1678    ,
			      17.3175    ,
			      17.4853    ,
			      17.6304    ,
			      17.8027    ,
			      18.0567    ,
			      18.3605    ,
			      18.6054    ,
			      18.7732    ,
			      18.9456    ,
			      19.1224    ,
			      19.2993    ,
			      19.4444    ,
			      19.5941    ,
			      19.7483    ,
			      19.8662    ,
			      20.0023    ,
			      20.1111    ,
			      20.229 ,
			      20.3333    ,
			      20.424 ,
			      20.5147    ,
			      20.5918    ,
			      20.6825    ,
			      20.7642    ,
			      20.8594 };

  double Kotera_max_y[42] = { -14.2157   ,
			      -14.1838   ,
			      -14.2051   ,
			      -14.2583   ,
			      -14.3222   ,
			      -14.4181   ,
			      -14.4874   ,
			      -14.5672   ,
			      -14.7217   ,
			      -14.9028   ,
			      -15.1425   ,
			      -15.3928   ,
			      -15.6698   ,
			      -15.8828   ,
			      -15.968    ,
			      -16.016    ,
			      -16    ,
			      -15.9734   ,
			      -15.9467   ,
			      -15.9627   ,
			      -16.0107   ,
			      -16.1438   ,
			      -16.4154   ,
			      -16.6551   ,
			      -16.8682   ,
			      -17.0759   ,
			      -17.273    ,
			      -17.4807   ,
			      -17.6937   ,
			      -17.9228   ,
			      -18.1625   ,
			      -18.3648   ,
			      -18.6152   ,
			      -18.8602   ,
			      -19.1265   ,
			      -19.3768   ,
			      -19.6218   ,
			      -19.8828   ,
			      -20.1385   ,
			      -20.49 ,
			      -20.7776   ,
			      -21.1718 };



  LogToLine(42, Kotera_max_x);
  LogToLine(42, Kotera_max_y);

  TGraph *g_Kotera_max = new TGraph( 42, Kotera_max_x, Kotera_max_y );


  //
  // Kotera shade is obtained from http://arxiv.org/pdf/1009.1382v2.pdf, Fig 9
  // upper bound is the pink dot-dashed line
  // lower bound is the lower bound of gray shaded area
  //
  //
  TGraph *g_Kotera_shade = new TGraph ( 42 + 37 ); // 42 for max and 37 for low
  for (int i=0; i<42; i++) {
    g_Kotera_shade->SetPoint(i, Kotera_max_x[i], Kotera_max_y[i]);
  }
  for (int i=0; i<37; i++) {
    g_Kotera_shade->SetPoint(i+42, Kotera_low_x[37-i-1], Kotera_low_y[37-i-1]);
  }

  // test shade area betwee ESS base and ESS strong
  TGraph *g_Kotera_shade_ver2 = new TGraph ( 42 + 37 + 2 ); // 42 for max and 37 for low
  for (int i=0; i<42; i++) {
    g_Kotera_shade_ver2->SetPoint(i, Kotera_max_x[i], Kotera_max_y[i]);
  }
  g_Kotera_shade_ver2->SetPoint(42, pow(10., 21.05), pow(10., -22) );
  g_Kotera_shade_ver2->SetPoint(43, pow(10., 19.88), pow(10., -22) );
  for (int i=0; i<37; i++) {
    g_Kotera_shade_ver2->SetPoint(i+44, Kotera_low_x[37-i-1], Kotera_low_y[37-i-1]);
  }


  /*

    double GRB_BP_ALL_x[11] = { 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5 };

    double GRB_BP_ALL_y[11] = { -8.52242, -8.52143, -9.00864, -9.39794, -9.71026, -9.88606, -9.96064, -11.5229, -12.0229, -12.5229, -13.0229 }; // this is log[E^2 dN] flux with unit in log [GeV cm^-2 sr^-1 s^-1]

    // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
    for (int i=0; i<11; i++) {
    GRB_BP_ALL_y[i] = GRB_BP_ALL_y[i] - (GRB_BP_ALL_x[i] - 9.);
    }


    LogToLine(11, GRB_BP_ALL_x);
    LogToLine(11, GRB_BP_ALL_y);

    TGraph *g_GRB_BP_ALL = new TGraph( 11, GRB_BP_ALL_x, GRB_BP_ALL_y );

  */




  double GRB_max_x[13] = { 11.0585 ,
			   11.7016    ,
			   12.8222    ,
			   13.7576    ,
			   14.3325    ,
			   14.8392    ,
			   15.2972    ,
			   15.9598    ,
			   16.6711    ,
			   17.2655    ,
			   18.0646    ,
			   18.7954    ,
			   19.9354    };

  double GRB_max_y[13] = { -9.92308  ,
			   -9.67692   ,
			   -9.1641    ,
			   -8.72308   ,
			   -8.45641   ,
			   -8.30256   ,
			   -8.2   ,
			   -8.05641   ,
			   -7.93333   ,
			   -8.05641   ,
			   -8.35385   ,
			   -8.65128   ,
			   -9.13333   };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<13; i++) {
    GRB_max_y[i] = GRB_max_y[i] - (GRB_max_x[i] - 9.);
  }

  double GRB_min_x[15] = { 11.0202   ,
			   12.5439    ,
			   13.2862    ,
			   13.5399    ,
			   13.9206    ,
			   14.2428    ,
			   14.4966    ,
			   14.8869    ,
			   15.2379    ,
			   15.4329    ,
			   15.7061    ,
			   15.8031    ,
			   16.1049    ,
			   16.28  ,
			   16.9415    };

  double GRB_min_y[15] = { -12.4062  ,
			   -10.9337   ,
			   -10.1975   ,
			   -10.1159   ,
			   -9.97302   ,
			   -9.7789    ,
			   -9.64611   ,
			   -9.5647    ,
			   -9.74945   ,
			   -9.85209   ,
			   -9.86268   ,
			   -10.3645   ,
			   -10.9178   ,
			   -11.3071   ,
			   -12.9769   };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<15; i++) {
    GRB_min_y[i] = GRB_min_y[i] - (GRB_min_x[i] - 9.);
  }

  LogToLine(13, GRB_max_x);
  LogToLine(13, GRB_max_y);
  LogToLine(15, GRB_min_x);
  LogToLine(15, GRB_min_y);

  TGraph *g_GRB_BATSE = new TGraph( 13, GRB_max_x, GRB_max_y );

  // test shade area betwee ESS base and ESS strong
  TGraph *g_GRB_shade = new TGraph ( 13 + 15 ); // 42 for max and 37 for low
  for (int i=0; i<13; i++) {
    g_GRB_shade->SetPoint(i, GRB_max_x[i], GRB_max_y[i]);
  }
  for (int i=0; i<15; i++) {
    g_GRB_shade->SetPoint(i+13, GRB_min_x[15-i-1], GRB_min_y[15-i-1]);
  }




  double GRB_IC_FC_x[13] = { 4.29727 ,   // this is in GeV already
			     4.67336    ,
			     4.849  ,
			     5.04783    ,
			     5.3001 ,
			     5.56767    ,
			     5.78115    ,
			     5.91063    ,
			     6.05538    ,
			     6.25354    ,
			     6.35269    ,
			     6.48978    ,
			     7.27744    };

  double GRB_IC_FC_y[13] = {-2.00018 , // unit of log[ E^2 F ] where GeV cm^-2 (no time or str there!!!)
			    -1.31443   ,
			    -1.22613   ,
			    -1.05247   ,
			    -0.857497  ,
			    -0.647294  ,
			    -0.629132  ,
			    -0.653601  ,
			    -0.671982  ,
			    -0.675152  ,
			    -0.656919  ,
			    -0.684442  ,
			    -2.002 };

  double normal = 1./ ( 24. * 60. * 60. * 4. * PI ); //one day in sec with 4 pi unit [1/s * str]

  // change unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<13; i++) {
    GRB_IC_FC_x[i] = GRB_IC_FC_x[i] + 9.;   // from GeV to eV
    GRB_IC_FC_y[i] = GRB_IC_FC_y[i] - (GRB_IC_FC_x[i] - 9.); // now in cm^-2
    //GRB_IC_FC_y[i] = GRB_IC_FC_y[i] * normal; // now in cm^-2 sec^-1 str^-1
  }

  LogToLine(13, GRB_IC_FC_x);
  LogToLine(13, GRB_IC_FC_y);

  // change unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<13; i++) {
    GRB_IC_FC_y[i] = GRB_IC_FC_y[i] * normal; // now in cm^-2 sec^-1 str^-1 (in linear)
    //std::cout<<"GRB_IC_FC["<<i<<"] : "<<log( GRB_IC_FC_y[i])<<std::endl;
  }


  TGraph *g_GRB_IC_FC = new TGraph( 13, GRB_IC_FC_x, GRB_IC_FC_y );


  double GRB_NFC_x[18] = { 4.89942   ,
			   5.03716    ,
			   5.20528    ,
			   5.33516    ,
			   5.59466    ,
			   5.90753    ,
			   6.10588    ,
			   6.2584 ,
			   6.4108 ,
			   6.58595    ,
			   6.73825    ,
			   6.86774    ,
			   6.99731    ,
			   7.12681    ,
			   7.26376    ,
			   7.38531    ,
			   7.51432    ,
			   7.63567    };

  double GRB_NFC_y[18] = { -1.99445  ,
			   -1.85124   ,
			   -1.73854   ,
			   -1.6563    ,
			   -1.565 ,
			   -1.46763   ,
			   -1.41897   ,
			   -1.39772   ,
			   -1.41001   ,
			   -1.44975   ,
			   -1.48642   ,
			   -1.50785   ,
			   -1.50792   ,
			   -1.52934   ,
			   -1.59345   ,
			   -1.69719   ,
			   -1.84665   ,
			   -2.00222   };

  //double normal = 1./ ( 24. * 60. * 60. * 4. * PI ); //one day in sec with 4 pi unit [1/s * str]

  // change unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<18; i++) {
    GRB_NFC_x[i] = GRB_NFC_x[i] + 9.;   // from GeV to eV
    GRB_NFC_y[i] = GRB_NFC_y[i] - (GRB_NFC_x[i] - 9.); // now in cm^-2
  }

  LogToLine(18, GRB_NFC_x);
  LogToLine(18, GRB_NFC_y);

  // change unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<18; i++) {
    GRB_NFC_y[i] = GRB_NFC_y[i] * normal; // now in cm^-2 sec^-1 str^-1 (in linear)
    //std::cout<<"GRB_IC_FC["<<i<<"] : "<<log( GRB_IC_FC_y[i])<<std::endl;
  }


  TGraph *g_GRB_NFC = new TGraph( 18, GRB_NFC_x, GRB_NFC_y );

  // test shade area for GRB IC-FC and NFC
  TGraph *g_GRB_shade2 = new TGraph ( 13 + 18 ); // 42 for max and 37 for low
  for (int i=0; i<13; i++) {
    g_GRB_shade2->SetPoint(i, GRB_IC_FC_x[i], GRB_IC_FC_y[i]);
  }
  for (int i=0; i<18; i++) {
    g_GRB_shade2->SetPoint(i+13, GRB_NFC_x[18-i-1], GRB_NFC_y[18-i-1]);
  }



  double GRB_WIND_MAX_x[21] = { 12.4094  ,
				13.1726    ,
				14.0468    ,
				14.9696    ,
				15.5872    ,
				16.3851    ,
				16.8361    ,
				17.2801    ,
				17.7728    ,
				18.1891    ,
				18.3695    ,
				18.5568    ,
				18.7441    ,
				18.9107    ,
				19.0772    ,
				19.2923    ,
				19.4657    ,
				19.6461    ,
				19.778 ,
				19.882 ,
				19.9722    };

  double GRB_WIND_MAX_y[21] = { -13.4734 ,
				-12.5503   ,
				-11.5207   ,
				-10.5385   ,
				-9.95858   ,
				-9.35503   ,
				-9.08284   ,
				-8.95266   ,
				-8.88166   ,
				-8.85799   ,
				-8.89349   ,
				-9.02367   ,
				-9.21302   ,
				-9.49704   ,
				-9.84024   ,
				-10.3609   ,
				-10.9172   ,
				-11.5917   ,
				-12.1598   ,
				-12.6568   ,
				-13.1657   };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<21; i++) {
    GRB_WIND_MAX_y[i] = GRB_WIND_MAX_y[i] - (GRB_WIND_MAX_x[i] - 9.);
  }

  LogToLine(21, GRB_WIND_MAX_x);
  LogToLine(21, GRB_WIND_MAX_y);

  TGraph *g_GRB_WIND_MAX = new TGraph( 21, GRB_WIND_MAX_x, GRB_WIND_MAX_y );


  double GRB_WIND_MIN_x[21] = {14.1023   ,
			       14.5117    ,
			       14.9003    ,
			       15.3027    ,
			       15.6357    ,
			       16.0173    ,
			       16.3781    ,
			       16.7112    ,
			       17.0095    ,
			       17.4536    ,
			       17.7728    ,
			       17.9879    ,
			       18.2723    ,
			       18.4874    ,
			       18.7095    ,
			       18.8968    ,
			       19.0286    ,
			       19.2229    ,
			       19.3894    ,
			       19.5212    ,
			       19.6184    };

  double GRB_WIND_MIN_y[21] = { -13.4852 ,
				-13.0237   ,
				-12.6213   ,
				-12.2426   ,
				-11.9349   ,
				-11.6391   ,
				-11.3669   ,
				-11.1775   ,
				-11.0355   ,
				-10.9172   ,
				-10.8817   ,
				-10.858    ,
				-10.8698   ,
				-10.9882   ,
				-11.1775   ,
				-11.4734   ,
				-11.7574   ,
				-12.1834   ,
				-12.6686   ,
				-13.1302   ,
				-13.497    };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<21; i++) {
    GRB_WIND_MIN_y[i] = GRB_WIND_MIN_y[i] - (GRB_WIND_MIN_x[i] - 9.);
  }

  LogToLine(21, GRB_WIND_MIN_x);
  LogToLine(21, GRB_WIND_MIN_y);


  // test shade area for GRB WIND shade
  TGraph *g_GRB_WIND_shade = new TGraph ( 21 + 21 ); // 42 for max and 37 for low
  for (int i=0; i<21; i++) {
    g_GRB_WIND_shade->SetPoint(i, GRB_WIND_MAX_x[i], GRB_WIND_MAX_y[i]);
  }
  for (int i=0; i<21; i++) {
    g_GRB_WIND_shade->SetPoint(i+21, GRB_WIND_MIN_x[21-i-1], GRB_WIND_MIN_y[21-i-1]);
  }



  double GRB_ISM_MAX_x[22] = {14.2134    ,
			      14.7962    ,
			      15.2402    ,
			      15.6219    ,
			      16.1006    ,
			      16.4961    ,
			      16.9263    ,
			      17.2454    ,
			      17.6062    ,
			      17.8838    ,
			      18.1683    ,
			      18.4666    ,
			      18.7233    ,
			      19.0078    ,
			      19.1674    ,
			      19.3478    ,
			      19.4727    ,
			      19.6114    ,
			      19.7433    ,
			      19.8404    ,
			      19.9167    ,
			      19.9792    };

  double GRB_ISM_MAX_y[22] = { -13.4852  ,
			       -12.787    ,
			       -12.2544   ,
			       -11.8166   ,
			       -11.2722   ,
			       -10.8462   ,
			       -10.432    ,
			       -10.1361   ,
			       -9.8639    ,
			       -9.69823   ,
			       -9.61539   ,
			       -9.55621   ,
			       -9.55621   ,
			       -9.60355   ,
			       -9.69823   ,
			       -9.87574   ,
			       -10.0769   ,
			       -10.4083   ,
			       -10.7988   ,
			       -11.1775   ,
			       -11.568    ,
			       -11.8876   };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<22; i++) {
    GRB_ISM_MAX_y[i] = GRB_ISM_MAX_y[i] - (GRB_ISM_MAX_x[i] - 9.);
  }

  LogToLine(22, GRB_ISM_MAX_x);
  LogToLine(22, GRB_ISM_MAX_y);



  double GRB_ISM_MIN_x[19] = {15.941 ,
			      16.2047    ,
			      16.4406    ,
			      16.8083    ,
			      17.1275    ,
			      17.4467    ,
			      17.6618    ,
			      17.9115    ,
			      18.1197    ,
			      18.4042    ,
			      18.6539    ,
			      18.8968    ,
			      19.1049    ,
			      19.2437    ,
			      19.4172    ,
			      19.5837    ,
			      19.7086    ,
			      19.8057    ,
			      19.8959    };

  double GRB_ISM_MIN_y[19] = { -13.4734  ,
			       -13.1893   ,
			       -12.929    ,
			       -12.5503   ,
			       -12.2544   ,
			       -12.0059   ,
			       -11.8521   ,
			       -11.7101   ,
			       -11.6154   ,
			       -11.5562   ,
			       -11.568    ,
			       -11.5917   ,
			       -11.6746   ,
			       -11.7811   ,
			       -11.9941   ,
			       -12.3373   ,
			       -12.716    ,
			       -13.0947   ,
			       -13.4852   };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<19; i++) {
    GRB_ISM_MIN_y[i] = GRB_ISM_MIN_y[i] - (GRB_ISM_MIN_x[i] - 9.);
  }

  LogToLine(19, GRB_ISM_MIN_x);
  LogToLine(19, GRB_ISM_MIN_y);



  // test shade area for GRB WIND shade
  TGraph *g_GRB_ISM_shade = new TGraph ( 22 + 19 ); // 42 for max and 37 for low
  for (int i=0; i<22; i++) {
    g_GRB_ISM_shade->SetPoint(i, GRB_ISM_MAX_x[i], GRB_ISM_MAX_y[i]);
  }
  for (int i=0; i<19; i++) {
    g_GRB_ISM_shade->SetPoint(i+22, GRB_ISM_MIN_x[19-i-1], GRB_ISM_MIN_y[19-i-1]);
  }

  // test shade area for GRB WIND shade
  TGraph *g_GRB_ISM_shade_ver2 = new TGraph ( 22 + 19 + 1 ); // 42 for max and 37 for low
  for (int i=0; i<22; i++) {
    g_GRB_ISM_shade_ver2->SetPoint(i, GRB_ISM_MAX_x[i], GRB_ISM_MAX_y[i]);
  }
  for (int i=0; i<19; i++) {
    g_GRB_ISM_shade_ver2->SetPoint(i+22, GRB_ISM_MIN_x[19-i-1], GRB_ISM_MIN_y[19-i-1]);
  }
  g_GRB_ISM_shade_ver2->SetPoint(19+22, pow(10.,14), pow(10.,-21.1));






  double GRB_WB99_x[4] = { 12.0234   ,
			   13.9922    ,
			   15.9844    ,
			   19.0078    };

  double GRB_WB99_y[4] = { -10.4975  ,
			   -8.51737   ,
			   -8.52605   ,
			   -11.5484   };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<4; i++) {
    GRB_WB99_y[i] = GRB_WB99_y[i] - (GRB_WB99_x[i] - 9.);
  }

  LogToLine(4, GRB_WB99_x);
  LogToLine(4, GRB_WB99_y);

  TGraph *g_GRB_WB99 = new TGraph( 4, GRB_WB99_x, GRB_WB99_y );



  double f_pi = 0.05;  // not sure at all!!!

  double GRB_WB_AG_x[11];
  double GRB_WB_AG_y[11];

  for (int i=0; i<11; i++) {

    GRB_WB_AG_x[i] = pow(10., 14. + i*0.5);

    if ( GRB_WB_AG_x[i] < pow(10., 17.) ) {
      std::cout<<"lower than 10^17"<<std::endl;
      GRB_WB_AG_y[i] = f_pi * pow(10.,-10.) / 0.1 * (GRB_WB_AG_x[i] / pow(10.,17.));  // unit GeV cm^-2 s^-1 sr^-1 (linear)
    }
    else {
      std::cout<<"higher than 10^17"<<std::endl;
      GRB_WB_AG_y[i] = f_pi * pow(10.,-10.) / 0.1 * pow( (GRB_WB_AG_x[i] / pow(10.,17.)), 0.5 );  // unit GeV cm^-2 s^-1 sr^-1 (linear)
    }

  }

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<11; i++) {
    GRB_WB_AG_y[i] = GRB_WB_AG_y[i] / (GRB_WB_AG_x[i] / pow(10.,9.));
  }

  TGraph *g_GRB_WB_AG = new TGraph( 11, GRB_WB_AG_x, GRB_WB_AG_y );





  double GRB_Hummer_SE_MAX_x[29] = {13.298   ,
				    13.4954    ,
				    13.7587    ,
				    13.9634    ,
				    14.0951    ,
				    14.234 ,
				    14.3876    ,
				    14.5996    ,
				    14.8629    ,
				    15.0676    ,
				    15.2578    ,
				    15.3821    ,
				    15.5064    ,
				    15.6307    ,
				    15.7477    ,
				    15.8501    ,
				    15.9671    ,
				    16.0695    ,
				    16.1865    ,
				    16.2962    ,
				    16.3912    ,
				    16.5155    ,
				    16.6252    ,
				    16.7569    ,
				    16.8665    ,
				    16.9835    ,
				    17.1298    ,
				    17.2907    ,
				    17.4808    };

  double GRB_Hummer_SE_MAX_y[29] = { -10.4152    ,
				     -10.1336   ,
				     -9.80144   ,
				     -9.57401   ,
				     -9.46209   ,
				     -9.37184   ,
				     -9.28881   ,
				     -9.21661   ,
				     -9.12996   ,
				     -9.07581   ,
				     -9.05054   ,
				     -9.05054   ,
				     -9.07942   ,
				     -9.1083    ,
				     -9.14079   ,
				     -9.15523   ,
				     -9.15523   ,
				     -9.15884   ,
				     -9.19856   ,
				     -9.25271   ,
				     -9.34296   ,
				     -9.48014   ,
				     -9.62816   ,
				     -9.81227   ,
				     -9.9639    ,
				     -10.1227   ,
				     -10.2744   ,
				     -10.3791   ,
				     -10.4946   };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<29; i++) {
    GRB_Hummer_SE_MAX_y[i] = GRB_Hummer_SE_MAX_y[i] - (GRB_Hummer_SE_MAX_x[i] - 9.);
  }

  LogToLine(29, GRB_Hummer_SE_MAX_x);
  LogToLine(29, GRB_Hummer_SE_MAX_y);



  double GRB_Hummer_SE_MIN_x[24] = { 13.5905 ,
				     13.7002    ,
				     13.8026    ,
				     13.9781    ,
				     14.117 ,
				     14.2998    ,
				     14.5046    ,
				     14.7459    ,
				     14.9799    ,
				     15.1993    ,
				     15.3382    ,
				     15.4771    ,
				     15.6088    ,
				     15.7404    ,
				     15.872 ,
				     16.0256    ,
				     16.1645    ,
				     16.2815    ,
				     16.3839    ,
				     16.4644    ,
				     16.6179    ,
				     16.7422    ,
				     16.8373    ,
				     16.9104    };

  double GRB_Hummer_SE_MIN_y[24] = { -10.4982    ,
				     -10.3466   ,
				     -10.2166   ,
				     -10.0433   ,
				     -9.90975   ,
				     -9.79061   ,
				     -9.70758   ,
				     -9.63177   ,
				     -9.56318   ,
				     -9.50903   ,
				     -9.50181   ,
				     -9.52347   ,
				     -9.56679   ,
				     -9.59567   ,
				     -9.61011   ,
				     -9.62094   ,
				     -9.63899   ,
				     -9.70758   ,
				     -9.79061   ,
				     -9.8917    ,
				     -10.0758   ,
				     -10.2599   ,
				     -10.4043   ,
				     -10.4946   };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<24; i++) {
    GRB_Hummer_SE_MIN_y[i] = GRB_Hummer_SE_MIN_y[i] - (GRB_Hummer_SE_MIN_x[i] - 9.);
  }

  LogToLine(24, GRB_Hummer_SE_MIN_x);
  LogToLine(24, GRB_Hummer_SE_MIN_y);


  // test shade area for GRB WIND shade
  TGraph *g_GRB_Hummer_shade = new TGraph ( 29 + 24 ); // 42 for max and 37 for low
  for (int i=0; i<29; i++) {
    g_GRB_Hummer_shade->SetPoint(i, GRB_Hummer_SE_MAX_x[i], GRB_Hummer_SE_MAX_y[i]);
  }
  for (int i=0; i<24; i++) {
    g_GRB_Hummer_shade->SetPoint(i+29, GRB_Hummer_SE_MIN_x[24-i-1], GRB_Hummer_SE_MIN_y[24-i-1]);
  }


  double GRB_Hummer_AU_MAX_x[40] = { 13.3272  ,
				     13.5101 ,
				     13.6417 ,
				     13.7733 ,
				     13.8976 ,
				     14.0219 ,
				     14.1536 ,
				     14.2925 ,
				     14.3949 ,
				     14.4973 ,
				     14.6069 ,
				     14.7166 ,
				     14.819  ,
				     14.936  ,
				     15.0165 ,
				     15.1335 ,
				     15.2505 ,
				     15.3748 ,
				     15.4698 ,
				     15.5503 ,
				     15.6307 ,
				     15.6892 ,
				     15.777  ,
				     15.9013 ,
				     16.011  ,
				     16.0914 ,
				     16.1865 ,
				     16.3035 ,
				     16.4205 ,
				     16.5375 ,
				     16.6618 ,
				     16.8007 ,
				     16.9031 ,
				     16.9835 ,
				     17.1152 ,
				     17.2907 ,
				     17.4004 ,
				     17.5686 ,
				     17.7441 ,
				     17.9049 };

  double GRB_Hummer_AU_MAX_y[40] = { -9.1769 ,
				     -9.00361   ,
				     -8.9278    ,
				     -8.87004   ,
				     -8.81588   ,
				     -8.75812   ,
				     -8.70397   ,
				     -8.63899   ,
				     -8.60289   ,
				     -8.57762   ,
				     -8.57401   ,
				     -8.59567   ,
				     -8.62094   ,
				     -8.65343   ,
				     -8.67148   ,
				     -8.66787   ,
				     -8.66065   ,
				     -8.66426   ,
				     -8.68953   ,
				     -8.74729   ,
				     -8.81227   ,
				     -8.87004   ,
				     -8.94585   ,
				     -9.02888   ,
				     -9.01444   ,
				     -9.00722   ,
				     -9.02527   ,
				     -9.06137   ,
				     -9.15884   ,
				     -9.28881   ,
				     -9.44043   ,
				     -9.59567   ,
				     -9.74368   ,
				     -9.87365   ,
				     -10    ,
				     -10.1119   ,
				     -10.1877   ,
				     -10.2491   ,
				     -10.3141   ,
				     -10.4188   };


  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<40; i++) {
    GRB_Hummer_AU_MAX_y[i] = GRB_Hummer_AU_MAX_y[i] - (GRB_Hummer_AU_MAX_x[i] - 9.);
  }

  LogToLine(40, GRB_Hummer_AU_MAX_x);
  LogToLine(40, GRB_Hummer_AU_MAX_y);


  double GRB_Hummer_AU_MIN_x[13] = { 14.4826  ,
				     14.6362 ,
				     14.7751 ,
				     14.8995 ,
				     15.0603 ,
				     15.2066 ,
				     15.3894 ,
				     15.6161 ,
				     15.7404 ,
				     15.9086 ,
				     16.1133 ,
				     16.2596 ,
				     16.5082 };

  double GRB_Hummer_AU_MIN_y[13] = { -10.4946 ,
				     -10.3791    ,
				     -10.2599    ,
				     -10.1733    ,
				     -10.1155    ,
				     -10.0614    ,
				     -10.0181    ,
				     -9.97834    ,
				     -10.0181    ,
				     -10.0686    ,
				     -10.1227    ,
				     -10.2166    ,
				     -10.4982    };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<13; i++) {
    GRB_Hummer_AU_MIN_y[i] = GRB_Hummer_AU_MIN_y[i] - (GRB_Hummer_AU_MIN_x[i] - 9.);
  }

  LogToLine(13, GRB_Hummer_AU_MIN_x);
  LogToLine(13, GRB_Hummer_AU_MIN_y);

  // test shade area for GRB WIND shade
  TGraph *g_GRB_Hummer_shade2 = new TGraph ( 40 + 13 ); // 42 for max and 37 for low
  for (int i=0; i<40; i++) {
    g_GRB_Hummer_shade2->SetPoint(i, GRB_Hummer_AU_MAX_x[i], GRB_Hummer_AU_MAX_y[i]);
  }
  for (int i=0; i<13; i++) {
    g_GRB_Hummer_shade2->SetPoint(i+40, GRB_Hummer_AU_MIN_x[13-i-1], GRB_Hummer_AU_MIN_y[13-i-1]);
  }





  /*
    double IC40_Hummer_x[20] = { 13.6022    ,
    13.6464 ,
    13.6464 ,
    13.6759 ,
    13.8004 ,
    13.8738 ,
    13.9618 ,
    14.0132 ,
    14.1305 ,
    14.2259 ,
    14.3432 ,
    14.4532 ,
    14.5632 ,
    14.6949 ,
    14.7827 ,
    14.9142 ,
    15.0092 ,
    15.1701 ,
    15.2944 ,
    15.3676 };

    double IC40_Hummer_y[20] = { -9.30352   ,
    -9.22053    ,
    -9.22053    ,
    -9.15557    ,
    -9.09792    ,
    -9.04023    ,
    -8.9681 ,
    -8.90317    ,
    -8.83107    ,
    -8.73729    ,
    -8.65076    ,
    -8.56782    ,
    -8.49572    ,
    -8.47779    ,
    -8.46704    ,
    -8.48521    ,
    -8.51056    ,
    -8.51793    ,
    -8.50721    ,
    -8.50006    };

    // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
    for (int i=0; i<20; i++) {
    IC40_Hummer_y[i] = IC40_Hummer_y[i] - (IC40_Hummer_x[i] - 9.);
    }

    LogToLine(20, IC40_Hummer_x);
    LogToLine(20, IC40_Hummer_y);


    TGraph *g_IC40_Hummer = new TGraph( 20, IC40_Hummer_x, IC40_Hummer_y );
  */

  double IC40_Hummer_x[10] = { 13.6022    ,
			       //13.6464 ,
			       13.6464 ,
			       //13.6759 ,
			       13.8004 ,
			       //13.8738 ,
			       13.9618 ,
			       //14.0132 ,
			       14.1305 ,
			       //14.2259 ,
			       //14.3432 ,
			       14.4532 ,
			       //14.5632 ,
			       14.6949 ,
			       //14.7827 ,
			       14.9142 ,
			       //15.0092 ,
			       15.1701 ,
			       //15.2944 ,
			       15.3676 };

  double IC40_Hummer_y[10] = { -9.30352   ,
			       //-9.22053    ,
			       -9.22053    ,
			       //-9.15557    ,
			       -9.09792    ,
			       //-9.04023    ,
			       -8.9681 ,
			       //-8.90317    ,
			       -8.83107    ,
			       //-8.73729    ,
			       //-8.65076    ,
			       -8.56782    ,
			       //-8.49572    ,
			       -8.47779    ,
			       //-8.46704    ,
			       -8.48521    ,
			       //-8.51056    ,
			       -8.51793    ,
			       //-8.50721    ,
			       -8.50006    };

  // change unit to log[ E dN ] with unit in log [cm^-2 sr^-2 s^-1]
  for (int i=0; i<10; i++) {
    IC40_Hummer_y[i] = IC40_Hummer_y[i] - (IC40_Hummer_x[i] - 9.);
  }

  LogToLine(10, IC40_Hummer_x);
  LogToLine(10, IC40_Hummer_y);


  TGraph *g_IC40_Hummer = new TGraph( 10, IC40_Hummer_x, IC40_Hummer_y );






  double IC_new_x[13] = {
    14,        14.5,        15,        15.5,        16,        16.5,        17,        17.5,        18,        18.5,        19,        19.5,        20
  };


  double IC_new_y[13] = {
    -9.83847,        -10.8223,        -12.6195,        -13.6823,        -14.4311,        -14.9294,        -15.3288,        -15.7242,        -16.0738,        -16.3775,        -16.64342,        -16.863,        -17.07565
  };



  LogToLine(13, IC_new_x);
  LogToLine(13, IC_new_y);

  //EF_to_E2F(13, IC_new_x, IC_new_y);


  TGraph *g_IC_new = new TGraph( 13, IC_new_x, IC_new_y );





  double IC_evts_HESE_x[3] = { 13.77889, 14.54418, 15.30157 };


  double IC_evts_HESE_y[3] = { -12.22865, -12.99501, -13.75347 };

  LogToLine(3, IC_evts_HESE_x);
  LogToLine(3, IC_evts_HESE_y);

  //EF_to_E2F(13, IC_new_x, IC_new_y);


  TGraph *g_IC_evts_HESE = new TGraph( 3, IC_evts_HESE_x, IC_evts_HESE_y );






  // WB limit from IC plot

  double WB98_3p2_x [45] = {
    13.97632,         14.09466,         14.19722,         14.29979,         14.38657,         14.48125,         14.56803,         14.6706,         14.78105,         14.8915,         15.00985,         15.12819,         15.24653,         15.37277,         15.499,         15.62523,         15.74357,         15.86981,         15.99604,         16.12227,         16.25639,         16.38263,         16.51675,         16.64298,         16.76921,         16.89545,         17.02957,         17.17158,         17.29781,         17.43193,         17.56606,         17.69229,         17.81063,         17.99209,         18.17355,         18.35501,         18.53647,         18.71793,         18.89939,         19.08085,         19.26231,         19.44377,         19.62523,         19.80669,         19.98815 };

  double WB98_3p2_y [45] = {
    -12.16196,         -12.2803,         -12.38286,         -12.48543,         -12.57221,         -12.66689,         -12.75367,         -12.85624,         -12.96669,         -13.07714,         -13.19549,         -13.31383,         -13.43217,         -13.55841,         -13.68464,         -13.81087,         -13.92921,         -14.05545,         -14.18168,         -14.30791,         -14.44203,         -14.56827,         -14.70239,         -14.82862,         -14.95485,         -15.08109,         -15.21521,         -15.35722,         -15.48345,         -15.61757,         -15.7517,         -15.87793,         -15.99627,         -16.17773,         -16.35919,         -16.54065,         -16.72211,         -16.90357,         -17.08503,         -17.26649,         -17.44795,         -17.62941,         -17.81087,         -17.99233,         -18.17379 };

  LogToLine(45, WB98_3p2_x);
  LogToLine(45, WB98_3p2_y);

  //EF_to_E2F(13, IC_new_x, IC_new_y);


  TGraph *g_WB98_3p2 = new TGraph( 45, WB98_3p2_x, WB98_3p2_y );



  const int n_ANITA_4 = 7;

  // TFile* fLimit = TFile::Open("limitNums.root");
  // TGraph* grLimit = (TGraph*) fLimit->Get("grLimit");

  // LogToLine(grLimit->GetN(), grLimit->GetX());
  Double_t ANITA_4_x[n_ANITA_4]       = {18, 18.5, 19, 19.5, 20, 20.5, 21};
  Double_t ANITA_4_y[n_ANITA_4];
  Double_t ANITA_all_y[n_ANITA_4];
  Double_t ANITA_4_y_low[n_ANITA_4];
  Double_t ANITA_4_y_high[n_ANITA_4];  

  Double_t ANITA_4_effArea[n_ANITA_4] = {  0.00194,       // E18     in km^2  
					 0.0376,      // E18.5		  
					 0.62948,      // E19    	  
					 3.47982,      // E19.5		  
					 14.65220,     // E20		  
					 47.48070,     // E20.5		  
					117.02800 };    // E21   

  double ANITA_4_anaEff  = 0.8;
  
  double N90 = 2.30; 
  double ANITA_4_livetime = 27.3*24*3600; // //27.3*24*3600; // 27.3 days


  Double_t ANITA_3_effArea[n_ANITA_4] = { 0.00160,      // E18     in km^2  
					0.03451,      // E18.5		  
					0.34288,      // E19    	  
					1.75735,      // E19.5		  
					7.14396,      // E20		  
					22.77670,     // E20.5		  
					58.72190};    // E21
  
  Double_t ANITA_2_effArea[n_ANITA_4] = { 0.00029 ,      // E18     in km^2  
					0.02121 ,      // E18.5		  
					0.22009 ,      // E19    	  
					1.27190 ,      // E19.5		  
					6.38188 ,      // E20		  
					18.45390 ,     // E20.5		  
					52.53270 };    // E21              
  
  double ANITA_all_effArea[n_ANITA_4];
  
  double ANITA_1_livetime = 17.4*24*3600.;
  double ANITA_2_livetime = 28.5*24*3600.;
  double ANITA_3_livetime = 17.4*24*3600.; // 17.4 days
  double ANITA_all_livetime = ANITA_1_livetime+ANITA_2_livetime+ANITA_3_livetime+ANITA_4_livetime;
  for (int ibin=0; ibin<n_ANITA_4; ibin++){
    ANITA_all_effArea[ibin] = (
			       (ANITA_1_livetime+ANITA_2_livetime)*ANITA_2_effArea[ibin] +
			       ANITA_3_livetime*ANITA_3_effArea[ibin] +
			       ANITA_4_livetime*ANITA_4_effArea[ibin]
			       )/(ANITA_all_livetime);
  }

  std::cout << "HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERE===========================================" << std::endl;

  for (int i=0; i<n_ANITA_4; i++){

    // convert km^2 in cm^2
    ANITA_4_effArea[i]*=1e10;
    ANITA_all_effArea[i]*=1e10;
    
    ANITA_4_y[i] = N90/(ANITA_4_effArea[i]*ANITA_4_livetime*ANITA_4_anaEff*TMath::Log(10.));
    
    ANITA_4_y[i] *= TMath::Power(10, ANITA_4_x[i])/(TMath::Power(10, ANITA_4_x[i]+0.25) - TMath::Power(10, ANITA_4_x[i]-0.25));
    
    // Divide by 4 
    ANITA_4_y[i] /= 4.;
    
    ANITA_all_y[i] = N90/(ANITA_all_effArea[i]*ANITA_all_livetime*ANITA_4_anaEff*TMath::Log(10.));

    
    ANITA_all_y[i] *= TMath::Power(10, ANITA_4_x[i])/(TMath::Power(10, ANITA_4_x[i]+0.25) - TMath::Power(10, ANITA_4_x[i]-0.25));
    
    // Divide by 4 
    ANITA_all_y[i] /= 4.;
    
    ANITA_4_y_low[i]  = ANITA_4_y[i]*0.5;
    ANITA_4_y_high[i] = ANITA_4_y[i]*1.5;
    
    //   ANITA_4_y[i]*= TMath::Power(10, ANITA_4_x[i])/(TMath::Power(10, ANITA_4_x[i]+0.25)-TMath::Power(10, ANITA_4_x[i]-0.25));
    std::cout << ANITA_4_x[i] << " " << ANITA_4_y[i] << std::endl;
  }
      
      
  //  double ANITA_4_y[n_ANITA_4] = {3.97188e-17, 9.1763e-19};
      LogToLine(n_ANITA_4, ANITA_4_x);
    // LogToLine(grLimit->GetN(), grLimit->GetY());
  // TGraph* g_ANITA_4 = grLimit;
  TGraph *g_ANITA_4 = new TGraph(n_ANITA_4, ANITA_4_x, ANITA_4_y);
  TGraph *g_ANITA_all = new TGraph(n_ANITA_4, ANITA_4_x, ANITA_all_y);
  
  // TGraphErrors *g_ANITA_4_shade = new TGraphErrors(n_ANITA_4);
  // for (int i=0; i<n_ANITA_4; i++){
  //   g_ANITA_4_shade->SetPoint(i,       ANITA_4_x[i], ANITA_4_y[i] );
  //   g_ANITA_4_shade->SetPointError(i,   0, ANITA_4_y[i]*0.1 );
  // }
  // g_ANITA_4_shade->SetFillStyle(3001);
  // g_ANITA_4_shade->SetFillColor(kMagenta);
  // g_ANITA_4_shade->SetLineColor(kMagenta);
  
  
  // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  // // for(int i=0; i < g_ANITA_4->GetN(); i++){
  // for(int i=g_ANITA_4->GetN()-1; i >= 0; i--){
  //   if((i%2)==0){
  //     g_ANITA_4->RemovePoint(i);
  //   }
  //   // std::cout << g_ANITA_4->GetX()[i] << "\t" << g_ANITA_4->GetY()[i] << std::endl;
  //   std::cout << g_ANITA_4->GetX()[i] << "\t" << g_ANITA_4->GetY()[i] << std::endl;
  //   // g_ANITA_4->GetY()[i]*=g_ANITA_4->GetX()[i];
  // }
  // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  double ANITA_x[8] = { 18.0068  ,
			18.5102    ,
			19.0091    ,
			19.5034    ,
			20.0068    ,
			20.5102    ,
			21.0045    ,
			21.5079 };

  double ANITA_y[8] = { -13.0599 ,
			-15.1372   ,
			-16.4048   ,
			-17.2197   ,
			-17.8961   ,
			-18.3063   ,
			-18.6844   ,
			-19.0573 };

  LogToLine(8, ANITA_x);
  LogToLine(8, ANITA_y);

  TGraph *g_ANITA = new TGraph( 8, ANITA_x, ANITA_y );



  // ANITA_II Erratum
  //

  double ANITA_erratum_x[11] = { 18.00692,
				 18.49481,
				 18.9827,
				 19.5017,
				 19.9896,
				 20.4983,
				 20.9862,
				 21.4948,
				 22.0035,
				 22.5017,
				 23 };

  double ANITA_erratum_y[11] = { -12.91039624,
				 -14.93512624,
				 -16.21946624,
				 -17.02028924,
				 -17.70023424,
				 -18.13842024,
				 -18.50105624,
				 -18.86369624,
				 -19.15078624,
				 -19.39254624,
				 -19.57385624 };

  LogToLine(11, ANITA_erratum_x);
  LogToLine(11, ANITA_erratum_y);

  // for (int i=0; i<11; i++) cout << ANITA_erratum_x[i] << " " << ANITA_erratum_y[i] << endl;
  TGraph *g_ANITA_erratum = new TGraph( 11, ANITA_erratum_x, ANITA_erratum_y );





  double Auger_x[14] = { 17.1406 ,
			 17.2766    ,
			 17.4036    ,
			 17.4989    ,
			 17.7347    ,
			 18.0113    ,
			 18.288 ,
			 18.5238    ,
			 18.8413    ,
			 19.2132    ,
			 19.6168    ,
			 19.9342    ,
			 20.2245    ,
			 20.4739 };

  double Auger_y[14] = { -12.5273    ,
			 -13.2783   ,
			 -13.9334   ,
			 -14.2423   ,
			 -14.5992   ,
			 -14.9241   ,
			 -15.1638   ,
			 -15.3076   ,
			 -15.4248   ,
			 -15.494    ,
			 -15.5313   ,
			 -15.5419   ,
			 -15.5473   ,
			 -15.5206 };

  LogToLine(14, Auger_x);
  LogToLine(14, Auger_y);

  TGraph *g_Auger = new TGraph( 14, Auger_x, Auger_y );




  double Auger11_x[7] = {
    17.0056,
    17.4993,
    18.0057,
    18.4999,
    19.0004,
    19.5009,
    20.0015 };

  double Auger11_y[7] = {
    -13.88807,
    -14.9759,
    -15.563,
    -15.873,
    -16.0552,
    -16.1735,
    -16.2065 };

  LogToLine(7, Auger11_x);
  LogToLine(7, Auger11_y);

  TGraph *g_Auger11 = new TGraph( 7, Auger11_x, Auger11_y );


  // new Auger result from http://arxiv.org/pdf/1304.1630v1.pdf
  // fig. 11, Earth skim (3.5 yr), differential points
  //

  double Auger13_x[7] = { 16.9953,        17.4905,            17.9913,            18.4981,            18.9985,            19.4866,            19.9869 }; // values in log(eV)

  double Auger13_y[7] = { -6.27033,        -6.87202,            -6.91867,            -6.75478,            -6.43779,            -6.05384,            -5.58375 }; // values in log( E^2F ) where E^2F is [GeV cm-2 s-1 sr-1]

  // change E^2F to EF
  for (int i=0; i<7; i++) {

    Auger13_y[i] = Auger13_y[i] - ( Auger13_x[i] - 9. );
  }

  LogToLine(7, Auger13_x);
  LogToLine(7, Auger13_y);

  // multiply factor of 3 to account all three flavors
  for (int i=0; i<7; i++) {
    Auger13_y[i] = Auger13_y[i]*3.;
  }

  double Auger13_y0[7] = { 0., 0., 0., 0., 0., 0., 0. };

  // get x bar range (each direction, 0.25 decade)
  double Auger13_xl[7];
  double Auger13_xh[7];
  for (int i=0; i<7; i++) {
    Auger13_xl[i] = Auger13_x[i] - (Auger13_x[i]/pow(10., 0.25));

    Auger13_xh[i] = (Auger13_x[i]*pow(10., 0.25)) - Auger13_x[i];
  }



  TGraphAsymmErrors *g_Auger13_bar = new TGraphAsymmErrors( 7, Auger13_x, Auger13_y, Auger13_xl, Auger13_xh, Auger13_y0, Auger13_y0 ); // only x range bar

  TGraph *g_Auger13 = new TGraph( 7, Auger13_x, Auger13_y );


  double Auger15_x[8] = { 16.70533, 17.21212, 17.71369, 18.22048, 18.72205, 19.22362, 19.73041, 20.22675 }; // values in log(eV)

  double Auger15_y[8] = {-6.82741, -7.44332, -7.66667, -7.58545, -7.32826, -6.98985, -6.60406, -6.17090 }; // values in log( E^2F ) where E^2F is [GeV cm-2 s-1 sr-1]

  for (int i=0; i<8; i++) {

    Auger15_y[i] = Auger15_y[i] - ( Auger15_x[i] - 9. );
  }

  LogToLine(7, Auger15_x);
  LogToLine(7, Auger15_y);

  // multiply factor of 3 to account all three flavors
  for (int i=0; i<8; i++) {
    Auger15_y[i] = Auger15_y[i]*3.;
  }
  TGraph *g_Auger15 = new TGraph( 7, Auger15_x, Auger15_y );

















  double Icecube_x[17] = { 14.7551   ,
			   14.932 ,
			   15.1769    ,
			   15.3628    ,
			   15.5941    ,
			   15.8707    ,
			   16.2426    ,
			   16.6327    ,
			   17.0227    ,
			   17.4127    ,
			   17.8481    ,
			   18.3288    ,
			   18.7732    ,
			   19.195 ,
			   19.5941    ,
			   19.8707    ,
			   20.0884 };

  double Icecube_y[17] = { -12.1491  ,
			   -12.5539   ,
			   -13.0652   ,
			   -13.3262   ,
			   -13.5819   ,
			   -13.8642   ,
			   -14.1624   ,
			   -14.4767   ,
			   -14.7856   ,
			   -15.0839   ,
			   -15.3768   ,
			   -15.6751   ,
			   -15.9414   ,
			   -16.1598   ,
			   -16.3089   ,
			   -16.3782   ,
			   -16.3515 };

  LogToLine(17, Icecube_x);
  LogToLine(17, Icecube_y);

  TGraph *g_Icecube = new TGraph( 17, Icecube_x, Icecube_y );


  double RICE_x[11] = { 16.07463 ,
			16.50204   ,
			16.83446   ,
			17.21438   ,
			17.64179   ,
			18.08277   ,
			18.60516   ,
			19.0393    ,
			19.46  ,
			20.0366    ,
			20.7083 };

  double RICE_y[11] = { -11.51402    ,
			-12.69143  ,
			-13.43294  ,
			-14.1462   ,
			-14.81603  ,
			-15.42368  ,
			-16.0294   ,
			-16.45597  ,
			-16.82364  ,
			-17.26387  ,
			-17.69315 };

  LogToLine(11, RICE_x);
  LogToLine(11, RICE_y);

  TGraph *g_RICE = new TGraph( 11, RICE_x, RICE_y );


  //double offset_peter = 0.36; // 0.36 offset is when a factor 2.3 is applied (numerator) for 90% CL
  double offset_peter = 0.; // 0. offset is when a factor 2.3 is not applied for 90% CL

  double ARA37_Peter_x[10] = { 16.003    ,
			       16.4796    ,
			       16.997 ,
			       17.4962    ,
			       18.0045    ,
			       18.4947    ,
			       19.003 ,
			       19.5023    ,
			       19.997 ,
			       20.4962};

  double ARA37_Peter_y[10] = { -14.221+offset_peter   ,
			       -15.494 +offset_peter   ,
			       -16.3515 +offset_peter  ,
			       -17.0333 +offset_peter  ,
			       -17.5126 +offset_peter  ,
			       -17.8748 +offset_peter  ,
			       -18.1571 +offset_peter  ,
			       -18.4181 +offset_peter  ,
			       -18.6525 +offset_peter  ,
			       -18.8868+offset_peter};

  LogToLine(10, ARA37_Peter_x);
  LogToLine(10, ARA37_Peter_y);

  TGraph *g_ARA37_Peter = new TGraph( 10, ARA37_Peter_x, ARA37_Peter_y );








  // factor of two for 4 pi -> 2 pi
  //double str_offset = 2.;
  double str_offset = 1.; // no need to apply offset
  //











  //double offset = 0.72; // 0.36 offset is when a factor 2.3 is applied (numerator) for 90% CL
  double offset = 0.36; // 0.36 offset is when a factor 2.3 is not applied for 90% CL
  //
  //double offset = 0.; // 0.36 offset is when a factor 2.3 is applied (denominator) for 90% CL // WRONG!!

  //double ARA37_x[5] = { 16, 17, 18, 19, 20};
  //double ARA37_x[9] = { 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20};
  double ARA37_x[8] = { 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20};

  //double ARA37_y[5] = { -15.22+offset, -16.71+offset, -17.65+offset, -18.25+offset, -18.71+offset };
  //double ARA37_y[9] = { -14.86, -15.495, -16.34, -16.84, -17.286, -17.654, -17.888, -18.174, -18.345 };
  double ARA37_y[8] = { -15.495, -16.34, -16.84, -17.286, -17.654, -17.888, -18.174, -18.345 };

  //LogToLine(9, ARA37_x);
  //LogToLine(9, ARA37_y);
  LogToLine(8, ARA37_x);
  LogToLine(8, ARA37_y);

  //TGraph *g_ARA37 = new TGraph( 9, ARA37_x, ARA37_y );
  TGraph *g_ARA37 = new TGraph( 8, ARA37_x, ARA37_y );

  double error_shift = 0.398; // factor 2.5 in log
  //double ARA37_up_y[9] = { -14.86+error_shift, -15.495+error_shift, -16.34+error_shift, -16.84+error_shift, -17.286+error_shift, -17.654+error_shift, -17.888+error_shift, -18.174+error_shift, -18.345+error_shift };
  double ARA37_up_y[8] = { -15.495+error_shift, -16.34+error_shift, -16.84+error_shift, -17.286+error_shift, -17.654+error_shift, -17.888+error_shift, -18.174+error_shift, -18.345+error_shift };

  //LogToLine(9, ARA37_up_y);
  LogToLine(8, ARA37_up_y);


  //TGraph *g_ARA37_up = new TGraph( 9, ARA37_x, ARA37_up_y );
  TGraph *g_ARA37_up = new TGraph( 8, ARA37_x, ARA37_up_y );



  //TGraph *g_ARA37_shade = new TGraph ( 9 + 9 );
  TGraph *g_ARA37_shade = new TGraph ( 8 + 8 );
  //for (int i=0; i<9; i++) {
  for (int i=0; i<8; i++) {
    g_ARA37_shade->SetPoint(i, ARA37_x[i], ARA37_up_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<8; i++) {
    g_ARA37_shade->SetPoint(i+8, ARA37_x[8-i-1], ARA37_y[8-i-1]);
  }



  // ARA37 3 yrs after fixing the bugs
  // from Tdomain result
  //
  double ARA37_3yr_Tdomain_debug_x[5] = { 17, 18, 19, 20, 21 };

  //double ARA37_3yr_debug_y[4] = { -17.1272274, -18.01233351, -18.59531273, -18.96846617 }; // before LPM effect factor
  double ARA37_3yr_Tdomain_debug_y[5] = { -15.69944978, -17.01478844, -17.86863938, -18.44779273, -18.86937923 }; // after LPM effect factor


  double ARA37_3yr_Fdomain_debug_y[5] = { -16.51058925, -17.51840278, -18.206287, -18.67688602, -19.00828299 }; // LPM effect on


  LogToLine(5, ARA37_3yr_Tdomain_debug_x);

  LogToLine(5, ARA37_3yr_Tdomain_debug_y);

  LogToLine(5, ARA37_3yr_Fdomain_debug_y);


  TGraph *g_ARA37_debug = new TGraph ( 5, ARA37_3yr_Tdomain_debug_x, ARA37_3yr_Tdomain_debug_y );
  //



  //TGraph *g_ARA37_shade = new TGraph ( 9 + 9 );
  TGraph *g_ARA37_shade_debug = new TGraph ( 5 + 5 );
  //for (int i=0; i<9; i++) {
  for (int i=0; i<5; i++) {
    g_ARA37_shade_debug->SetPoint(i, ARA37_3yr_Tdomain_debug_x[i], ARA37_3yr_Tdomain_debug_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<5; i++) {
    g_ARA37_shade_debug->SetPoint(i+5, ARA37_3yr_Tdomain_debug_x[5-i-1], ARA37_3yr_Fdomain_debug_y[5-i-1]);
  }






  // ARA37 3 yrs after fixing the bugs and use same POSNU_RADIUS as old set (bugged version) -> oldR
  // from Tdomain result
  //
  double ARA37_3yr_Tdomain_debug_oldR_x[3] = { 17, 19, 21 };

  double ARA37_3yr_Tdomain_debug_oldR_y[3] = { -15.68481778, -17.87660923, -18.84290393 }; // after LPM effect factor

  double ARA37_3yr_Fdomain_debug_oldR_y[3] = { -16.44011024, -18.11373247, -18.87821451 }; // LPM effect on


  LogToLine(3, ARA37_3yr_Tdomain_debug_oldR_x);

  LogToLine(3, ARA37_3yr_Tdomain_debug_oldR_y);

  LogToLine(3, ARA37_3yr_Fdomain_debug_oldR_y);


  TGraph *g_ARA37_debug_oldR = new TGraph ( 3, ARA37_3yr_Tdomain_debug_oldR_x, ARA37_3yr_Tdomain_debug_oldR_y );
  //



  //TGraph *g_ARA37_shade = new TGraph ( 9 + 9 );
  TGraph *g_ARA37_shade_debug_oldR = new TGraph ( 3 + 3 );
  //for (int i=0; i<9; i++) {
  for (int i=0; i<3; i++) {
    g_ARA37_shade_debug_oldR->SetPoint(i, ARA37_3yr_Tdomain_debug_oldR_x[i], ARA37_3yr_Tdomain_debug_oldR_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<3; i++) {
    g_ARA37_shade_debug_oldR->SetPoint(i+3, ARA37_3yr_Tdomain_debug_oldR_x[3-i-1], ARA37_3yr_Fdomain_debug_oldR_y[3-i-1]);
  }





  // ARA37 3 yrs after fixing the bugs and use same POSNU_RADIUS as new set (bugged version), don't apply phase response from the detector
  // from Tdomain result
  //
  double ARA37_3yr_Tdomain_debug_NoPhase_x[3] = { 17, 19, 21 };

  double ARA37_3yr_Tdomain_debug_NoPhase_y[3] = { -15.83861331, -17.90028506, -18.87336119 }; // after LPM effect factor

  double ARA37_3yr_Fdomain_debug_NoPhase_y[5] = { -16.51058925, -18.206287, -19.00828299 }; // LPM effect on


  LogToLine(3, ARA37_3yr_Tdomain_debug_NoPhase_x);

  LogToLine(3, ARA37_3yr_Tdomain_debug_NoPhase_y);

  LogToLine(3, ARA37_3yr_Fdomain_debug_NoPhase_y);


  TGraph *g_ARA37_debug_NoPhase = new TGraph ( 3, ARA37_3yr_Tdomain_debug_NoPhase_x, ARA37_3yr_Tdomain_debug_NoPhase_y );
  //



  //TGraph *g_ARA37_shade = new TGraph ( 9 + 9 );
  TGraph *g_ARA37_shade_debug_NoPhase = new TGraph ( 3 + 3 );
  //for (int i=0; i<9; i++) {
  for (int i=0; i<3; i++) {
    g_ARA37_shade_debug_NoPhase->SetPoint(i, ARA37_3yr_Tdomain_debug_NoPhase_x[i], ARA37_3yr_Tdomain_debug_NoPhase_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<3; i++) {
    g_ARA37_shade_debug_NoPhase->SetPoint(i+3, ARA37_3yr_Tdomain_debug_NoPhase_x[3-i-1], ARA37_3yr_Fdomain_debug_NoPhase_y[3-i-1]);
  }




  // test oldRF (Fdomain) mode, with additional 10deg off cut applied (like Tdomain)
  //
  double ARA37_3yr_Fdomain_debug_10degCut_y[3] = { -16.44117777, -18.16807101, -18.9478815 }; // LPM effect on

  LogToLine(3, ARA37_3yr_Fdomain_debug_10degCut_y);

  TGraph *g_ARA37_shade_debug_10degCut = new TGraph ( 3 + 3 );
  //for (int i=0; i<9; i++) {
  for (int i=0; i<3; i++) {
    g_ARA37_shade_debug_10degCut->SetPoint(i, ARA37_3yr_Tdomain_debug_NoPhase_x[i], ARA37_3yr_Tdomain_debug_NoPhase_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<3; i++) {
    g_ARA37_shade_debug_10degCut->SetPoint(i+3, ARA37_3yr_Tdomain_debug_NoPhase_x[3-i-1], ARA37_3yr_Fdomain_debug_10degCut_y[3-i-1]);
  }





  // test increased RAYSOL_RANGE from 5km to 10km
  //
  double ARA37_3yr_RAYSOL10km_x[2] = { 17, 21 };

  double ARA37_3yr_Fdomain_debug_RAYSOL10km_y[2] = { -16.45116502, -19.03873193 }; // LPM effect on

  double ARA37_3yr_Tdomain_debug_RAYSOL10km_y[2] = { -15.83238865, -18.92906607 }; // LPM effect accounted


  LogToLine(2, ARA37_3yr_RAYSOL10km_x);
  LogToLine(2, ARA37_3yr_Fdomain_debug_RAYSOL10km_y);
  LogToLine(2, ARA37_3yr_Tdomain_debug_RAYSOL10km_y);

  TGraph *g_ARA37_shade_debug_RAYSOL10km = new TGraph ( 2 + 2 );
  //for (int i=0; i<9; i++) {
  for (int i=0; i<2; i++) {
    g_ARA37_shade_debug_RAYSOL10km->SetPoint(i, ARA37_3yr_RAYSOL10km_x[i], ARA37_3yr_Tdomain_debug_RAYSOL10km_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<2; i++) {
    g_ARA37_shade_debug_RAYSOL10km->SetPoint(i+2, ARA37_3yr_RAYSOL10km_x[2-i-1], ARA37_3yr_Fdomain_debug_RAYSOL10km_y[2-i-1]);
  }







  // ARA37 3 yrs after fixing the bugs and use same POSNU_RADIUS as new set (bugged version), don't apply phase response from the detector
  // use same shower length for all energies (same as 10^16 eV shower)
  // from Tdomain result
  //
  double ARA37_3yr_Tdomain_debug_SameShowerL_y[3] = { -15.87789628, -17.93215395, -18.92906607 }; // after LPM effect factor


  LogToLine(3, ARA37_3yr_Tdomain_debug_SameShowerL_y);



  //TGraph *g_ARA37_shade = new TGraph ( 9 + 9 );
  TGraph *g_ARA37_shade_debug_SameShowerL = new TGraph ( 3 + 3 );
  //for (int i=0; i<9; i++) {
  for (int i=0; i<3; i++) {
    g_ARA37_shade_debug_SameShowerL->SetPoint(i, ARA37_3yr_Tdomain_debug_NoPhase_x[i], ARA37_3yr_Tdomain_debug_SameShowerL_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<3; i++) {
    g_ARA37_shade_debug_SameShowerL->SetPoint(i+3, ARA37_3yr_Tdomain_debug_NoPhase_x[3-i-1], ARA37_3yr_Fdomain_debug_NoPhase_y[3-i-1]);
  }




  // test increased Vsat (to 10^8)
  // NoPhase mode in Tdomain,
  // Raysol range default 5km
  // no offcone cut to Fdomain
  // PosnuR new (10km, 13km)
  //
  double ARA37_3yr_Vsat_x[2] = { 17, 21 };

  double ARA37_3yr_Fdomain_debug_Vsat_y[2] = { -16.44743474, -18.96290223 }; // LPM effect on

  double ARA37_3yr_Tdomain_debug_Vsat_y[2] = { -15.81833614, -18.88307494 }; // LPM effect accounted

  LogToLine(2, ARA37_3yr_Vsat_x);
  LogToLine(2, ARA37_3yr_Fdomain_debug_Vsat_y);
  LogToLine(2, ARA37_3yr_Tdomain_debug_Vsat_y);

  TGraph *g_ARA37_shade_debug_Vsat = new TGraph ( 2 + 2 );
  //for (int i=0; i<9; i++) {
  for (int i=0; i<2; i++) {
    g_ARA37_shade_debug_Vsat->SetPoint(i, ARA37_3yr_Vsat_x[i], ARA37_3yr_Tdomain_debug_Vsat_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<2; i++) {
    g_ARA37_shade_debug_Vsat->SetPoint(i+2, ARA37_3yr_Vsat_x[2-i-1], ARA37_3yr_Fdomain_debug_Vsat_y[2-i-1]);
  }






  // test new MakeArrays function for Fdomain
  // test increased Vsat (to 10^8)
  // NoPhase mode in Tdomain,
  // Raysol range 10km
  // no offcone cut to Fdomain
  // PosnuR new (10km, 13km)
  //
  double ARA37_3yr_MakeArraySat_x[2] = { 17, 21 };

  double ARA37_3yr_Fdomain_debug_MakeArraySat_y[2] = { -16.09143751, -19.02142348 }; // LPM effect on

  double ARA37_3yr_Tdomain_debug_MakeArraySat_y[2] = { -15.81253952, -18.94510683 }; // LPM effect accounted

  LogToLine(2, ARA37_3yr_MakeArraySat_x);
  LogToLine(2, ARA37_3yr_Fdomain_debug_MakeArraySat_y);
  LogToLine(2, ARA37_3yr_Tdomain_debug_MakeArraySat_y);

  TGraph *g_ARA37_shade_debug_MakeArraySat = new TGraph ( 2 + 2 );
  //for (int i=0; i<9; i++) {
  for (int i=0; i<2; i++) {
    g_ARA37_shade_debug_MakeArraySat->SetPoint(i, ARA37_3yr_MakeArraySat_x[i], ARA37_3yr_Tdomain_debug_MakeArraySat_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<2; i++) {
    g_ARA37_shade_debug_MakeArraySat->SetPoint(i+2, ARA37_3yr_MakeArraySat_x[2-i-1], ARA37_3yr_Fdomain_debug_MakeArraySat_y[2-i-1]);
  }





  // same above but
  // with Phase in Tdomain,
  //
  /*
    double ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_y[2] = { -15.68733137, -18.94056872 }; // LPM effect accounted

    LogToLine(2, ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_y);

    TGraph *g_ARA37_shade_debug_MakeArraySat_wPhase = new TGraph ( 2 + 2 );
    //for (int i=0; i<9; i++) {
    for (int i=0; i<2; i++) {
    g_ARA37_shade_debug_MakeArraySat_wPhase->SetPoint(i, ARA37_3yr_MakeArraySat_x[i], ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_y[i]);
    }
    //for (int i=0; i<9; i++) {
    for (int i=0; i<2; i++) {
    g_ARA37_shade_debug_MakeArraySat_wPhase->SetPoint(i+2, ARA37_3yr_MakeArraySat_x[2-i-1], ARA37_3yr_Fdomain_debug_MakeArraySat_y[2-i-1]);
    }

  */






  // ARA37 3 yrs after fixing the bugs, new MakeArraysforFFT, normalizations, Saturations
  // new results for both oldRF, newRF modes
  //
  double ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_x[5] = { 17, 18, 19, 20, 21 };

  double ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_y[5] = { -15.7089, -16.9754, -17.8609, -18.4642, -18.9438 }; // after LPM effect factor

  double ARA37_3yr_Fdomain_debug_MakeArraySat_wPhase_y[5] = { -15.92446568, -17.15433689, -18.0054038, -18.57129149, -19.02351568 }; // LPM effect on


  LogToLine(5, ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_x);

  LogToLine(5, ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_y);

  LogToLine(5, ARA37_3yr_Fdomain_debug_MakeArraySat_wPhase_y);


  // apply str correction factor from 4pi to 2pi
  for (int bin=0; bin<5; bin++) {

    ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_y[bin] *= str_offset;

    ARA37_3yr_Fdomain_debug_MakeArraySat_wPhase_y[bin] *= str_offset;
  }



  TGraph *g_ARA37_shade_debug_MakeArraySat_wPhase = new TGraph ( 5 + 5 );
  for (int i=0; i<5; i++) {
    g_ARA37_shade_debug_MakeArraySat_wPhase->SetPoint(i, ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_x[i], ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<5; i++) {
    g_ARA37_shade_debug_MakeArraySat_wPhase->SetPoint(i+5, ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_x[5-i-1], ARA37_3yr_Fdomain_debug_MakeArraySat_wPhase_y[5-i-1]);
  }



  for (int i=0; i<5; i++) {
    //pint difference between two models
    //

    cout<<"at "<<ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_x[i]<<" eV, from the default "<< ARA37_3yr_Fdomain_debug_MakeArraySat_wPhase_y[i] /ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_y[i] <<" factor from Tdomain mode (Fdomain / Tdomain)"<<endl;

  }



  for (int i=0; i<5; i++) {
    //pint difference between two models
    //

    cout<<"at "<<ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_x[i]<<" eV, from the default "<< ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_y[i]/ARA37_3yr_Fdomain_debug_MakeArraySat_wPhase_y[i] <<" factor from Tdomain mode (Tdomain / Fdomain)"<<endl;

  }




  // ARA37 3 yrs after fixing the bugs, new MakeArraysforFFT, normalizations, Saturations
  // new results for both oldRF, newRF modes
  // Tdomain wo Phase
  //
  double ARA37_3yr_Tdomain_debug_MakeArraySat_woPhase_x[5] = { 17, 18, 19, 20, 21 };

  double ARA37_3yr_Tdomain_debug_MakeArraySat_woPhase_y[5] = { -15.8044, -17.0794, -17.9099, -18.4704, -18.9579 }; // after LPM effect factor

  double ARA37_3yr_Fdomain_debug_MakeArraySat_woPhase_y[5] = { -15.92446568, -17.15433689, -18.0054038, -18.57129149, -19.02351568 }; // LPM effect on


  LogToLine(5, ARA37_3yr_Tdomain_debug_MakeArraySat_woPhase_x);

  LogToLine(5, ARA37_3yr_Tdomain_debug_MakeArraySat_woPhase_y);

  LogToLine(5, ARA37_3yr_Fdomain_debug_MakeArraySat_woPhase_y);


  TGraph *g_ARA37_shade_debug_MakeArraySat_woPhase = new TGraph ( 5 + 5 );
  for (int i=0; i<5; i++) {
    g_ARA37_shade_debug_MakeArraySat_woPhase->SetPoint(i, ARA37_3yr_Tdomain_debug_MakeArraySat_woPhase_x[i], ARA37_3yr_Tdomain_debug_MakeArraySat_woPhase_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<5; i++) {
    g_ARA37_shade_debug_MakeArraySat_woPhase->SetPoint(i+5, ARA37_3yr_Tdomain_debug_MakeArraySat_woPhase_x[5-i-1], ARA37_3yr_Fdomain_debug_MakeArraySat_woPhase_y[5-i-1]);
  }






  // test ARA37 Fdomain wo secondaries
  //
  double ARA37_3yr_debug_MakeArraySat_woSecondary_x[5] = { 17, 18, 19, 20, 21 };

  double ARA37_3yr_Tdomain_debug_MakeArraySat_woSecondary_y[5] = { -15.6320, -16.9561, -17.7814, -18.3732, -18.8908 }; // LPM effect accounted

  double ARA37_3yr_Fdomain_debug_MakeArraySat_woSecondary_y[5] = { -15.88310477, -17.09541769, -17.93644085, -18.51426185, -18.95488316 }; // LPM effect on

  LogToLine(5, ARA37_3yr_debug_MakeArraySat_woSecondary_x);

  LogToLine(5, ARA37_3yr_Tdomain_debug_MakeArraySat_woSecondary_y);

  LogToLine(5, ARA37_3yr_Fdomain_debug_MakeArraySat_woSecondary_y);

  //TGraph *g_ARA37_Fdomain_debug_MakeArraySat_woSecondary = new TGraph ( 5, ARA37_3yr_Fdomain_debug_MakeArraySat_woSecondary_x, ARA37_3yr_Fdomain_debug_MakeArraySat_woSecondary_y );
  TGraph *g_ARA37_shade_debug_MakeArraySat_woSecondary = new TGraph ( 5 + 5 + 1 ); // last connection
  for (int i=0; i<5; i++) {
    g_ARA37_shade_debug_MakeArraySat_woSecondary->SetPoint(i, ARA37_3yr_debug_MakeArraySat_woSecondary_x[i], ARA37_3yr_Tdomain_debug_MakeArraySat_woSecondary_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<5; i++) {
    g_ARA37_shade_debug_MakeArraySat_woSecondary->SetPoint(i+5, ARA37_3yr_debug_MakeArraySat_woSecondary_x[5-i-1], ARA37_3yr_Fdomain_debug_MakeArraySat_woSecondary_y[5-i-1]);
  }
  g_ARA37_shade_debug_MakeArraySat_woSecondary->SetPoint(5+5, ARA37_3yr_debug_MakeArraySat_woSecondary_x[0], ARA37_3yr_Tdomain_debug_MakeArraySat_woSecondary_y[0]);





  // try to apply expected analysis efficiency based on Testbed with expected improvements
  //
  double DeepStationEfficiencies[5] = { 0.252411, 0.388351, 0.408001, 0.402681, 0.28912 };

  // apply expected efficiency to ARA37
  double ARA37_3yr_Tdomain_ExpAnalysis_y[5];
  double ARA37_3yr_Fdomain_ExpAnalysis_y[5];

  for (int bin=0; bin<5; bin++) {

    ARA37_3yr_Tdomain_ExpAnalysis_y[bin] = ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_y[bin] / DeepStationEfficiencies[bin];

    ARA37_3yr_Fdomain_ExpAnalysis_y[bin] = ARA37_3yr_Fdomain_debug_MakeArraySat_wPhase_y[bin] / DeepStationEfficiencies[bin];
  }

  TGraph *g_ARA37_shade_ExpAnalysis = new TGraph ( 5 + 5 );
  for (int i=0; i<5; i++) {
    g_ARA37_shade_ExpAnalysis->SetPoint(i, ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_x[i], ARA37_3yr_Tdomain_ExpAnalysis_y[i]);
  }
  //for (int i=0; i<9; i++) {
  for (int i=0; i<5; i++) {
    g_ARA37_shade_ExpAnalysis->SetPoint(i+5, ARA37_3yr_Tdomain_debug_MakeArraySat_wPhase_x[5-i-1], ARA37_3yr_Fdomain_ExpAnalysis_y[5-i-1]);
  }


















  //
  // Thomas' ARA2 & 3 limit
  //
  double ARA2_3_Thomas_x[10] = { 16., 16.5, 17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5 };

  double ARA2_3_Thomas_y_limit[10] = { 1.69910786982e-05, 4.37725904955e-06, 2.06553689551e-06, 1.34756618659e-06, 1.22484233758e-06, 1.4107730099e-06, 1.82651912547e-06, 3.13593754551e-06, 4.70454551863e-06,  9.17992277445e-06 };

  double ARA2_3_Thomas_y_lower[10] = { 5.90402431425e-06, 2.25190476981e-06, 1.33235426426e-06, 1.01398457667e-06, 9.37209085932e-07, 9.99657179314e-07, 1.17517807396e-06, 1.88080289613e-06, 2.46068175412e-06, 4.13387361751e-06 };

  double ARA2_3_Thomas_y_upper[10] = { 5.91455791925e-05, 1.23336135464e-05, 4.82400145283e-06, 2.75573835673e-06, 2.29144913147e-06, 2.49764437533e-06, 3.33882646265e-06, 6.18896662598e-06, 1.02199648477e-05, 2.30232912e-05 };

  double ARA2_3_Thomas_y_limit_20150430[10] = { 0.00038408633228, 2.18791782261e-05, 7.10104734165e-06, 4.15985851743e-06, 2.98380028661e-06, 2.59657510311e-06, 2.87760250081e-06, 4.34606544278e-06, 6.00080292597e-06,  1.10441829088e-05 };

  double ARA2_3_Thomas_y_lower_20150430[10] = { 0.000133461511469, 1.12558624585e-05, 4.58046076395e-06, 3.13011147042e-06, 2.28310587692e-06, 1.83990261029e-06, 1.85144262514e-06, 2.60658650019e-06, 3.13868071028e-06, 4.97338130999e-06 };

  double ARA2_3_Thomas_y_upper_20150430[10] = { 0.00133699625469 , 6.16480144076e-05, 1.65842899089e-05, 8.50680418453e-06 , 5.58212789145e-06 , 4.59699835188e-06, 5.26017781293e-06, 8.57722884763e-06, 1.30359021331e-05, 2.76988647314e-05 };


  LogToLine(10, ARA2_3_Thomas_x);

  // change from E^2*F(E) in GeV/(cm^2 s^1 sr^1) to E*F(E) in 1/(cm^2 s^1 sr^1)
  for (int i=0; i<10; i++) {
    ARA2_3_Thomas_y_limit[i] = ARA2_3_Thomas_y_limit[i] / (ARA2_3_Thomas_x[i]/1.e9);

    ARA2_3_Thomas_y_lower[i] = ARA2_3_Thomas_y_lower[i] / (ARA2_3_Thomas_x[i]/1.e9);

    ARA2_3_Thomas_y_upper[i] = ARA2_3_Thomas_y_upper[i] / (ARA2_3_Thomas_x[i]/1.e9);

    ARA2_3_Thomas_y_limit_20150430[i] = ARA2_3_Thomas_y_limit_20150430[i] / (ARA2_3_Thomas_x[i]/1.e9);

    ARA2_3_Thomas_y_lower_20150430[i] = ARA2_3_Thomas_y_lower_20150430[i] / (ARA2_3_Thomas_x[i]/1.e9);

    ARA2_3_Thomas_y_upper_20150430[i] = ARA2_3_Thomas_y_upper_20150430[i] / (ARA2_3_Thomas_x[i]/1.e9);


  }


  TGraph *g_ARA2_3_Thomas_points = new TGraph ( 10, ARA2_3_Thomas_x, ARA2_3_Thomas_y_limit );

  TGraph *g_ARA2_3_Thomas_points_low = new TGraph ( 10, ARA2_3_Thomas_x, ARA2_3_Thomas_y_lower );

  TGraph *g_ARA2_3_Thomas_points_20150430 = new TGraph ( 10, ARA2_3_Thomas_x, ARA2_3_Thomas_y_limit_20150430 );

  TGraph *g_ARA2_3_Thomas_points_low_20150430 = new TGraph ( 10, ARA2_3_Thomas_x, ARA2_3_Thomas_y_lower_20150430 );



  TGraph *g_ARA2_3_Thomas_shade = new TGraph ( 10 + 10 );
  for (int i=0; i<10; i++) {
    g_ARA2_3_Thomas_shade->SetPoint(i, ARA2_3_Thomas_x[i], ARA2_3_Thomas_y_lower[i]);
  }
  for (int i=0; i<10; i++) {
    g_ARA2_3_Thomas_shade->SetPoint(i+10, ARA2_3_Thomas_x[10-i-1], ARA2_3_Thomas_y_upper[10-i-1]);
  }

  TGraph *g_ARA2_3_Thomas_shade_20150430 = new TGraph ( 10 + 10 );
  for (int i=0; i<10; i++) {
    g_ARA2_3_Thomas_shade_20150430->SetPoint(i, ARA2_3_Thomas_x[i], ARA2_3_Thomas_y_lower_20150430[i]);
  }
  for (int i=0; i<10; i++) {
    g_ARA2_3_Thomas_shade_20150430->SetPoint(i+10, ARA2_3_Thomas_x[10-i-1], ARA2_3_Thomas_y_upper_20150430[10-i-1]);
  }


  // error bar ver
  //
  //double ARA2_3_Thomas_dx = 2.5;
  double ARA2_3_Thomas_dx = pow(10.,0.5)/1.8; // half of half decade

  double ARA2_3_Thomas_h_x[10];
  double ARA2_3_Thomas_l_x[10];
  double ARA2_3_Thomas_h_y[10];
  double ARA2_3_Thomas_l_y[10];

  for (int i=0; i<10; i++) {

    ARA2_3_Thomas_h_x[i] = ARA2_3_Thomas_x[i]*ARA2_3_Thomas_dx - ARA2_3_Thomas_x[i];
    ARA2_3_Thomas_l_x[i] = ARA2_3_Thomas_x[i] - ARA2_3_Thomas_x[i]/ARA2_3_Thomas_dx;


    ARA2_3_Thomas_h_y[i] = ARA2_3_Thomas_y_upper[i] - ARA2_3_Thomas_y_limit[i];
    ARA2_3_Thomas_l_y[i] = ARA2_3_Thomas_y_limit[i] - ARA2_3_Thomas_y_lower[i];
  }

  TGraphAsymmErrors *g_ARA2_3_Thomas_error = new TGraphAsymmErrors( 10, ARA2_3_Thomas_x, ARA2_3_Thomas_y_limit, ARA2_3_Thomas_l_x, ARA2_3_Thomas_h_x, ARA2_3_Thomas_l_y, ARA2_3_Thomas_h_y );

  TGraph *g_ARA2_3_Thomas = new TGraph( 10, ARA2_3_Thomas_x, ARA2_3_Thomas_y_limit);
  TGraph *g_ARA2_3_Thomas_20150430 = new TGraph( 10, ARA2_3_Thomas_x, ARA2_3_Thomas_y_limit_20150430);






















  // make Veff from old vs. new results for ARA37
  //
  //
  double E_Veff_we_new[5] = { 17, 18, 19, 20, 21 };

  double Veff_we_new_Fdomain[5] = { 22.1668, 168.86, 576.55, 1074.68, 1608.88 };

  TGraph *gVeff_new_Fdomain = new TGraph( 5, E_Veff_we_new, Veff_we_new_Fdomain );


  double Veff_we_new_Tdomain[5] = { 15.5007, 144.435, 574.76, 1179.47, 1681.37 };

  double LPM_effect_factor[5] = { 0.8689, 0.7719, 0.7183, 0.712, 0.796 };

  for (int n=0; n<5; n++) {
    Veff_we_new_Tdomain[n] *= LPM_effect_factor[n];
  }

  TGraph *gVeff_new_Tdomain = new TGraph( 5, E_Veff_we_new, Veff_we_new_Tdomain );



  double E_Veff_we_old_Fdomain[5] = { 16, 17, 18, 19, 20 };

  double Veff_we_old_Fdomain[5] = { 5.22266, 60.3621, 257.911, 584.753, 908.165 };

  TGraph *gVeff_old_Fdomain = new TGraph( 5, E_Veff_we_old_Fdomain, Veff_we_old_Fdomain );



  // old Tdomain, short range, no LPM accounted?
  double E_Veff_we_old_Tdomain[4] = { 17.5, 18, 19, 20 };

  double Veff_we_old_Tdomain[4] = { 32.8477, 87.6905, 246.852, 375.004 };

  TGraph *gVeff_old_Tdomain = new TGraph( 4, E_Veff_we_old_Tdomain, Veff_we_old_Tdomain );


  // Veff plots!
  TCanvas *cVeff_compare = new TCanvas("cVeff_compare","",200,10,800,600);
  //cVeff_compare->Divide(1,2);

  // cd 1 for Fdomain
  cVeff_compare->cd();
  cVeff_compare->cd()->SetLogy();
  cVeff_compare->cd()->SetGrid();

  gVeff_new_Tdomain->SetTitle("Effective volume for ARA37 (km^{3} sr)");
  gVeff_new_Tdomain->GetHistogram()->SetXTitle("log E");
  gVeff_new_Tdomain->GetHistogram()->SetYTitle("km^3 sr");

  gVeff_new_Tdomain->GetHistogram()->SetMaximum(1.e4);
  gVeff_new_Tdomain->GetHistogram()->SetMinimum(1.e0);
  //
  //gVeff_ARA37[domain]->SetLineStyle(8);
  gVeff_new_Tdomain->SetLineWidth(2);
  gVeff_new_Tdomain->Draw("al");

  gVeff_new_Fdomain->SetLineColor(kBlue);
  //gVeff_ARA1[domain]->SetLineStyle(2);
  gVeff_new_Fdomain->SetLineWidth(2);
  gVeff_new_Fdomain->Draw("l");

  gVeff_old_Tdomain->SetLineColor(kRed);
  //gVeff_TestBed[domain]->SetLineStyle(2);
  gVeff_old_Tdomain->SetLineWidth(2);
  gVeff_old_Tdomain->Draw("l");

  gVeff_old_Fdomain->SetLineColor(kGreen+2);
  //gVeff_TestBed[domain]->SetLineStyle(2);
  gVeff_old_Fdomain->SetLineWidth(2);
  gVeff_old_Fdomain->Draw("l");


  TLegend *Leg_Veff_1 = new TLegend(0.95, 0.15, 0.5,0.45);
  Leg_Veff_1-> AddEntry(gVeff_new_Tdomain, "newRF, after fixes", "l");
  Leg_Veff_1-> AddEntry(gVeff_new_Fdomain, "oldRF, after fixes", "l");
  Leg_Veff_1-> AddEntry(gVeff_old_Tdomain, "newRF, before fixes", "l");
  Leg_Veff_1-> AddEntry(gVeff_old_Fdomain, "oldRF, before fixes", "l");
  //Leg_Veff_1-> AddEntry(gVeff_TB_AS_100m, "TestBed AraSim phase 100m deep", "l");
  Leg_Veff_1 -> Draw();

  cVeff_compare->Print("test_Veff_compare.pdf");

  delete cVeff_compare;


































  double ARA3_x[5] = { 16, 17, 18, 19, 20};

  double ARA3_y[5] = { -14.09+offset, -15.60+offset, -16.56+offset, -17.20+offset, -17.67+offset };

  LogToLine(5, ARA3_x);
  LogToLine(5, ARA3_y);

  TGraph *g_ARA3 = new TGraph( 5, ARA3_x, ARA3_y );

  double ARA3_1yr_x[5] = { 16, 17, 18, 19, 20};

  double ARA3_1yr_y[5] = { -13.61+offset, -15.12+offset, -16.09+offset, -16.72+offset, -17.20+offset };

  LogToLine(5, ARA3_1yr_x);
  LogToLine(5, ARA3_1yr_y);

  TGraph *g_ARA3_1yr = new TGraph( 5, ARA3_1yr_x, ARA3_1yr_y );





  double ARA6_2016_x[4] = { 17, 18, 19, 20};

  double ARA6_2016_y[4] = { -15.29, -16.331, -17.036, -17.562};

  LogToLine(4, ARA6_2016_x);
  LogToLine(4, ARA6_2016_y);


  TGraph *g_ARA6_2016 = new TGraph( 4, ARA6_2016_x, ARA6_2016_y );



  double ARA10_2017_x[4] = { 17, 18, 19, 20};

  double ARA10_2017_y[4] = { -15.57, -16.613, -17.321, -17.844};

  LogToLine(4, ARA10_2017_x);
  LogToLine(4, ARA10_2017_y);


  TGraph *g_ARA10_2017 = new TGraph( 4, ARA10_2017_x, ARA10_2017_y );





  // testbed with cwdB 10, peakcut 7.2 applied
  //
  //double offsetTB = 0.; // 0 offset is when factor 2.3 is applied
  //double offsetTB = -0.36; // 0 offset is when factor 2.3 is applied
  double offsetTB = 0.36; // 0 offset is when factor 2.3 is applied (can't remember I need to apply +? -? offset to TestBed)

  double TestBed_3yr_x[7] = { 17., 17.5, 18., 18.5, 19., 19.5, 20. };

  // with opt cuts (CW dB 10, peakcut 7.2) not done correctly!
  //double TestBed_3yr_y[7] = { -13.87+offsetTB, -14.44+offsetTB, -14.87+offsetTB, -15.18+offsetTB, -15.43+offsetTB, -15.55+offsetTB, -15.71+offsetTB };
  // without cuts only trigger
  double TestBed_3yr_y[7] = { -14.20+offsetTB, -14.71+offsetTB, -15.09+offsetTB, -15.39+offsetTB, -15.69+offsetTB, -15.87+offsetTB, -16.15+offsetTB };

  LogToLine(7, TestBed_3yr_x);
  LogToLine(7, TestBed_3yr_y);

  TGraph *g_TestBed_3yr = new TGraph( 7, TestBed_3yr_x, TestBed_3yr_y );



  // TestBed '11 from Davez
  //
  double TestBed_Dave_x[6] = { 16.6176, 16.8607, 17.1207, 18.2206, 19.2074, 20.9303};

  double TestBed_Dave_y[6] = { -12.485, -13.3952, -13.897, -15.2914, -15.9952, -17.4976};

  LogToLine(6, TestBed_Dave_x);
  LogToLine(6, TestBed_Dave_y);

  TGraph *g_TestBed_Dave11 = new TGraph( 6, TestBed_Dave_x, TestBed_Dave_y );





  // testbed trig lvl sensitivity from phase version AraSim (5 months, same as Dave above)
  //
  //double TestBed_phase_x[7] = { 17., 17.5, 18., 18.5, 19., 19.5, 20. };

  // below 5 months, 3yrs wrong!
  // for 5 months
  //double TestBed_phase_y[7] = { -12.91609888, -13.43843395, -13.93516574, -14.34885047, -14.73255535, -15.02403108, -15.27025796 };
  // for 3 yrs
  //double TestBed_phase_y[7] = { -13.76997084, -14.29230592, -14.7890377, -15.20272243, -15.58642732, -15.87790305, -16.12412993 };
  //
  //LogToLine(7, TestBed_phase_x);
  //LogToLine(7, TestBed_phase_y);

  //TGraph *g_TestBed_phase = new TGraph( 7, TestBed_phase_x, TestBed_phase_y );
  //
  // New TetBed value
  //
  double TestBed_phase_x[9] = { 17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5, 21. };

  // for 3 months
  //double TestBed_phase_y[9] = { -12.69145337, -13.21312467, -13.71018384, -14.12354119, -14.50724607, -14.7987218, -15.04494868, -15.25151774, -15.4176158};

  // for 5 months
  double TestBed_phase_y[9] = { -12.91364942, -13.43532072, -13.93237989, -14.34573723, -14.72944212, -15.02091785, -15.26714473, -15.47371379, -15.63981185};

  // for 3 yrs
  //double TestBed_phase_y[9] = { -13.77063462, -14.29230592, -14.78936509, -15.20272243, -15.58642732, -15.87790305, -16.12412993, -16.33069899, -16.49679705};

  LogToLine(9, TestBed_phase_x);
  LogToLine(9, TestBed_phase_y);

  TGraph *g_TestBed_phase = new TGraph( 9, TestBed_phase_x, TestBed_phase_y );






  double TestBed_limit_2012_75days_x[9] = { 17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5, 21. };

  // for 75 livetime
  //double TestBed_phase_y[9] = { -12.69145337, -13.21312467, -13.71018384, -14.12354119, -14.50724607, -14.7987218, -15.04494868, -15.25151774, -15.4176158};

  // for 75 livetime
  double TestBed_limit_2012_75days_y[9] = { -11.77881854, -12.43674036, -13.00393343, -13.4369922, -13.73268694, -13.97149584, -14.11804774, -14.20274014, -14.24987852 };

  // for 3 yrs
  //double TestBed_phase_y[9] = { -13.77063462, -14.29230592, -14.78936509, -15.20272243, -15.58642732, -15.87790305, -16.12412993, -16.33069899, -16.49679705};

  LogToLine(9, TestBed_limit_2012_75days_x);
  LogToLine(9, TestBed_limit_2012_75days_y);

  TGraph *g_TestBed_limit_2012_75days = new TGraph( 9, TestBed_limit_2012_75days_x, TestBed_limit_2012_75days_y );



  // trig level
  double TestBed_limit_2012_75days_trig_y[9] = { -12.60520791, -13.12702076, -13.62649855, -14.0374033, -14.42121004, -14.71238293, -14.95963247, -15.16458373, -15.33020796 };

  LogToLine(9, TestBed_limit_2012_75days_trig_y);

  TGraph *g_TestBed_limit_2012_75days_trig = new TGraph( 9, TestBed_limit_2012_75days_x, TestBed_limit_2012_75days_trig_y );





  double TestBed_limit_20112012_210days_x[9] = { 17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5, 21. };

  // for 210 livetime
  double TestBed_limit_20112012_210days_y[9] = { -12.132, -12.783, -13.375, -13.823, -14.126, -14.367, -14.515, -14.603, -14.668 };

  // for 3 yrs
  //double TestBed_phase_y[9] = { -13.77063462, -14.29230592, -14.78936509, -15.20272243, -15.58642732, -15.87790305, -16.12412993, -16.33069899, -16.49679705};

  LogToLine(9, TestBed_limit_20112012_210days_x);
  LogToLine(9, TestBed_limit_20112012_210days_y);

  TGraph *g_TestBed_limit_20112012_210days = new TGraph( 9, TestBed_limit_20112012_210days_x, TestBed_limit_20112012_210days_y );







  double TestBed_limit_20112012_75_210days_x[9] = { 17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5, 21. };

  // for 210 livetime
  double TestBed_limit_20112012_75_210days_y[9] = { -12.29, -12.944, -13.529, -13.9722, -14.275, -14.514, -14.661, -14.749, -14.809 };

  // for 3 yrs
  //double TestBed_phase_y[9] = { -13.77063462, -14.29230592, -14.78936509, -15.20272243, -15.58642732, -15.87790305, -16.12412993, -16.33069899, -16.49679705};

  LogToLine(9, TestBed_limit_20112012_75_210days_x);
  LogToLine(9, TestBed_limit_20112012_75_210days_y);

  TGraph *g_TestBed_limit_20112012_75_210days = new TGraph( 9, TestBed_limit_20112012_75_210days_x, TestBed_limit_20112012_75_210days_y );



  // for 224 livetime (0.86 factor need)
  // 021115 we changed the definition of the livetime based on the paper reviewer's comment
  // so we moved noisy time and good baseline cuts to analysis (instead of putting it in livetime cut)
  // under the new definition of the cuts, livetime become 382 days
  // and another new livetime definition (exclude calpulser timing cut), livetime become 415 days
  // even though livetime has changed, overall limit is unchanged as change in livetime will get compensated by new analysis efficiency
  //
  double TestBed_limit_20112012_224days_y[9];
  for (int bin=0; bin<9; bin++) {

    TestBed_limit_20112012_224days_y[bin] = TestBed_limit_20112012_75_210days_y[bin] * (1./0.86); // Aeff * T reduced by factor 0.86
  }//



  TGraph *g_TestBed_limit_20112012_224days = new TGraph( 9, TestBed_limit_20112012_75_210days_x, TestBed_limit_20112012_224days_y );



  // for 224 livetime (and fixed AraSim result for 2nd analysis only)
  // for 224 livetime (and fixed AraSim result for 2nd analysis only) and applied LPM effect factors
  //
  //double TestBed_limit_20112012_224days_debug_x[4] = { 17, 18, 19, 20 };
  //double TestBed_limit_20112012_224days_debug_x[5] = { 17, 18, 19, 20, 21 };
  double TestBed_limit_20112012_224days_debug_x[9] = { 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21 };

  //double TestBed_limit_20112012_224days_debug_y[4] = { -12.3187, -13.6781, -14.4922, -14.9693 };
  //double TestBed_limit_20112012_224days_debug_y[5] = { -12.35229767, -13.70330109, -14.50192003, -14.92969621, -15.26244302 };
  //double TestBed_limit_20112012_224days_debug_y[9] = { -12.29126747, -12.98935833, -13.59086213, -14.00967102, -14.35822589, -14.66114254, -14.86142219, -15.02522553, -15.16335609 };
  double TestBed_limit_20112012_224days_debug_y[9] = { -12.29126747, -12.98935833, -13.59086213, -14.00967102, -14.37362444, -14.68501972, -14.88814973, -15.04107829, -15.18006601 }; // after fix bugs and without missing peakrms files

  //double TestBed_limit_20112012_224days_debug_y[9] = { -12.29126747, -12.98935833, -13.59086213, -14.00967102, -14.35822589, -14.66114254, -14.78217621, -14.99295223, -15.10754106 }; // test modified LPM effect at high energies

  LogToLine(9, TestBed_limit_20112012_224days_debug_x);
  LogToLine(9, TestBed_limit_20112012_224days_debug_y);

  // apply str correction factor from 4pi to 2pi
  for (int bin=0; bin<9; bin++) {
    TestBed_limit_20112012_224days_debug_y[bin] *= str_offset;
  }



  //TGraph *g_TestBed_limit_20112012_224days_debug = new TGraph( 4, TestBed_limit_20112012_224days_debug_x, TestBed_limit_20112012_224days_debug_y );
  //TGraph *g_TestBed_limit_20112012_224days_debug = new TGraph( 5, TestBed_limit_20112012_224days_debug_x, TestBed_limit_20112012_224days_debug_y );
  TGraph *g_TestBed_limit_20112012_224days_debug = new TGraph( 9, TestBed_limit_20112012_224days_debug_x, TestBed_limit_20112012_224days_debug_y );


  double ARA2_limit_010616_x[10]={16., 16.5, 17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5   };
  double ARA2_limit_010616_y[10]={-10.72872, -12.36658, -13.35999, -14.16004, -14.80544, -15.32194, -15.79976, -16.16159, -16.47185, -16.71767 };

  LogToLine(10,ARA2_limit_010616_x);
  LogToLine(10,ARA2_limit_010616_y);

  TGraph *g_ARA2_limit_010616 = new TGraph( 10, ARA2_limit_010616_x, ARA2_limit_010616_y  );




  double ARIANNA_HRA3_limit_paper_x[8]={17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5 };
  double ARIANNA_HRA3_limit_paper_y[8]={-4.42842, -4.67521, -4.72756, -4.64530, -4.50321, -4.33120, -4.12179, -3.90491};

  for (int bin=0;bin<8;bin++) {
    ARIANNA_HRA3_limit_paper_y[bin]=ARIANNA_HRA3_limit_paper_y[bin]-(ARIANNA_HRA3_limit_paper_x[bin]-9.);
  }

  LogToLine(8,ARIANNA_HRA3_limit_paper_x);
  LogToLine(8,ARIANNA_HRA3_limit_paper_y);

  TGraph *g_ARIANNA_HRA3_limit_paper = new TGraph( 8, ARIANNA_HRA3_limit_paper_x, ARIANNA_HRA3_limit_paper_y  );










  // Kotera flux limit
  //
  // assume Kotera flux and use min factor f
  //
  double TestBed_2012_Optimized_KoteraFlux_x[9] = { 17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5, 21. };

  double TestBed_2012_Optimized_KoteraFlux_y[9] = { -8.0148, -7.4483, -7.1141, -7.0519, -7.1365, -7.2788, -7.611, -8.3405, -9.7552 }; // this is in GeV/cm2 sr s unit

  // change unit to /cm2 sr s
  for ( int i=0; i<9; i++) {

    TestBed_2012_Optimized_KoteraFlux_y[i] = TestBed_2012_Optimized_KoteraFlux_y[i] - (TestBed_2012_Optimized_KoteraFlux_x[i]-9.);
  }

  double factor_f = 907.; // got new Uexp, S_up files with new Clust_2, Grad cuts, also livetime : 75 days
  //double factor_f = 792.; // previous cut result (from prev Clust_2, Grad cut result, the only thing changed = livetime from 3 months to 75 days )

  double TestBed_limit_2012_75days_Kotera_y[9];
  for ( int i=0; i<9; i++) {
    TestBed_limit_2012_75days_Kotera_y[i] = TestBed_2012_Optimized_KoteraFlux_y[i] + log10( factor_f );
    //TestBed_limit_2012_75days_Kotera_y[i] = TestBed_2012_Optimized_KoteraFlux_y[i];
  }

  LogToLine(9, TestBed_2012_Optimized_KoteraFlux_x );
  LogToLine(9, TestBed_limit_2012_75days_Kotera_y );

  TGraph *g_TestBed_limit_2012_75days_Kotera = new TGraph( 9, TestBed_2012_Optimized_KoteraFlux_x, TestBed_limit_2012_75days_Kotera_y );






  // from 2015 paper
  const unsigned int nPoints=12;
  const double e[nPoints]={65161.306028,137109.233631,288498.544499,607044.529198,1277313.412688,2687660.419882,5655243.623724,11899487.080692,25038318.807265,52684406.012047,110855950.761230,233257670.521449};
  double elow[nPoints]={41983.333910,88339.278147,185879.189115,391117.899883,822971.158513,1731655.666862,3643665.172901,7666822.074546,16132152.086837,33944485.527694,71424326.496400,150287575.026061};
  double ehigh[nPoints]={88339.278147,185879.189115,391117.899883,822971.158513,1731655.666862,3643665.172901,7666822.074546,16132152.086837,33944485.527694,71424326.496400,150287575.026061,316227766.016838};
  double flux[nPoints]={2.041836,1.220496,1.495319,0.000000,0.798685,0.974741,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000};
  double fluxlow[nPoints]={0.909091,0.606061,0.909091,0.000000,0.303030,0.303030,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000};
  double fluxhigh[nPoints]={3.636364,2.121212,2.121212,0.303030,1.818182,2.121212,0.303030,0.909091,1.212121,1.818182,3.333333,5.757576};







  // IceCube 3 years HESE limit
  //
  double only_mu_to_all = 3.; // only mu flavor to all flavors factor
  //double only_mu_to_all = 1.; // leave as it is
  //


  double IC3yrHESE_lowarrow = 2.;


  // below are pure data point values from Albrecht
  vector<double> IC3yrHESE_x;
  "entering this part.\n";
  vector<double> IC3yrHESE_xl;
  vector<double> IC3yrHESE_xh;
  vector<double> IC3yrHESE_y;
  vector<double> IC3yrHESE_yl;

  vector<double> IC3yrHESE_yh;
  vector<double> IC3yrHESE_y0;

  int whichpub=1;  //0 for IC 3 yr paper, 1 for 2015


  const int IC3yrHESE_bin = 12;

  if (whichpub==0) {
    IC3yrHESE_x.push_back(65161.30603);
    IC3yrHESE_x.push_back(  137109.2336);
    IC3yrHESE_x.push_back(       288498.5445);
    IC3yrHESE_x.push_back(      607044.5292);
    IC3yrHESE_x.push_back(     1277313.413);
    IC3yrHESE_x.push_back(     2687660.42);
    IC3yrHESE_x.push_back(      5655243.624);
    IC3yrHESE_x.push_back(      11899487.08);
    IC3yrHESE_x.push_back(    25038318.81);
    IC3yrHESE_x.push_back(   52684406.01);
    IC3yrHESE_x.push_back(     110855950.8);
    IC3yrHESE_x.push_back(    233257670.5);


    IC3yrHESE_xl.push_back(41983.33391);
    IC3yrHESE_xl.push_back( 88339.27815);
    IC3yrHESE_xl.push_back(      185879.1891);
    IC3yrHESE_xl.push_back(  391117.8999);
    IC3yrHESE_xl.push_back(   822971.1585);
    IC3yrHESE_xl.push_back(   1731655.667);
    IC3yrHESE_xl.push_back(  3643665.173);
    IC3yrHESE_xl.push_back( 7666822.075);
    IC3yrHESE_xl.push_back( 16132152.09);
    IC3yrHESE_xl.push_back( 33944485.53);
    IC3yrHESE_xl.push_back(  71424326.5);
    IC3yrHESE_xl.push_back(  150287575);


    IC3yrHESE_xh.push_back(88339.27815);
    IC3yrHESE_xh.push_back(185879.1891);
    IC3yrHESE_xh.push_back(       391117.8999);
    IC3yrHESE_xh.push_back(       822971.1585);
    IC3yrHESE_xh.push_back(       1731655.667);
    IC3yrHESE_xh.push_back(   3643665.173);
    IC3yrHESE_xh.push_back(     7666822.075);
    IC3yrHESE_xh.push_back(    16132152.09);
    IC3yrHESE_xh.push_back(     33944485.53);
    IC3yrHESE_xh.push_back(    71424326.5);
    IC3yrHESE_xh.push_back(    150287575);
    IC3yrHESE_xh.push_back(  316227766 );


    IC3yrHESE_y.push_back(2.041751);
    IC3yrHESE_y.push_back(1.22059);
    IC3yrHESE_y.push_back(      1.49531);
    IC3yrHESE_y.push_back(    0.336603);
    IC3yrHESE_y.push_back(      0.798582);
    IC3yrHESE_y.push_back(    0.974725);
    IC3yrHESE_y.push_back(      0.246002);
    IC3yrHESE_y.push_back(     0.879309);
    IC3yrHESE_y.push_back(      1.360814);
    IC3yrHESE_y.push_back(      1.866319);
    IC3yrHESE_y.push_back(    3.285333);
    IC3yrHESE_y.push_back(   5.622056 );



    /*
      double IC3yrHESE_yl[IC3yrHESE_bin] = { 0.776708,
      0.560106,
      0.893709,
      0.336603 - IC3yrHESE_lowarrow,
      0.222002,
      0.310203,
      0.246002 - IC3yrHESE_lowarrow,
      0.879309 - IC3yrHESE_lowarrow,
      1.360814 - IC3yrHESE_lowarrow,
      1.866319 - IC3yrHESE_lowarrow,
      3.285333 - IC3yrHESE_lowarrow,
      5.622056 - IC3yrHESE_lowarrow };
    */

    IC3yrHESE_yl.push_back(0.776708);
    IC3yrHESE_yl.push_back(   0.560106);
    IC3yrHESE_yl.push_back(     0.893709);
    IC3yrHESE_yl.push_back(   0.);
    IC3yrHESE_yl.push_back(   0.222002);
    IC3yrHESE_yl.push_back(    0.310203);
    IC3yrHESE_yl.push_back(    0.);
    IC3yrHESE_yl.push_back(    0.);
    IC3yrHESE_yl.push_back(   0.);
    IC3yrHESE_yl.push_back(   0.);
    IC3yrHESE_yl.push_back(   0.);
    IC3yrHESE_yl.push_back(  0. );


    IC3yrHESE_yh.push_back(3.729637);
    IC3yrHESE_yh.push_back(2.074821);
    IC3yrHESE_yh.push_back(     2.265923);
    IC3yrHESE_yh.push_back(   0.336603);
    IC3yrHESE_yh.push_back( 1.696517);
    IC3yrHESE_yh.push_back( 2.148021);
    IC3yrHESE_yh.push_back( 0.246002);
    IC3yrHESE_yh.push_back( 0.879309);
    IC3yrHESE_yh.push_back( 1.360814);
    IC3yrHESE_yh.push_back(1.866319);
    IC3yrHESE_yh.push_back( 3.285333);
    IC3yrHESE_yh.push_back( 5.622056 );


    for (int bin=0;bin<IC3yrHESE_bin;bin++) {
      IC3yrHESE_y0.push_back(0.);
    }
  }
  else if (whichpub==1) {
    for (int bin=0;bin<nPoints;bin++) {
      IC3yrHESE_x.push_back(e[bin]);
      IC3yrHESE_xl.push_back(elow[bin]);
      IC3yrHESE_xh.push_back(ehigh[bin]);
      IC3yrHESE_y.push_back(flux[bin]);
      IC3yrHESE_yl.push_back(fluxlow[bin]);
      IC3yrHESE_yh.push_back(fluxhigh[bin]);
      IC3yrHESE_y0.push_back(0.);
    }

  }



  int yhl_count = 0;
  vector <double> IC3yrHESE_hl_x;
  vector <double> IC3yrHESE_hl_y;
  vector <double> IC3yrHESE_hl_yh;
  vector <double> IC3yrHESE_hl_yl;
  vector <double> IC3yrHESE_hl_0;

  int yh_count = 0;
  vector <double> IC3yrHESE_h_x;
  vector <double> IC3yrHESE_h_y;
  vector <double> IC3yrHESE_h_yl;
  vector <double> IC3yrHESE_h_0;
  double tmp=0.;

  // now calculate the correct values
  //
  int N=0;
  if (whichpub==0)
    N=IC3yrHESE_bin;
  else if (whichpub==1)
    N=nPoints;

  for (int bin=0; bin<N; bin++) {

    // cout<<"before, IC 3yr HESE at "<<IC3yrHESE_x[bin]<<", EF : "<<IC3yrHESE_y[bin]<<endl;

    // cout << "tmp is " << tmp << "\n";
    IC3yrHESE_x.at(bin)= log10(IC3yrHESE_x[bin]) + 9.; // in eV, log

    IC3yrHESE_xl.at(bin) = log10(IC3yrHESE_xl[bin]) + 9.; // in eV, log

    IC3yrHESE_xh.at(bin) = log10(IC3yrHESE_xh[bin]) + 9. ; // in eV, log
    // cout<<"after, IC 3yr HESE at "<<IC3yrHESE_x[bin]<<", EF : "<<IC3yrHESE_y[bin]<<endl;


    // to linear
    IC3yrHESE_x.at(bin) = pow(10., IC3yrHESE_x[bin]);
    IC3yrHESE_xl.at(bin) = pow(10., IC3yrHESE_xl[bin]);
    IC3yrHESE_xh.at(bin) = pow(10., IC3yrHESE_xh[bin]);



    // get difference
    IC3yrHESE_xl.at(bin) = IC3yrHESE_x[bin] - IC3yrHESE_xl[bin];
    IC3yrHESE_xh.at(bin) = IC3yrHESE_xh[bin] - IC3yrHESE_x[bin];

    if (whichpub==1) {
      if (IC3yrHESE_y[bin]==0)
	IC3yrHESE_y.at(bin)=IC3yrHESE_yh[bin];

    }
    // calculate proper flux limit value (from x10^-8 in E^2F [GeV cm-2 s-1 sr-1] to EF)
    IC3yrHESE_y.at(bin) = IC3yrHESE_y[bin]*1.e-8; // now in linear E^2F [GeV cm-2 s-1 sr-1]
    // cout << "I'm here1.\n";
    IC3yrHESE_yl.at(bin) = IC3yrHESE_yl[bin]*1.e-8; // now in linear E^2F [GeV cm-2 s-1 sr-1]
    // cout << "I'm here2.\n";
    // cout << "yh is " <<  IC3yrHESE_yh[bin] << "\n";
    IC3yrHESE_yh.at(bin) = IC3yrHESE_yh[bin]*1.e-8; // now in linear E^2F [GeV cm-2 s-1 sr-1]
    // cout<<"IC 3yr HESE at "<<IC3yrHESE_x[bin]<<", E^2F : "<<IC3yrHESE_y[bin]<<endl;



    IC3yrHESE_y.at(bin) = IC3yrHESE_y[bin] / (IC3yrHESE_x[bin]/1.e9); // now in linear EF [cm-2 s-1 sr-1]
    IC3yrHESE_yl.at(bin) = IC3yrHESE_yl[bin] / (IC3yrHESE_x[bin]/1.e9); // now in linear EF [cm-2 s-1 sr-1]
    IC3yrHESE_yh.at(bin) = IC3yrHESE_yh[bin] / (IC3yrHESE_x[bin]/1.e9); // now in linear EF [cm-2 s-1 sr-1]


    // get difference
    if ( IC3yrHESE_yl[bin]== 0. ) {
      IC3yrHESE_yl.at(bin) = IC3yrHESE_y[bin] - IC3yrHESE_y[bin]/IC3yrHESE_lowarrow;
    }
    else {
      IC3yrHESE_yl.at(bin) = IC3yrHESE_y[bin] - IC3yrHESE_yl[bin];
    }
    IC3yrHESE_yh.at(bin) = IC3yrHESE_yh[bin] - IC3yrHESE_y[bin];


    // apply flavor factor?
    IC3yrHESE_y.at(bin) = IC3yrHESE_y[bin]*only_mu_to_all; // flavor factor
    IC3yrHESE_yl.at(bin) = IC3yrHESE_yl[bin]*only_mu_to_all; // flavor factor
    IC3yrHESE_yh.at(bin) = IC3yrHESE_yh[bin]*only_mu_to_all; // flavor factor




    // find only H bins
    if ( IC3yrHESE_yl[bin]== IC3yrHESE_y[bin] - IC3yrHESE_y[bin]/IC3yrHESE_lowarrow ) {
      //if ( IC3yrHESE_yl[bin]==0) {

      yh_count ++;

      IC3yrHESE_h_x.push_back( IC3yrHESE_x[bin] );
      IC3yrHESE_h_y.push_back( IC3yrHESE_y[bin] );
      IC3yrHESE_h_yl.push_back( IC3yrHESE_yl[bin] );

      IC3yrHESE_h_0.push_back( 0. );
    }
    else {

      yhl_count ++;

      IC3yrHESE_hl_x.push_back( IC3yrHESE_x[bin] );
      IC3yrHESE_hl_y.push_back( IC3yrHESE_y[bin] );
      IC3yrHESE_hl_yh.push_back( IC3yrHESE_yh[bin] );
      IC3yrHESE_hl_yl.push_back( IC3yrHESE_yl[bin] );

      IC3yrHESE_hl_0.push_back( 0. );
    }


  }





  /*
    IC3yrHESE_x[bin] = log10(IC3yrHESE_x[bin]) + 9.; // in eV, log
    IC3yrHESE_xl[bin] = log10(IC3yrHESE_xl[bin]) + 9.; // in eV, log
    IC3yrHESE_xh[bin] = log10(IC3yrHESE_xh[bin]) + 9.; // in eV, log

    IC3yrHESE_y[bin] = IC3yrHESE_y[bin]*only_mu_to_all; // flavor factor
    IC3yrHESE_yl[bin] = IC3yrHESE_yl[bin]*only_mu_to_all; // flavor factor
    IC3yrHESE_yh[bin] = IC3yrHESE_yh[bin]*only_mu_to_all; // flavor factor
  */


  //
  TGraphAsymmErrors *g_IC_3yr_HESE_X = new TGraphAsymmErrors(IC3yrHESE_x.size(),&(IC3yrHESE_x[0]), &(IC3yrHESE_y[0]), &(IC3yrHESE_xl[0]), &(IC3yrHESE_xh[0]), &(IC3yrHESE_y0[0]), &(IC3yrHESE_y0[0]) ); // only x range bar

  TGraphAsymmErrors *g_IC_3yr_HESE_HL = new TGraphAsymmErrors(IC3yrHESE_hl_x.size(), &(IC3yrHESE_hl_x[0]), &(IC3yrHESE_hl_y[0]), &(IC3yrHESE_hl_0[0]), &(IC3yrHESE_hl_0[0]), &(IC3yrHESE_hl_yl[0]), &(IC3yrHESE_hl_yh[0]) );

  TGraphAsymmErrors *g_IC_3yr_HESE_H = new TGraphAsymmErrors(IC3yrHESE_h_x.size(), &(IC3yrHESE_h_x[0]), &(IC3yrHESE_h_y[0]), &(IC3yrHESE_h_0[0]), &(IC3yrHESE_h_0[0]), &(IC3yrHESE_h_yl[0]), &(IC3yrHESE_h_0[0]) );










  double EVA_3fly_x[5] = { 1.e17, 3.16e17, 1.e18, 3.16e18, 1.e19};

  double EVA_3fly_y[5] = { 1.66337E-12, 1.08234E-15, 3.35229E-17, 6.09518E-18, 1.56944E-18 };

  /*
    LogToLine(5, EVA_3fly_x);
    LogToLine(5, EVA_3fly_y);
  */

  TGraph *g_EVA_3fly = new TGraph( 5, EVA_3fly_x, EVA_3fly_y );

  //E [10^20 eV]
  double NuMoon2014_x[31] = {
    479.74,
    480.06,
    490.45,
    501.00,
    511.72,
    512.18,
    523.17,
    556.75,
    592.35,
    630.33,
    670.68,
    713.69,
    807.57,
    932.51,
    1076.71,
    1295.23,
    1622.61,
    2032.62,
    2707.08,
    3390.09,
    4606.79,
    6259.04,
    8858.11,
    11790.21,
    16345.61,
    21754.81,
    30781.09,
    40131.35,
    54508.11,
    68215.54,
    81965.86
  };



  //E^2 dN/dE/dA/dOmega/dt [GeV/cm^2/sr/s]
  double NuMoon2014_y[31] = {
    4.02E-05,
    3.35E-05,
    2.56E-05,
    2.03E-05,
    1.66E-05,
    1.29E-05,
    1.04E-05,
    8.07E-06,
    6.71E-06,
    5.31E-06,
    4.35E-06,
    3.44E-06,
    2.63E-06,
    2.12E-06,
    1.73E-06,
    1.37E-06,
    1.18E-06,
    1.03E-06,
    9.20E-07,
    8.75E-07,
    8.47E-07,
    8.62E-07,
    8.92E-07,
    9.39E-07,
    1.02E-06,
    1.09E-06,
    1.21E-06,
    1.36E-06,
    1.51E-06,
    1.72E-06,
    1.87E-06
  };

  for (int ipoints=0;ipoints<31;ipoints++) {
    NuMoon2014_x[ipoints]*=1.e20;
    NuMoon2014_y[ipoints]=NuMoon2014_y[ipoints]/NuMoon2014_x[ipoints]*1.e9;
  }

  TGraph *g_NuMoon2014 = new TGraph( 31, NuMoon2014_x, NuMoon2014_y );




  // start of the plot I want


  //TCanvas *cConst = new TCanvas("cConst","A Simple Graph Example",200,10,1000,1000);
  //TCanvas *cConst_2 = new TCanvas("cConst_2","A Simple Graph Example",200,10,1000,1200);
  //TCanvas *cConst_2 = new TCanvas("cConst_2","A Simple Graph Example",200,10,1200,900); // wider
  TCanvas *cConst_2 = new TCanvas("cConst_2","A Simple Graph Example",200,10,1400,1400); // wider

  cConst_2->cd();
  cConst_2->SetLogy();
  cConst_2->SetLogx();
  //g_ESS_base->SetTitle("Sensitivity");
  g_Kotera_shade->SetTitle("");
  g_Kotera_shade->GetHistogram()->SetXTitle("E (eV)");
  g_Kotera_shade->GetHistogram()->SetYTitle("E dN/dE dA d#Omega dt (cm^{-2} sr ^{-1} s^{-1} )");
  //g_vhvh->GetHistogram()->SetMaximum(100);
  //g_vhvh->GetHistogram()->SetMinimum(0.1);
  //g_Kotera_shade->GetHistogram()->SetMaximum(1.e-11);
  //g_Kotera_shade->GetHistogram()->SetMaximum(3.e-12);
  g_Kotera_shade->GetHistogram()->SetMaximum(3.e-11);
  //g_Kotera_shade->GetHistogram()->SetMinimum(1.e-20);
  //g_Kotera_shade->GetHistogram()->SetMaximum(1.e-13);
  //g_Kotera_shade->GetHistogram()->SetMinimum(1.e-20);
  g_Kotera_shade->GetHistogram()->SetMinimum(1.e-21);
  //g_Kotera_shade->GetHistogram()->SetMinimum(1.e-21);
  //g_Kotera_shade->GetHistogram()->SetMinimum(1.e-22);

  g_Kotera_shade->GetXaxis()->SetLimits(3.2e14,1.e24);
  //g_Kotera_shade->GetXaxis()->SetLimits(9.9e14,1.e22); // zoom little bit
  g_Kotera_shade->GetHistogram()->SetTitleSize  ( 0.04,"X");
  g_Kotera_shade->GetHistogram()->SetLabelOffset( 0.006,"X");
  g_Kotera_shade->GetHistogram()->SetLabelSize( 0.04,"X");
  g_Kotera_shade->GetHistogram()->SetLabelSize( 0.04,"Y");
  g_Kotera_shade->GetHistogram()->SetLabelOffset( 0.007,"Y");
  g_Kotera_shade->GetHistogram()->SetTitleSize  ( 0.04,"Y");
  //g_Kotera_shade->GetHistogram()->SetTitleOffset( 2.0,"Y");
  g_Kotera_shade->GetHistogram()->SetTitleOffset( 1.8,"Y");
  g_Kotera_shade->GetHistogram()->SetTitleOffset( 1.3,"X");
  gPad->SetLeftMargin(0.15);
  //g_Kotera_shade->SetLineColor(39);
  //g_Kotera_shade->SetLineWidth(3);
  //g_Kotera_shade->Draw("al");
  // ESS shade
  //g_Kotera_shade->SetFillStyle(3001);
  g_Kotera_shade->SetFillStyle(1001);
  //g_Kotera_shade->SetFillStyle(3004);
  //g_Kotera_shade->SetFillColor(16);
  //g_Kotera_shade->SetFillColor(13);
  g_Kotera_shade->SetFillColor(15);
  g_Kotera_shade->SetLineColor(0);
  //g_Kotera_shade->SetLineColor(1);
  //g_Kotera_shade->SetLineStyle(3);
  //g_Kotera_shade->SetLineWidth(5);
  g_Kotera_shade->Draw("af");
  //  g_Kotera_shade->Draw("afaxis");
  gPad->RedrawAxis();

  // g_ARIANNA_1296_shade->SetFillStyle(1001);
  //  g_ARIANNA_1296_shade->SetFillColor(15);
  //  g_ARIANNA_1296_shade->SetLineColor(0);
  //  g_ARIANNA_1296_shade->Draw("f");


  /*
    g_ESS_shade_ver2->SetTitle("");
    g_ESS_shade_ver2->GetHistogram()->SetXTitle("E (eV)");
    g_ESS_shade_ver2->GetHistogram()->SetYTitle("E F(E) (cm^{-2} s^{-1} sr ^{-1} )");
    //g_vhvh->GetHistogram()->SetMaximum(100);
    //g_vhvh->GetHistogram()->SetMinimum(0.1);
    g_ESS_shade_ver2->GetHistogram()->SetMaximum(1.e-11);
    //g_ESS_shade->GetHistogram()->SetMinimum(1.e-20);
    //g_ESS_shade->GetHistogram()->SetMaximum(1.e-13);
    //g_ESS_shade->GetHistogram()->SetMinimum(1.e-21);
    g_ESS_shade_ver2->GetHistogram()->SetMinimum(1.e-22);

    g_ESS_shade_ver2->GetXaxis()->SetLimits(3.2e14,1.e22);
    //g_ESS_shade->GetXaxis()->SetLimits(9.9e14,1.e22); // zoom little bit
    g_ESS_shade_ver2->GetHistogram()->SetTitleSize  ( 0.04,"X");
    g_ESS_shade_ver2->GetHistogram()->SetLabelOffset( 0.006,"X");
    g_ESS_shade_ver2->GetHistogram()->SetLabelSize( 0.04,"X");
    g_ESS_shade_ver2->GetHistogram()->SetLabelSize( 0.04,"Y");
    g_ESS_shade_ver2->GetHistogram()->SetLabelOffset( 0.007,"Y");
    g_ESS_shade_ver2->GetHistogram()->SetTitleSize  ( 0.04,"Y");
    g_ESS_shade_ver2->GetHistogram()->SetTitleOffset( 2.0,"Y");
    //g_ESS_shade->SetLineColor(39);
    //g_ESS_shade->SetLineWidth(3);
    //g_ESS_shade->Draw("al");
    g_ESS_shade_ver2->SetFillStyle(3001);
    //g_ESS_shade->SetFillStyle(3004);
    //g_ESS_shade->SetFillColor(16);
    g_ESS_shade_ver2->SetFillColor(13);
    g_ESS_shade_ver2->SetLineColor(0);
    //g_ESS_shade->SetLineColor(1);
    //g_ESS_shade->SetLineStyle(3);
    //g_ESS_shade->SetLineWidth(5);
    //g_ESS_shade->Draw("f");
    g_ESS_shade_ver2->Draw("af");
  */






  /*
  //g_ESS_strong->SetLineStyle(7);
  //g_ESS_strong->SetLineColor(8);
  g_ESS_strong->SetLineColor(1);
  g_ESS_strong->SetLineWidth(1);
  g_ESS_strong->Draw("l");
  g_ESS_base->SetLineColor(1);
  g_ESS_base->SetLineWidth(1);
  g_ESS_base->Draw("l");
  */



  /*
    g_Ave07Femix->SetLineColor(28);
    //g_Ave07Femix->SetLineStyle(3);
    g_Ave07Femix->SetLineWidth(3);
    g_Ave07Femix->Draw("l");
  */


  //g_Ahlers->SetLineStyle(4);
  g_Ahlers->SetLineColor(43);
  g_Ahlers->SetLineWidth(3);
  g_Ahlers->Draw("l");



  /*
  // Kotera shade
  g_Kotera_shade->SetFillStyle(3002);
  //g_Kotera_shade->SetFillStyle(3005);
  //g_Kotera_shade->SetFillColor(16);
  //g_Kotera_shade->SetFillColor(13);
  g_Kotera_shade->SetFillColor(15);
  g_Kotera_shade->SetLineColor(0);
  g_Kotera_shade->Draw("f");
  */


  /*
  // Kotera shade_ver2
  g_Kotera_shade_ver2->SetFillStyle(3002);
  //g_Kotera_shade->SetFillStyle(3005);
  //g_Kotera_shade->SetFillColor(16);
  //g_Kotera_shade->SetFillColor(13);
  g_Kotera_shade_ver2->SetFillColor(15);
  g_Kotera_shade_ver2->SetLineColor(0);
  g_Kotera_shade_ver2->Draw("f");
  */


  /*
    g_Kotera_low->SetLineColor(1);
    g_Kotera_low->SetLineWidth(1);
    g_Kotera_low->Draw("l");
    g_Kotera_max->SetLineColor(1);
    g_Kotera_max->SetLineWidth(1);
    g_Kotera_max->Draw("l");
  */



  /*
  //g_Kotera_low->SetLineStyle(8);
  g_Kotera_low->SetLineColor(46);
  g_Kotera_low->SetLineWidth(3);
  g_Kotera_low->Draw("l");
  */

  /*
  //g_Kotera_mid->SetLineStyle(8);
  g_Kotera_mid->SetLineColor(34);
  g_Kotera_mid->SetLineWidth(3);
  g_Kotera_mid->Draw("l");
  */


  /*
  //g_Kotera_max->SetLineStyle(8);
  g_Kotera_max->SetLineColor(40);
  g_Kotera_max->SetLineWidth(3);
  g_Kotera_max->Draw("l");
  */


  /*
    g_GRB_BP_ALL->SetLineStyle(8);
    //g_GRB_BP_ALL->SetLineColor(8);
    g_GRB_BP_ALL->SetLineWidth(3);
    g_GRB_BP_ALL->Draw("l");
  */

  /*
  // GRB shade
  g_GRB_shade->SetFillStyle(3013);
  //g_Kotera_shade->SetFillStyle(3005);
  //g_Kotera_shade->SetFillColor(16);
  //g_Kotera_shade->SetFillColor(13);
  g_GRB_shade->SetFillColor(15);
  g_GRB_shade->SetLineColor(0);
  g_GRB_shade->Draw("f");
  */

  /*
    g_GRB_BATSE->SetLineStyle(3);
    g_GRB_BATSE->SetLineColor(kMagenta+2);
    g_GRB_BATSE->SetLineWidth(7);
    g_GRB_BATSE->Draw("l");
  */


  //g_GRB_IC_FC->Draw("l");

  //g_GRB_NFC->Draw("l");

  /*  wrong one!!!
  // GRB IC_FC and NFC shade
  //g_GRB_shade2->SetFillStyle(3013);
  g_GRB_shade2->SetFillStyle(3144);
  //g_Kotera_shade->SetFillStyle(3005);
  //g_GRB_shade2->SetFillColor(13);
  g_GRB_shade2->SetFillColor(kGreen+3);
  g_GRB_shade2->SetLineColor(0);
  g_GRB_shade2->Draw("f");
  */

  /*
  // GRB from Hummer
  //g_GRB_shade2->SetFillStyle(3013);
  g_GRB_Hummer_shade->SetFillStyle(3144);
  g_GRB_Hummer_shade->SetFillColor(kGreen+3);
  g_GRB_Hummer_shade->SetLineColor(0);
  g_GRB_Hummer_shade->Draw("f");
  */

  // GRB from Hummer AU
  /*
  //g_GRB_shade2->SetFillStyle(3013);
  g_GRB_Hummer_shade2->SetFillStyle(3144);
  g_GRB_Hummer_shade2->SetFillColor(kGreen+3);
  g_GRB_Hummer_shade2->SetLineColor(0);
  g_GRB_Hummer_shade2->Draw("f");
  */


  /*
  // GRB WIND shade
  //g_GRB_shade2->SetFillStyle(3013);
  //g_GRB_WIND_shade->SetFillStyle(3007);
  g_GRB_WIND_shade->SetFillStyle(3144);
  //g_Kotera_shade->SetFillStyle(3005);
  //g_GRB_WIND_shade->SetFillColor(13);
  g_GRB_WIND_shade->SetFillColor(kMagenta+2);
  //g_GRB_WIND_shade->SetFillColor(kRed);
  g_GRB_WIND_shade->SetLineColor(0);
  g_GRB_WIND_shade->Draw("f");


  // GRB ISM shade
  //g_GRB_shade2->SetFillStyle(3013);
  //g_GRB_ISM_shade->SetFillStyle(3006);
  g_GRB_ISM_shade->SetFillStyle(3144);
  //g_Kotera_shade->SetFillStyle(3005);
  //g_GRB_ISM_shade->SetFillColor(13);
  g_GRB_ISM_shade->SetFillColor(kMagenta-7);
  //g_GRB_ISM_shade->SetFillColor(kRed);
  g_GRB_ISM_shade->SetLineColor(0);
  g_GRB_ISM_shade->Draw("f");
  */

  /*
    g_GRB_WIND_MAX->SetLineStyle(3);
    g_GRB_WIND_MAX->SetLineColor(kMagenta-7);
    g_GRB_WIND_MAX->SetLineWidth(7);
    g_GRB_WIND_MAX->Draw("l");
  */



  /*
  // GRB ISM shade_ver2
  //g_GRB_shade2->SetFillStyle(3013);
  g_GRB_ISM_shade_ver2->SetFillStyle(3006);
  //g_Kotera_shade->SetFillStyle(3005);
  //g_GRB_ISM_shade->SetFillColor(13);
  g_GRB_ISM_shade_ver2->SetFillColor(kMagenta+2);
  g_GRB_ISM_shade_ver2->SetLineColor(0);
  g_GRB_ISM_shade_ver2->Draw("f");
  */


  /*
    g_GRB_WB99->SetLineStyle(3);
    g_GRB_WB99->SetLineColor(kMagenta+2);
    g_GRB_WB99->SetLineWidth(7);
    g_GRB_WB99->Draw("l");
  */


  /*
    g_GRB_WB_AG->SetLineStyle(2);
    //g_GRB_WB_AG->SetLineColor(kMagenta-7);
    g_GRB_WB_AG->SetLineColor(kRed+1);
    g_GRB_WB_AG->SetLineWidth(7);
    g_GRB_WB_AG->Draw("l");
  */



  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  //g_ARA37_shade->SetFillColor(kBlue+2);
  //g_ARA37_shade->SetFillColor(kRed+2);
  g_ARA37_shade->SetFillColor(kBlue+3);
  g_ARA37_shade->SetLineColor(0);
  g_ARA37_shade->Draw("f");
  */




  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_debug->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA37_shade_debug->SetFillColor(kBlue+2);
  g_ARA37_shade_debug->SetLineColor(0);
  g_ARA37_shade_debug->Draw("f");
  */



  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_debug_oldR->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA37_shade_debug_oldR->SetFillColor(kBlue+2);
  g_ARA37_shade_debug_oldR->SetLineColor(0);
  g_ARA37_shade_debug_oldR->Draw("f");
  */



  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_debug_NoPhase->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA37_shade_debug_NoPhase->SetFillColor(kBlue+2);
  g_ARA37_shade_debug_NoPhase->SetLineColor(0);
  g_ARA37_shade_debug_NoPhase->Draw("f");
  */



  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_debug_10degCut->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA37_shade_debug_10degCut->SetFillColor(kBlue+2);
  g_ARA37_shade_debug_10degCut->SetLineColor(0);
  g_ARA37_shade_debug_10degCut->Draw("f");
  */


  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_debug_RAYSOL10km->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA37_shade_debug_RAYSOL10km->SetFillColor(kBlue+2);
  g_ARA37_shade_debug_RAYSOL10km->SetLineColor(0);
  g_ARA37_shade_debug_RAYSOL10km->Draw("f");
  */



  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_debug_SameShowerL->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA37_shade_debug_SameShowerL->SetFillColor(kBlue+2);
  g_ARA37_shade_debug_SameShowerL->SetLineColor(0);
  g_ARA37_shade_debug_SameShowerL->Draw("f");
  */



  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_debug_Vsat->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA37_shade_debug_Vsat->SetFillColor(kBlue+2);
  g_ARA37_shade_debug_Vsat->SetLineColor(0);
  g_ARA37_shade_debug_Vsat->Draw("f");
  */


  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_debug_MakeArraySat->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA37_shade_debug_MakeArraySat->SetFillColor(kBlue+2);
  g_ARA37_shade_debug_MakeArraySat->SetLineColor(0);
  g_ARA37_shade_debug_MakeArraySat->Draw("f");
  */


  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_debug_MakeArraySat_wPhase->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA37_shade_debug_MakeArraySat_wPhase->SetFillColor(kBlue+2);
  //g_ARA37_shade_debug_MakeArraySat_wPhase->SetFillColor(kRed+2);
  g_ARA37_shade_debug_MakeArraySat_wPhase->SetLineColor(0);
  g_ARA37_shade_debug_MakeArraySat_wPhase->Draw("f");
  */

  /*
  //g_ARA37_shade_ExpAnalysis
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_ExpAnalysis->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA37_shade_ExpAnalysis->SetFillColor(kGreen+3);
  //g_ARA37_shade_ExpAnalysis->SetFillColor(kRed+2);
  g_ARA37_shade_ExpAnalysis->SetLineColor(0);
  g_ARA37_shade_ExpAnalysis->Draw("f");
  */




  // Thomas ARA2_3 limit band test
  /*
  // g_ARA2_3_Thomas_shade
  g_ARA2_3_Thomas_shade->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA2_3_Thomas_shade->SetFillColor(kViolet);
  //g_ARA2_3_Thomas_shade->SetFillColor(kRed+2);
  g_ARA2_3_Thomas_shade->SetLineColor(kViolet+2);
  //g_ARA2_3_Thomas_shade->Draw("f");


  // g_ARA2_3_Thomas_shade
  g_ARA2_3_Thomas_shade_20150430->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  g_ARA2_3_Thomas_shade_20150430->SetFillColor(kViolet);
  //g_ARA2_3_Thomas_shade->SetFillColor(kRed+2);
  g_ARA2_3_Thomas_shade_20150430->SetLineWidth(3);
  g_ARA2_3_Thomas_shade_20150430->SetLineColor(kViolet+2);
  //    g_ARA2_3_Thomas_shade_20150430->Draw("f");

  g_ARA2_3_Thomas->SetLineWidth(3);
  g_ARA2_3_Thomas->SetLineColor(kViolet+2);
  //    g_ARA2_3_Thomas->Draw("l");

  g_ARA2_3_Thomas_20150430->SetLineColor(kViolet+2);
  g_ARA2_3_Thomas_20150430->SetLineWidth(3);
  //    g_ARA2_3_Thomas_20150430->Draw("l");
  */
  // points test
  // g_ARA2_3_Thomas_points
  // this moved to below to draw it over other experiments
  /*
    g_ARA2_3_Thomas_points->SetLineWidth(5);
    //g_TestBed_limit_2012_75days_Kotera->SetLineColor(kBlue);
    g_ARA2_3_Thomas_points->SetLineColor(kGreen+2);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerStyle(26);
    g_ARA2_3_Thomas_points->SetMarkerStyle(22);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerColor(12);
    g_ARA2_3_Thomas_points->SetMarkerColor(kGreen+2);
    g_ARA2_3_Thomas_points->SetMarkerSize(3);
    g_ARA2_3_Thomas_points->Draw("lp");
  */



  // g_ARA2_limit_010616->SetLineWidth(5);
  //  //g_ARA2_limit_010616->SetLineColor(kBlue);
  //  g_ARA2_limit_010616->SetLineColor(kGray+2);
  //  //g_ARA2_limit_010616->SetMarkerStyle(26);
  //  g_ARA2_limit_010616->SetMarkerStyle(20);
  //  //g_ARA2_limit_010616->SetMarkerColor(12);
  //  g_ARA2_limit_010616->SetMarkerColor(kGray+2);
  //  g_ARA2_limit_010616->SetMarkerSize(3);
  //  g_ARA2_limit_010616->Draw("lp");



  // g_ARIANNA_HRA3_limit_paper->SetLineWidth(5);
  //  //g_ARIANNA_HRA3_limit_paper->SetLineColor(kBlue);
  //  g_ARIANNA_HRA3_limit_paper->SetLineColor(kGray+2);
  //  //g_ARIANNA_HRA3_limit_paper->SetMarkerStyle(26);
  //  g_ARIANNA_HRA3_limit_paper->SetMarkerStyle(22);
  //  //g_ARIANNA_HRA3_limit_paper->SetMarkerColor(12);
  //  g_ARIANNA_HRA3_limit_paper->SetMarkerColor(kGray+2);
  //  g_ARIANNA_HRA3_limit_paper->SetMarkerSize(3);
  //  g_ARIANNA_HRA3_limit_paper->Draw("lp");




  /*
  //g_ARA37_shade->SetFillStyle(3001);
  //g_ARA37_shade_debug_MakeArraySat_woSecondary->SetFillStyle(1001);
  //g_ARA37_shade_debug_MakeArraySat_woSecondary->SetFillStyle(3001);
  //g_ARA37_shade->SetFillColor(kBlue);
  //g_ARA37_shade_debug_MakeArraySat_wPhase->SetFillColor(kBlue+2);
  //g_ARA37_shade_debug_MakeArraySat_woSecondary->SetFillColor(kGreen+2);
  g_ARA37_shade_debug_MakeArraySat_woSecondary->SetLineColor(1);
  //g_ARA37_shade_debug_MakeArraySat_woSecondary->Draw("f");
  g_ARA37_shade_debug_MakeArraySat_woSecondary->SetLineWidth(3);
  g_ARA37_shade_debug_MakeArraySat_woSecondary->Draw("l");
  */




  /*
    g_ARA37_Fdomain_debug_MakeArraySat_woSecondary->SetLineColor(kGreen+2);
    g_ARA37_Fdomain_debug_MakeArraySat_woSecondary->SetLineWidth(3);
    g_ARA37_Fdomain_debug_MakeArraySat_woSecondary->Draw("l");
  */




  /*
  //g_ARA37_shade->SetFillStyle(3001);
  g_ARA37_shade_debug_MakeArraySat_woPhase->SetFillStyle(1001);
  //g_ARA37_shade->SetFillColor(kBlue);
  //g_ARA37_shade_debug_MakeArraySat_woPhase->SetFillColor(kBlue+2);
  g_ARA37_shade_debug_MakeArraySat_woPhase->SetFillColor(kRed+2);
  g_ARA37_shade_debug_MakeArraySat_woPhase->SetLineColor(0);
  g_ARA37_shade_debug_MakeArraySat_woPhase->Draw("f");
  */









  /*
    g_ARA37_up->SetLineStyle(7);
    g_ARA37_up->SetLineWidth(5);
    g_ARA37_up->SetLineColor(kBlue+1);
    g_ARA37_up->SetMarkerStyle(26);
    g_ARA37_up->SetMarkerSize(3);
    g_ARA37_up->SetMarkerColor(kBlue+1);
    g_ARA37_up->Draw("lp");
  */




  /*
  //g_Icecube->SetLineColor(9);
  g_IC40_Hummer->SetLineWidth(3);
  //g_IC40_Hummer->SetMarkerStyle(21);
  //g_IC40_Hummer->SetMarkerStyle(25);
  g_IC40_Hummer->SetMarkerStyle(20);
  g_IC40_Hummer->SetMarkerSize(2.2);
  //g_Icecube->SetLineStyle(3);
  g_IC40_Hummer->Draw("lp");
  */





  /*
    g_ARA37_Peter -> SetLineColor(kBlue);
    g_ARA37_Peter -> SetLineStyle(8);
    g_ARA37_Peter->SetLineWidth(5);
    g_ARA37_Peter->SetMarkerStyle(23);
    //g_ARA37_Peter->SetMarkerStyle(26);
    g_ARA37_Peter->SetMarkerColor(kBlue);
    g_ARA37_Peter->SetMarkerSize(3);
    g_ARA37_Peter->Draw("lp");
  */


  /*
    g_ARA3->SetLineStyle(7);
    g_ARA3->SetLineWidth(5);
    //g_ARA3->SetLineColor(12);
    g_ARA3->SetMarkerStyle(26);
    g_ARA3->SetMarkerSize(3);
    g_ARA3->Draw("lp");
  */


  /*
    g_ARA3_1yr->SetLineStyle(9);
    //g_ARA3_1yr->SetLineStyle(3);
    g_ARA3_1yr->SetLineWidth(5);
    g_ARA3_1yr->SetLineColor(12);
    g_ARA3_1yr->SetMarkerStyle(26);
    g_ARA3_1yr->SetMarkerColor(12);
    g_ARA3_1yr->SetMarkerSize(3);
    g_ARA3_1yr->Draw("lp");
  */


  /*
    g_ARA6_2016->SetLineStyle(7);
    g_ARA6_2016->SetLineWidth(5);
    g_ARA6_2016->SetLineColor(4);
    g_ARA6_2016->SetMarkerStyle(26);
    g_ARA6_2016->SetMarkerSize(3);
    g_ARA6_2016->Draw("lp");
  */



  /*
    g_ARA10_2017->SetLineStyle(9);
    //g_ARA10_2017->SetLineStyle(3);
    g_ARA10_2017->SetLineWidth(5);
    g_ARA10_2017->SetLineColor(kMagenta+3);
    g_ARA10_2017->SetMarkerStyle(26);
    g_ARA10_2017->SetMarkerColor(12);
    g_ARA10_2017->SetMarkerSize(3);
    g_ARA10_2017->Draw("lp");
  */





  /*
  //g_TestBed_3yr->SetLineStyle(9);
  g_TestBed_3yr->SetLineWidth(4);
  g_TestBed_3yr->SetLineColor(kRed);
  //g_TestBed_3yr->SetMarkerStyle(26);
  g_TestBed_3yr->SetMarkerStyle(23);
  g_TestBed_3yr->SetMarkerColor(kBlue);
  g_TestBed_3yr->SetMarkerSize(3);
  g_TestBed_3yr->Draw("lp");
  */



  /*
  //g_ANITA->SetLineColor(36);
  g_ANITA->SetLineWidth(3);
  g_ANITA->SetMarkerStyle(24);
  g_ANITA->SetMarkerSize(2.2);
  //g_ANITA->SetLineStyle(10);
  g_ANITA->SetLineColor(kGray+1);
  g_ANITA->SetMarkerColor(kGray+1);
  g_ANITA->Draw("lp");
  */





  g_NuMoon2014->SetLineWidth(3);
  g_NuMoon2014->SetMarkerStyle(21);
  g_NuMoon2014->SetMarkerSize(0.);
  g_NuMoon2014->SetMarkerColor(kOrange+2);
  //g_NuMoon2014->SetLineStyle(5);
  g_NuMoon2014->SetLineColor(kOrange+2);
  g_NuMoon2014->Draw("lp");



  /*
  //g_Auger->SetLineColor(14);
  g_Auger->SetLineWidth(3);
  g_Auger->SetMarkerStyle(34);
  g_Auger->SetMarkerSize(2.8);
  //g_Auger->SetLineStyle(9);
  g_Auger->SetLineColor(kGray+1);
  g_Auger->SetMarkerColor(kGray+1);
  g_Auger->Draw("lp");

  //g_Auger->SetLineColor(14);
  g_Auger11->SetLineWidth(3);
  g_Auger11->SetMarkerStyle(34);
  g_Auger11->SetMarkerSize(2.8);
  //g_Auger->SetLineStyle(9);
  g_Auger11->SetLineColor(kGray+1);
  g_Auger11->SetMarkerColor(kGray+1);
  g_Auger11->Draw("lp");
  */




  /*
  //g_Icecube->SetLineColor(9);
  g_Icecube->SetLineWidth(3);
  g_Icecube->SetMarkerStyle(21);
  g_Icecube->SetMarkerSize(2.2);
  //g_Icecube->SetLineStyle(3);
  g_Icecube->Draw("lp");
  */

  


  //g_IC_new->SetLineColor(9);
  g_IC_new->SetLineWidth(3);
  g_IC_new->SetMarkerStyle(21);
  g_IC_new->SetMarkerSize(2.2);
  //g_IC_new->SetLineStyle(3);
  g_IC_new->SetLineColor(kGreen+2);
  g_IC_new->SetMarkerColor(kGreen+2);
  // g_IC_new->Draw("lp");



  g_RICE->SetLineWidth(3);
  g_RICE->SetMarkerStyle(33);
  g_RICE->SetMarkerSize(3.5);
  //g_Icecube->SetLineStyle(3);
  g_RICE->SetLineColor(kBlue+2);
  g_RICE->SetMarkerColor(kBlue+2);
  // g_RICE->Draw("lp");






  /*
  // updated Auger13
  g_Auger13->SetLineWidth(3);
  g_Auger13->SetMarkerStyle(34);
  g_Auger13->SetMarkerSize(2.8);
  g_Auger13->SetLineColor(kViolet+1);
  g_Auger13->SetMarkerColor(kViolet+1);
  g_Auger13->Draw("lp");
  g_Auger13->Draw("p"); // for error bar case
  */


  // updated Auger15
  g_Auger15->SetLineWidth(3);
  g_Auger15->SetMarkerStyle(34);
  g_Auger15->SetMarkerSize(2.8);
  g_Auger15->SetLineColor(kViolet+2);
  g_Auger15->SetMarkerColor(kViolet+2);
  g_Auger15->Draw("lp");
  g_Auger15->Draw("p"); // for error bar case



  //g_ANITA->SetLineColor(36);
  g_ANITA_erratum->SetLineWidth(3);
  g_ANITA_erratum->SetMarkerStyle(20);
  g_ANITA_erratum->SetMarkerSize(2.2);
  //g_ANITA->SetLineStyle(10);
  g_ANITA_erratum->SetLineColor(kRed+2);
  g_ANITA_erratum->SetMarkerColor(kRed+2);
  g_ANITA_erratum->Draw("lp");



  //  g_ANITA->SetLineColor(36);
  g_ANITA_4->SetLineColor(kBlue);
  g_ANITA_4->SetLineWidth(3);
  g_ANITA_4->SetLineStyle(2);
  g_ANITA_4->Draw("l");
  //g_ANITA_4_shade->Draw("l");


  g_ANITA_all->SetLineColor(kMagenta);
  g_ANITA_all->SetLineWidth(3);
  g_ANITA_all->SetLineStyle(2);
  g_ANITA_all->Draw("l");
  
  // g_ANITA_4->SetLineWidth(3);
  // g_ANITA_4->SetMarkerStyle(8);
  // g_ANITA_4->SetMarkerSize(2.2);
  // //g_ANITA->SetLineStyle(10);
  // g_ANITA_4->SetLineColor(kRed+2);
  // g_ANITA_4->SetMarkerColor(kRed+2);
  // g_ANITA_4->Draw("lp");
  



  /*
    g_IC_evts_HESE->SetLineColor(kGray+2);
    g_IC_evts_HESE->SetLineWidth(3);
    g_IC_evts_HESE->SetMarkerStyle(25);
    g_IC_evts_HESE->SetMarkerSize(2.2);
    //g_IC_evts_HESE->SetLineStyle(3);
    g_IC_evts_HESE->Draw("lp");
  */



  /*
  // add WB98 x 3/2 limit
  //
  g_WB98_3p2->SetLineColor(kGray+2);
  g_WB98_3p2->SetLineWidth(3);
  g_WB98_3p2->SetLineStyle(5);
  g_WB98_3p2->Draw("l");
  */






  /*
    g_EVA_3fly->SetLineWidth(5);
    g_EVA_3fly->SetMarkerStyle(8);
    g_EVA_3fly->SetMarkerSize(3);
    g_EVA_3fly->SetMarkerColor(4);
    g_EVA_3fly->SetLineStyle(2);
    g_EVA_3fly->SetLineColor(4);
    g_EVA_3fly->Draw("lp");
  */








  /*
    g_ARA37->SetLineStyle(7);
    g_ARA37->SetLineWidth(5);
    g_ARA37->SetMarkerStyle(22);
    //g_ARA37->SetMarkerStyle(26);
    //g_ARA37->SetMarkerColor(kRed);
    g_ARA37->SetLineColor(kBlue);
    g_ARA37->SetMarkerColor(kBlue);
    g_ARA37->SetMarkerSize(3);
    g_ARA37->Draw("lp");
  */



  /*
    g_TestBed_Dave11->SetLineStyle(9);
    //g_TestBed_Dave11->SetLineStyle(3);
    g_TestBed_Dave11->SetLineWidth(5);
    //g_TestBed_Dave11->SetLineColor(kBlue);
    g_TestBed_Dave11->SetLineColor(kRed);
    //g_TestBed_Dave11->SetMarkerStyle(26);
    g_TestBed_Dave11->SetMarkerStyle(22);
    //g_TestBed_Dave11->SetMarkerColor(12);
    g_TestBed_Dave11->SetMarkerColor(kRed);
    g_TestBed_Dave11->SetMarkerSize(3);
    g_TestBed_Dave11->Draw("lp");
  */



  /*
    g_TestBed_phase->SetLineStyle(9);
    //g_TestBed_phase->SetLineStyle(3);
    g_TestBed_phase->SetLineWidth(5);
    //g_TestBed_phase->SetLineColor(kBlue);
    g_TestBed_phase->SetLineColor(kBlue);
    //g_TestBed_phase->SetMarkerStyle(26);
    g_TestBed_phase->SetMarkerStyle(22);
    //g_TestBed_phase->SetMarkerColor(12);
    g_TestBed_phase->SetMarkerColor(kBlue);
    g_TestBed_phase->SetMarkerSize(3);
    g_TestBed_phase->Draw("lp");
  */


  /*
    g_TestBed_limit_2012_75days_trig->SetLineStyle(9);
    //g_TestBed_limit_2012_75days_trig->SetLineStyle(3);
    g_TestBed_limit_2012_75days_trig->SetLineWidth(5);
    //g_TestBed_limit_2012_75days_trig->SetLineColor(kBlue);
    g_TestBed_limit_2012_75days_trig->SetLineColor(kRed);
    //g_TestBed_limit_2012_75days_trig->SetMarkerStyle(26);
    g_TestBed_limit_2012_75days_trig->SetMarkerStyle(22);
    //g_TestBed_limit_2012_75days_trig->SetMarkerColor(12);
    g_TestBed_limit_2012_75days_trig->SetMarkerColor(kRed);
    g_TestBed_limit_2012_75days_trig->SetMarkerSize(3);
    g_TestBed_limit_2012_75days_trig->Draw("lp");
  */



  /*
    g_TestBed_limit_2012_75days->SetLineStyle(9);
    //g_TestBed_limit_2012_75days->SetLineStyle(3);
    g_TestBed_limit_2012_75days->SetLineWidth(5);
    //g_TestBed_limit_2012_75days->SetLineColor(kBlue);
    g_TestBed_limit_2012_75days->SetLineColor(kBlue);
    //g_TestBed_limit_2012_75days->SetMarkerStyle(26);
    g_TestBed_limit_2012_75days->SetMarkerStyle(21);
    //g_TestBed_limit_2012_75days->SetMarkerColor(12);
    g_TestBed_limit_2012_75days->SetMarkerColor(kBlue);
    g_TestBed_limit_2012_75days->SetMarkerSize(3);
    g_TestBed_limit_2012_75days->Draw("lp");
  */



  /*
    g_TestBed_limit_2012_75days_Kotera->SetLineStyle(9);
    //g_TestBed_limit_2012_75days_Kotera->SetLineStyle(3);
    g_TestBed_limit_2012_75days_Kotera->SetLineWidth(5);
    //g_TestBed_limit_2012_75days_Kotera->SetLineColor(kBlue);
    g_TestBed_limit_2012_75days_Kotera->SetLineColor(kGreen+2);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerStyle(26);
    g_TestBed_limit_2012_75days_Kotera->SetMarkerStyle(22);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerColor(12);
    g_TestBed_limit_2012_75days_Kotera->SetMarkerColor(kGreen+2);
    g_TestBed_limit_2012_75days_Kotera->SetMarkerSize(3);
    g_TestBed_limit_2012_75days_Kotera->Draw("lp");
  */





  /*
    g_TestBed_limit_20112012_210days->SetLineStyle(9);
    //g_TestBed_limit_2012_75days_Kotera->SetLineStyle(3);
    g_TestBed_limit_20112012_210days->SetLineWidth(5);
    //g_TestBed_limit_2012_75days_Kotera->SetLineColor(kBlue);
    g_TestBed_limit_20112012_210days->SetLineColor(kGreen+2);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerStyle(26);
    g_TestBed_limit_20112012_210days->SetMarkerStyle(22);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerColor(12);
    g_TestBed_limit_20112012_210days->SetMarkerColor(kGreen+2);
    g_TestBed_limit_20112012_210days->SetMarkerSize(3);
    g_TestBed_limit_20112012_210days->Draw("lp");
  */




  /*
  //g_TestBed_limit_20112012_75_210days->SetLineStyle(9);
  //g_TestBed_limit_2012_75days_Kotera->SetLineStyle(3);
  g_TestBed_limit_20112012_75_210days->SetLineWidth(5);
  //g_TestBed_limit_2012_75days_Kotera->SetLineColor(kBlue);
  g_TestBed_limit_20112012_75_210days->SetLineColor(kRed+2);
  //g_TestBed_limit_2012_75days_Kotera->SetMarkerStyle(26);
  g_TestBed_limit_20112012_75_210days->SetMarkerStyle(23);
  //g_TestBed_limit_2012_75days_Kotera->SetMarkerColor(12);
  g_TestBed_limit_20112012_75_210days->SetMarkerColor(kRed+2);
  g_TestBed_limit_20112012_75_210days->SetMarkerSize(3);
  g_TestBed_limit_20112012_75_210days->Draw("lp");
  */




  /*
    g_TestBed_limit_20112012_224days->SetLineWidth(5);
    //g_TestBed_limit_2012_75days_Kotera->SetLineColor(kBlue);
    //g_TestBed_limit_20112012_224days->SetLineColor(kRed+2);
    g_TestBed_limit_20112012_224days->SetLineColor(kBlue+3);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerStyle(26);
    g_TestBed_limit_20112012_224days->SetMarkerStyle(23);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerColor(12);
    g_TestBed_limit_20112012_224days->SetMarkerColor(kBlue+3);
    g_TestBed_limit_20112012_224days->SetMarkerSize(3);
    g_TestBed_limit_20112012_224days->Draw("lp");
  */



  /*
    g_TestBed_limit_20112012_224days_debug->SetLineWidth(5);
    //g_TestBed_limit_2012_75days_Kotera->SetLineColor(kBlue);
    g_TestBed_limit_20112012_224days_debug->SetLineColor(kRed+2);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerStyle(26);
    g_TestBed_limit_20112012_224days_debug->SetMarkerStyle(23);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerColor(12);
    g_TestBed_limit_20112012_224days_debug->SetMarkerColor(kRed+2);
    g_TestBed_limit_20112012_224days_debug->SetMarkerSize(3);
    g_TestBed_limit_20112012_224days_debug->Draw("lp");
  */







  /*
    g_ARA2_3_Thomas_points->SetLineWidth(5);
    //g_TestBed_limit_2012_75days_Kotera->SetLineColor(kBlue);
    g_ARA2_3_Thomas_points->SetLineColor(kGreen+2);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerStyle(26);
    g_ARA2_3_Thomas_points->SetMarkerStyle(22);
    //g_TestBed_limit_2012_75days_Kotera->SetMarkerColor(12);
    g_ARA2_3_Thomas_points->SetMarkerColor(kGreen+2);
    g_ARA2_3_Thomas_points->SetMarkerSize(3);
    g_ARA2_3_Thomas_points->Draw("lp");
  */
  /*
  // Thomas error bar
  // g_ARA2_3_Thomas_error
  g_ARA2_3_Thomas_error->SetLineWidth(3);
  g_ARA2_3_Thomas_error->SetMarkerStyle(1);
  g_ARA2_3_Thomas_error->SetMarkerSize(2);
  g_ARA2_3_Thomas_error->SetMarkerColor(kGreen+2);
  g_ARA2_3_Thomas_error->SetLineColor(kGreen+2);
  //g_ARA2_3_Thomas_error->Draw("p");

  //g_ARA2_3_Thomas_points_low->Draw("pl");
  */






  gStyle->SetEndErrorSize(7);

  g_IC_3yr_HESE_X->SetLineWidth(2);
  g_IC_3yr_HESE_X->SetMarkerStyle(1);
  g_IC_3yr_HESE_X->SetMarkerSize(2);
  //g_IC_3yr_HESE_X->SetLineColor(kGray+1);
  g_IC_3yr_HESE_X->Draw("p");

  g_IC_3yr_HESE_HL->SetLineWidth(2);
  g_IC_3yr_HESE_HL->SetMarkerStyle(1);
  g_IC_3yr_HESE_HL->SetMarkerSize(2);
  //g_IC_3yr_HESE_HL->SetLineColor(kGray+1);
  g_IC_3yr_HESE_HL->Draw("p");

  g_IC_3yr_HESE_H->SetLineWidth(2);
  g_IC_3yr_HESE_H->SetMarkerStyle(1);
  g_IC_3yr_HESE_H->SetMarkerSize(2.5);
  //g_IC_3yr_HESE_H->SetLineColor(kGray+1);
  //g_IC_3yr_HESE_H->Draw("p");
  g_IC_3yr_HESE_H->Draw(">");


  // add preliminary text
  //  Just add text in keynote!

  TText *tt = new TText();
  tt->SetTextAlign(5);
  tt->SetTextSize(0.025);

  //tt->SetTextAngle(45);
  //tt->DrawTextNDC(0.64,0.6,"Connolly \& Vieregg 2016");
  // tt->DrawTextNDC(0.28,0.163,"Connolly \& Vieregg 2016");


  /*
    TText *t1 = new TText();
    t1->SetTextFont(62);
    t1->SetTextColor(1);   // 4
    t1->SetTextAlpha(1);   // 4
    t1->SetTextAlign(5);
    t1->SetTextSize(0.06);
    t1->SetTextAngle(45);
    t1->DrawTextNDC(0.38,0.92," Preliminary");
  */



  TLegend *Leg_Const_2 = new TLegend(0.48, 0.6, 0.75, 0.89); // for NOT little bit zoom in plot // L, B, R,
  //Leg_Const_2 -> AddEntry(g_ESS_base, "ESS base", "l");
  //Leg_Const_2 -> AddEntry(g_ARA3, "ARA3 (3yrs)", "lp");
  //Leg_Const_2 -> AddEntry(g_ARA6_2016, "ARA6 (2016)", "lp");
  //Leg_Const_2 -> AddEntry(g_ARA10_2017, "ARA10 (2017)", "lp");
  //Leg_Const_2 -> AddEntry(g_ARA3_1yr, "ARA3 (1yr) AraSim", "lp");
  //Leg_Const_2 -> AddEntry(g_TestBed_3yr, "TestBed (3yr) AraSim", "lp");
  //Leg_Const_2 -> AddEntry(g_ARA2_limit_010616, "ARA2 (10 x 0.75 mo)", "lp");
  //Leg_Const_2 -> AddEntry(g_ARIANNA_HRA3_limit_paper, "ARIANNA HRA-3 (3 x 0.58 yrs)", "lp");
  //Leg_Const_2 -> AddEntry(g_ARA37_Peter, "ARA37 (3yrs) UH", "lp");
  //Leg_Const_2 -> AddEntry(g_EVA_3fly, "EVA (3flights)", "lp");
  // Leg_Const_2 -> AddEntry(g_IC_new, "IceCube '12 (333.5 days)", "lp");
  Leg_Const_2 -> AddEntry(g_IC_3yr_HESE_X, "IceCube '15 (659.5 days)", "l");
  // Leg_Const_2 -> AddEntry(g_RICE, "RICE '11 (4.8 years)", "lp");
  Leg_Const_2 -> AddEntry(g_Auger15, "Auger '15 (6.5 yrs)", "lp");
  //Leg_Const_2 -> AddEntry(g_ANITA, "ANITA II", "lp");
  Leg_Const_2 -> AddEntry(g_ANITA_erratum, "ANITA 2 '10 (28.5 days)", "lp");
  Leg_Const_2 -> AddEntry(g_NuMoon2014, "NuMoon '10 (47.6 hrs) ", "lp");
  Leg_Const_2 -> AddEntry(g_ANITA_4, "ANITA-4 sensitivity (27.3 days)", "lp");
  Leg_Const_2 -> AddEntry(g_ANITA_all, "ANITA 1-4 sensitivity (99.3 days)", "l");//27.3 days)", "lf");

  //Leg_Const_2_2 -> AddEntry(g_Auger, "Auger '09", "lp");
  //Leg_Const_2 -> AddEntry(g_Auger11, "Auger '11", "lp");
  //Leg_Const_2 -> AddEntry(g_Auger13, "Auger '13", "lp");
  //Leg_Const_2 -> AddEntry(g_Auger13, "#splitline{Auger '13 (3.5yrs)}{Earth-skim}", "lp");
  //Leg_Const_2 -> AddEntry(g_Auger13, "Auger '13 (3.5yrs)", "lp");
  //Leg_Const_2 -> AddEntry(g_Auger13, "Auger '13 (3.5 yrs)", "lp");
  //Leg_Const_2_2 -> AddEntry(g_Icecube, "IceCube40", "lp");
  //Leg_Const_2 -> AddEntry(g_IC_new, "IceCube '12 (2yrs)", "lp");

  //Leg_Const_2 -> AddEntry(g_IC_3yr_HESE_X, "IceCube 3yrs HESE", "l");
  //Leg_Const_2 -> AddEntry(g_IC_3yr_HESE_X, "IceCube (HESE 3 yrs)", "l");

  //Leg_Const_2 -> AddEntry(g_IC_evts_HESE, "IC79+86 Starting Events", "lp");
  //Leg_Const_2_2 -> AddEntry(g_IC40_Hummer, "IC40", "lp");






  //Leg_Const_2 -> AddEntry(g_TestBed_limit_20112012_224days, "#splitline{ARA TestBed arXiv}{2011-2012 (224 days)}", "lp");
  //Leg_Const_2 -> AddEntry(g_TestBed_limit_20112012_224days_debug, "#splitline{ARA TestBed new}{2011-2012 (224 days)}", "lp");
  //Leg_Const_2 -> AddEntry(g_TestBed_limit_20112012_224days, "ARA TestBed arXiv", "lp");
  //Leg_Const_2 -> AddEntry(g_TestBed_limit_20112012_224days_debug, "ARA TestBed new", "lp");

  //Leg_Const_2 -> AddEntry(g_Kotera_shade, "GZK, Kotera '10", "f");
  //Leg_Const_2 -> AddEntry(g_ESS_shade, "GZK, ESS '01", "f");
  //Leg_Const_2 -> AddEntry(g_WB98_3p2, "W&B '98 x 3/2", "l");

  /*
    Leg_Const_2 -> AddEntry(g_ANITA, "ANITA II", "lp");
    Leg_Const_2 -> AddEntry(g_Auger, "Auger", "lp");
    Leg_Const_2 -> AddEntry(g_Icecube, "IC40, Diff. Limit", "lp");
    Leg_Const_2 -> AddEntry(g_IC40_Hummer, "IC40", "lp");
    Leg_Const_2 -> AddEntry(g_RICE, "RICE '11", "lp");
  */

  Leg_Const_2 -> SetBorderSize(0);
  Leg_Const_2 -> SetTextSize(0.025);
  Leg_Const_2 -> SetTextFont(42);
  Leg_Const_2 -> SetFillColor(0);
  Leg_Const_2 -> SetMargin(0.32);
  Leg_Const_2 -> Draw();


  //TLegend *Leg_Const2_2 = new TLegend(0.14, 0.15, 0.5, 0.4);
  //TLegend *Leg_Const2_2 = new TLegend(0.14, 0.14, 0.5, 0.32);
  TLegend *Leg_Const2_2 = new TLegend(0.16, 0.12, 0.43, 0.254);
  // TLegend *Leg_Const2_2 = new TLegend(0.16, 0.2, 0.48, 0.334);
  //TLegend *Leg_Const2_2 = new TLegend(0.14, 0.15, 0.5, 0.27);
  //TLegend *Leg_Const2_2 = new TLegend(0.14, 0.12, 0.5, 0.4);
  //TLegend *Leg_Const2_2 = new TLegend(0.21, 0.179, 0.425, 0.42);
  //TLegend *Leg_Const2_2 = new TLegend(0.19, 0.179, 0.415, 0.42);
  //TLegend *Leg_Const2_2 = new TLegend(0.19, 0.179, 0.415, 0.37);
  //TLegend *Leg_Const2_2 = new TLegend(0.19, 0.179, 0.425, 0.35);
  //TLegend *Leg_Const2_2 = new TLegend(0.205, 0.18, 0.4, 0.42); // org (110812)
  //
  //Leg_Const2_2 -> AddEntry(g_ANITA, "ANITA II", "lp");
  //Leg_Const2_2 -> AddEntry(g_Auger, "Auger '09", "lp");
  //Leg_Const2_2 -> AddEntry(g_Auger11, "Auger '11", "lp");
  //Leg_Const2_2 -> AddEntry(g_Icecube, "IceCube40", "lp");
  //Leg_Const2_2 -> AddEntry(g_IC_new, "IceCube2012 (2yrs)", "lp");
  //Leg_Const2_2 -> AddEntry(g_IC40_Hummer, "IC40", "lp");
  //Leg_Const2_2 -> AddEntry(g_RICE, "RICE '11", "lp");

  //Leg_Const2_2 -> AddEntry(g_TestBed_3yr, "TestBed (3yr) AraSim", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_Dave11, "TestBed (5 months)", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_phase, "TestBed AS phase (3yrs)", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_phase, "TestBed AS phase (5months)", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_limit_2012_75days_trig, "TestBed 2012 trig (75days)", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_limit_2012_75days, "TestBed 2012 limit (75days)", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_limit_20112012_210days, "TestBed 20112012 limit (210days)", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_limit_20112012_75_210days, "#splitline{ARA TestBed}{2011-2012 limit (285 days)}", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_limit_2012_75days_Kotera, "TestBed 2012 limit (75days) Kotera", "lp");
  //Leg_Const2_2 -> AddEntry(g_ARA37_up, "ARA37 (3yrs)", "lp");
  //Leg_Const_2 -> AddEntry(g_ARA37, "ARA37 (3yrs)", "lp");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade, "ARA37-3yr arXiv", "f");
  //
  //Leg_Const2_2 -> AddEntry(g_TestBed_limit_20112012_224days_debug, "ARA TestBed new", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_limit_20112012_224days_debug, "#splitline{ARA TestBed}{2011-2012 (224 days)}", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_limit_20112012_224days_debug, "#splitline{ARA TestBed}{2011-2012 (382 days)}", "lp");
  //Leg_Const2_2 -> AddEntry(g_TestBed_limit_20112012_224days_debug, "#splitline{ARA TestBed}{2011-2012 (415 days)}", "lp");

  //Leg_Const2_2 -> AddEntry(g_ARA2_3_Thomas_shade, "ARA A2, A3 10 months", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA2_3_Thomas_points, "ARA A2, A3 10 months", "lp");
  //Leg_Const2_2 -> AddEntry(g_ARA2_3_Thomas_shade_20150430, "ARA A2, A3 10 months", "f");
  //
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug, "ARA37 (3yrs) debug", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_oldR, "ARA37-3yr debug oldR", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_NoPhase, "ARA37-3yr debug NoPhase", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_10degCut, "ARA37-3yr debug 10degCut", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_RAYSOL10km, "ARA37-3yr debug RAYSOL10km", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_SameShowerL, "ARA37-3yr debug SameShowerL", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_Vsat, "ARA37-3yr debug Vsat", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_MakeArraySat, "ARA37-3yr debug MakeArraySat", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_MakeArraySat_wPhase, "ARA37-3yr new", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_MakeArraySat_wPhase, "ARA37 (3yr)", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_MakeArraySat_wPhase, "ARA37 (3 yrs)", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_MakeArraySat_wPhase, "ARA37 (3 yrs) Trigger (simulated)", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_ExpAnalysis, "ARA37 (3 yrs) Analysis (simulated)", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_MakeArraySat_woPhase, "ARA37-3yr new", "f");

  //Leg_Const2_2 -> AddEntry(g_ARA37_Fdomain_debug_MakeArraySat_woSecondary, "ARA37-3yr woSecondary oldRF", "lp");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_MakeArraySat_woSecondary, "ARA37-3yr new woSecondary", "f");
  //Leg_Const2_2 -> AddEntry(g_ARA37_shade_debug_MakeArraySat_woSecondary, "ARA37-3yr new woSecondary", "l");

  //Leg_Const2_2 -> AddEntry(g_ARA37_debug, "ARA37 (3yrs) debug", "f");

  Leg_Const2_2 -> AddEntry(g_Kotera_shade, "GZK, Kotera '10", "f");

  //Leg_Const2_2 -> AddEntry(g_Ave07Femix, "Ave '07 Fe mix", "l");
  //Leg_Const2_2 -> AddEntry(g_Kotera_max, "Kotera '10 max", "l");
  //Leg_Const2_2 -> AddEntry(g_Kotera_mid, "Kotera '10 mid", "l");
  //Leg_Const2_2 -> AddEntry(g_Kotera_shade, "Kotera '10 range", "f");
  //Leg_Const2_2 -> AddEntry(g_Kotera_shade, "GZK, Kotera '10", "f");
  //Leg_Const2_2 -> AddEntry(g_Kotera_shade, "GZK, KAO '10", "f");
  //Leg_Const2_2 -> AddEntry(g_Kotera_shade_ver2, "GZK, KAO '10", "f");
  //Leg_Const2_2 -> AddEntry(g_Kotera_low, "Kotera '10 low", "l");
  Leg_Const2_2 -> AddEntry(g_Ahlers, "Ahlers '11, E_{min}=10^{18.5} eV", "l");
  //Leg_Const2_2 -> AddEntry(g_Ahlers, "Ahlers '10", "l");
  //Leg_Const2_2 -> AddEntry(g_ESS_strong, "ESS '01 strong", "l");
  //Leg_Const2_2 -> AddEntry(g_ESS_shade, "ESS '01 range", "f");
  //Leg_Const2_2 -> AddEntry(g_ESS_shade, "ESS '01 range", "f");
  //Leg_Const2_2 -> AddEntry(g_ESS_shade, "GZK, ESS '01", "f");
  //Leg_Const2_2 -> AddEntry(g_ESS_shade_ver2, "GZK, ESS '01", "f");
  //Leg_Const2_2 -> AddEntry(g_ESS_base, "ESS '01 base", "l");
  //Leg_Const2_2 -> AddEntry(g_GRB_IC_FC, "GRB IC-FC", "l");
  //Leg_Const2_2 -> AddEntry(g_GRB_NFC, "GRB NFC", "l");
  //Leg_Const2_2 -> AddEntry(g_GRB_shade2, "GRB IC-FC, NFC range", "f");
  //Leg_Const2_2 -> AddEntry(g_GRB_shade2, "GRB Prompt", "f");
  //Leg_Const2_2 -> AddEntry(g_GRB_Hummer_shade, "GRB Prompt", "f");
  //Leg_Const2_2 -> AddEntry(g_GRB_Hummer_shade2, "GRB Prompt", "f");
  //Leg_Const2_2 -> AddEntry(g_GRB_WIND_shade, "GRB WIND range", "f");
  //Leg_Const2_2 -> AddEntry(g_GRB_WIND_shade, "GRB AG (WIND), Murase '12", "f");
  //Leg_Const2_2 -> AddEntry(g_GRB_shade, "GRB BATSE, Swift, Konus range", "f");
  //Leg_Const2_2 -> AddEntry(g_GRB_BATSE, "GRB BATSE sample", "l");
  //Leg_Const2_2 -> AddEntry(g_GRB_ISM_shade, "GRB ISM range", "f");
  //Leg_Const2_2 -> AddEntry(g_GRB_ISM_shade, "GRB AG (ISM), Murase '12", "f");
  //Leg_Const2_2 -> AddEntry(g_GRB_ISM_shade_ver2, "GRB AG (ISM), Murase '12", "f");
  //Leg_Const2_2 -> AddEntry(g_GRB_WB99, "GRB WB '99", "l");
  //Leg_Const2_2 -> AddEntry(g_GRB_WIND_MAX, "GRB AG, Murase '12", "l");
  //Leg_Const2_2 -> AddEntry(g_GRB_WB_AG, "GRB AG WB '00", "l");

  Leg_Const2_2 -> SetBorderSize(0);
  Leg_Const2_2 -> SetFillColor(0);
  Leg_Const2_2 -> SetTextFont(42);
  Leg_Const2_2 -> SetTextSize(0.03);
  Leg_Const2_2 -> Draw();

  cConst_2->Print("test_Sensitivity_20170801.png");
  cConst_2->Print("test_Sensitivity_20170801.pdf");
  cConst_2->Print("test_Sensitivity_20170801.eps");












}





TStyle* RootStyle() {


  TStyle *RootStyle = new TStyle("Root-Style","The Perfect Style for Plots ;-)");

#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                          // save the global style reference

  gStyle = RootStyle;
#endif
  // otherwise you need to call TROOT::SetStyle("Root-Style")

  // Paper size

  RootStyle->SetPaperSize(TStyle::kUSLetter);

  // Canvas

  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderSize(10);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH      (600);
  RootStyle->SetCanvasDefW      (600);
  RootStyle->SetCanvasDefX      (10);
  RootStyle->SetCanvasDefY      (10);

  // Pads

  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderSize  (10);
  RootStyle->SetPadBorderMode  (0);
  //  RootStyle->SetPadBottomMargin(0.13);
  /*
    RootStyle->SetPadBottomMargin(0.16);
    RootStyle->SetPadTopMargin   (0.08);
    RootStyle->SetPadLeftMargin  (0.18);
    RootStyle->SetPadRightMargin (0.05);
  */
  RootStyle->SetPadGridX       (0);
  RootStyle->SetPadGridY       (0);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);

  // for wide plot
  RootStyle->SetPadBottomMargin(0.115);
  RootStyle->SetPadTopMargin   (0.05);
  RootStyle->SetPadLeftMargin  (0.12);
  RootStyle->SetPadRightMargin (0.04);

  // Frames

  RootStyle->SetFrameFillStyle ( 0);
  RootStyle->SetFrameFillColor ( 0);
  RootStyle->SetFrameLineColor ( 1);
  RootStyle->SetFrameLineStyle ( 0);
  RootStyle->SetFrameLineWidth ( 2);
  RootStyle->SetFrameBorderSize(10);
  RootStyle->SetFrameBorderMode( 0);


  // Histograms

  RootStyle->SetHistFillColor(0);
  RootStyle->SetHistFillStyle(1);
  RootStyle->SetHistLineColor(1);
  RootStyle->SetHistLineStyle(0);
  RootStyle->SetHistLineWidth(2);

  // Functions

  RootStyle->SetFuncColor(1);
  RootStyle->SetFuncStyle(0);
  RootStyle->SetFuncWidth(1);

  //Legends

  //  RootStyle->SetStatBorderSize(0);
  //  RootStyle->SetStatFill(0);
  RootStyle->SetStatFont      (42);
  //  RootStyle->SetOptStat       (111111);
  //  RootStyle->SetOptStat       (0);
  //RootStyle->SetStatColor     (0);
  RootStyle->SetStatColor     (2);
  //  RootStyle->SetStatX         (0.93);
  //  RootStyle->SetStatY         (0.90);
  //RootStyle->SetStatFontSize  (0.2);
  RootStyle->SetStatFontSize  (2);
  //  RootStyle->SetStatW         (0.2);
  //  RootStyle->SetStatH         (0.15);

  // Labels, Ticks, and Titles

  RootStyle->SetTickLength ( 0.015,"X");
  //RootStyle->SetTitleSize  ( 0.04,"X");
  RootStyle->SetTitleSize  ( 0.1,"X");
  RootStyle->SetTitleOffset( 1.6,"X");
  RootStyle->SetTitleBorderSize(0);
  RootStyle->SetTitleFillColor(0);

  //  RootStyle->SetTitleFontSize((double)0.5,"X");
  RootStyle->SetLabelOffset( 0.015,"X");
  //RootStyle->SetLabelSize  ( 0.050,"X");
  RootStyle->SetLabelSize  ( 0.040,"X");
  RootStyle->SetLabelFont  ( 42   ,"X");

  RootStyle->SetTickLength ( 0.015,"Y");
  RootStyle->SetTitleSize  ( 0.06,"Y");
  RootStyle->SetTitleSize  ( 0.05,"X");

  RootStyle->SetTitleOffset( 1.500,"Y");
  RootStyle->SetLabelOffset( 0.015,"Y");
  //RootStyle->SetLabelSize  ( 0.050,"Y");
  RootStyle->SetLabelSize  ( 0.040,"Y");
  RootStyle->SetLabelFont  ( 42   ,"Y");

  RootStyle->SetTitleFont  (42,"XY");
  RootStyle->SetTitleColor  (1);




  // Options

  RootStyle->SetOptFit     (1);

  RootStyle->SetMarkerStyle(20);
  RootStyle->SetMarkerSize(0.4);

  //  cout << ">> Style initialized with the Root Style!" << endl;
  //  cout << ">> " << modified << endl << endl;
  return RootStyle;
}


void LogToLine(int N, double *Data) {

  for (int i=0; i<N; i++) {
    Data[i] = pow( 10., Data[i] );
  }
}


void ApplyFactorArray( double *output, double *input, double *factor, int size ) {

  for (int i=0; i<size; i++) {
    output[i] = input[i] * factor[i];
  }
}




// get average factor difference between arrays
double GetAvgFactor( double *EArray1, double *FluxArray1, int bin1, double *EArray2, double *FluxArray2, int bin2 ) {

  double AvgFactor = 0.;
  double SumBin = 0.;

  for (int bin=0; bin<bin1; bin++) {

    for (int binloop=0; binloop<bin2; binloop++) {


      if ( EArray1[bin] == EArray2[binloop] ) {

	AvgFactor += FluxArray1[bin] / FluxArray2[binloop];

	SumBin ++;

	//break;
      }
    }

    //cout<<"loop bin "<<bin<<endl;

  }

  AvgFactor = AvgFactor / SumBin;

  return AvgFactor;

}



// get RMS factor difference between arrays
double GetRMSFactor( double *EArray1, double *FluxArray1, int bin1, double *EArray2, double *FluxArray2, int bin2, double AvgFactor ) {


  double RMSFactor = 0.;
  double SumBin = 0.;

  for (int bin=0; bin<bin1; bin++) {

    for (int binloop=0; binloop<bin2; binloop++) {


      if ( EArray1[bin] == EArray2[binloop] ) {

	RMSFactor += pow( AvgFactor - (FluxArray1[bin] / FluxArray2[binloop]) , 2 );

	SumBin ++;

	//break;
      }
    }

    //cout<<"loop bin "<<bin<<endl;

  }

  RMSFactor = RMSFactor / SumBin;

  return sqrt(RMSFactor);

}
