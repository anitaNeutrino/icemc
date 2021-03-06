#include "Acceptances.h"

void LogToLine(int N, double *Data);
TGraph *getKoteraShade();
TGraph *getANITA2erratum();
TGraph *getAhlers();
TGraphAsymmErrors *getIceCube();

TGraph *getLimitOldFormula(double N90, double effArea[n_ANITA], double eff[n_ANITA], double livetime);

TGraph *getLimitNewFormula(double N90, double effArea[n_ANITA], double eff[n_ANITA], double livetime);

TGraph* getCombinedLimitNoDelta(double denom[n_ANITA], double N90);

TGraph* getCombinedLimit(double denom[n_ANITA], double N90);

TGraph *auger2015();
TGraph *auger2017();
TGraph *icecube();
TGraph *icecube2018();

void anita4paperLimit(){

  string outname = "LimitTemp_ANITA4";

  //   if (!gROOT->GetClass("TFeldmanCousins")) gSystem->Load("libPhysics");

  TFeldmanCousins f;

  // calculate either the upper or lower limit for 10 observed
  // events with an estimated background of 3.  The calculation of
  // either upper or lower limit will return that limit and fill
  // data members with both the upper and lower limit for you.
  Double_t Nobserved   = 1.0;
  Double_t Nbackground = 0.644;
//  Double_t Nbackground = 0.34;

  Double_t NobsAll  = 4.0;
  Double_t NbkgAll  = (Nbackground+0.7+0.98+1.1);

  Double_t A3ul = 3.471 ;
  Double_t A4ul = f.CalculateUpperLimit(Nobserved, Nbackground);
//  printf("a4 ul = %g\n", A4ul);
  //   Double_t ll = f.GetLowerLimit();


  
  Double_t ulA123 = f.CalculateUpperLimit(3, 0.7+0.98+1.1);
  Double_t ulAll = f.CalculateUpperLimit(NobsAll, NbkgAll);

  cout << "ANITA 1-4 Upper Limit = " <<  ulAll << endl;
     
  LogToLine(n_ANITA, ANITA_x);
  TGraph *g_Kotera_shade = getKoteraShade();  
  
  Double_t ANITA_3_effArea[n_ANITA];             
  Double_t ANITA_4_effArea[n_ANITA];             
  Double_t ANITA_3_effAreaUp[n_ANITA];
  Double_t ANITA_3_effAreaLow[n_ANITA];
  Double_t ANITA_3_effAreaReno[n_ANITA];
  Double_t ANITA_3_geomAverage[n_ANITA];
  Double_t ANITA_4_geomAverage[n_ANITA];
  
  for (int i=0; i<n_ANITA; i++){
    ANITA_3_effArea[i]     = ANITA_3_effVol[i]/intLength_CONNOLLY_nuCC[i]; 
    ANITA_4_effArea[i]     = ANITA_4_effVol[i]/intLength_CONNOLLY_nuCC[i]; 
    // ANITA_3_effAreaUp[i]   = ANITA_3_effVol[i]/intLength_CONNOLLY_nuCCup[i];
    // ANITA_3_effAreaLow[i]  = ANITA_3_effVol[i]/intLength_CONNOLLY_nuCClow[i]; 
    ANITA_3_effAreaReno[i] = ANITA_3_effVol[i]/intLength_RENO[i];
    ANITA_2_effArea_Peter[i]  = ANITA_2_effArea_published[i]*ANITA_2_effArea_published[i]/ANITA_2_effArea_icemc2010[i];
    ANITA_3_geomAverage[i] = TMath::Sqrt(ANITA_3_effArea[i]*ANITA_2_effArea_Peter[i]*ANITA_3_effArea[i]/ANITA_2_effArea_icemc2010[i]);
    ANITA_4_geomAverage[i] = TMath::Sqrt(ANITA_4_effArea[i]*ANITA_2_effArea_Peter[i]*ANITA_4_effArea[i]/ANITA_2_effArea_icemc2010[i]);
    // cout  << "ANITA-3 geom  " << ANITA_3_geomAverage[i] << " " << endl;
    ANITA_3_effAreaUp[i]   = ANITA_3_geomAverage[i]*intLength_CONNOLLY_nuCC[i]/intLength_CONNOLLY_nuCCup[i];
    ANITA_3_effAreaLow[i]  = ANITA_3_geomAverage[i]*intLength_CONNOLLY_nuCC[i]/intLength_CONNOLLY_nuCClow[i];
    cout <<  ANITA_4_geomAverage[i] << endl;
  }
  

  double N90A3 = A3ul;
  
  cout << "ANITA3" << endl;
  TGraph *g_ANITA_3_combined = getLimitOldFormula(N90A3, ANITA_3_geomAverage, ANITA_3_eff, ANITA_3_livetime);
  g_ANITA_3_combined->SetLineStyle(1);
  g_ANITA_3_combined->SetLineColor(kRed+2);

  double N90A4 = A4ul;
  cout << "ANITA4" << endl;
  TGraph *g_ANITA_4_combined = getLimitOldFormula(N90A4, ANITA_4_geomAverage, ANITA_4_eff, ANITA_4_livetime);
  g_ANITA_4_combined->SetLineStyle(1);
  g_ANITA_4_combined->SetLineColor(kBlack);
  // g_ANITA_4_combined->SetLineStyle(2);


  Double_t ANITA_2_effArea[n_ANITA];
  for (int i=0; i<n_ANITA; i++){
    ANITA_2_effArea[i]    = ANITA_2_effVol[i]/intLength_CONNOLLY_nuCC[i];
    // cout << ANITA_2_effVol[i]/intLength_RENO[i] << endl;
  }


  double N90A2_sens = 2.3;
  double N90A2_lim  = 3.39;
  
  
  double ANITA_123_denom[n_ANITA];
  double ANITA_1234_denom[n_ANITA];
  for (int ibin=0; ibin<n_ANITA; ibin++){
    ANITA_123_denom[ibin] = (
			     ANITA_1_livetime*ANITA_1_effArea[ibin]*ANITA_1_eff[ibin] +
			     ANITA_2_livetime*ANITA_2_effArea[ibin]*ANITA_2_eff[ibin] +
			     ANITA_3_livetime*ANITA_3_geomAverage[ibin]*ANITA_3_eff[ibin] 
			     );
    ANITA_1234_denom[ibin] = (
			     ANITA_1_livetime*ANITA_1_effArea[ibin]*ANITA_1_eff[ibin] +
			     ANITA_2_livetime*ANITA_2_effArea[ibin]*ANITA_2_eff[ibin] +
			     ANITA_3_livetime*ANITA_3_geomAverage[ibin]*ANITA_3_eff[ibin] +
			     ANITA_4_livetime*ANITA_4_geomAverage[ibin]*ANITA_4_eff[ibin] 
			     );
  }


  cout << "Combined limit ANITA 1-3" << endl;
  TGraph *g_ANITA_123 = getCombinedLimit(ANITA_123_denom, ulA123);
  g_ANITA_123->SetLineColor(kCyan);

  TGraph *g_ANITA_1234 = getCombinedLimit(ANITA_1234_denom, ulAll);
  g_ANITA_1234->SetLineColor(kCyan);
  // g_ANITA_1234->SetLineStyle(2);

  // cout << "Combined limit no delta" << endl;
  // TGraph *gtemp = getCombinedLimitNoDelta(ANITA_123_denom);
  // gtemp->SetLineColor(kRed);
  
  TCanvas *cConst_2 = new TCanvas("cConst_2","A Simple Graph Example",800,800); // wider

  cConst_2->cd();
  cConst_2->SetLogy();
  cConst_2->SetLogx();
  //  cConst_2->SetGrid();
  g_Kotera_shade->SetTitle(";E (eV);E dN/dE dA d#Omega dt (cm^{-2} sr ^{-1} s^{-1} )");
  g_Kotera_shade->GetHistogram()->SetMaximum(1.e-12);
  g_Kotera_shade->GetHistogram()->SetMinimum(1.e-19);
  g_Kotera_shade->GetXaxis()->SetLimits(1e16,1.e21); // zoom little bit
  g_Kotera_shade->GetHistogram()->SetTitleSize  ( 0.04,"X");
  g_Kotera_shade->GetHistogram()->SetLabelOffset( 0.006,"X");
  g_Kotera_shade->GetHistogram()->SetLabelSize( 0.04,"X");
  g_Kotera_shade->GetHistogram()->SetLabelSize( 0.04,"Y");
  g_Kotera_shade->GetHistogram()->SetLabelOffset( 0.007,"Y");
  g_Kotera_shade->GetHistogram()->SetTitleSize  ( 0.04,"Y");
  g_Kotera_shade->GetHistogram()->SetTitleOffset( 1.8,"Y");
  g_Kotera_shade->GetHistogram()->SetTitleOffset( 1.3,"X");
  gPad->SetLeftMargin(0.15);
  g_Kotera_shade->SetFillStyle(1001);
  g_Kotera_shade->SetFillColor(15);
  g_Kotera_shade->SetLineColor(0);
  g_Kotera_shade->Draw("af");
  gPad->RedrawAxis();

  TGraph *g_Ahlers=getAhlers();
  g_Ahlers->Draw("l");

  
  TGraph *g_ANITA_2_erratum = getANITA2erratum();
  // g_ANITA_2_erratum->Draw("l");
  // gtemp->Draw("l");
  // g_ANITA_3_shade->Draw("f");
  g_ANITA_3_combined->Draw("l");
  // g_ANITA_3up->Draw("l");
  // g_ANITA_3low->Draw("l");
  // g_ANITA_3Reno->Draw("l");
  //  g_ANITA_123->Draw("l");
  g_ANITA_1234->Draw("l");
  g_ANITA_4_combined->Draw("l");

  TGraph *gAuger   = auger2017();
  TGraph *gIcecube = icecube2018();


  gAuger->SetLineColor(kRed);
  gAuger->SetLineStyle(2);
  gAuger->Draw("l");

  gIcecube->SetLineColor(kBlue);
  gIcecube->SetLineStyle(3);
  gIcecube->Draw("l");
  
  TLegend *leg = new TLegend(0.5, 0.7, 0.89, 0.89);
  leg->AddEntry(gAuger,   "Auger 2017", "l");
  leg->AddEntry(gIcecube, "IceCube 2018", "l");
  leg->AddEntry(g_ANITA_3_combined,    "ANITA-III",  "l" );
  // leg->AddEntry(g_ANITA_2_erratum, "ANITA-II erratum", "l");
  // leg->AddEntry(gtemp, "ANITA I-III arxiv", "l");
  // // leg->AddEntry(g_ANITA_3,    "#sigma Connolly et al, Nominal", "l");
  // // leg->AddEntry(g_ANITA_3up,  "#sigma Connolly et al, Upper bound",   "l");
  // // leg->AddEntry(g_ANITA_3low, "#sigma Connolly et al, Lower bound",   "l");
  // // leg->AddEntry(g_ANITA_3Reno,        "#sigma Reno et al",     "l");
  // leg->AddEntry(g_ANITA_123,       "ANITA I-III",  "l" );
  leg->AddEntry(g_ANITA_4_combined,       "ANITA IV",  "l" );
  leg->AddEntry(g_ANITA_1234,       "ANITA I-IV",  "l" );
  leg->Draw();

  
  TLegend *Leg_Const2_2 = new TLegend(0.16, 0.12, 0.43, 0.254);
  Leg_Const2_2 -> AddEntry(g_Kotera_shade, "GZK, Kotera '10", "f");  
  Leg_Const2_2 -> AddEntry(g_Ahlers, "Ahlers '12, E_{min}=10^{18.5} eV", "l");
  Leg_Const2_2 -> SetBorderSize(0);
  Leg_Const2_2 -> SetFillColor(0);
  Leg_Const2_2 -> Draw();

  //  Preliminary();
  
  cConst_2->Print(Form("%s.png", outname.c_str()));
  cConst_2->Print(Form("%s.pdf", outname.c_str()));
  cConst_2->Print(Form("%s.C", outname.c_str()));
}


TGraph *getKoteraShade(){



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
  

  TGraph *g_Kotera_shade = new TGraph ( 42 + 37 ); // 42 for max and 37 for low
  for (int i=0; i<42; i++) {
    g_Kotera_shade->SetPoint(i, Kotera_max_x[i], Kotera_max_y[i]);
  }
  for (int i=0; i<37; i++) {
    g_Kotera_shade->SetPoint(i+42, Kotera_low_x[37-i-1], Kotera_low_y[37-i-1]);
  }
  
  return g_Kotera_shade;
}


TGraph *getAhlers(){
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
  g_Ahlers->SetLineColor(43);
  g_Ahlers->SetLineWidth(3);
  
  return g_Ahlers;
}

TGraph *getANITA2erratum(){

  
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

  // for (int i=0; i<11; i++) cout << "ANITA-2 erratum : " << ANITA_erratum_x[i] << " " << ANITA_erratum_y[i] << endl;
  TGraph *g_ANITA_erratum = new TGraph( 11, ANITA_erratum_x, ANITA_erratum_y );

  g_ANITA_erratum->SetLineWidth(3);
  g_ANITA_erratum->SetMarkerStyle(20);
  g_ANITA_erratum->SetMarkerSize(2.2);
  g_ANITA_erratum->SetLineColor(kRed+2);
  g_ANITA_erratum->SetMarkerColor(kRed+2);

  
  return g_ANITA_erratum;

}


TGraph* getLimitOldFormula(double N90, double effArea[n_ANITA], double eff[n_ANITA], double livetime){

  Double_t ANITA_4_y[n_ANITA];

  for (int i=0; i<n_ANITA; i++){
    
    ANITA_4_y[i] = N90/(effArea[i]*1e10*livetime*eff[i]);

    double exponent = TMath::Log10(ANITA_x[i]);
    
    //  ANITA_4_y[i] *= TMath::Power(10, exponent)/(TMath::Power(10, exponent+0.25) - TMath::Power(10, exponent-0.25));
    
    // Divide by 4 
    ANITA_4_y[i] /= 4.;

    
    printf("%8.3e %8.3e \n", ANITA_x[i], ANITA_4_y[i]);
    printf("eff area %8.3e %8.3e \n", ANITA_x[i], effArea[i]);
  }
      
      
  TGraph *g_ANITA_4 = new TGraph(n_ANITA, ANITA_x, ANITA_4_y);


  g_ANITA_4->SetLineWidth(3);
  g_ANITA_4->SetLineStyle(2);
  
  return g_ANITA_4;
}

TGraph* getCombinedLimit(double denom[n_ANITA], double N90){

  Double_t ANITA_4_y[n_ANITA];

  for (int i=0; i<n_ANITA; i++){
    
    ANITA_4_y[i] = N90/(denom[i]*1e10);

    ANITA_4_y[i] /= 4.;
    printf("%8.3e %8.3e \n", ANITA_x[i], ANITA_4_y[i]);

  }
      
      
  TGraph *g_ANITA_4 = new TGraph(n_ANITA, ANITA_x, ANITA_4_y);


  g_ANITA_4->SetLineWidth(3);
  g_ANITA_4->SetLineStyle(1);
  
  return g_ANITA_4;
}

TGraph* getLimitNewFormula(double N90, double effArea[n_ANITA], double eff[n_ANITA], double livetime){

  Double_t ANITA_4_y[n_ANITA];
  
  for (int i=0; i<n_ANITA; i++){
    
    ANITA_4_y[i] = N90/(effArea[i]*1e10*livetime*eff[i]*TMath::Log(10.));

    double exponent = TMath::Log10(ANITA_x[i]);
    
    ANITA_4_y[i] *= TMath::Power(10, exponent)/(TMath::Power(10, exponent+0.25) - TMath::Power(10, exponent-0.25));

    
    printf("%8.3e %8.3e \n", ANITA_x[i], ANITA_4_y[i]);
  }
      
      
  TGraph *g_ANITA_4 = new TGraph(n_ANITA, ANITA_x, ANITA_4_y);


  g_ANITA_4->SetLineWidth(3);
  g_ANITA_4->SetLineStyle(2);
  
  return g_ANITA_4;
}



TGraph* getCombinedLimitNoDelta(double denom[n_ANITA], double N90){

  Double_t ANITA_4_y[n_ANITA];
  
  for (int i=0; i<n_ANITA; i++){
    
    ANITA_4_y[i] = N90/(denom[i]*1e10*TMath::Log(10.));

    double exponent = TMath::Log10(ANITA_x[i]);
    
    ANITA_4_y[i] *= TMath::Power(10, exponent)/(TMath::Power(10, exponent+0.25) - TMath::Power(10, exponent-0.25));
    
    printf("%8.3e %8.3e \n", ANITA_x[i], ANITA_4_y[i]);
  }
      
      
  TGraph *g_ANITA_4 = new TGraph(n_ANITA, ANITA_x, ANITA_4_y);


  g_ANITA_4->SetLineWidth(3);
  g_ANITA_4->SetLineStyle(1);
  
  return g_ANITA_4;
}

void LogToLine(int N, double *Data) {

  for (int i=0; i<N; i++) {
    Data[i] = pow( 10., Data[i] );
  }
}



TGraphAsymmErrors *getIceCube(){

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

  
  TGraphAsymmErrors *g_IC_3yr_HESE_X = new TGraphAsymmErrors(IC3yrHESE_x.size(),&(IC3yrHESE_x[0]), &(IC3yrHESE_y[0]), &(IC3yrHESE_xl[0]), &(IC3yrHESE_xh[0]), &(IC3yrHESE_y0[0]), &(IC3yrHESE_y0[0]) ); // only x range bar

  TGraphAsymmErrors *g_IC_3yr_HESE_HL = new TGraphAsymmErrors(IC3yrHESE_hl_x.size(), &(IC3yrHESE_hl_x[0]), &(IC3yrHESE_hl_y[0]), &(IC3yrHESE_hl_0[0]), &(IC3yrHESE_hl_0[0]), &(IC3yrHESE_hl_yl[0]), &(IC3yrHESE_hl_yh[0]) );

  TGraphAsymmErrors *g_IC_3yr_HESE_H = new TGraphAsymmErrors(IC3yrHESE_h_x.size(), &(IC3yrHESE_h_x[0]), &(IC3yrHESE_h_y[0]), &(IC3yrHESE_h_0[0]), &(IC3yrHESE_h_0[0]), &(IC3yrHESE_h_yl[0]), &(IC3yrHESE_h_0[0]) );


  g_IC_3yr_HESE_X->SetMarkerColor(kBlack);
  g_IC_3yr_HESE_X->SetLineWidth(2);
  g_IC_3yr_HESE_X->SetMarkerStyle(1);
  g_IC_3yr_HESE_X->SetMarkerSize(2);
  
  g_IC_3yr_HESE_H->SetMarkerColor(kBlack);
  g_IC_3yr_HESE_H->SetLineWidth(2);
  g_IC_3yr_HESE_H->SetMarkerStyle(1);
  g_IC_3yr_HESE_H->SetMarkerSize(2);
  
  g_IC_3yr_HESE_HL->SetMarkerColor(kBlack);
  g_IC_3yr_HESE_HL->SetLineWidth(2);
  g_IC_3yr_HESE_HL->SetMarkerStyle(1);
  g_IC_3yr_HESE_HL->SetMarkerSize(2);
  
  return g_IC_3yr_HESE_X;

}


TGraph *auger2015(){

  for (int i=0; i<8; i++) {

    Auger15_y[i] = Auger15_y[i] - ( Auger15_x[i] - 9. );
  }
 
  LogToLine(8, Auger15_x);
  LogToLine(8, Auger15_y);

  // multiply factor of 3 to account all three flavors
  for (int i=0; i<8; i++) {
    Auger15_y[i] = (Auger15_y[i]*3.);
  }
  TGraph *g_Auger15 = new TGraph( 8 , Auger15_x, Auger15_y );
  g_Auger15->SetLineWidth(3);
  
  return g_Auger15;

}


TGraph *icecube(){

  LogToLine(11, IceCube2017x);
  LogToLine(11, IceCube2017y);
  
  for (int i=0; i<12; i++) {
    IceCube2017y[i]*=2;
  }
  
  TGraph *g_Icecube = new TGraph( 11, IceCube2017x, IceCube2017y );

  g_Icecube->SetLineWidth(3);
  
  return g_Icecube;
}

TGraph *icecube2018(){
  TGraph* g = new TGraph("icecube_2018.dat");
  g->SetLineWidth(3);

  return g;
}

TGraph *auger2017(){
  TGraph* g = new TGraph("auger_2017.dat");
  g->SetLineWidth(3);

  return g;
}
