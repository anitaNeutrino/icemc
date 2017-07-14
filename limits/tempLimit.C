const int n_ANITA = 7;
Double_t ANITA_4_x[n_ANITA]       = {18, 18.5, 19, 19.5, 20, 20.5, 21};
void LogToLine(int N, double *Data);
TGraph *getKoteraShade();
TGraph *getANITA2erratum();
TGraph *getAhlers();
TGraph *getLimit(double effArea[n_ANITA], double eff[n_ANITA], double livetime);


void tempLimit(){

  LogToLine(n_ANITA, ANITA_4_x);
  TGraph *g_Kotera_shade = getKoteraShade();

  TGraph *g_ANITA_2_erratum = getANITA2erratum();
  
  Double_t ANITA_4_effArea[n_ANITA] = { 0.00153,      // E18     in km^2  
					0.04991,      // E18.5		  
					0.53408,      // E19    	  
					2.98348,      // E19.5		  
					13.02480,     // E20		  
					41.35640,     // E20.5		  
					107.672 };    // E21               
  
  double ANITA_4_eff[n_ANITA] = { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};

  double ANITA_4_livetime = 27.3*24*3600.; //27.3*24*3600; // 27.3 days
  
  TGraph *g_ANITA_4 = getLimit(ANITA_4_effArea, ANITA_4_eff, ANITA_4_livetime);

  Double_t ANITA_3_effArea[n_ANITA] = { 0.00160,      // E18     in km^2  
					0.03451,      // E18.5		  
					0.34288,      // E19    	  
					1.75735,      // E19.5		  
					7.14396,      // E20		  
					22.77670,     // E20.5		  
					58.72190};    // E21              


  double ANITA_3_eff[n_ANITA] = { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};

  double ANITA_3_livetime = 17.4*24*3600.; // 17.4 days
  
  TGraph *g_ANITA_3 = getLimit(ANITA_3_effArea, ANITA_3_eff, ANITA_3_livetime);
  g_ANITA_3->SetLineColor(kBlue);

  
  Double_t ANITA_2_effArea[n_ANITA] = { 0.00029 ,      // E18     in km^2  
					0.02121 ,      // E18.5		  
					0.22009 ,      // E19    	  
					1.27190 ,      // E19.5		  
					6.38188 ,      // E20		  
					18.45390 ,     // E20.5		  
					52.53270 };    // E21              

  double ANITA_2_eff[n_ANITA] = { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};

  double ANITA_2_livetime = 28.5*24*3600.; // 17.4 days
  
  TGraph *g_ANITA_2 = getLimit(ANITA_2_effArea, ANITA_2_eff, ANITA_2_livetime);
  g_ANITA_2->SetLineColor(kRed);


  
  TCanvas *cConst_2 = new TCanvas("cConst_2","A Simple Graph Example",200,10,1400,1400); // wider

  cConst_2->cd();
  cConst_2->SetLogy();
  cConst_2->SetLogx();
  //g_ESS_base->SetTitle("Sensitivity");
  g_Kotera_shade->SetTitle(";E (eV);E dN/dE dA d#Omega dt (cm^{-2} sr ^{-1} s^{-1} )");
  // g_Kotera_shade->GetHistogram()->SetXTitle("E (eV)");
  // g_Kotera_shade->GetHistogram()->SetYTitle("E dN/dE dA d#Omega dt (cm^{-2} sr ^{-1} s^{-1} )");
  //g_vhvh->GetHistogram()->SetMaximum(100);
  //g_vhvh->GetHistogram()->SetMinimum(0.1);
  //g_Kotera_shade->GetHistogram()->SetMaximum(1.e-11);
  g_Kotera_shade->GetHistogram()->SetMaximum(1.e-12);
  //g_Kotera_shade->GetHistogram()->SetMaximum(3.e-11);
  //g_Kotera_shade->GetHistogram()->SetMinimum(1.e-20);
  //g_Kotera_shade->GetHistogram()->SetMaximum(1.e-13);
  //g_Kotera_shade->GetHistogram()->SetMinimum(1.e-20);
  g_Kotera_shade->GetHistogram()->SetMinimum(1.e-19);
  //g_Kotera_shade->GetHistogram()->SetMinimum(1.e-21);
  //g_Kotera_shade->GetHistogram()->SetMinimum(1.e-22);

  //  g_Kotera_shade->GetXaxis()->SetLimits(3.2e14,1.e24);
  g_Kotera_shade->GetXaxis()->SetLimits(1e17,1.e22); // zoom little bit
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

  TGraph *g_Ahlers=getAhlers();
  g_Ahlers->Draw("l");
  
  g_ANITA_2_erratum->Draw("l");

  g_ANITA_2->Draw("l");
  g_ANITA_3->Draw("l");
  g_ANITA_4->Draw("l");


  TLegend *leg = new TLegend(0.6, 0.6, 0.89, 0.89);
  leg->AddEntry(g_ANITA_2_erratum, "A2 erratum",                  "lp" );
  leg->AddEntry(g_ANITA_2,         "A2 icemc #epsilon_{ANA}=0.8",  "l" );
  leg->AddEntry(g_ANITA_3,         "A3 icemc #epsilon_{ANA}=0.8",  "l" );
  leg->AddEntry(g_ANITA_4,         "A4 icemc #epsilon_{ANA}=0.8",  "l" );
  leg->Draw();
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

  // for (int i=0; i<11; i++) cout << ANITA_erratum_x[i] << " " << ANITA_erratum_y[i] << endl;
  TGraph *g_ANITA_erratum = new TGraph( 11, ANITA_erratum_x, ANITA_erratum_y );

  g_ANITA_erratum->SetLineWidth(3);
  g_ANITA_erratum->SetMarkerStyle(20);
  g_ANITA_erratum->SetMarkerSize(2.2);
  g_ANITA_erratum->SetLineColor(kRed+2);
  g_ANITA_erratum->SetMarkerColor(kRed+2);

  
  return g_ANITA_erratum;

}


TGraph* getLimit(double effArea[n_ANITA], double eff[n_ANITA], double livetime){

  Double_t ANITA_4_y[n_ANITA];
  
  double N90 = 2.30; 

  for (int i=0; i<n_ANITA; i++){

    // convert km^2 in cm^2
    effArea[i]*=1e10;
    
    ANITA_4_y[i] = N90/(effArea[i]*livetime*eff[i]);

    //    ANITA_4_y[i] *= TMath::Power(10, ANITA_4_x[i])/(TMath::Power(10, ANITA_4_x[i]+0.25) - TMath::Power(10, ANITA_4_x[i]-0.25));
    
    // Divide by 4 
    ANITA_4_y[i] /= 4.;

    std::cout << ANITA_4_x[i] << " " << ANITA_4_y[i] << std::endl;
  }
      
      
  TGraph *g_ANITA_4 = new TGraph(n_ANITA, ANITA_4_x, ANITA_4_y);


  g_ANITA_4->SetLineWidth(3);
  g_ANITA_4->SetLineStyle(2);
  
  return g_ANITA_4;
}

void LogToLine(int N, double *Data) {

  for (int i=0; i<N; i++) {
    Data[i] = pow( 10., Data[i] );
  }
}
