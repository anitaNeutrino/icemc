const int n_ANITA = 7;
Double_t ANITA_4_x[n_ANITA]       = {18, 18.5, 19, 19.5, 20, 20.5, 21};

void tempAcceptance(){

 Double_t ANITA_2_effArea[n_ANITA] = { 0.00029 ,      // E18     in km^2  
					0.02121 ,      // E18.5		  
					0.22009 ,      // E19    	  
					1.27190 ,      // E19.5		  
					6.38188 ,      // E20		  
					18.45390 ,     // E20.5		  
					52.53270 };    // E21              

  Double_t ANITA_2_effArea_published[n_ANITA] = {0.00043,
						 0.05000,
						 0.92000,
						 6.60,
						 36.00,
						 108.00,
						 259.00};

  Double_t ANITA_2_effVolume_old[n_ANITA] = { 0.5,      // E18
					      63.1341 ,    // E18.5
					      633.598 ,    // E19
					      3118.19 ,    // E19.5
					      9414.39 ,    // E20
					      20304   ,    // E20.5
					      38143   };   // E21

  Double_t ANITA_2_effVolume_new[n_ANITA] = { ,    // E18
					      ,    // E18.5
					      207.579,    // E19
					      1333.15,    // E19.5
					      3624.32,    // E20
					      15310.8,    // E20.5
					      26217.4   };   // E21
   

  Double_t intLength_CONNOLLY_nuNC[n_ANITA] = { 3826.5,
						2617.44,
						1817.98,
						1279.9,
						912.003,
						 656.936,
						 477.867}; // E21

  
  Double_t intLength_CONNOLLY_nuCC[n_ANITA] = { 1544.84,
						1065.44,
						745.54,
						528.431,
						378.864,
						 274.445,
						 200.67}; // E21

  
  Double_t intLength_CONNOLLY_nubarNC[n_ANITA] = { 3868.47,
						   2647.45,
						   1840.89,
						   1297.89,
						   926.252,
						   668.193 ,
						   486.701 }; // E21
  
  Double_t intLength_CONNOLLY_nubarCC[n_ANITA] = { 1671.17,
						   1151.12,
						   805.022,
						   570.46,
						   408.962,
						   296.22,
						   216.548}; // E21
  
  Double_t intLength_RENO[n_ANITA] = { 1131.34,
				       793.952,
				       557.179,
				       391.016,
				       274.407,
				       192.573 ,
				       135.143 }; // E21

  TGraph *g_intLength_RENO             = new TGraph(n_ANITA, ANITA_4_x, intLength_RENO);
  TGraph *g_intLength_CONNOLLY_nuNC    = new TGraph(n_ANITA, ANITA_4_x, intLength_CONNOLLY_nuNC);
  TGraph *g_intLength_CONNOLLY_nuCC    = new TGraph(n_ANITA, ANITA_4_x, intLength_CONNOLLY_nuCC);
  TGraph *g_intLength_CONNOLLY_nubarNC = new TGraph(n_ANITA, ANITA_4_x, intLength_CONNOLLY_nubarNC);
  TGraph *g_intLength_CONNOLLY_nubarCC = new TGraph(n_ANITA, ANITA_4_x, intLength_CONNOLLY_nubarCC);

  g_intLength_RENO             ->SetLineColor(kBlack);
  g_intLength_CONNOLLY_nuNC    ->SetLineColor(kBlue);
  g_intLength_CONNOLLY_nuCC    ->SetLineColor(kRed);
  g_intLength_CONNOLLY_nubarNC ->SetLineColor(kGreen);
  g_intLength_CONNOLLY_nubarCC ->SetLineColor(kViolet);
  
  g_intLength_RENO             ->SetLineWidth(2);
  g_intLength_CONNOLLY_nuNC    ->SetLineWidth(2);
  g_intLength_CONNOLLY_nuCC    ->SetLineWidth(2);
  g_intLength_CONNOLLY_nubarNC ->SetLineWidth(2);
  g_intLength_CONNOLLY_nubarCC ->SetLineWidth(2);
  
  TCanvas *c1 = new TCanvas("c1");
  g_intLength_RENO->SetTitle("Interaction Lengths for different cross-section parametrizations;Log10(E_{#nu});Interaction Length [km]");
  g_intLength_RENO->GetYaxis()->SetRangeUser(0, 4000);
  g_intLength_RENO->Draw("Al");
  g_intLength_CONNOLLY_nuNC->Draw("l");
  g_intLength_CONNOLLY_nuCC->Draw("l");
  g_intLength_CONNOLLY_nubarNC->Draw("l");
  g_intLength_CONNOLLY_nubarCC->Draw("l");

  TLegend *leg1 = new TLegend (0.5, 0.65, 0.89, 0.89);
  leg1->AddEntry(g_intLength_RENO,                "Reno et all, #nu/#bar{#nu} CC",  "l");
  leg1->AddEntry(g_intLength_CONNOLLY_nuNC,       "Connolly et all, #nu NC",        "l");
  leg1->AddEntry(g_intLength_CONNOLLY_nuCC,       "Connolly et all, #nu CC",        "l");
  leg1->AddEntry(g_intLength_CONNOLLY_nubarNC,    "Connolly et all, #bar{#nu} NC",  "l");
  leg1->AddEntry(g_intLength_CONNOLLY_nubarCC,    "Connolly et all, #bar{#nu} CC",  "l");
  leg1->Draw();

  c1->Print("Compare_intLengths.png");
  c1->Print("Compare_intLengths.pdf");
  c1->Print("Compare_intLengths.eps");
  c1->Print("Compare_intLengths.root");

  
  Double_t ANITA_2_effArea_old[n_ANITA];
  Double_t ANITA_2_effArea_nowRENO[n_ANITA];
  Double_t ANITA_2_effArea_nowCONNOLLY_nuNC[n_ANITA];
  Double_t ANITA_2_effArea_nowCONNOLLY_nuCC[n_ANITA];
  Double_t ANITA_2_effArea_nowCONNOLLY_nubarNC[n_ANITA];
  Double_t ANITA_2_effArea_nowCONNOLLY_nubarCC[n_ANITA];
  for (int i=0; i<n_ANITA; i++){
    ANITA_2_effArea_old[i]                 = ANITA_2_effVolume_old[i]/intLength_RENO[i]; 
    ANITA_2_effArea_nowRENO[i]             = ANITA_2_effArea[i]*intLength_CONNOLLY_nubarNC[i]/intLength_RENO[i];
    ANITA_2_effArea_nowCONNOLLY_nuNC[i]    = ANITA_2_effArea[i]*intLength_CONNOLLY_nubarNC[i]/intLength_CONNOLLY_nuNC[i]; 
    ANITA_2_effArea_nowCONNOLLY_nuCC[i]    = ANITA_2_effArea[i]*intLength_CONNOLLY_nubarNC[i]/intLength_CONNOLLY_nuCC[i]; 
    ANITA_2_effArea_nowCONNOLLY_nubarNC[i] = ANITA_2_effArea[i]; 
    ANITA_2_effArea_nowCONNOLLY_nubarCC[i] = ANITA_2_effArea[i]*intLength_CONNOLLY_nubarNC[i]/intLength_CONNOLLY_nubarCC[i]; 
  }
  
  TGraph *g_A2icemc    = new TGraph (n_ANITA, ANITA_4_x, ANITA_2_effArea);
  TGraph *g_A2pub      = new TGraph (n_ANITA, ANITA_4_x, ANITA_2_effArea_published);
  TGraph *g_A2icemcold = new TGraph (n_ANITA, ANITA_4_x, ANITA_2_effArea_old);

  TGraph *g_A2_icemc_nowRENO             = new TGraph(n_ANITA, ANITA_4_x, ANITA_2_effArea_nowRENO);
  TGraph *g_A2_icemc_nowCONNOLLY_nuNC    = new TGraph(n_ANITA, ANITA_4_x, ANITA_2_effArea_nowCONNOLLY_nuNC);
  TGraph *g_A2_icemc_nowCONNOLLY_nuCC    = new TGraph(n_ANITA, ANITA_4_x, ANITA_2_effArea_nowCONNOLLY_nuCC);
  TGraph *g_A2_icemc_nowCONNOLLY_nubarNC = new TGraph(n_ANITA, ANITA_4_x, ANITA_2_effArea_nowCONNOLLY_nubarNC);
  TGraph *g_A2_icemc_nowCONNOLLY_nubarCC = new TGraph(n_ANITA, ANITA_4_x, ANITA_2_effArea_nowCONNOLLY_nubarCC);


  g_A2_icemc_nowRENO             ->SetLineColor(kBlack);
  g_A2_icemc_nowCONNOLLY_nuNC    ->SetLineColor(kBlue);
  g_A2_icemc_nowCONNOLLY_nuCC    ->SetLineColor(kRed);
  g_A2_icemc_nowCONNOLLY_nubarNC ->SetLineColor(kGreen);
  g_A2_icemc_nowCONNOLLY_nubarCC ->SetLineColor(kViolet);
  
  g_A2_icemc_nowRENO             ->SetLineStyle(2);
  g_A2_icemc_nowCONNOLLY_nuNC    ->SetLineStyle(2);
  g_A2_icemc_nowCONNOLLY_nuCC    ->SetLineStyle(2);
  g_A2_icemc_nowCONNOLLY_nubarNC ->SetLineStyle(2);
  g_A2_icemc_nowCONNOLLY_nubarCC ->SetLineStyle(2);
  
  g_A2_icemc_nowRENO             ->SetLineWidth(2);
  g_A2_icemc_nowCONNOLLY_nuNC    ->SetLineWidth(2);
  g_A2_icemc_nowCONNOLLY_nuCC    ->SetLineWidth(2);
  g_A2_icemc_nowCONNOLLY_nubarNC ->SetLineWidth(2);
  g_A2_icemc_nowCONNOLLY_nubarCC ->SetLineWidth(2);
  g_A2pub                        ->SetLineWidth(2);
  
  TCanvas *c2 = new TCanvas("c2");
  c2->SetLogy();
  g_A2_icemc_nowRENO->SetTitle("ANITA-2 acceptance for different cross-section parametrizations;Log10(E_{#nu});Acceptance [km^{2} str]");
  g_A2_icemc_nowRENO->GetYaxis()->SetRangeUser(0, 4000);
  g_A2_icemc_nowRENO->Draw("Al");
  g_A2_icemc_nowCONNOLLY_nuNC->Draw("l");
  g_A2_icemc_nowCONNOLLY_nuCC->Draw("l");
  g_A2_icemc_nowCONNOLLY_nubarNC->Draw("l");
  g_A2_icemc_nowCONNOLLY_nubarCC->Draw("l");
  g_A2pub->Draw("l");
  
  TLegend *leg2 = new TLegend (0.11, 0.65, 0.45, 0.89);
  leg2->AddEntry(g_A2_icemc_nowRENO,                "Reno et all, #nu/#bar{#nu} CC",  "l");
  leg2->AddEntry(g_A2_icemc_nowCONNOLLY_nuNC,       "Connolly et all, #nu NC",        "l");
  leg2->AddEntry(g_A2_icemc_nowCONNOLLY_nuCC,       "Connolly et all, #nu CC",        "l");
  leg2->AddEntry(g_A2_icemc_nowCONNOLLY_nubarNC,    "Connolly et all, #bar{#nu} NC",  "l");
  leg2->AddEntry(g_A2_icemc_nowCONNOLLY_nubarCC,    "Connolly et all, #bar{#nu} CC",  "l");
  leg2->AddEntry(g_A2pub,       "ANITA-2 published",      "l");
  leg2->Draw();

  
  c2->Print("Compare_A2_acceptance_xsec.png");
  c2->Print("Compare_A2_acceptance_xsec.pdf");
  c2->Print("Compare_A2_acceptance_xsec.C");
  c2->Print("Compare_A2_acceptance_xsec.root");
  
  TCanvas *c = new TCanvas("c");
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();
  pad1->SetLogy();

  g_A2pub->SetTitle("ANITA-2 acceptance;;Acceptance [km^{2} str]");
  g_A2icemc->SetLineColor(kBlue);
  g_A2icemcold->SetLineColor(kRed);
  g_A2pub->SetLineColor(kBlack);
  g_A2pub->GetYaxis()->SetTitleSize(20);
  g_A2pub->GetYaxis()->SetTitleFont(43);
  g_A2pub->GetYaxis()->SetTitleOffset(0.9);
  g_A2pub->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  g_A2pub->GetYaxis()->SetLabelSize(15);
  
  g_A2pub->Draw("Al");
  g_A2icemc->Draw("l");
  g_A2icemcold->Draw("l");
  
  TLegend *leg = new TLegend (0.5, 0.11, 0.89, 0.3);
  leg->AddEntry(g_A2pub,       "ANITA-2 published",      "l");
  leg->AddEntry(g_A2icemcold,  "ANITA-2 icemc 2011",     "l");
  leg->AddEntry(g_A2icemc,     "ANITA-2 icemc NOW",      "l");
  leg->Draw();

  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();

  double ratio1[n_ANITA];
  for (int i=0; i<n_ANITA; i++) ratio1[i] = ANITA_2_effArea_published[i]/ANITA_2_effArea_old[i];
  double ratio2[n_ANITA];
  for (int i=0; i<n_ANITA; i++) ratio2[i] = ANITA_2_effArea_old[i]/ANITA_2_effArea[i];

  TGraph *gratio1 = new TGraph(n_ANITA, ANITA_4_x, ratio1);
  TGraph *gratio2 = new TGraph(n_ANITA, ANITA_4_x, ratio2);

  gratio1->SetLineStyle(2);
  gratio2->SetLineStyle(2);

  gratio1->SetLineColor(kRed);
  gratio2->SetLineColor(kViolet);

  gratio1->GetYaxis()->SetRangeUser(0, 7);
  gratio1->SetTitle(";Energy exponent;Ratio");
  gratio1->GetYaxis()->SetTitleSize(20);
  gratio1->GetYaxis()->SetTitleFont(43);
  gratio1->GetYaxis()->SetTitleOffset(0.9);
  gratio1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  gratio1->GetYaxis()->SetLabelSize(15);
  
  // X axis ratio plot settings
  gratio1->GetXaxis()->SetTitleSize(0.14);
  gratio1->GetXaxis()->SetTitleOffset(0.7);
  gratio1->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  gratio1->GetXaxis()->SetLabelSize(15);
   
  gratio1->Draw("Al");
  gratio2->Draw("l");
  
  TLegend *leg3 = new TLegend (0.11, 0.75, 0.4, 0.99);
  leg3->AddEntry(gratio1,       "ANITA-2 published / icemc 2011", "l");
  leg3->AddEntry(gratio2,       "icemc 2011 / icemc NOW ", "l");
  leg3->Draw();
  
  c->Print("Compare_A2_acceptance.png");
  c->Print("Compare_A2_acceptance.pdf");
  c->Print("Compare_A2_acceptance.C");
  c->Print("Compare_A2_acceptance.root");
}
