
void plotGraphs(const char* fileName, const double yThresh){

  TFile* f = TFile::Open(fileName);

  int numGraph = 0;
  // f->GetListOfKeys()->Print();
  TIter next(f->GetListOfKeys());
  TCanvas* c = NULL;

  TKey* k = (TKey*) next();
  TGraph* grFirst = NULL;
  while(k){
    TGraph* gr = (TGraph*) f->Get(k->GetName());

    bool reasonableSize = false;
    for(int pt=0; pt < gr->GetN(); pt++){
      if(fabs(gr->GetY()[pt]) > yThresh){
	reasonableSize = true;
	break;
      }
    }

    if(reasonableSize){
      TString opt = numGraph > 0 ? "lsame" : "al";
      if(!numGraph){
	grFirst = gr;
	c = new TCanvas();
      }
      numGraph++;
      
      gr->SetLineColor(numGraph);
      gr->SetFillColor(0);
      gr->SetTitle(gr->GetName());
      gr->Draw(opt);
    }

    k = (TKey*) next();
  }
  if(c && grFirst){
    c->BuildLegend();
    grFirst->SetTitle(TString::Format("There are %d graphs with fabs(y) > %.2e", numGraph, yThresh));
  }
}

void plotQuick(){

  TFile* f1 = TFile::Open("testRX.root");
  TGraph* gr1_f = (TGraph*) f1->Get("gr_PSD_RX_0_2_afterGain");
  gr1_f->SetLineColor(kBlack);
  gr1_f->SetFillColor(0);
  gr1_f->SetTitle("ANITA");
  gr1_f->Draw("al");
  
  TFile* f2 = TFile::Open("oldChanTrigger.root");
  TGraph* gr2_f = (TGraph*) f2->Get("grF_pol0_ant2_iband4_aftergain");
  gr2_f->SetLineColor(kRed);
  gr2_f->SetFillColor(0);  
  gr2_f->SetTitle("ChanTrigger");
  gr2_f->Draw("lsame");

  gPad->BuildLegend();

  std::cout << gr1_f->GetN() << "\t" << gr2_f->GetN() << "\t" << std::endl;
}


void plotCompare(){
  TFile* f1 = TFile::Open("oldChanTrigger.root");
  TFile* f2 = TFile::Open("fSeaveysDebug.root");

  auto c = new TCanvas();
  c->Divide(16, 3);

  const int markerStyle = 2;
  for(int i=0; i < 48; i++){
    c->cd(i+1);
    // TString name1 = TString::Format("gr_pol0_ant%d_iband4_aftergain", i);
    TString name1 = TString::Format("gr_pol1_ant%d_iband4_aftergain", i);
    TGraph* gr1 = (TGraph*) f1->Get(name1);
    gr1->Draw("alp");
    gr1->SetMarkerStyle(markerStyle);


    // TString name2 = TString::Format("grV_after_%d", i);
    TString name2 = TString::Format("grH_after_%d", i);
    TGraph* gr2 = (TGraph*) f2->Get(name2);
    gr2->SetLineColor(kRed);
    gr2->SetMarkerStyle(markerStyle);
    gr2->SetMarkerColor(kRed);
    gr2->Draw("lpsame");

    int i1 = TMath::LocMax(gr1->GetN(), gr1->GetY());
    int i2 = TMath::LocMax(gr2->GetN(), gr2->GetY());
    double tMax1 = gr1->GetX()[i1];
    double tMax2 = gr2->GetX()[i2];

    // std::cout << gr1->GetX()[1]-gr1->GetX()[0] << "\t";
    // std::cout << gr2->GetX()[1]-gr2->GetX()[0] << std::endl;

    double p1 = 0;
    double dt1 = gr1->GetX()[1] - gr1->GetX()[0];

    double p2 = 0;
    double dt2 = gr1->GetX()[1] - gr1->GetX()[0];

    for(int j=0; j < gr1->GetN(); j++){
      gr1->GetX()[j] -= tMax1;
      p1 += gr1->GetY()[j]*gr1->GetY()[j]*dt1;
    }
    for(int j=0; j < gr2->GetN(); j++){
      gr2->GetX()[j] -= tMax2;
      p2 += gr2->GetY()[j]*gr2->GetY()[j]*dt2;
    }

    if(p1 > 0){
      std::cout << (p1 > 0 ? (p1 - p2)/p1 : 0) << std::endl;

      new TCanvas();
      gr1->Draw("alp");
      gr2->Draw("lpsame");

      // for(int k=0; k < gr2->GetN(); k++){
      // 	double t = gr2->GetX()[k];
      // 	double y2 = gr2->GetY()[k];
      // 	double y1 = gr1->Eval(t);

      // 	std::cout << y2 / y1 << std::endl;
      // }
    }
  }
}



void plotRX(){
  // plotGraphs("fTest.root", 680);
  // plotGraphs("testRX.root", 6.5e-6);
  // plotGraphs("testRX.root", 6.5e-6);
  // plotGraphs("oldChanTrigger.root", 6.5e-6);
  //  plotQuick();
  plotCompare();
}



