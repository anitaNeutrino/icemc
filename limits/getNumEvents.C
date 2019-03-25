#include "Acceptances.h"

Double_t ANITA_livetime = ANITA_3_livetime;

Double_t ANITA_effArea[n_ANITA];
Double_t ANITA_efficiency[n_ANITA];

double deltaLog=0.01;
double Emin = 18.;
double Emax = 21.;
int nBins = floor( (Emax - Emin) / deltaLog);
TGraph *gFluence;
TGraph *gEffArea;
TH1D *hRate;
TH1D *hNumObs;

Double_t EmaxModel=30;


void GetFlux(string name);
string GetFluxFromNumber(int EXPONENT);
void GetEffArea();
void GetRate();
void GetNumObs();

void getNumEvents(){

  int expMin = 30;
  int expMax = 250;

  // int expMin = 200;
  // int expMax = expMin+1;

  bool debug=true;
  TCanvas *c;

  if (debug){
    c = new TCanvas("c");
  }
  
  //  for (int iexp=expMin;  iexp<expMax; iexp++){
  for (int iexp=101;  iexp<=101; iexp++){
    
    string fluxname = GetFluxFromNumber(iexp);
    if (fluxname=="") continue;
    GetFlux(fluxname);
    GetEffArea();
    GetRate();
    GetNumObs();

    if (debug){
      c->SetLogy();
      gFluence->Draw("Al");
      Simulation();
      c->Print("Fluence.png");
      gEffArea->Draw("Al");
      Simulation();
      c->Print("Acceptance.png");
      c->SetLogy(0);
      hRate->Draw("histo");
      Simulation();
      c->Print("Rate.png");
      hNumObs->Draw("histo");
      Simulation();
      c->Print("NumObs.png");
    }
    
    //  cout << "Expected number of events is " << hNumObs->Integral() << endl;
    char buffer[100];
    snprintf(buffer, sizeof(buffer), "%s:", fluxname.c_str());
    double integral, error;
    integral = hNumObs->IntegralAndError(1, nBins, error, "");
    if (integral>1e-4){
      printf("%-*s %8.5f +/- %8.5f \n", 25, buffer, integral, error);
    } else {
      printf("%-*s %3.2e +/- %3.2e \n", 25, buffer, integral, error);
    }
    if (!debug){
      delete hRate;
      delete hNumObs;
    }
  }
}


void GetNumObs(){
  hNumObs = new TH1D("hNumObs", "", nBins, Emin, Emax);

  double tempRate   = 0;
  double tempEnergy = 0;
  for (int i=1; i<nBins+1; i++){
    tempEnergy = Emin + deltaLog*(i+0.5);
    if (tempEnergy>EmaxModel) break;
    tempRate = hRate->GetBinContent(i);
    //    cout << i << " " << tempEnergy << " " << tempRate << endl;
    hNumObs->Fill(tempEnergy, tempRate*ANITA_livetime);
  }

  hNumObs->SetTitle(";log_{10} #left(#frac{E_{#nu}}{eV}#right);Number per Bin");
}

void GetRate(){
  hRate = new TH1D("hRate", "", nBins, Emin, Emax);

  double tempRate   = 0;
  double tempEnergy = 0;
  for (int i=0; i<nBins+1; i++){
    tempEnergy = Emin + deltaLog*(i+0.5);
    if (tempEnergy>EmaxModel) break;
    tempRate = TMath::Log(10.)*deltaLog*gEffArea->Eval(tempEnergy)*gFluence->Eval(tempEnergy);
    // cout << i << " " << tempEnergy << " " << tempRate << endl;
    hRate->Fill(tempEnergy, tempRate);
  }

  hRate->SetTitle(";log_{10} #left(#frac{E_{#nu}}{eV}#right);Rate per bin, s^{-1}");
}

void GetEffArea(){
  
  for (int i=0; i<n_ANITA; i++){
    
    ANITA_effArea[i] = ANITA_3_effVol[i]/intLength_CONNOLLY_nuCC[i]; 
    ANITA_effArea[i] = TMath::Sqrt(ANITA_effArea[i]*ANITA_2_effArea_Peter[i]*ANITA_effArea[i]/ANITA_2_effArea_icemc2010[i])*ANITA_3_eff[i];
    // ANITA_effArea[i] = ANITA_1_effArea[i]*ANITA_1_eff[i];
    cout << ANITA_effArea[i] << endl;
    ANITA_effArea[i] *= 1e10; // from km^2 to cm^2
  }
  
  gEffArea = new TGraph (n_ANITA, ANITA_x, ANITA_effArea);
  gEffArea->SetTitle(";log_{10} #left(#frac{E_{#nu}}{eV}#right);[A#Omega], cm^{2} str");
}


string GetFluxFromNumber(int EXPONENT){

  if (EXPONENT<113) {
    switch (EXPONENT)
      {
      case 32:  // ESS3
	return ("essfig9.dat");
	break;
      case 33:  // ESS4
	return ("essbaseline.dat");
	break;
      case 34:  // ESS5
	return ("ess_n0.dat");
	break;
      case 40:
	return ("ahlers.dat");
	break;
      case 41:
	return ("ahlers2012.dat");
	break;
      case 50:
	return ("allard.dat");
	break;
      case 60:
	return ("ave_max.dat");
	break;
      case 61:
	return ("ave_min.dat");
	break;
      case 70:
	return ("kkss_envo.dat");
	break;
      case 80:
	return ("gzk_peter.dat");
	break;
      case 90:
	return ("waxgzk.dat");
	break;
      case 100:
	return ("e-2.dat");
	break;
      case 101:
	return ("e-2_icecube.dat");
	break;
      case 110:
	return ("yuksel_grb.dat");
	break;
      case 111:
	return ("yuksel_qso.dat");
	break;
      case 112:
	return ("yuksel_sfh.dat");
	break;

      default:
	return "";
      }
  
  } else if (EXPONENT==200)   // Iron model
    {
      return ("Ave2005_Fe_Emax21.0.dat");
    }

  else if (EXPONENT==201)
    {
      return ("Ave2005_Fe_Emax21.5.dat");
    }
  
  else if (EXPONENT==202)
    {
      return ("Ave2005_Fe_Emax22.0.dat");
    }

  else if (EXPONENT==203)
    {
      return ("Ave2005_Fe_hi_evo.dat");
    }

  else if (EXPONENT==204)
    {
      return ("Ave2005_Fe_low_evo.dat");
    }

  else if (EXPONENT==210)
    {
      return ("Stanev2008_heavy.dat");
    }

  else if (EXPONENT==220)
    {
      return ("Kotera2010_Fe_only.dat");
    }

  else if (EXPONENT==221)
    {
      return ("Kotera2010_Fe_rich.dat");
    }

  else if (EXPONENT==222)
    {
      return ("Kotera2010_mix_max.dat");
    }

  else if (EXPONENT==223)
    {
      return ("Kotera2010_mix_min.dat");
    }
  
  else if (EXPONENT==224)
    {
      return ("Kotera2010_proton.dat");
    }
  // Recent papers, optimistic flux spectra
  else if (EXPONENT==300)
    {
      return ("heinze2018C.dat");
    }
  else if (EXPONENT==301)
    {
      return ("vliet2018C.dat");
    }
  else if (EXPONENT==302)
    {
      return ("biehlTDE2018C.dat");
    }
  else if (EXPONENT==303)
    {
      return ("fang2017C.dat");
    }
  else if (EXPONENT==304)
    {
      return ("fangPulsar2016C.dat");
    }
  else if (EXPONENT==305)
    {
      return ("bonciolo2019C.dat");
    }
  else if (EXPONENT==306)
    {
      return ("bonciolo2019cosmoMuC.dat");
    }
  



  return "";

}


void GetFlux(string filename){

  ifstream influx(("../fluxes/"+filename).c_str());
  int NLINES;
  influx >> NLINES;   // Read how much lines in the file.
  // cout<<"We are using "<<filename.c_str()<<" as the flux data."<<endl;
  // cout<<"total lines in the file are "<<NLINES<<endl;
 
  double E2dNdEdAdt[100]; //E2dNdE GeV cm^-2 s^-1 sr^-1
  double EdNdEdAdt[100]; //E2dNdE GeV cm^-2 s^-1 sr^-1
  double energy[100];
  
  for(int i=0;i<NLINES;i++) {
    influx>>energy[i]>>E2dNdEdAdt[i];
  }
    
  for(int i=0;i<NLINES;i++) {
    //        EdNdEdAdt[i]=pow(10.,(flux[i]+9.-energy[i]));
    EdNdEdAdt[i] = E2dNdEdAdt[i] + 9. - energy[i];  // change from GeV to eV and E2dN -> EdN
    EdNdEdAdt[i] = TMath::Power(10,  EdNdEdAdt[i]);

  }

  EmaxModel = energy[NLINES-1];
  
  gFluence = new TGraph(NLINES, energy, EdNdEdAdt);
  gFluence->SetTitle(";log_{10} #left(#frac{E_{#nu}}{eV}#right);EF(E), cm^{-2} s^{-1} str^{-1}");
}
