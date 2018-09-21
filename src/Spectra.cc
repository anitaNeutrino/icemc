#include "Spectra.h"
#include "Tools.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"

#include "Report.h"
#include "Settings.h"

double icemc::Source::Spectra::pickNeutrinoEnergy() {
  if(fSettings->USEDARTBOARD){
    return GetNuEnergy();
  }
  else {
    return GetCDFEnergy();
  }
}

icemc::Source::Spectra::Spectra(const Settings* settings) : fSettings(settings){
  
  EXPONENT=fSettings->EXPONENT;
  E_bin = 12;

  double Emuons[E_bin]; // E dN/dE/dA/dt for neutrinos that are produced as muon neutrinos or muon antineutrinos.
  double Eelectrons[E_bin];// E dN/dE/dA/dt for neutrinos that are produced as electron neutrinos or muon antineutrinos.
  
  for (int i=0;i<E_bin;i++) {
    const double EMin = 18.0;
    const double EMax = 22.0;
    double dE = (EMax-EMin)/(E_bin-1);
    energy[i]=EMin+((double)i)*dE; // initial E_bin = 50. energy = 18 ~ 22 eV
    Emuons[i]=-30.;
    Eelectrons[i]=-30.;
  } //for

  // end of initialization!!
  if (EXPONENT==1)  // dNdEdAdt ~ E^-1
    {
      E_bin = 12;
      for (int i=0;i<E_bin;i++) {   // set energy, EdNdEdAdt and E2dNdEdAdt
	energy[i] = 16.+((double)i)/2.;   // in log, in eV
	EdNdEdAdt[i] = 0.;    // in log
	E2dNdEdAdt[i] = EdNdEdAdt[i] + (energy[i] - 9.);    // in log, in GeV
      }
    }

  else if (EXPONENT==2) // dNdEdAdt ~ E^-2
    {
      E_bin = 12;
      for (int i=0;i<E_bin;i++) {   // set energy, EdNdEdAdt and E2dNdEdAdt
	energy[i] = 16.+((double)i)/2.;   // in log, in eV
	E2dNdEdAdt[i] = 0.;    // in log, in GeV
	EdNdEdAdt[i] = E2dNdEdAdt[i] - (energy[i] - 9.);    // in log
      }
    }

  else if (EXPONENT==3) // dNdEdAdt ~ E^-3
    {
      E_bin = 12;
      for (int i=0;i<E_bin;i++) {   // set energy, EdNdEdAdt and E2dNdEdAdt
	energy[i] = 16.+((double)i)/2.;   // in log, in eV
	E2dNdEdAdt[i] = -(energy[i] - 9.);    // in log, in GeV
	EdNdEdAdt[i] = E2dNdEdAdt[i] - (energy[i] - 9.);    // in log
      }
    }

  else if (EXPONENT==4) // dNdEdAdt ~ E^-4
    {
      E_bin = 12;
      for (int i=0;i<E_bin;i++) {   // set energy, EdNdEdAdt and E2dNdEdAdt
	energy[i] = 16.+((double)i)/2.;   // in log, in eV
	E2dNdEdAdt[i] = -2. * (energy[i] - 9.);    // in log, in GeV
	EdNdEdAdt[i] = E2dNdEdAdt[i] - (energy[i] - 9.);    // in log
      }
    }

  else if (EXPONENT==30) // ESS baseline model. Used to be EXPONENT "0"
    {
      E_bin = 9;

      // electron component of Figure 4 of ES&S
      // astro-ph/0101216
      Eelectrons[0]=-17.2; // 16.
      Eelectrons[1]=-17.35; // 16.5
      Eelectrons[2]=-17.2; // 17.
      Eelectrons[3]=-17.1; // 17.5
      Eelectrons[4]=-17.2; // 18.
      Eelectrons[5]=-17.5; // 18.5
      Eelectrons[6]=-18.0; // 19
      Eelectrons[7]=-18.5; // 19.5
      Eelectrons[8]=-19.4; // 20.
      Eelectrons[9]=-30.; // 20.5 punt------ not using above here
      Eelectrons[10]=-30.; // 21.0 punt
      Eelectrons[11]=-30.; // 21.5 punt
  
      // lower curve of Figure 9 of ES&S
      // astro-ph/0101216
      Emuons[0]=-17.1;  //16.
      Emuons[1]=-16.6;  //16.5
      Emuons[2]=-16.3;  //17.
      Emuons[3]=-16.2; // 17.5
      Emuons[4]=-16.4; // 18.
      Emuons[5]=-16.7; // 18.5
      Emuons[6]=-17.3; // 19
      Emuons[7]=-17.95; // 19.5
      Emuons[8]=-18.85; // 20.
      Emuons[9]=-19.9; // 20.5 punt------not using above here
      Emuons[10]=-30.; // 21.0 punt
      Emuons[11]=-30.; // 21.5 punt
  
      for (int i=0;i<E_bin;i++) {
        energy[i] = 16.+((double)i)/2.;   // in log, in eV
        EdNdEdAdt[i] = log10( pow(10., Eelectrons[i]) + pow(10., Emuons[i]) );   // in log
        // I know that it have to change this to non-log later but I want to make them same way.
        E2dNdEdAdt[i] = EdNdEdAdt[i] + (energy[i] - 9.);    // in log, in GeV
      }
    }

  else if (EXPONENT==31) // ESS-cosmological constant. Use Nu_el : Nu_mu ratio. used to be EXPONENT "5"
    {
      E_bin = 9;
      double emuratio[E_bin];

      // electron component of Figure 4 of ES&S
      // astro-ph/0101216
      Eelectrons[0]=-17.2; // 16.
      Eelectrons[1]=-17.35; // 16.5
      Eelectrons[2]=-17.2; // 17.
      Eelectrons[3]=-17.1; // 17.5
      Eelectrons[4]=-17.2; // 18.
      Eelectrons[5]=-17.5; // 18.5
      Eelectrons[6]=-18.0; // 19
      Eelectrons[7]=-18.5; // 19.5
      Eelectrons[8]=-19.4; // 20.
      Eelectrons[9]=-30.; // 20.5 punt------ not using above here
      Eelectrons[10]=-30.; // 21.0 punt
      Eelectrons[11]=-30.; // 21.5 punt
    
      // muon component of Figure 4 of ES&S
      // astro-ph/0101216
      Emuons[0]=-17.8; // 16.
      Emuons[1]=-17.4; // 16.5
      Emuons[2]=-17.; // 17.
      Emuons[3]=-16.75; // 17.5
      Emuons[4]=-16.9; // 18.
      Emuons[5]=-17.2; // 18.5
      Emuons[6]=-17.7; // 19
      Emuons[7]=-18.3; // 19.5
      Emuons[8]=-19.1; // 20.
      Emuons[9]=-30.; // 20.5 punt
      Emuons[10]=-30.; // 21.0 punt
      Emuons[11]=-30.; // 21.5 punt
    
      for(int i=0;i<E_bin;i++)
	emuratio[i]=Eelectrons[i]/Emuons[i];
    
      // upper curve in Figure 9 of ES&S
      // astro-ph/0101216
      Emuons[0]=-16.85;  //16.
      Emuons[1]=-16.4;  //16.5
      Emuons[2]=-16.05;  //17.
      Emuons[3]=-16.; // 17.5
      Emuons[4]=-16.15; // 18.
      Emuons[5]=-16.5; // 18.5
      Emuons[6]=-17.1; // 19
      Emuons[7]=-17.7; // 19.5
      Emuons[8]=-18.65; // 20.
      Emuons[9]=-19.75; // 20.5 punt
      Emuons[10]=-30.; // 21.0 punt
      Emuons[11]=-30.; // 21.5 punt
    
    
      for(int i=0;i<E_bin;i++) {
	Eelectrons[i]=Emuons[i]*emuratio[i];
	std::cout << "Eelectrons, Emuons are " << Eelectrons[i] << " " << Emuons[i] << "\n";
      }
  
      for (int i=0;i<E_bin;i++) {
	energy[i] = 16.+((double)i)/2.;   // in log, in eV
	EdNdEdAdt[i]=log10( pow(10.,Eelectrons[i]) + pow(10.,Emuons[i]) );  // in log.
	E2dNdEdAdt[i] = EdNdEdAdt[i] + (energy[i] - 9.);    // in log, in GeV
	std::cout << "EdNdEdAdt are " << EdNdEdAdt[i] << "\n";
      }
    
    } // end if ESS-cosmological constant

  else if (EXPONENT>31 && EXPONENT<200)  // use digitized flux from different models
    {
      switch (EXPONENT)
	{
	case 32:  // ESS3
	  GetFlux("essfig9.dat");
	  break;
	case 33:  // ESS4
	  GetFlux("essbaseline.dat");
	  break;
	case 34:  // ESS5
	  GetFlux("ess_n0.dat");
	  break;
	case 40:
	  GetFlux("ahlers.dat");
	  break;
	case 41:
	  GetFlux("ahlers2012.dat");
	  break;
	case 50:
	  GetFlux("allard.dat");
	  break;
	case 60:
	  GetFlux("ave_max.dat");
	  break;
	case 61:
	  GetFlux("ave_min.dat");
	  break;
	case 70:
	  GetFlux("kkss_envo.dat");
	  break;
	case 80:
	  GetFlux("gzk_peter.dat");
	  break;
	case 90:
	  GetFlux("waxgzk.dat");
	  break;
	case 100:
	  GetFlux("e-2.dat");
	  break;
	case 101:
	  GetFlux("e-2_icecube.dat");
	  break;
	case 110:
	  GetFlux("yuksel_grb.dat");
	  break;
	case 111:
	  GetFlux("yuksel_qso.dat");
	  break;
	case 112:
	  GetFlux("yuksel_sfh.dat");
	  break;
	  //          case 100:
	  //              GetFlux("berezinsky_saturate.dat");
	  //              break;
	default:
	  std::cout<<"Error: Wrong input of EXPONENT!"<<std::endl;
	}
    }

  else if (EXPONENT==200)   // Iron model
    {
      GetFlux("Ave2005_Fe_Emax21.0.dat");
    }

  else if (EXPONENT==201)
    {
      GetFlux("Ave2005_Fe_Emax21.5.dat");
    }
  
  else if (EXPONENT==202)
    {
      GetFlux("Ave2005_Fe_Emax22.0.dat");
    }

  else if (EXPONENT==203)
    {
      GetFlux("Ave2005_Fe_hi_evo.dat");
    }

  else if (EXPONENT==204)
    {
      GetFlux("Ave2005_Fe_low_evo.dat");
    }

  else if (EXPONENT==210)
    {
      GetFlux("Stanev2008_heavy.dat");
    }

  else if (EXPONENT==220)
    {
      GetFlux("Kotera2010_Fe_only.dat");
    }

  else if (EXPONENT==221)
    {
      GetFlux("Kotera2010_Fe_rich.dat");
    }

  else if (EXPONENT==222)
    {
      GetFlux("Kotera2010_mix_max.dat");
    }

  else if (EXPONENT==223)
    {
      GetFlux("Kotera2010_mix_min.dat");
    }
  else if (EXPONENT==224)
    {
      GetFlux("Kotera2010_proton.dat");
    }
  

  if(EXPONENT <6 || EXPONENT >29){//is spectrum
    GetCDF();
  }
  //
  // End of selecting the Model!!!
  //

  /*
    if (EXPONENT<=10. || EXPONENT>=100.) {
    gspectrum[(int)EXPONENT]=new TGraph(12,energy,EdNdEdAdt);
    //    maxflux=gspectrum[(int)EXPONENT]->GetMaximum();
    maxflux=Tools::dMax(EdNdEdAdt,12);

    for (int i=0;i<12;i++) {
    E2dNdEdAdt[i]=log10(EdNdEdAdt[i])+energy[i]-9.;
    }
    
    spectrum=new TSpline5("spectrum",energy,EdNdEdAdt,12);
    
    }
    else
    // find the max so we can normalise it
    maxflux=Tools::dMax(EdNdEdAdt,12);
  */


  // From log to linear!!  
  for (int i=0;i<E_bin;i++) {
    EdNdEdAdt[i] = pow(10, EdNdEdAdt[i]);     //to linear 
    E2dNdEdAdt[i] = pow(10, E2dNdEdAdt[i]);   //to linear
    //      E2dNdEdAdt[i]=log10(EdNdEdAdt[i])+energy[i]-9.;
  }

  gEdNdEdAdt = new TGraph(E_bin, energy, EdNdEdAdt);
  gE2dNdEdAdt = new TGraph(E_bin, energy, E2dNdEdAdt);

  sEdNdEdAdt = new TSpline3("sEdNdEdAdt", gEdNdEdAdt);
  sE2dNdEdAdt = new TSpline3("sE2dNdEdAdt", gE2dNdEdAdt);

  maxflux = Tools::dMax(EdNdEdAdt,E_bin);
}


double  icemc::Source::Spectra::GetNuEnergy() {
  double thisenergy=16.; // arbitrary initialisation
  double thisflux=2.; // initialise higher than max
  double max=1.;
  int energybin=0; // arbitrary initialisation
  double maxenergy=Tools::dMax(energy,E_bin);
  double minenergy=Tools::dMin(energy,E_bin);
  // this uses the dartboard approach
  //cout << "minenergy, maxenergy are " << minenergy << " " << maxenergy << "\n";
  while(thisflux>max) {
    // pick an energy  
    // thisenergy=Rand3.Rndm()*(maxenergy-minenergy)+minenergy; // pick energy at random between the highest and lowest
    thisenergy = pickUniform(minenergy, maxenergy);
    // the energy array is actually filled with energy exponents 
    // and thisenergy starts from 0 so it has an offset
    
    energybin=Tools::Getifreq(thisenergy,minenergy,maxenergy,E_bin);
    //max=gspectrum[(int)EXPONENT]->Eval(thisenergy,0,"S")/maxflux; // this is the maximum the normalized flux can be in this bin, always less than 1
    max=EdNdEdAdt[energybin]/maxflux;
    thisflux=pickUniform(); //Rand3.Rndm(); // pick the flux at random between 0 and 1, if it's less than max it's a winner
  } //while
  return pow(10.,thisenergy);
} //Pick Neutrino Energy



void icemc::Source::Spectra::GetCDF(){//set up CDF and inverse CDF;
  std::cout<<"in CDF " << EXPONENT << "\t" << E_bin << std::endl;;
  for(auto e :  energy){
    std::cout << e << ",  ";
  }
  std::cout << std::endl;
  double y_val=0.;
  double E_min = 18;//energy[0];
 
  double E_max = energy[E_bin-1];
  if(E_max > 21) E_max=21;
  double step_size =.25;//in logE
  int n =(int) floor((E_max-E_min)/step_size);

  std::cout<< n << "\t" << E_bin << "\t" << E_max << "\t" << E_min << "\t" << step_size << std::endl;

  std::vector<double> E(n);
  std::vector<double> N(n);
  double E_tmp=0.;
  double integral=0.;

  TGraph *hEdNdE = new TGraph(E_bin,energy,EdNdEdAdt);
  //cout<<"E_min, Max, n are "<<E_min<<" "<<E_max<<" "<<n<<"\n";
  for(int i=0;i<n;i++){
    E_tmp = E_min+i*step_size;
    y_val=hEdNdE->Eval(E_tmp,0,"s");//get interpolated value
    
    integral +=y_val*step_size;//integrate in log space
    
    E.at(i)=E_tmp;
    N.at(i)=integral;
  }
  
  for(int i=0;i<n;i++){
    N.at(i)=N.at(i)/integral;    
  }

  CDF = new TGraph(n,&E[0],&N[0]);
  inverse_CDF = new TGraph(n,&N[0],&E[0]);
}

double icemc::Source::Spectra::GetCDFEnergy(){//get Energy from 'CDF'

  if(!CDF){
    GetCDF(); // do initialization if needed
  }

  double ran = gRandom->Rndm();
 
  double thisenergy=0.;

  thisenergy = inverse_CDF->Eval(ran,0,"S");
  //cout<<"ran is "<<ran<<" thisenergy is "<<thisenergy<<"\n";

  while(thisenergy <18){//redundant?
    //cout<<"thisenergy was "<<thisenergy<<" ran was "<<ran<<"\n";
    ran = gRandom->Rndm();
    thisenergy = inverse_CDF->Eval(ran,0,"S");
    //cout<<"ran is "<<ran<<" thisenergy is "<<thisenergy<<"\n";
  }
  
  return pow(10.,thisenergy);

  
}

inline void icemc::Source::Spectra::GetFlux(std::string filename)
{
  
  const std::string ICEMC_SRC_DIR=std::getenv("ICEMC_SRC_DIR");
  std::ifstream influx((ICEMC_SRC_DIR+"/fluxes/"+filename).c_str());
  int NLINES;
  influx >> NLINES;   // Read how much lines in the file.
  std::cout<<"We are using "<<filename.c_str()<<" as the flux data."<<std::endl;
  std::cout<<"total lines in the file are "<<NLINES<<std::endl;
  E_bin = NLINES;
    
  //    double flux[NLINES];//E2dNdE GeV cm^-2 s^-1 sr^-1
    
  for(int i=0;i<NLINES;i++) {
    influx>>energy[i]>>E2dNdEdAdt[i];
  }
    
  for(int i=0;i<E_bin;i++) {
    //        EdNdEdAdt[i]=pow(10.,(flux[i]+9.-energy[i]));
    EdNdEdAdt[i] = E2dNdEdAdt[i] + 9. - energy[i];  // change from GeV to eV and E2dN -> EdN
  }
  // maxflux=Tools::dMax(EdNdEdAdt,12);
  // cout<<maxflux<<endl;
}


TGraph *icemc::Source::Spectra::GetGEdNdEdAdt() {
  return gEdNdEdAdt;
}


TGraph *icemc::Source::Spectra::GetGE2dNdEdAdt() {
  return gE2dNdEdAdt;
}


TSpline3 *icemc::Source::Spectra::GetSEdNdEdAdt() {
  return sEdNdEdAdt;
}


TSpline3 *icemc::Source::Spectra::GetSE2dNdEdAdt() {
  return sE2dNdEdAdt;
}


double *icemc::Source::Spectra::Getenergy() {
  return energy;
}


double *icemc::Source::Spectra::GetEdNdEdAdt() {
  return EdNdEdAdt;
}


double *icemc::Source::Spectra::GetE2dNdEdAdt() {
  return E2dNdEdAdt;
}


double icemc::Source::Spectra::GetEdNdEdAdt(double E_val) {
  double tmp_Get;
  if (E_val < energy[0]) {
    std::cout<<"Energy value is smaller than the energy boundary!\n";
    std::cout<<"Energy value is replaced to minimum value of energy bound : "<<energy[0]<<"\n";
    tmp_Get = sEdNdEdAdt->Eval(energy[0]);
  }
  else if (E_val > energy[E_bin-1]) {
    std::cout<<"Energy value is bigger than the energy boundary!\n";
    std::cout<<"Energy value is replaced to maximum value of energy bound : "<<energy[E_bin-1]<<"\n";
    tmp_Get = sEdNdEdAdt->Eval(energy[E_bin-1]);
  }
  else {
    tmp_Get = sEdNdEdAdt->Eval(E_val);
  }
  return tmp_Get;
}


double icemc::Source::Spectra::GetE2dNdEdAdt(double E_val) {
  double tmp_Get;
  if (E_val < energy[0]) {
    std::cout<<"Energy value is smaller than the energy boundary!\n";
    std::cout<<"Energy value is replaced to minimum value of energy bound : "<<energy[0]<<"\n";
    tmp_Get = sE2dNdEdAdt->Eval(energy[0]);
  }
  else if (E_val > energy[E_bin-1]) {
    std::cout<<"Energy value is bigger than the energy boundary!\n";
    std::cout<<"Energy value is replaced to maximum value of energy bound : "<<energy[E_bin-1]<<"\n";
    tmp_Get = sE2dNdEdAdt->Eval(energy[E_bin-1]);
  }
  else {
    tmp_Get = sE2dNdEdAdt->Eval(E_val);
  }
  return tmp_Get;
}


double icemc::Source::Spectra::Getmaxflux() {
  return maxflux;
}


int icemc::Source::Spectra::GetE_bin() {
  return E_bin;
}


int icemc::Source::Spectra::IsSpectrum() {
  int out;
  if (EXPONENT<=10||EXPONENT>=30) {
    out = 1;
  }
  else {
    out = 0;
  }
  return out;
}


int icemc::Source::Spectra::IsMonoenergetic() {
  int out;
  if (EXPONENT>10&&EXPONENT<30) {
    out = 1;
  }
  else {
    out = 0;
  }
  return out;
}


void icemc::Source::Spectra::savePlots2(const TString& fileName){
  TCanvas *ctest1 = new TCanvas("ctest1", "", 880, 800);

  GetGEdNdEdAdt()->Draw("al");
  GetGEdNdEdAdt()->GetHistogram()->SetMinimum(-0.2*Getmaxflux());
  GetSEdNdEdAdt()->SetLineColor(2);
  GetSEdNdEdAdt()->Draw("l same");
 
  ctest1->Print(fileName);
  delete ctest1;

}

void icemc::Source::Spectra::savePlots(const TString& fileNameNoSuffix){

  TCanvas *ctemp = new TCanvas(fileNameNoSuffix);

  TFile *out = new TFile(fileNameNoSuffix + ".root", "recreate");
  TGraph *g1 = GetGEdNdEdAdt();
  g1->Draw("Al");
  ctemp->Print(fileNameNoSuffix + "1.png");
  int n = g1->GetN();
  double *x = g1->GetX();
  double x2[20];
  for (int i=0;i<n;i++) x2[i] = TMath::Power(10., x[i]);
  TGraph *g2 = new TGraph(n, x2, g1->GetY());
  g2->Draw("Al");
  ctemp->SetLogy();  
  ctemp->SetLogx();
  std::cout << g2->Integral() << " " << g2->Integral(1, n) << std::endl;  
  ctemp->Print(fileNameNoSuffix + "2.png");
  g1->Write();
  g2->Write();

  out->Close();
  delete out;
  delete ctemp;
}
