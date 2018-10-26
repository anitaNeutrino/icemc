#include "TVector3.h"
#include "TRandom3.h"
#include "Geoid.h"
#include "ConnollyEtAl2011.h"
#include "ShowerModel.h"
#include "Antarctica.h"
#include "Tools.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "TH1F.h"
#include "Constants.h"
#include "Settings.h"
#include "TTreeIndex.h"
#include "TChain.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h" 
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
#include "TauModel.h"
#include "Interaction.h"

ClassImp(icemc::TauModel);

using std::cout;
using std::stringstream;
using std::setprecision;
using std::accumulate;
using std::max_element;
using std::partial_sum;
using std::max;

icemc::TauModel::TauModel() {
  /**For Total Tau Survival probability equation
   //n.b. not in SI units.
   ////from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
   //Equation 16  &  used in Equation 30. 
   */

  //////////////////////////|Units////|Description////////////////////////
  B0 = 1.2*pow(10.,-7);  ///| m^2/kg  |
  B1 = 0.16*pow(10.,-7); ///| m^2/kg  | }parameterization using a logarithmic dependence on energy for B,
  E0.set(pow(10.,19), Energy::Unit::eV);      ///| eV      | the tau elecromagnetic energy loss parameter.
  mT.set(1.777E9, Energy::Unit::GeV);          ///| eV      |Mass of Tau
  cT = 0.00008693;       ///| m       |Tau Decay length (86.93 microMeters)
                         ///|         |
  Mn = 1.672622E-24;     ///| g       |nucleon/ proton mass in grams,also equal to 0.938 GeV. 
  A = 1.;                ///| none    |constant that sets the total probability to unity
		
  /// these last two constants from Connolly Calc 2011, used in d_dzPsurvNu().
  mydensityvector.clear();
  myavgdensityvector.clear();
  myenergyvector.clear();
  myPsurvvector.clear();
  etaufarray.clear();
  PDFarray.clear();
	
}//

//-------------------------------------------------------------

/**
   GetTauWeight is the function that will calculate the probability that a tau neutrino will interact along its path through the earth,and the tau will survive the rest of the journey and decay in the ice. This probability is calculated for final energies from 10^15.5 to the energy of the neutrino.
*/
// double icemc::TauModel::GetTauWeight(ConnollyEtAl2011 *primary1, const Settings *settings1, const Antarctica *antarctica1,Interaction *interaction1, double pnu, int nu_nubar, double& ptauf, int& crust_entered){ // 1 or 0
// double icemc::TauModel::GetTauWeight(ConnollyEtAl2011 *primary1, const Settings *settings1, const Antarctica *antarctica1,Interaction *interaction1, double pnu, Neutrino::L leptonNumber, double& ptauf, int& crust_entered){ // 1 or 0
// double icemc::TauModel::GetTauWeight(ConnollyEtAl2011 *primary1, const Settings *settings1, const Antarctica *antarctica1,Interaction *interaction1, Energy pnu, Neutrino::L leptonNumber, Energy& ptauf, int& crust_entered){ // 1 or 0
double icemc::TauModel::GetTauWeight(CrossSectionModel *primary1, const Settings *settings1, const Antarctica *antarctica1,Interaction *interaction1, Energy pnu, Neutrino::L leptonNumber, Energy& ptauf, int& crust_entered){ // 1 or 0    
			      // int& mantle_entered, // 1 or 0
			      // int& core_entered){//add secondaries?

  TVector3 chord3;///vector from earth_in to "interaction point". Sets the path direction
  TVector3 nchord;///normalized chord3
  
  
  /** Bring in useful variables from other classes */
  Neutrino::Interaction::Current current = interaction1->current;
  const Geoid::Position earth_in = interaction1->r_in;
  double TauWeight = 0;
  const Geoid::Position r_enterice = interaction1->r_enterice;
  const Geoid::Position nuexitice = interaction1->nuexitice;
  int inu=4;
  //cout<<"inu is "<<inu<<"\n";


 ///Find the chord, its length and its unit vector.
  chord3 = nuexitice - earth_in;///posnu-earth_in; 
  double Distance=chord3.Mag();
  nchord = chord3 * (1./Distance);///normalized chord3
  
  
  Energy Etaui,Etau_final;
  double y=0;
  double yweight;
  
  double  zdistance=0; //m
  TauWeight=0;///total tau survival probability.
  double TauWeight_tmp=0;
  double prob_at_z=0; 
  
  
  double Prob_Nu_Surv;
    
  double tau_surv;
  double sigma = 0;

  primary1->getSigma(pnu,leptonNumber,current);
  double len_int_kgm2 = CrossSectionModel::getInteractionLength(sigma);
  
  double step=TMath::Min(len_int_kgm2/densities[1]/10,25.0); ///how big is the step size
  
  ///set up stuff to be used later.
  Geoid::Position posnunow;
  double avgdensity=0;
  
  // double Etau_now;//=Etau_final;
  Energy Emin(1E15, Energy::Unit::eV);
 
  int i=0;
 
  double taudensity=0;
  double dEtauidEtauf=0;
   
  int etauint =0;
  double startingz=0;
  
  double totaltaudistance=0;
  int totalnusteps=0;
  TVector3 nchord1;
  double enter_ice_mag = r_enterice.Distance(earth_in);///where the particle enters the ice.
  
  mydensityvector.clear();
  myavgdensityvector.clear();
  this->GetDensityVectors(antarctica1,interaction1,nchord,step,Distance,totalnusteps, crust_entered);///Get the density vectors this function needs
  
  
  for(double logEtau_final=log10(Emin.in(Energy::Unit::eV));logEtau_final<=log10(pnu.in(Energy::Unit::eV));logEtau_final+=.01){//integral over energy (in log space?)
    Etau_final.set(pow(10,logEtau_final), Energy::Unit::eV);
    i=0;
    int totalsteps=0;
    TauWeight_tmp=0;
    // Etau_now=Etau_final;
    double gamma = Etau_final/mT;
   
    //calculate the initial energy needed at the step so the tau will end at the correct final energy
    
    this->GetEnergyVector(Etau_final,step,totalnusteps,totalsteps,totaltaudistance, pnu);///set energy vector
    
    this->GetTauSurvVector(step,totalsteps);///set tau surv vector
    
    startingz = Distance-totaltaudistance;///Set the starting position for tau. First possible interaction point for neutrino
    if(inu==7857){
       //cout<<"total nu steps is "<<totalnusteps<<"\n";
       //cout<<"density_TVector3 is "<<mydensityvector[mydensityvector.size()-1]<<"\n";
       //  cout<<"Distance is "<<Distance<<"\n";
       //cout<<"total nu steps is "<<totalnusteps<<"\n";
       //cout<<"vector size is "<<myenergyvector.size()<<"\n";
       //cout<<"totaltaudistance is "<<totaltaudistance<<"\n";
       //cout<<"crust, mantle, core are "<<crust_entered<<","<<mantle_entered<<","<<core_entered<<"\n";
       //cout<<"*********************************************************************************************** \n";
    }
    //////////////////Integral over distance/////////////////////////////////////////
    for(zdistance = startingz; zdistance<=Distance; zdistance +=step){
      int nustep = (int)(zdistance/step);
      int tauenergystep=totalsteps-i;///step number to get the correct initial energy at that point;
           
      nchord1 =  zdistance*nchord;
      posnunow = earth_in + nchord1; ///vector pointing to the step we are currently on.
      
      avgdensity = myavgdensityvector[nustep];///average density neutrino has seen to that step
      taudensity=mydensityvector[nustep];///density of the current step
      Etaui=myenergyvector[tauenergystep];///Energy vector is filled backwards (final to initial), pull out correct energy
      tau_surv = myPsurvvector[i];///chance tau survives from creation to ice
      y=1.-Etaui/pnu;///inelasticity
      
      if(zdistance >=enter_ice_mag){
        tau_surv=1;///if neutrino interacts in ice, the tau will already be in ice.
      }
      
      if (Etaui>=pnu || Etaui<Etau_final || Etaui!=Etaui){//to catch anything that might make it through. Precaution
        prob_at_z =0;
      }
      else{
        // yweight=primary1->Getyweight(pnu,y,leptonNumber,current);
        yweight=0; ///@todo FIXME
	
        Prob_Nu_Surv = avgdensity/len_int_kgm2*exp(-zdistance*avgdensity*1./len_int_kgm2);///prob neutrino survives to this point
        dEtauidEtauf = exp(B1*taudensity*(Distance-zdistance))*(Etaui/Etau_final);///how the initial tau energy is related to final
        prob_at_z = Prob_Nu_Surv*yweight*1./pnu.in(Energy::Unit::eV)*dEtauidEtauf*step*tau_surv;///probability that tau survives to ice
	/*	if(i==100 && Etaui > 1E18){
	  cout<<"Prob_Nu_Surv is "<<Prob_Nu_Surv<<"\n";
	  cout<<"avgdensity is "<<avgdensity<<"\n";
	  cout<<"len_int_kgm2 is "<<len_int_kgm2<<"\n";
	  cout<<"zdistance is "<<zdistance<<"\n";
	  cout<<"yweight is "<<yweight<<"\n";
	  cout<<"dE is "<<dEtauidEtauf<<"\n";
	  cout<<"prob at z is "<<prob_at_z<<"\n";
	  cout<<"TauWeight_tmp is "<<TauWeight_tmp<<"\n";
	  
	  }*/
        if(zdistance<=enter_ice_mag){
          prob_at_z*=1.-exp(-1.*(r_enterice.Distance(nuexitice)/(gamma*cT)));///chance to decay in ice
        }
        else if(zdistance>enter_ice_mag){
          prob_at_z*=1.-exp(-1.*(posnunow.Distance(nuexitice)/(gamma*cT)));///chance to decay in ice if already in ice
        }
        TauWeight_tmp+=prob_at_z;///total prob for neutrino to make tau and for tau to decay in ice at this energy
	//cout<<"prob at z is "<<prob_at_z<<"\n";

        i++;
      }//etaui<pnu
    }//zdist loop
     
     
    TauWeight_tmp*=log(10)*Etau_final.in(Energy::Unit::eV);///dP/dE ->dP/d(log10(E))   
    TauWeight_tmp*=.01; ///dP/dlog10(E) * dlog10(E) 
    
    TauWeight+=TauWeight_tmp;
    // cout<<"TauWeight_tmp is "<<TauWeight_tmp<<" at Etau_final is "<<Etau_final<<"\n";
    //prob_at_z_vector.push_back(row);
    etaufarray.push_back(Etau_final);
    PDFarray.push_back(TauWeight);
    etauint++;
    
  }//Etau_final loop
  
   // }//xloop
   ///We have the weight. Now use a PDF to find the final energy of the tau.///
  double xrandom= TauWeight*gRandom->Rndm();
 
  //cout<<"TauWeight is "<<TauWeight<<" xrandom is "<<xrandom<<"\n";
  for(size_t loopthrough =0; loopthrough<=PDFarray.size();loopthrough++){
    if(xrandom >PDFarray[loopthrough] && xrandom <PDFarray[loopthrough+1]){
      ptauf = etaufarray[loopthrough];
      break;
    }//if xrandom
    
  }//loopthrough
   //cout<<"ptauf is "<<ptauf<<"\n";
  weight_tau_prob = TauWeight;
  return 1;
} //GetTauWeight


/**
Get Density Vectors sets two density vectors. One has the density at each step along the path, the other has an average density from the starting point to the current step.
 */

void icemc::TauModel::GetDensityVectors(const Antarctica *antarctica1,Interaction *interaction1, TVector3 nchord, double step, double Distance, int &totalnusteps,int &crust_entered){
   
    TVector3 nchord1;
    double avgdensity =0;//initilize average density.
    double density_total=0;//initilize running sum
    double density_now=0;//density at this step
    Geoid::Position posnunow;
    Geoid::Position postaunow;
    // double altitude_tau;
    double lat_tau;
    const Geoid::Position earth_in = interaction1->r_in;
    std::ofstream myNewFile;
    
   
  for(double taudistance=0;taudistance<=Distance;taudistance+=step){
     nchord1=taudistance*nchord;
     postaunow=earth_in+nchord1;
     lat_tau=postaunow.Latitude();
     // altitude_tau = postaunow.Mag()-Geoid::getGeoidRadius(lat_tau);
     density_now=antarctica1->Density(postaunow,crust_entered);
     mydensityvector.push_back(density_now);///filled with density at that step
     
     density_total+=density_now;
     avgdensity=density_total/(totalnusteps+1);///avgdensity
     myavgdensityvector.push_back(avgdensity);///avgdensity up to that step for neutrino
     totalnusteps++;///total number of steps neutrino will take through earth
     
  }
  
}//Get Density Vectors

/**
   Get Energy TVector3 sets the energy of tau particle at every step along the path. It starts from the final energy and works back towards the nuetrino interaction point.
 */
void icemc::TauModel::GetEnergyVector(Energy Etau_final, double step,int totalnusteps, int &totalsteps, double &totaltaudistance, Energy pnu){
  
  myenergyvector.clear();
  myenergyvector.push_back(Etau_final);
  Energy Etau_now=Etau_final;
  double density_now;
  std::ofstream myNewFile_1;
 
  // double pnu;//SET UP EVENT CLASS
  ///calculate the initial energy needed at the step so the tau will end at the correct final energy
  for(int energysteps=totalnusteps-1;energysteps>0;energysteps--){///start at nuexit, work backwards
    density_now=mydensityvector[energysteps];///get the density at this step

    double dE_eV = B0 + B1*log(Etau_now/E0)*density_now*Etau_now.in(Energy::Unit::eV)*step;///E_i=E_last +dE/dz*dz
    Etau_now = Etau_now + Energy(dE_eV, Energy::Unit::eV);
        
    if(Etau_now <=pnu){///as long as the energy is less than the neutrino energy
      totaltaudistance=(totalnusteps-energysteps)*step;///distance tau can travel
      totalsteps++;///number of steps tau can take
      myenergyvector.push_back(Etau_now);///fill energy vector
    }//Etau_now<=pnu
    else if(Etau_now >pnu){///Initial energy cannot go above pnu
      break;
    }
   
  } //energy steps


}//energy vector

/**
   Get Tau Surv Vector gets a vector that is filled with the probability the tau will survive from neutrino interaction point to the ice.
 */

void icemc::TauModel::GetTauSurvVector(double step, int totalsteps){
  myPsurvvector.clear();
  Energy Etau_now;
  Energy nothing(0, Energy::Unit::eV);
  double tau_surv;
  
  for(int k1=totalsteps;k1>=0;k1--){///tau surv vector
    tau_surv=1;
    for(int k2=k1-1;k2>=0;k2--){
      Etau_now=myenergyvector[k2];///energy vector starts at the endpoint and goes to earth_in;
      if (Etau_now > nothing)
        tau_surv=tau_surv*(1-step*mT/(cT*Etau_now));///calculate chance for tau to survive
      else
        tau_surv = 0;///is the tau runs out of energy, it wont make it to the ice to decay
    }
    myPsurvvector.push_back(tau_surv);///fill the vector
  }//tau surv

}//Tau surv vector
