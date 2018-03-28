////////////////////////////////////////////////////////////////////////////////////////////////
//class Taumodel:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef TAU_MODEL_H
#define TAU_MODEL_H

#include <algorithm>
#include <numeric>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#ifndef __CINT__
#include "Primaries.h"
#include "Settings.h"
#endif

#include "vector.hh"

class TH1F;

namespace icemc{
  class Primaries;
  class Interaction;
  class IceModel;

  class Taumodel: public TObject
  {

  private: //stuff other programs arent allowed to touch
 	
    /**For Tau Weight & Survival probability equations.
       n.b. not in SI units.
       from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
       Equation 16  &  used in Equation 30. 
    */
    double B0,B1,E0;/*!<parameterization using a logarithmic dependence on energy>*/ 
    ///for B, the tau elecromagnetic energy loss parameter. 
    double mT;/*!<mass of the Tau in Gev> */
    double cT;/*!<Tau Decay length in cm> */
    // double p;/*!<Density of Standard Rock. g/cm^3> */
  
    ///Used in Connolly Calc 2011.(d_dzPsurvNu())
    //p, the Density of Standard Rock. g/cm^3
    double A; /*!<constant that sets the total probability to unity> */
    double Mn;/*! < nucleon/ proton mass in grams,also equal to 0.938 GeV.> */
  
    std::vector<double> mydensityvector;/*!<Filled with density at each step> */
  
    std::vector<double> myavgdensityvector;/*!<Filled with avg density from earth-in to interaction point.> */

    std::vector<double> myenergyvector;/*!<vector to hold initial energies.> */
    std::vector<double> myPsurvvector;/*!<vector to hold the chance tau would survive from that step.> */
    std::vector<double> etaufarray;/*! <filled with all the possible tau final energies.> */
    std::vector<double> PDFarray;/*! <Filled with the total weight(a running sum) up to the current final energy.> */
  
    /** \brief  Get Density Vectors sets two density vectors. One has the density at each step along the path, the other has an average density from the starting point to the current step.
     */
    void GetDensityVectors(IceModel *antarctica1, Interaction *interaction1, Vector nchord, double step, double Distance, int &totalnusteps, int &crust_entered);
    /** \brief   Get Energy Vector sets the energy of tau particle at every step along the path. It starts from the final energy and works back towards the nuetrino interaction point.
     */
    void GetEnergyVector(double Etau_final, double step,int totalnusteps, int &totalsteps, double &totaltaudistance, double pnu);
  
    /** \brief  Get Tau Surv Vector gets a vector that is filled with the probability the tau will survive from neutrino interaction point to the ice.
     */
    void GetTauSurvVector(double step, int totalsteps);
  public:
    Taumodel();
    double ptauf; ///< Final energy of the tau.
    double weight_tau_prob; ///< Weight for tau neutrino to interact, create a tau, tau survives and decays in the ice.
    int inu; 
    double weight_nu_prob;
 
    /** \brief GetTauWeight is the function that will calculate the probability that a tau neutrino will interact along its path through the earth,and the tau will survive the rest of the journey and decay in the ice. This probability is calculated for final energies from 10^15.5 to the energy of the neutrino.
     */

    double GetTauWeight(Primaries *primary1, const Settings *settings1, IceModel*antarctica1, Interaction *interaction1, double pnu, int nu_nubar, double& ptauf, int& crust_entered); // 1 or 0
    // int& mantle_entered, // 1 or 0
    // int& core_entered);//include secondaries?
 
    ClassDef(Taumodel,1);
  }; //class Taumodel

}
#endif //TAU_MODEL_H
