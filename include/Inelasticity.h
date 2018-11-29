#ifndef ICEMC_INELASTICITY_H
#define ICEMC_INELASTICITY_H

#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TF3.h"
#include "TH2D.h"

#include "Geoid.h"
#include "TVector3.h"
#include "Neutrino.h"
#include "RNG.h"
#include "Interaction.h"

namespace icemc {

  // class Settings;

  /**
   * @class YGenerator
   * @brief Inelasticity distributions: stores parameterizations and picks inelasticities
   */
  class YGenerator {
  public:
    ///@todo weight?
    virtual double pickY(Energy energy, Neutrino::L leptonNumber, Interaction::Current current) = 0;
  };

  class ConstantY : public YGenerator {
    double fY;
  public:
    ConstantY(double y) : fY(y){
      if(fY <= 0 || fY > 1){
	fY = 0.2;
      }
    }
    double pickY(Energy energy, Neutrino::L leptonNumber, Interaction::Current current) override {
      return fY;
    }
  };
  

  /**
   * @namespace GhandiEtAl
   * @brief hep-ph/9512364 (the curves are not in their later article)
   * 
   * There is also a slow energy dependence.
   */

  namespace GhandiEtAl {
    const double R1 = 0.36787944;			///< 1/e, used for Gandhi et al.
    const double R2 = 0.63212056;			///< 1-R1, used for Gandhi et al.
    
    class YGenerator : public icemc::YGenerator, public RNG {
    public:
      virtual ~YGenerator(){};
      
      /** 
       * @brief THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM Ghandhi,Reno,Quigg,Sarcevic hep-ph/9512364
       * 
       * The curves are not in their later article.
       * There is also a slow energy dependence.
       * 
       * @return a parameterized y value
       */
      virtual double pickY(Energy energy, Neutrino::L leptonNumber, Interaction::Current current) override;            
    };    
  };


}




#endif
