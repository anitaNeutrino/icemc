#ifndef CONNOLLY_ET_AL_2011_H
#define CONNOLLY_ET_AL_2011_H

#include <iostream>
#include "CrossSectionModel.h"
#include "RNG.h"

#include "TF1.h"
#include "TF2.h"
#include "TF3.h"

#include "Geoid.h"
#include "TVector3.h"
#include "Neutrino.h"
#include "Inelasticity.h"

// std::vector;
namespace icemc {  
  
  class Interaction;
  class Settings;

  /**
   * @class ConnollyEtAl2011
   * @brief To compute Cross sections given in https://arxiv.org/pdf/1102.0691.pdf
   */
  class ConnollyEtAl2011 : public CrossSectionModel, public RNG {

    // static constexpr int NSIGMAS=2;///< number of possible cross section models
    ///< 0=Gandhi et al.
    ///< 1=Connolly et al. 2011
    // double mine[NSIGMAS]; ///< minimum energy for cross section parametrizations, in eV
    // double maxe[NSIGMAS]; ///< maximum energy for cross section parametrizations, in eV
    
  public:
    double pickY(double pnu, Neutrino::L leptonNumber, Neutrino::Current currentint);///<pick inelasticity y according to chosen model    
    double Getyweight(double pnu,double y,Neutrino::L leptonNumber,Neutrino::Current currentint);///< in case you choose y from a flat distribution, this is the weight you should give it according to Connolly et al. (2011)

    ConnollyEtAl2011(const Settings* settings); ///< Constructor 
    
    /// Neutrino-nucleon cross-sections using model chosen
    virtual double getSigma(double pnu, Neutrino::L leptonNumber, Neutrino::Current currentint) const override;

  protected:
    const Settings* fSettings;
    
  private:
    
    Y m_myY; ///< defined in Inelasticity.h
    
    typedef std::map<std::pair<Neutrino::L, Neutrino::Current>, double> DoubleLCC; // represent a double for nu/nubar, neutral/charged currents
    std::array<double, 4> A_low; ///< Table V of Connolly et al. for use in Eq. 16.
    DoubleLCC A0_high; ///< Table V of Connolly et al. for use in Eq. 16.
    DoubleLCC A1_high; ///< Table V of Connolly et al. for use in Eq. 16.
    DoubleLCC A2_high; ///< Table V of Connolly et al. for use in Eq. 16.
    DoubleLCC A3_high; ///< Table V of Connolly et al. for use in Eq. 16.
    double b0;         ///<  Eq. 17 of Connolly et al.
    double b1;         ///<  Eq. 17 of Connolly et al.
    
    // TF1* m_fy[2][2];
    

    std::map<std::pair<Neutrino::L, Neutrino::Current>, TF1> fSigma;
    
    DoubleLCC c0;      ///< Table V of Connolly et al. for Eq. 7
    DoubleLCC c1;      ///< Table V of Connolly et al. for Eq. 7
    DoubleLCC c2;      ///< Table V of Connolly et al. for Eq. 7
    DoubleLCC c3;      ///< Table V of Connolly et al. for Eq. 7
    DoubleLCC c4;      ///< Table V of Connolly et al. for Eq. 7
    
  };///<ConnollyEtAl2011
  
}



#endif
