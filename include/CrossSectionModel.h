#ifndef CROSS_SECTION_MODEL_H
#define CROSS_SECTION_MODEL_H

#include "Neutrino.h"
#include "Constants.h"

class TCanvas;

namespace icemc {

  class Settings;

  /**
   * @class CrossSectionModel
   * @brief Pure virtual class to model neutrino cross section
   */

  class CrossSectionModel {
  public:
    virtual ~CrossSectionModel(){;}
    // virtual double getSigma(double energy_eV, Neutrino::L leptonNumber, Neutrino::Interaction::Current current) const = 0;
    virtual double getSigma(Energy energy, Neutrino::L leptonNumber, Neutrino::Interaction::Current current) const = 0;

    inline static double getInteractionLength(double sigma){
      return constants::M_NUCL/sigma; // kg/m^2
    }
    
    TCanvas* plotSigma(Neutrino::L l, Neutrino::Interaction::Current c, int nSamples = 2000) const;

    inline bool validEnergy(const Energy& energy) const {
      return energy >= fMinEnergy && energy <= fMaxEnergy;
    }
    
  protected:
    Energy fMinEnergy;
    Energy fMaxEnergy;
  };




  

  class MHReno : public CrossSectionModel {
  public:
    MHReno(const Settings* settings)
      : fSettings(settings)
    {
      fMinEnergy = Energy(1.2E15, Energy::Unit::eV);
      fMaxEnergy = Energy(1.E21, Energy::Unit::eV);
    }
    virtual ~MHReno(){;}

    virtual double getSigma(Energy energy, Neutrino::L leptonNumber, Neutrino::Interaction::Current current) const override;
  private:
    const Settings* fSettings;
  };
  
}

#endif //CROSS_SECTION_MODEL_H
