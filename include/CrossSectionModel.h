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

    virtual double getSigma(double energy_eV, Neutrino::L leptonNumber, Neutrino::Current current) const = 0;

    inline static double getInteractionLength(double sigma){
      return constants::M_NUCL/sigma; // kg/m^2
    }
    
    TCanvas* plotSigma(Neutrino::L l, Neutrino::Current c, int nSamples = 2000) const;

    inline bool validEnergy(double energy_eV) const {
      return energy_eV >= fMinEnergy_eV && energy_eV <= fMaxEnergy_eV;
    }
    
  protected:
    double fMinEnergy_eV = 0;
    double fMaxEnergy_eV = 0;
  };


  class Reno : public CrossSectionModel {
  public:
    Reno(const Settings* settings) : fSettings(settings){
      fMinEnergy_eV = 1.2E15;
      fMaxEnergy_eV = 1.E21;
    }

    virtual double getSigma(double energy_eV, Neutrino::L leptonNumber, Neutrino::Current current) const override;
  private:
    const Settings* fSettings;
  };
  
}

#endif //CROSS_SECTION_MODEL_H
