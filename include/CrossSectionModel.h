#ifndef CROSS_SECTION_MODEL_H
#define CROSS_SECTION_MODEL_H

#include "Neutrino.h"
#include "Interaction.h"
#include "Constants.h"

class TCanvas;

namespace icemc {

  class Settings;

  /**
   * @class CrossSectionModel
   * @brief Pure virtual class to model neutrino-nucleon cross section
   */

  class CrossSectionModel {
  public:
    CrossSectionModel(const Settings* settings) : fSettings(settings) {;}
    virtual ~CrossSectionModel(){;}
    virtual double getSigma(Energy energy, Neutrino::L leptonNumber, Interaction::Current current) const = 0;

    inline static double getInteractionLength(double sigma){
      return constants::M_NUCL/sigma; // kg/m^2
    }
    
    TCanvas* plotSigma(Neutrino::L l, Interaction::Current c, int nSamples = 2000) const;

    inline bool validEnergy(const Energy& energy) const {
      return energy >= fMinEnergy && energy <= fMaxEnergy;
    }
    
  protected:
    const Settings* fSettings = nullptr;
    Energy fMinEnergy;
    Energy fMaxEnergy;
  };





  
  namespace MHReno {

    class CrossSectionModel : public icemc::CrossSectionModel {
    public:
      CrossSectionModel(const Settings* settings) : icemc::CrossSectionModel(settings)
      {
	fMinEnergy = Energy(1.2e15, Energy::Unit::eV);
	fMaxEnergy = Energy(1.0e21, Energy::Unit::eV);
      }
      virtual ~CrossSectionModel(){;}

      virtual double getSigma(Energy energy, Neutrino::L leptonNumber, Interaction::Current current) const override;
    };    
  }


  
}

#endif //CROSS_SECTION_MODEL_H
