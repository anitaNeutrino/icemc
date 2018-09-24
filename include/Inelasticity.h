
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

namespace icemc {

  class Interaction;
  class Primaries;
  class IceModel;
  class Counting;
  class Settings;


  /**
   * @class Y
   * @brief Inelasticity distributions: stores parametrizations and picks inelasticities
   */
  class Y :  public RNG {  

  private:

    std::unique_ptr<TF1> ffrac;          ///< This is the fraction of the distribution in the low y region given by Equation 18. 
    // std::unique_ptr<TF1> fC1_high[2][2]; ///< parameterization of parameter C1 in the high y region according to Equation 16
    std::map<std::pair<Neutrino::L, Neutrino::Current>, TF1> fC1_high; ///< parameterization of parameter C1 in the high y region according to Equation 16    
    std::unique_ptr<TF1> fC1_low;        ///< parameterization of parameter C1 in the low y region according to Equation 16.
    std::unique_ptr<TF1> fC2;            ///< parameterization of parameter C2 in the low y region according to Equation 17.
    std::unique_ptr<TF3> fy0_low;        ///< For picking inelasticity in low y region according to Equation 14.
    std::unique_ptr<TF2> fy0_high;       ///< For picking inelasticity in high y region according to Equation 15.

    double pickYConnollyetal2011(Neutrino::L leptonNumber,Neutrino::Current CURRENT,double e);	///< pick an inelasticity using recipe in Connolly et al. (2011)
    // L =0: nubar, NU=1: nu
    // CURRENT=0: CC, CURRENT-1: NC

    static constexpr double ymin_low  = 0.00002;		///< Minimum y in low-y region, Connolly et al.
    static constexpr double ymax_low  = 0.001;			///< Maximum y in low-y region, Connolly et al.
    static constexpr double ymin_high = 0.001;			///< Minimum y in high-y region, Connolly et al.
    static constexpr double ymax_high = 1.;			///< Maximum y in high-y region, Connolly et al.
    static constexpr double dy_low    = 0.00002;		///< y step in low region, Connolly et al.
    static constexpr double dy_high   = 0.001;			///< y step in high region, Connolly et al.

    /** 
     * @brief THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM Ghandhi,Reno,Quigg,Sarcevic hep-ph/9512364
     * 
     * The curves are not in their later article.
     * There is also a slow energy dependence.
     * 
     * @return a parameterized y value
     */
    double pickYGandhietal();

    static constexpr double R1 = 0.36787944;			///< 1/e, used for Gandhi et al.
    static constexpr double R2 = 0.63212056;			///< 1-R1, used for Gandhi et al.

    
  public:
    Y();

    /** 
     * Pick inelasticity y according to chosen model
     * 
     * @param settings1 are the icemc settings 
     * @param pnu 
     * @param nu_nubar 
     * @param currentint 
     * 
     * @return 
     */
    double pickY(const Settings *settings1, double pnu, Neutrino::L leptonNumber, Neutrino::Current currentint);

    /** 
     * @brief If you want to choose y from a flat distribution this is the weight it should have according to Connolly et al. (2011)
     * 
     * @param pnu 
     * @param y 
     * @param nu_nubar 
     * @param currentint 
     * 
     * @return 
     */
    // double Getyweight(double pnu, double y, int nu_nubar, Neutrino::Current currentint);
    double Getyweight(double pnu, double y, Neutrino::L leptonNumber, Neutrino::Current currentint);    
    
  };//Y



}




#endif
