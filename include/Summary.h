#ifndef ICEMC_SUMMARY
#define ICEMC_SUMMARY

#include "TMath.h"
#include "Event.h"

namespace icemc{

  /**
   * @class FlavorSummary
   * @brief For given flavorSummarize results from run, including weighted events passing and effective volume
   */  
  class FlavorSummary{
    /**
     * @class Flavor
     * @brief enum for neutrino flavour, compares to Neutrino::Flavor other than 'all'
     */

    
  public:
    enum class Flavor : int {
      all = 0,
      e = 1,
      mu = 2,
      tau = 3
    };    
    FlavorSummary(Flavor f){ flavor = f; };

    /**
     * If an event has passed the trigger, add it to the summary 
     *
     * @param event -- simulated event that passed trigger
     */    
    void addEvent(const EventSummary& event, bool passed=false);

    /**
     * After the simulation is done, do the final calculations 
     *
     * @param n -- number of neutrinos that were simulated of this flavor
     * @param iceVolume -- volume of ice in Antarctica, in km^3
     */    
    void summarize(double iceVolume);
    void reportSummary();
    
    
    //private:

    Flavor flavor;
    int nTotal = 0; // number of simulated neutrinos
    int nPass = 0; // number of passing neutrinos
    double nWeighted = 0; // weighted number of passing neutrinos
    double length = 0; // interaction length for neutrinos, averaged
    double effectiveVolume = -1; // km^3 str
    double effectiveArea = -1; // km^2 str

    /**
     * Calculates the effective volume of the 'detector' at the given energy, stores it in effectiveVolume
     *
     * @param iceVolume -- volume of ice in Antarctica, in km^3
     */    
    void calculateEffectiveVolume(double iceVolume);

    /**
     * Calculates the effective area of the 'detector' at the given energy, stores it in effectiveArea
     */    
    void calculateEffectiveArea();

  };

  /**
   * @class Summary
   * @brief Summarize results from run, including weighted events passing and effective volume
   */
    class Summary : public FlavorSummary {

    public:
      Summary(double e, double vol) : FlavorSummary(FlavorSummary::Flavor::all) {
	exponent = e;
	iceVolume = vol;
      }

      void addEvent(const EventSummary& event, bool passed=false);
      void summarize(); // how many of each we ended up simulating

      double exponent; // simulated energy exponent, 0 if not monoenergetic
      double iceVolume; // m^3 in all of antarctica    

      //private:
    
      FlavorSummary eSummary = FlavorSummary(FlavorSummary::Flavor::e); // electon neutrinos summary
      FlavorSummary muSummary = FlavorSummary(FlavorSummary::Flavor::mu); // mu neutrinos summary
      FlavorSummary tauSummary = FlavorSummary(FlavorSummary::Flavor::tau); // tau neutrinos summary
    
      ClassDef(Summary, 0);
    };

}

#endif //ICEMC_SUMMARY
