#ifndef ICEMC_RNG_H
#define ICEMC_RNG_H

#include <memory>
#include <map>
#include <string>

class TRandom3;

namespace icemc {
   
  /**
   * @class RNG
   * @brief A random number generator with a fairly deterministic seed setting function.
   * 
   * This class wraps a TRandom3 object.
   * An interface is provided to reset all seeds with a single function call.
   * @see icemc::RNG::newSeeds(int run, UInt_t eventNumber);
   * The class is designed to be inhierited
   * The seed depends on:
   * - the name of the class, provided by typeid::name
   * - the name ID, #fID, which counts how many class objects with the same name there are.
   *   - This fID is assigned the first time the pick* functions are first called.
   * - the run and eventNumber given to RNG::newSeeds(int run, UInt_t eventNumber)
   */
  
  class RNG {
  public:

    RNG();
    virtual ~RNG();
    
    /** 
     * Call this function to update the seeds for =every= RNG.
     * 
     * @param run the run
     * @param eventNumber the eventNumber
     */
    static void newSeeds(int run, unsigned int eventNumber);
    
    double pickUniform(double x1=1);
    double pickUniform(double x1, double x2);    
    double pickGaus(double mean=0, double sigma=1);
    int pickPoisson(double mean);

  private:

    void updateSeed();
    std::unique_ptr<TRandom3> fRandom;
    int fRun = 0;
    int fEventNumber = 0;
    int fID = -1;

    static int gRun; // global run
    static unsigned int gEventNumber; // global eventNumber
    static std::map<std::string, unsigned int> fNameCount; // for counting instances with the same class name;

  };

}



#endif // ICEMC_INTEGRATOR_H
