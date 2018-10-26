#include "Constants.h"
#include "TVector3.h"
#include "Geoid.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include "Settings.h"
#include "Crust2.h"
#include "Antarctica.h"
#include "ConnollyEtAl2011.h"
#include "Report.h"
#include "RayTracer.h"

#include <cmath>

#include "TH2D.h"
#include "TCanvas.h"

#include "Inelasticity.h"



//! THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM Ghandhi,Reno,Quigg,Sarcevic  hep-ph/9512364 (the curves are not in their later article).  There is also a slow energy dependence.
double icemc::GhandiEtAl::YGenerator::pickY(Energy energy, Neutrino::L leptonNumber, Neutrino::Interaction::Current current) {
  (void) energy;
  (void) leptonNumber;
  (void) current;
    
  // generate according to Ghandi fig. 6 
  // adjust exponent until looks like the curve
  //  and has right mean.
  //  (Note this is not the fcn, but the inverse of the integral...)
  double rnd = pickUniform(); //gRandom->Rndm(1); // (0,1)
  //  cout << "R1, R2, rnd are " << R1 << " " << R2 << " " << rnd << "\n";
  double x = pow(-log(R1+rnd*R2),2.5);
  return x;
}
