#include "IcemcStatistics.h"
#include <iostream>

void icemc::statistics::getPoissonError(int n, double& poissonErrorPlus, double& poissonErrorMinus){

  const int np = 21;
  static const double poissonerror_minus[np] = {0.-0.00, 1.-0.37, 2.-0.74, 3.-1.10, 4.-2.34,
						5.-2.75, 6.-3.82, 7.-4.25, 8.-5.30, 9.-6.33,
						10.-6.78, 11.-7.81, 12.-8.83, 13.-9.28,
						14.-10.30, 15.-11.32, 16.-12.33, 17.-12.79,
						18.-13.81, 19.-14.82, 20.-15.83};
  static const double poissonerror_plus[np] = {1.29-0., 2.75-1., 4.25-2., 5.30-3., 6.78-4.,
					       7.81-5., 9.28-6., 10.30-7., 11.32-8., 12.79-9.,
					       13.81-10., 14.82-11., 16.29-12., 17.30-13.,
					       18.32-14., 19.32-15., 20.80-16., 21.81-17.,
					       22.82-18., 23.82-19., 25.30-20.};  
  if(n > -1 && n < np){
    poissonErrorPlus = poissonerror_plus[n];
    poissonErrorMinus = poissonerror_minus[n];
  }
  else{
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", don't have poisson errors for " << n << " events." << std::endl;
    poissonErrorPlus = -1;
    poissonErrorMinus = -1;
  }
}
