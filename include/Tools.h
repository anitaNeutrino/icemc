////////////////////////////////////////////////////////////////////////////////////////////////
//namespace Tools:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TOOLS_H_
#define TOOLS_H_

// c++ libraries thingies
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>

// ROOT
#include "TObject.h"

class TH1;
class TGraph;
class TH2F;
class TRandom3;
class TTree;
class TStyle;

namespace icemc{

  //! Functions to make life easier.  Many of these probably exist other places.
  namespace Tools {
    // double dMax(const double*,int);
    int Getifreq(double freq,double freq_low,double freq_high,int n);
    void NormalTimeOrdering(const int n,double *volts);
    void reverseTimeOrdering(const int n,double *bitsin,double *bitsout);
    void reverseTimeOrdering(const int n,int *bitsin,int *bitsout);
    void ShiftLeft(double *x,const int n,int ishift);  
    void ShiftRight(double *x,const int n,int ishift);  
    double GetFWHM(TH1 *h1);
    int findIndex(const double *freqlab,double freq,int npoints,double min,double max); // this is the same thing but more general
    TGraph *getInterpolatedGraph(TGraph *grIn, Double_t deltaT);
    
    template<class V> double max(const V& v){
      return v.size() > 0 ? *std::max_element(v.begin(), v.end()) : 0;
    }

    template<class V> double min(const V& v){
      return v.size() > 0 ? *std::min_element(v.begin(), v.end()) : 0;
    }

    template<class V> double RMS(const V& v){
      double mean = std::accumulate(v.begin(), v.end(), 0.0)/v.size();

      // First param in the lambda is the init value for the 0th elem, and the accumulated value for the nth element.
      // https://stackoverflow.com/questions/29685003/how-can-i-use-stdaccumulate-and-a-lambda-to-calculate-a-mean      
      auto accumulate_square = [](double x, double y){return x + y*y;};
      double mean_of_square = std::accumulate(v.begin(), v.end(), 0.0, accumulate_square)/v.size();
      double rms = mean_of_square - mean*mean;
      return rms;
    }

    template<class V1, class V2> double calculateSNR(const V1& v1, const V2& v2){
      double p2p = max(v1) - min(v1);
      double rms = RMS(v2);
      return p2p/(2*rms);
    }
    
    template <class T, class U> void vector_element_convert(const std::vector<T>& input, std::vector<U>& output){
      output.clear();
      output.reserve(input.size());
      for (unsigned int index = 0; index < input.size(); index++){
	output.push_back(U (input[index]));
      }
    }

    template <class T, class U> void nested_vector_element_convert(const std::vector< std::vector<T> >& input, std::vector< std::vector<U> >& output){
      output.clear();
      output.reserve(input.size());
      for (unsigned int index = 0; index < input.size(); index++){
	std::vector<U> temp;
	vector_element_convert<T,U>(input[index], temp);
	output.push_back(temp);
      }
    }
  };
}
#endif
