////////////////////////////////////////////////////////////////////////////////////////////////
//namespace Tools:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TOOLS_H_
#define TOOLS_H_

#include "TRandom3.h" 
#include <vector>

// c++ libraries thingies
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>

// ROOT
#include "TGraph.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

//using std::string;
//using std::vector;

using namespace std;

class TSpline5;
class TH1;
class TGraph;
class TH2F;
class TRandom3;
class TTree;

//! Functions to make life easier.  Many of these probably exist other places.
namespace Tools {
    double dMax(double,double);
    double dMax(const double*,int);
    double dvMax(const vector<double>);
    int findLonLatPair(vector<int>,vector<int>,int,int);
    double dsMax(TSpline5 *sp);
    double dMin(const double*,int);
    double dMinNotZero(const double*,int);
    double dMin(double,double);
    double getMaxMagnitude(vector<double> v);
    int Getifreq(double freq,double freq_low,double freq_high,int n);
    void InterpolateReal(double* array, const unsigned n);
    void InterpolateComplex(double *array, const unsigned n);

    void four1(double *data, const int isign,int nsize);
    void realft(double *data, const int isign, int nsize);

    void SWAP(double &a, double &b);// swaps two numbers
    void NormalTimeOrdering(const int n,double *volts);
    void reverseTimeOrdering(const int n,double *bitsin,double *bitsout);
    void reverseTimeOrdering(const int n,int *bitsin,int *bitsout);

    void ShiftLeft(double *x,const int n,int ishift);  
    void ShiftRight(double *x,const int n,int ishift);  
    double GetFWHM(TH1 *h1);
    void  MakeGraph(int index, int n,double *time,double *volts,TGraph *&mygraph,TH2F *&h2, double scalex,double scaley,string xaxistitle,string yaxistitle);

    double dDot(double*,double*,int);
    void dCross(double*,double*,double*);
    double dSquare(double*);
    double Step(double x);
    double dGetTheta(double*);
    double dGetPhi(double*);
    int WhichIsMax(double *x,int n);
    int WhichIsMin(double *x,int n);

    double dSum(double*,int);
    int iSum(int*,int);
    void Print(double*,int);
    void Print(int*,int);
    void Zero(double *anarray,int n);
    void Zero(int *anarray,int n);
    int NonZero(double *anarray,int n); // how many bins are nonzero
    void GetNumbersAsStringArray(ifstream& fin, ofstream& fout,vector<string>& vnumbers, int nelements);
    void GetNext2NumbersAsString(ifstream& fin,ofstream& fout,string& number1,string& number2, string& stherest);
    void GetNextNumberAsString(ifstream& fin,ofstream& fout,string& number);

    int findIndex(double *freqlab,double freq,int npoints,double min,double max); // this is the same thing but more general
    void get_random_rician(double signal, double signal_phase, double sigma, double& amplitude, double &phase);
    void get_circular_bivariate_normal_random_variable(double& rand_gauss_a, double& rand_gauss_b);
    int round(double number);

    double AbbyPhiCalc(double x_abby, double y_abby);

    TGraph *getInterpolatedGraph(TGraph *grIn, Double_t deltaT);

    double calculateSNR(double justSig[512], double justNoise[512]);
    
    template <class T, class U> void vector_element_convert(const vector<T>& input, vector<U>& output){
    output.clear();
        for (unsigned int index = 0; index < input.size(); index++){
            output.push_back(U (input[index]));
        }
    }

    template <class T, class U> void nested_vector_element_convert(const vector< vector<T> >& input, vector< vector<U> >& output){
    output.clear();
        for (unsigned int index = 0; index < input.size(); index++){
            vector<U> temp;
            vector_element_convert<T,U>(input[index], temp);
            output.push_back(temp);
        }
    }
};
#endif
