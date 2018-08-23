// -*- C++ -*-.
/**************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Give this class a TTree plus a branch name to serve as an 'x-axis', normally this is a time.
             FancyTTreeInterpolator generates TGraphs of any TTree branch variable as a function of your x-axis.
             That allows us to do two slightly clever things.
             The first slightly clever thing is that it sorts the data so it is x-axis (time) ordered.
             The second slightly clever thing is that it can then interpolate the data to values between entries.
             This interpolation is equivalent to TGraph::Eval.
**************************************************************************************************************/

#ifndef FANCYTTREEINTERPOLATOR_H
#define FANCYTTREEINTERPOLATOR_H

#include "TObject.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TMath.h"
#include "TAxis.h"
#include "TCanvas.h"

#include <iostream>
#include <exception>
#include <stdexcept>

namespace icemc
{

  /**
   * @class FancyTTreeInterpolator
   * @brief A class to interpolate sparse, but continuous data in a TTree. 
   * 
   * This was developed with adu5Pat variables like heading in mind.
   */
  class FancyTTreeInterpolator{

  public:
    FancyTTreeInterpolator(TTree* t, TString xAxisText);
    ~FancyTTreeInterpolator();

    void add(TString yAxisText);
    void add(TString yAxisText, TString cut);
    void add(TString yAxisText, Double_t wrapValue);
    void add(TString yAxisText, TString cut, Double_t wrapValue);
    Double_t interp(TString yAxisText, Double_t xAxisValue);
    TGraph* get(TString yAxisText);
    TGraph* makeSortedTGraph(TString yAxisText);
    TGraph* makeSortedTGraph(TString yAxisText, TString cutString);
    TGraph* makeSortedTGraph(TString yAxisText, Double_t wrapValue);
    TGraph* makeSortedTGraph(TString yAxisText, TString cutString, Double_t wrapValue);

    TTree* fTree; //!< TTree with which the intepolater was initialized.
    TString fXAxisText; //!< Branch name with which the intepolater was initialized.
    std::map<TString,TGraph*> fStringToGraph; //!< Internally stored TGraphs, accessed by TTree branch name.
    std::map<TString, Double_t> fStringToWrapValue; //!< Internally stored wrapValues, accessed by TTree branch name.
    Double_t fXmin; //!< Stored x-axis lower limit.
    Double_t fXmax; //!< Stored x-axis lower limit.
  
  };
}

#endif //FANCYTTREESORTERANDINTERPOLATOR_H
