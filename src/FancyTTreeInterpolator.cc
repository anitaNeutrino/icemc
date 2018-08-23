#include "FancyTTreeInterpolator.h"





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * @param t is a TTree
 * @param xAxisText is the name of a time-like branch in the TTree. realTime would be a good choice.
 *
 * Does not presume to own t.
 * If t s deleted out from under it, this class will fall over.
 */
icemc::FancyTTreeInterpolator::FancyTTreeInterpolator(TTree* t, TString xAxisText){
  /* Do use this constructor */
  fXAxisText = xAxisText;
  fTree = t;
  TGraph* gr = makeSortedTGraph(xAxisText + ":" + xAxisText);
  //  TGraph* gr = makeSortedTGraph(xAxisText + ":" + xAxisText);
  fXmin = gr->GetX()[0];
  fXmax = gr->GetX()[gr->GetN()-1];

  // TFile* fout = new TFile("ftti.root", "recreate");
  // gr->Write();
  // fout->Close();

  delete gr;

  /* For error reporting, want precision of unsigned int ~10 digits*/
  std::cerr.precision(10);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Desonstructor
 *
 * Deletes interally stored TGraphs.
 */
icemc::FancyTTreeInterpolator::~FancyTTreeInterpolator(){
  std::map<TString, TGraph*> ::iterator i;
  for(i=fStringToGraph.begin(); i!=fStringToGraph.end(); ++i){
    delete (*i).second;
    (*i).second = NULL;
  }  
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Makes a sorted TGraph from fTree->Draw() with no cuts.
 *
 * @param drawText is the text to pass to fTree::Draw();
 * @returns a pointer to the created, sorted TGraph.
 */
TGraph* icemc::FancyTTreeInterpolator::makeSortedTGraph(TString drawText){
  return makeSortedTGraph(drawText, "", 0);
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Makes a sorted TGraph from fTree->Draw() with cuts.
 *
 * @param drawText is the text to pass to fTree->Draw();
 * @param cutString defines the cuts to padd to fTree->Draw();
 * @returns a pointer to the created, sorted TGraph.
 */
TGraph* icemc::FancyTTreeInterpolator::makeSortedTGraph(TString drawText, TString cutString){
  return makeSortedTGraph(drawText, cutString, 0);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Makes a sorted TGraph from fTree->Draw() with no cuts.
 *
 * @param drawText is the text to pass to fTree->Draw();
 * @param wrapValue value to unwrap around. e.g. 360 for rotations in degrees.
 * @returns a pointer to the created, sorted TGraph.
 * 
 * For example if you want to interpolate between 359 and 1 degrees, these will be unwrapped as  359, 361.
 * These values can then be interpolated between correctly.
 */
TGraph* icemc::FancyTTreeInterpolator::makeSortedTGraph(TString drawText, Double_t wrapHigh, Double_t wrapLow){
  return makeSortedTGraph(drawText, "", wrapHigh, wrapLow);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Makes a sorted TGraph from fTree->Draw() with cuts.
 *
 * @param drawText is the text to pass to fTree->Draw();
 * @param cutString defines the cuts to padd to fTree->Draw();
 * @param wrapValue value to unwrap around. e.g. 360 for rotations in degrees.
 * @returns a pointer to the created, sorted TGraph.
 * 
 * For example if you want to interpolate between 359 and 1 degrees, these will be unwrapped as  259, 361.
 * These values can then be interpolated between correctly.
 */
TGraph* icemc::FancyTTreeInterpolator::makeSortedTGraph(TString drawText, TString cutString, Double_t wrapHigh, Double_t wrapLow){
  /*
    This function:
       Uses TTree::Draw to select the data.
       Sorts it using TMath::Sort.
       Returns the graph.
       Note: This does not add it to the std::map of graphs that this class uses for storage.
             To create a graph and add it to the std::map use icemc::FancyTTreeInterpolator::add.
  */
  
  // Draw
  const Int_t nEntries = fTree->Draw(drawText, cutString, "goff"); // "goff" means graphics off

  // Sort
  std::vector<Int_t> sortedIndices(nEntries);
  TMath::Sort(nEntries, fTree->GetV2(), &sortedIndices.front(), kFALSE);
  std::vector<Double_t> newX(nEntries);
  std::vector<Double_t> newY(nEntries);
  
  for(int i=0; i<nEntries; i++){
    newX.at(i) = fTree->GetV2()[sortedIndices.at(i)];
    newY.at(i) = fTree->GetV1()[sortedIndices.at(i)];
  }

  // Unwrap data here so interpolation works smoothly
  // will unwrap when getting entries from graph
  if(wrapHigh != 0){
    for(int i=1; i<nEntries; i++){
      // If y[i] >> y[i-1] => then we went below zero to the wrap value
      // can only modify y[i], so we should subtract wrapValue enough times 
      // that the graph becomes smooth
      while (newY.at(i) - newY.at(i-1) > (wrapHigh - wrapLow)/2){
	newY.at(i) -= (wrapHigh - wrapLow);
      }

      // If y[i] << y[i-1] => then we went add zero to the wrap value
      // can only modify y[i], so we should add wrapValue enough times 
      // that the graph becomes smooth
      while (newY.at(i) - newY.at(i-1) < -(wrapHigh - wrapLow)/2){
	newY.at(i) += (wrapHigh - wrapLow);
      }    

    }   
  }

  // sorted TGraph
  TGraph* gr(new TGraph(nEntries,&newX.front(), &newY.front()));
  gr->SetTitle(drawText + ", " + cutString);
  gr->GetXaxis()->SetTitle(fXAxisText);
  gr->GetYaxis()->SetTitle(drawText);

  return gr;

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Adds a TGraph to the interally stored TGraphs.
 *
 * @param yAxisText is defines the y-axis of the TGraph, drawn against fXaxisText
 */
void icemc::FancyTTreeInterpolator::add(TString yAxisText){
  add(yAxisText, "", 0);
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Adds a TGraph to the interally stored TGraphs.
 *
 * @param yAxisText is defines the y-axis of the TGraph, drawn against fXaxisText
 * @param cutString defines the cuts to padd to fTree->Draw();
 */
void icemc::FancyTTreeInterpolator::add(TString yAxisText, TString cutString){
  add(yAxisText, cutString, 0);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Adds a TGraph to the interally stored TGraphs.
 *
 * @param yAxisText is defines the y-axis of the TGraph, drawn against fXaxisText
 * @param wrapValue value to unwrap around. e.g. 360 for rotations in degrees.
 */
void icemc::FancyTTreeInterpolator::add(TString yAxisText, Double_t wrapHigh, Double_t wrapLow){
  add(yAxisText, "", wrapHigh, wrapLow);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Adds a TGraph to the interally stored TGraphs.
 *
 * @param yAxisText is defines the y-axis of the TGraph, drawn against fXaxisText
 * @param cutString defines the cuts to padd to fTree->Draw();
 * @param wrapValue value to unwrap around. e.g. 360 for rotations in degrees.
 */
void icemc::FancyTTreeInterpolator::add(TString yAxisText, TString cutString, Double_t wrapHigh, Double_t wrapLow){
  TString drawText = yAxisText + ":" + fXAxisText;
  TGraph* gr = makeSortedTGraph(drawText, cutString, wrapHigh, wrapLow);
  fStringToGraph[yAxisText] = gr;
  fStringToWrapValues[yAxisText] = std::make_pair(wrapHigh, wrapLow);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Access the interally stored TGraph via the yAxisText
 *
 * @param yAxisText is defines the y-axis of the TGraph, drawn against fXaxisText
 * @return A pointer to accessed TGraph, returns NULL if no match is found.
 */
TGraph* icemc::FancyTTreeInterpolator::get(TString yAxisText){
  /* Use this to access a graph you've made with the add function */

  if(fStringToGraph.count(yAxisText)==0){
    std::cerr << "Can't find TGraph in FancyTTreeInterpolator created with text " + yAxisText << std::endl;

    // Hack for c++98... in 2014
    throw std::invalid_argument("Can't find TGraph in FancyTTreeInterpolator created with text " + std::string(yAxisText.Data()));
  }
  else{
    return fStringToGraph.find(yAxisText)->second;
  }
  return NULL;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get interpolated yAxisText variable at xAxisValue (time).
 *
 * @param yAxisText is defines the y-axis of the TGraph, drawn against fXaxisText
 * @param xAxisValue is the value of the fXaxisText branch to interpolate the yAxisText branch variable at.
 * @return the interpolated value.
 */
Double_t icemc::FancyTTreeInterpolator::interp(TString yAxisText, Double_t xAxisValue){
  /* This function handles the getting of values from the appropriate TGraph */

  TGraph* gr = get(yAxisText);

  if(xAxisValue >= fXmin && xAxisValue <= fXmax){
    Double_t tempVal = gr->Eval(xAxisValue);

    // do the unwrapping if required
    std::pair<Double_t, Double_t> wrap = fStringToWrapValues.find(yAxisText)->second;
    Double_t wrapHigh = wrap.first;
    Double_t wrapLow = wrap.second;
    
    if(wrapHigh != 0){
      while (tempVal < wrapLow){
        tempVal += (wrapHigh - wrapLow);
      }
      while (tempVal >= wrapHigh){
        tempVal -= (wrapHigh - wrapLow);
      }
    }
    return tempVal;
  }
  else{
    std::cerr << "Value " << xAxisValue << " lies outside range of " << fXAxisText << std::endl;
    std::cerr << "xMin = " << fXmin << ", xMax = " << fXmax << std::endl;
    throw std::domain_error("");
  }
  return 0;

}
