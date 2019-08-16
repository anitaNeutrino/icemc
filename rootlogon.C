{


  gSystem->Load("${ANITA_UTIL_INSTALL_DIR}/lib/libRootFftwWrapper.so") && gSystem->Load("$ANITA_UTIL_INSTALL_DIR/lib/libRootFftwWrapper.dylib");
  gSystem->Load("${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaEvent.so") && gSystem->Load("$ANITA_UTIL_INSTALL_DIR/lib/libAnitaEvent.dylib"); 
  gSystem->Load("${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaCorrelator.so") && gSystem->Load("$ANITA_UTIL_INSTALL_DIR/lib/libAnitaCorrelator.dylib"); 
  gSystem->Load("${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaAnalysis.so") && gSystem->Load("$ANITA_UTIL_INSTALL_DIR/lib/libAnitaAnalysis.dylib"); 
  gSystem->Load("${ANITA_UTIL_INSTALL_DIR}/lib/libUCorrelator.so") && gSystem->Load("$ANITA_UTIL_INSTALL_DIR/lib/libUCorrelator.dylib"); 
//  gSystem->Load("build/libUCorrelator.so"); 


//  const Int_t NRGBs = 5;
//  const Int_t NCont = 255;
//  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
//  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
//  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
//  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
//  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(255);
  gStyle->SetPalette(kTemperatureMap);
}
