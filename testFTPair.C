void makeSineWave(std::vector<double>& x, int n, double dt, double f){

  x.clear();
  x.reserve(n);
  // double f = 1/period;
  double omega = TMath::TwoPi()*f;
  for(int i=0; i < n; i++){
    double t = dt*i;
    x.push_back(sin(omega*t));
  }
}



std::vector<double> makeGaussWave(){

  std::vector<double> x;
  const int n = 256;
  const double dt = 1e-9/2.6;
  x.clear();
  x.reserve(n);

  double t = 0;
  const double tMax=30e-9;
  const double sigmaMax = 3e-9;
  const double A = 100; //mV
  for(int i=0; i < n; i++){

    double y = A*exp(-(t - tMax)*(t - tMax)/(2*sigmaMax*sigmaMax));
    x.push_back(y);
    t += dt;
  }
  return x;
}


std::vector<double> makeDelta(int maxSamp = 10){

  std::vector<double> x;
  const int n = 256; //pow(2, 12);
  x.reserve(n);
  for(int i=0; i < n; i++){
    double y = i == maxSamp ? 1 : 0;
    x.push_back(y);
  }
  return x;
}


void testFTPair(){

  std::vector<EColor> cols {kBlack, kRed, kOrange, kYellow, kGreen, kCyan, kBlue, kViolet, kMagenta};

  const double dt = 0.5;
  const double maxTime = 10;
  
  icemc::FTPair p(makeDelta(TMath::Nint(maxTime/dt)), dt);
  auto c1 = new TCanvas();
  for(int i=0; i < cols.size(); i++){
    double delaySec = -10;
    p.applyConstantGroupDelay(delaySec);
    auto gr = new TGraph();
    *gr = p.getTimeDomain();

    TString opt = i == 0 ? "al" : "lsame";
    gr->SetLineColor(cols.at(i));
    gr->Draw(opt);
  }
}
