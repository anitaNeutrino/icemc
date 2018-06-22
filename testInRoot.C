void makeSineWave(std::vector<double>& x, int n, double dt, double f){

  x.clear();
  x.reserve(n);
  // double f = 1/period;
  double omega = TMath::TwoPi()*f;
  for(int i=0; i < n; i++){
    double t = dt*i;
    x.push_back(sin(omega*t));
  }
  // for(int i=0; i < n; i++){
  //   if((i%2)==0){
  //     x.push_back(sin(omega*dt*i));
  //   }
  //   else{
  //     x.push_back(0);
  //   }
  // }
}








void testInRoot(){

  gSystem->AddIncludePath("include/");
  gSystem->Load("../../build/components/icemc/libAnitaIceMC.so");
  gROOT->ProcessLine("#include \"FTPair.h\"");

}
