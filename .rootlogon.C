{
  gSystem->AddIncludePath("include/");
  gSystem->Load("../../build/components/icemc/libAnitaIceMC.so");
  gROOT->ProcessLine("#include \"FTPair.h\"");
}
