void supernovaA3()
{

  TFile of("supernovaA3.root","RECREATE"); 

  std::string name;
  double ra; 
  double dec; 
  std::string type; 
  int tdiscovery;
  int texplode;

  TTree * snTree = new TTree("sn","sn");

  snTree->Branch("name",&name);
  snTree->Branch("ra",&ra);
  snTree->Branch("dec",&dec);
  snTree->Branch("type",&type);
  snTree->Branch("tdiscovery",&tdiscovery);
  snTree->Branch("texplode",&texplode);

  FILE * f = fopen("supernovaA3.tsv","r"); 
  char buf[128]; 
  while (!feof(f))
  {
    fscanf(f,"%s", buf); 
    name = buf; 

    int d,min; 
    float fsec;
    fscanf(f,"%d:%d:%f", &d,&min,&fsec); 
    ra = d + min/60. + fsec/3600.; 
    fscanf(f,"%d:%d:%f", &d,&min,&fsec); 
    dec = abs(d) + min/60. + fsec/3600.; 
    if (d < 0) dec*=-1; 
    fscanf(f,"%s", buf); 
    type = buf; 
    int year,month,day,hour,sec; 

    fscanf(f,"\t\"%04d-%02d-%02d %02d:%02d:%02d\"", &year,&month,&day,&hour,&min,&sec); 
    printf("\t\"%04d-%02d-%02d %02d:%02d:%02d\"\n", year,month,day,hour,min,sec); 
    TTimeStamp ts_discovery(year,month,day,hour,min,sec); 
    tdiscovery = ts_discovery.GetSec(); 
    fscanf(f,"\t\"%04d-%02d-%02d %02d:%02d:%02d\"\n", &year,&month,&day,&hour,&min,&sec); 
    printf("\t\"%04d-%02d-%02d %02d:%02d:%02d\"\n", year,month,day,hour,min,sec); 
    TTimeStamp ts_explode(year,month,day,hour,min,sec); 
    texplode = ts_explode.GetSec(); 
    printf("%s %g %g %s %d %d\n", name.c_str(), ra, dec, type.c_str(), tdiscovery, texplode); 
        
    snTree->Fill(); 
  }

  snTree->Write(); 

}
