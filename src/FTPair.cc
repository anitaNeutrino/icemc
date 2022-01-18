#include "FTPair.h"
#include "TFile.h"

/** 
 * Print the state of the updater variables
 * 
 * @param funcName is supposed to be the __PRETTY_FUNCTION__ macro
 * @param fTime should be fNeedToUpdateTimeDomain
 * @param fFreq should be fNeedToUpdateFreqDomain
 */
void PRINT_STATE(const char* funcName, bool fTime,  bool fFreq){
  std::cerr << "In " << (funcName)
	    << ": fNeedToUpdateTimeDomain = " << fTime << ", "
	    << "fNeedToUpdateFreqDomain = " << fFreq  << std::endl;
}

/**
 * Print the state of the updater variables
 */
#define PRINT_STATE_IF_DEBUG()						               \
  if(fDebug){                                                                          \
    PRINT_STATE(__PRETTY_FUNCTION__, fNeedToUpdateTimeDomain, fNeedToUpdateFreqDomain);\
  }






icemc::FTPair::FTPair()
  : fTimeDomainGraph(),
    fFreqDomain(),
    fNeedToUpdateTimeDomain(false),
    fNeedToUpdateFreqDomain(true),
    fDebug(false),
    fDoNormalTimeDomainOrdering(false)
{
  PRINT_STATE_IF_DEBUG();

  // we need to put some kind of state in here as
  // otherwise things can go very wrong.
  const int n = 8;
  const double t0 = 0;
  const double dt = 1e-9;
  for(int i=0; i < n; i++){
    fTimeDomainGraph.SetPoint(i, t0 + dt*i, 0);
  }

  zeroPadTimeDomainLengthToPowerOf2();
}


icemc::FTPair::FTPair(int n, const double* timeDomainAmplitude,  double dt, double t0)
  : fTimeDomainGraph(n, timeDomainAmplitude, timeDomainAmplitude),
    fFreqDomain(),
    fNeedToUpdateTimeDomain(false),
    fNeedToUpdateFreqDomain(true),
    fDebug(false),
    fDoNormalTimeDomainOrdering(false)
{
  PRINT_STATE_IF_DEBUG();
  for(int i=0; i < n; i++){
    fTimeDomainGraph.GetX()[i] = t0 + dt*i;
  }

  zeroPadTimeDomainLengthToPowerOf2();
}



icemc::FTPair::FTPair(const std::vector<double>& timeDomainAmplitude, double dt, double t0)
  : fTimeDomainGraph(timeDomainAmplitude.size(), &timeDomainAmplitude[0], &timeDomainAmplitude[0]),
    fFreqDomain(),
    fNeedToUpdateTimeDomain(false),
    fNeedToUpdateFreqDomain(true),
    fDebug(false),
    fDoNormalTimeDomainOrdering(false)
{
  PRINT_STATE_IF_DEBUG();
  
  for(int i=0; i < timeDomainAmplitude.size(); i++){
    fTimeDomainGraph.GetX()[i] = t0 + dt*i;
  }
  zeroPadTimeDomainLengthToPowerOf2();
}



icemc::FTPair::FTPair(const TGraph& grTimeDomain)
  : fTimeDomainGraph(grTimeDomain),
    fFreqDomain(),
    fNeedToUpdateTimeDomain(false),
    fNeedToUpdateFreqDomain(true),
    fDebug(false),
    fDoNormalTimeDomainOrdering(false)
{
  PRINT_STATE_IF_DEBUG();
  zeroPadTimeDomainLengthToPowerOf2();
}








icemc::FTPair::FTPair(const std::vector<std::complex<double> >& freqDomainPhasors, double df, bool doNormalTimeOrdering, double t0)
  : fTimeDomainGraph(),
    fFreqDomain(freqDomainPhasors),
    fNeedToUpdateTimeDomain(true),
    fNeedToUpdateFreqDomain(false),
    fDebug(false),
    fDoNormalTimeDomainOrdering(doNormalTimeOrdering)    
{
  PRINT_STATE_IF_DEBUG();
  int nf = zeroPadFreqDomainSoTimeDomainLengthIsPowerOf2(df);
  double dt = getDeltaT(nf, df);
  int nt = getNumTimes(nf);
  fTimeDomainGraph.Set(nt);
  for(int i=0; i < nt; i++){
    fTimeDomainGraph.GetX()[i] = t0 + dt*i;
  }
}


icemc::FTPair::FTPair(int nf, const std::complex<double>* freqDomainPhasors, double df, bool doNormalTimeOrdering, double t0)
  : fTimeDomainGraph(),
    fFreqDomain(freqDomainPhasors, freqDomainPhasors+nf),
    fNeedToUpdateTimeDomain(true),
    fNeedToUpdateFreqDomain(false),
    fDebug(false),
    fDoNormalTimeDomainOrdering(doNormalTimeOrdering)
{
  PRINT_STATE_IF_DEBUG();
  nf = zeroPadFreqDomainSoTimeDomainLengthIsPowerOf2(df);
  double dt = getDeltaT(nf, df);
  int nt = getNumTimes(nf);
  fTimeDomainGraph.Set(nt);
  for(int i=0; i < nt; i++){
    fTimeDomainGraph.GetX()[i] = t0 + dt*i;
  }
}







const TGraph& icemc::FTPair::getTimeDomain() const{
  PRINT_STATE_IF_DEBUG();
  maybeUpdateTimeDomain();
  return fTimeDomainGraph;
}

TGraph& icemc::FTPair::changeTimeDomain() {
  PRINT_STATE_IF_DEBUG();
  maybeUpdateTimeDomain();
  fNeedToUpdateFreqDomain = true;
  return fTimeDomainGraph;
}

const std::vector<std::complex<double> >& icemc::FTPair::getFreqDomain() const {
  PRINT_STATE_IF_DEBUG();
  maybeUpdateFreqDomain();
  return fFreqDomain;
}

std::vector<std::complex<double> >& icemc::FTPair::changeFreqDomain() {
  PRINT_STATE_IF_DEBUG();
  maybeUpdateFreqDomain();
  fNeedToUpdateTimeDomain = true;
  return fFreqDomain;
}


void icemc::FTPair::forceUpdateTimeDomain() const {
  PRINT_STATE_IF_DEBUG();
  fNeedToUpdateTimeDomain = true;
  maybeUpdateTimeDomain();
};

void icemc::FTPair::forceUpdateFreqDomain() const {
  PRINT_STATE_IF_DEBUG();
  fNeedToUpdateFreqDomain = true;
  maybeUpdateFreqDomain();
};



TGraph icemc::FTPair::makePowerSpectralDensityGraph() const {
  TGraph grPowerSpectrum = makePowerSpectrumGraph();
  double df = grPowerSpectrum.GetX()[1] - grPowerSpectrum.GetX()[0];
  for(int j=0; j < grPowerSpectrum.GetN(); j++){
    grPowerSpectrum.GetY()[j]/=df;
  }
  return grPowerSpectrum;
}


TGraph icemc::FTPair::makePowerSpectrumGraph() const {

  PRINT_STATE_IF_DEBUG();
  maybeUpdateFreqDomain();

  double df = getDeltaF();
  const TGraph& grT = getTimeDomain();
  const double dt = grT.GetX()[1] - grT.GetX()[0];
  const int N = grT.GetN();

  // the factor of two comes from using the numerical methods FFT
  // @todo if you implement FFTW change this!
  double scaleFactor = 2.0*dt/N;

  TGraph grF(fFreqDomain.size());

  for(int j=0; j < grF.GetN(); j++){
    grF.GetY()[j] = scaleFactor*std::norm(fFreqDomain[j]);
    grF.GetX()[j] = df*j;
  }

  if(fDebug){
    // let's check we've got a normalization consistent with Parseval's Theorem
    // http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf

    double sumSqTime = 0;
    // double dt = grT.GetX()[1] - grT.GetX()[0];
    for(int i=0; i < grT.GetN(); i++){
      sumSqTime += grT.GetY()[i]*grT.GetY()[i]*dt;
    }

    double sumSqFreq = 0;
    for(int j=0; j < grF.GetN(); j++){
      // std::norm gives the square of the magnitude,
      // so don't square here
      sumSqFreq += grF.GetY()[j];
    }

    const double fractionalDifferenceInPowers = (sumSqTime - sumSqFreq)/sumSqTime;
    const double tolerance = 1e-10;
    if(TMath::Abs(fractionalDifferenceInPowers) > tolerance){
      std::cerr << "Debug info in "  << __PRETTY_FUNCTION__
		<< ", fractional difference in powers is " << (sumSqTime - sumSqFreq)/sumSqTime
		<< ", time domain power = " << sumSqTime
		<< ", freq domain power = " << sumSqFreq << std::endl;
    }
  }

  return grF;
}





void icemc::FTPair::doNormalTimeDomainOrdering() const {
  PRINT_STATE_IF_DEBUG();
  
  TGraph grOld = getTimeDomain();
  TGraph& grNew = fTimeDomainGraph;

  const double t0 = grOld.GetX()[0];
  const double dt = grOld.GetX()[1] - t0;
  const int offset = grOld.GetN()/2;
  const double timeOffset = -offset*dt;
  for(int i=0; i < offset; i++){
    grNew.GetX()[i] += timeOffset;
    grNew.GetY()[i] = grOld.GetY()[i+offset];
  }
  for(int i=offset; i < grNew.GetN(); i++){
    grNew.GetX()[i] += timeOffset;
    grNew.GetY()[i] = grOld.GetY()[i-offset];
  }
  fNeedToUpdateFreqDomain = true;
}





void icemc::FTPair::applyConstantGroupDelay(double delaySeconds, bool circularDelay){
  PRINT_STATE_IF_DEBUG();

  if(delaySeconds==0) {
    // then we don't do anything
    return;
  }

  if(circularDelay){
    double freqHz = 0;
    double dfHz = getDeltaF();
    for(auto& c : changeFreqDomain()){
      // for each frequency bin, a phase delay of pi radians delays that frequency by 1./f
      const double phaseShift = TMath::TwoPi()*delaySeconds*freqHz;
      std::complex<double> phasor = std::polar(1.0,  phaseShift);
      c *= phasor;
      freqHz += dfHz;
    }
  }
  else{
    // For a non-circular delay, the strategy is this:
    // 1. construct a longer graph,
    // 2. apply the group delay,
    // 3. copy the portion of the extended graph corresponding to the original back
    // 4. profit

    const int n = getNumTimes();
    const double T = n*getDeltaT();
    if(TMath::Abs(delaySeconds) >= T){
      if(fDebug){
	std::cerr << "Requested non-circular group delay of " << delaySeconds
		  << " but whole graph is only " << T << " long. "
		  << " This results in an empty graph!" << std::endl;
      }
      TGraph& grChange = changeTimeDomain();
      for(int i=0; i < grChange.GetN(); i++){
	grChange.GetY()[i] = 0;
      }
      return;
    }

    const int padFactor = 2;
    const int n2 = padFactor*n;
    TGraph gr2(n2); // padded graph
    const TGraph& gr = getTimeDomain();
    double t0 = gr.GetX()[0];
    const double dt = gr.GetX()[1] - t0;
    if(delaySeconds > 0){
      // shifting forwards in time, so pad at the back
      for(int i=0; i < n; i++){
	gr2.SetPoint(i, t0 + i*dt, gr.GetY()[i]);
      }
      for(int i=n; i < n2; i++){
	gr2.SetPoint(i, t0 + i*dt, 0);
      }
    }
    else { // pad at the front
      t0 = t0 - T;
      for(int i=0; i < n; i++){
	gr2.SetPoint(i, t0 + i*dt, 0);
      }
      for(int i=n; i < n2; i++){
	gr2.SetPoint(i, t0 + i*dt, gr.GetY()[i-n]);
      }
    }

    FTPair extended(gr2);
    extended.applyConstantGroupDelay(delaySeconds, true);
    const TGraph& grDelayed = extended.getTimeDomain();
    int offset = delaySeconds > 0 ? 0 : n;

    TGraph& grChange = changeTimeDomain();
    for(int i=0; i < n; i++){
      grChange.GetY()[i] = grDelayed.GetY()[i + offset];
    }
  }
}





int icemc::FTPair::zeroPadFreqDomainSoTimeDomainLengthIsPowerOf2(double df) const {
  PRINT_STATE_IF_DEBUG();
  
  // assume freq domain vector is not dirty, but everything else is...

  int nt = FTPair::getNumTimes(fFreqDomain.size()); 
  if(!isPowerOf2(nt)){
    int newNt = nextPowerOf2(nt);


    int newNf = FTPair::getNumFreqs(newNt);
    if(fDebug){
      std::cerr << "Info in " << __PRETTY_FUNCTION__ << ": "
		<< "Padding freq domain: "
		<< fFreqDomain.size() << "->" << newNf << std::endl;
    }

    while(fFreqDomain.size() < newNf){
      fFreqDomain.push_back(0);
    }

    /**
     * OK, here's an FFT normalization subtlety.
     * Without any normalization, a forward FFT then an inverse FFT
     * would give you the original signal scaled up by some factor.
     * For the numerical methods implementation, this factor is N/2.
     * This is accouted for as a scaleFactor=2/N in maybeUpdateTimeDomain(),
     * and scaleFactor=2/N in makePowerSpectrumGraph.
     *
     * This is all fine if you don't do any padding in the frequency domain.
     * However, if you have a time domain signal, FFT for the set of complex
     * frequencies and zero pad the frequencies, then do an inverse FFT you
     * need to account for the fact that your inverse FFT will scale down your
     * signal by a different amount than it was scaled up by the forward FFT.
     *
     * i.e. I'm treating time domain power as definitive, and you only ever
     * see the frequency domain values scaled by some implied scale factor.
     *
     * This is a long winded way of saying that the complex frequencies
     * need to get scaled as you pad, which is done here.
     */
    double scaleFactorInPadding = double(newNt)/nt;
    for(auto& c : fFreqDomain){
      c *= scaleFactorInPadding;
    }
    
    // make sure there are at least 2 points in the time graph...
    
    if(fTimeDomainGraph.GetN() < 2){
      fTimeDomainGraph.Set(2);
    }
    // so we can encode deltaT as the difference between the first two points

    double dt = FTPair::getDeltaT(newNf, df);
    fTimeDomainGraph.GetX()[1] = fTimeDomainGraph.GetX()[0] + dt;
  }
  
  return fFreqDomain.size();
}


 
int icemc::FTPair::zeroPadTimeDomainLengthToPowerOf2() const {
  PRINT_STATE_IF_DEBUG();
  
  // asssume unpadded time domain is good, everything else is dirty
  
  int n  = fTimeDomainGraph.GetN();
  if(isPowerOf2(n)){
    return n;
  }

  int n2 = nextPowerOf2(n);

  if(fDebug){
    std::cerr << "Info in " << __PRETTY_FUNCTION__ << ": " 
	      << "Padding time domain waveform to next power of 2: "
	      << n << "->" << n2 << std::endl;
  }
  double t0 = fTimeDomainGraph.GetX()[0];
  double dt = fTimeDomainGraph.GetX()[1] - t0;
  for(int i=n; i < n2; i++){
    fTimeDomainGraph.SetPoint(i, t0 + dt*i, 0);
  }

  return n2;
}





void icemc::FTPair::maybeUpdateFreqDomain() const {
  PRINT_STATE_IF_DEBUG();
  
  if(fNeedToUpdateFreqDomain){

    // time domain has been mucked around with
    // ... so it could be any length, let's sort that out
    zeroPadTimeDomainLengthToPowerOf2();


    std::vector<double> temp(fTimeDomainGraph.GetY(), fTimeDomainGraph.GetY() + fTimeDomainGraph.GetN());
    // temp.reserve(fTimeDomainGraph.GetN());

    // // Forward and then inv FFT scales output by N/2
    // // here we take that out by scaling down just before the forward
    // int n = fTimeDomainGraph.GetN();
    // double scaleFactor = 2./n; // numerical recipes scaling
    // for(int i=0; i < n; i++){
    //   temp.push_back(fTimeDomainGraph.GetY()[i]*scaleFactor);
    // }

    realft(&temp[0], 1,  temp.size());
    fFreqDomain.clear();
    int nf = getNumFreqs(fTimeDomainGraph.GetN());
    fFreqDomain.reserve(nf);

    // unpack the 0th bin...
    fFreqDomain.push_back(std::complex<double>(temp.at(0), 0));
    int j2=2;
    for(int j=1; j < nf-1; j++){
      fFreqDomain.push_back(std::complex<double>(temp.at(j2), temp.at(j2+1)));
      j2 += 2;
    }
    // ... and the nyquist bin
    fFreqDomain.push_back(std::complex<double>(temp.at(1), 0));
    fNeedToUpdateFreqDomain = false;
  }
}



void icemc::FTPair::maybeUpdateTimeDomain() const {
  PRINT_STATE_IF_DEBUG();

  if(fNeedToUpdateTimeDomain){
    // if(fDebug){
    //   std::cerr << __PRETTY_FUNCTION__ << ", updating time domain!"  << std::endl;
    // }

    int n = fTimeDomainGraph.GetN();
    double dt = fTimeDomainGraph.GetX()[1] - fTimeDomainGraph.GetX()[0];
    double df = getDeltaF(n, dt);

    // frequency domain is dirty! it could be of any length.
    // let's make sure it's an accepatable length for a dft.
    int nf = zeroPadFreqDomainSoTimeDomainLengthIsPowerOf2(df);
    int nNew = getNumTimes(nf); // this the the required power of 2

    // Here we follow the conventional place for the normalization as used in icemc.
    // Forward and then inv FFT scales output by N/2
    // here we take that out by scaling down just before the inverse FT
    double scaleFactor = 2./nNew; // numerical recipes scaling

    std::vector<double> temp;
    temp.reserve(nNew);
    temp.push_back(fFreqDomain.at(0   ).real()*scaleFactor); // DC offset
    temp.push_back(fFreqDomain.at(nf-1).real()*scaleFactor); // nqyuist

    for(int j=1; j < nf-1; j++){
      temp.push_back(fFreqDomain.at(j).real()*scaleFactor);
      temp.push_back(fFreqDomain.at(j).imag()*scaleFactor);
    }
    realft(&temp[0], -1, temp.size());

    fTimeDomainGraph.Set(nNew);
    double t0 = fTimeDomainGraph.GetX()[0];
    double dtNew = fTimeDomainGraph.GetX()[1] - t0;
    for(int i=0;  i < nNew; i++){
      fTimeDomainGraph.SetPoint(i, t0 + i*dtNew,  temp.at(i));
    }

    fNeedToUpdateTimeDomain = false;
    
    if(fDoNormalTimeDomainOrdering){
      doNormalTimeDomainOrdering();
      // Unset the flag since it probably only makes sense to do this on the first inverse FT
      fDoNormalTimeDomainOrdering = false; 
    }
  }
}


/** 
 * @brief Swap the values of a and b
 * 
 * Used in these home-rolled FFTs.
 * 
 * @param a will take the value of b
 * @param b will take the value of a
 */
inline void SWAP(double &a, double &b){
  double dum = a;
  a = b;
  b = dum;
}

/** 
 * These functions seem to have been taken from the numerical recipes handbook (google it).
 * The conventions seem to be different from the FFTW we use elsewhere in ANITA.
 * NOTE: nsize MUST be a power of 2 otherwise this function will corrupt program memory!!!
 * 
 * @param data is the array to transform
 * @param isign is the direction, +1 is forward FT, -1 is inverse FT
 * @param nsize is the the length of data: MUST BE A POWER OF 2!
 */
void four1(double *data, const int isign, int nsize) {

  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
    
  int nn=nsize/2;
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j-1],data[i-1]);
      SWAP(data[j],data[i]);
    }
    m=nn;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j-1]-wi*data[j];
	tempi=wr*data[j]+wi*data[j-1];
	data[j-1]=data[i-1]-tempr;
	data[j]=data[i]-tempi;
	data[i-1] += tempr;
	data[i] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}


/** 
 * These functions seem to have been taken from the numerical recipes handbook (google it).
 * The conventions seem to be different from the FFTW we use elsewhere in ANITA.
 * NOTE: nsize MUST be a power of 2 otherwise this function will corrupt program memory!!!
 * 
 * @param data is the array to transform
 * @param isign is the direction, +1 is forward FT, -1 is inverse FT
 * @param nsize is the the length of data: MUST BE A POWER OF 2!
 */

void icemc::FTPair::realft(double *data, const int isign, int nsize){
  int i, i1, i2, i3, i4;
  double c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;
  theta=3.141592653589793238/(nsize>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,1,nsize);
  }
  else {
    c2 = 0.5;
    theta = -theta;
  }
  wtemp = sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  wr = 1.0+wpr;
  wi = wpi;
  for (i=1;i<(nsize>>2);i++) {
    i2=1+(i1=i+i);
    i4=1+(i3=nsize-i1);
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r= -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4]= -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[0] = (h1r=data[0])+data[1];
    data[1] = h1r-data[1];
  }
  else {
    data[0]=c1*((h1r=data[0])+data[1]);
    data[1]=c1*(h1r-data[1]);
    four1(data,-1,nsize);
  }
}



void icemc::FTPair::dump(const char* fileName) const {
  
  TFile* f= new TFile(fileName, "recreate");

  TGraph grTime = getTimeDomain();
  grTime.SetName("grTimeDomain");
  grTime.SetTitle("Time domain");
  grTime.Write();


  TGraph grFreq = makePowerSpectralDensityGraph();
  grFreq.SetName("grFreqDomain");
  grFreq.SetTitle("Power spectral density");
  grFreq.Write();
  
  f->Write();
  f->Close();
  delete f;  
}
