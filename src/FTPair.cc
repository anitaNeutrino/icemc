#include "FTPair.h"

/** 
 * Print the state of the updater variables
 * 
 * @param funcName must be cout-able, is supposed to be the __PRETTY_FUNCTION__ macro
 */

#define PRINT_STATE_IF_DEBUG(funcName)                                               \
  if(fDebug){                                                                        \
    std::cerr << "In " << (funcName)						     \
              << ": fMustUpdateTimeDomain = " << fNeedToUpdateTimeDomain << ", "     \
              << "fMustUpdateFreqDomain = " << fNeedToUpdateFreqDomain << std::endl; \
   }



icemc::FTPair::FTPair(int n, const double* timeDomainAmplitude,  double dt, double t0)
  : fTimeDomainGraph(n, timeDomainAmplitude, timeDomainAmplitude),
    fFreqDomain(),
    fNeedToUpdateTimeDomain(false),
    fNeedToUpdateFreqDomain(true),
    fDebug(false)
{
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
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
    fDebug(false)
{
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  
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
    fDebug(false)
{
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__) 
  zeroPadTimeDomainLengthToPowerOf2();
}








icemc::FTPair::FTPair(const std::vector<std::complex<double> >& freqDomainPhasors, double df)
  : fTimeDomainGraph(),
    fFreqDomain(freqDomainPhasors),
    fNeedToUpdateTimeDomain(true),
    fNeedToUpdateFreqDomain(false),
    fDebug(false)
{
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  zeroPadFreqDomainSoTimeDomainLengthIsPowerOf2(df);
}



icemc::FTPair::FTPair(int nf, const std::complex<double>* freqDomainPhasors, double df)
  : fTimeDomainGraph(),
    fFreqDomain(freqDomainPhasors, freqDomainPhasors+nf),
    fNeedToUpdateTimeDomain(true),
    fNeedToUpdateFreqDomain(false),
    fDebug(false)
{
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  zeroPadFreqDomainSoTimeDomainLengthIsPowerOf2(df);
}







const TGraph& icemc::FTPair::getTimeDomain() const{
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  maybeUpdateTimeDomain();
  return fTimeDomainGraph;
}

TGraph& icemc::FTPair::changeTimeDomain() {
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  maybeUpdateTimeDomain();
  fNeedToUpdateFreqDomain = true;
  return fTimeDomainGraph;
}

const std::vector<std::complex<double> >& icemc::FTPair::getFreqDomain() const {
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  maybeUpdateFreqDomain();
  return fFreqDomain;
}

std::vector<std::complex<double> >& icemc::FTPair::changeFreqDomain() {
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  maybeUpdateFreqDomain();
  fNeedToUpdateTimeDomain = true;
  return fFreqDomain;
}


void icemc::FTPair::forceUpdateTimeDomain() const {
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  fNeedToUpdateTimeDomain = true;
  maybeUpdateTimeDomain();
};

void icemc::FTPair::forceUpdateFreqDomain() const {
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  fNeedToUpdateFreqDomain = true;
  maybeUpdateFreqDomain();
};

TGraph icemc::FTPair::makePowerSpectrumGraph() const {
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  maybeUpdateFreqDomain();
  double df = getDeltaF();

  TGraph gr(fFreqDomain.size());
  for(int j=0; j < gr.GetN(); j++){
    gr.GetY()[j] = std::norm(fFreqDomain[j]);
    gr.GetX()[j] = df*j;
  }
  return gr;
}














int icemc::FTPair::zeroPadFreqDomainSoTimeDomainLengthIsPowerOf2(double df) const {
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  
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
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  
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
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)
  
  if(fNeedToUpdateFreqDomain){

    // time domain has been mucked around with
    // ... so it could be any length, let's sort that out
    zeroPadTimeDomainLengthToPowerOf2();

    
    std::vector<double> temp;
    temp.reserve(fTimeDomainGraph.GetN());

    // Forward and then inv FFT scales output by N/2
    // here we take that out by scaling down just before the forward
    int n = fTimeDomainGraph.GetN();
    double scaleFactor = 2./n; // numerical recipes scaling
    for(int i=0; i < n; i++){
      temp.push_back(fTimeDomainGraph.GetY()[i]*scaleFactor);
    }

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
  PRINT_STATE_IF_DEBUG(__PRETTY_FUNCTION__)

  if(fNeedToUpdateTimeDomain){
    if(fDebug){
      std::cerr << __PRETTY_FUNCTION__ << ", updating time domain!"  << std::endl;
    }

    int n = fTimeDomainGraph.GetN();
    double dt = fTimeDomainGraph.GetX()[1] - fTimeDomainGraph.GetX()[0];
    double df = getDeltaF(n, dt);

    // frequency domain is dirty! it could be of any length.
    // let's make sure it's an accepatable length for a dft.
    int nf = zeroPadFreqDomainSoTimeDomainLengthIsPowerOf2(df);
    int nNew = getNumTimes(nf); // this the the required power of 2

    std::vector<double> temp;
    temp.reserve(nNew);
    temp.push_back(fFreqDomain.at(0).real()); // DC offset
    temp.push_back(fFreqDomain.at(nf-1).real()); // nqyuist
    for(int j=1; j < nf-1; j++){
      temp.push_back(fFreqDomain.at(j).real());
      temp.push_back(fFreqDomain.at(j).imag());
    }
    realft(&temp[0], -1, temp.size());


    double t0 = fTimeDomainGraph.GetX()[0];
    double dtNew = fTimeDomainGraph.GetX()[1] - t0;
    for(int i=0;  i < nNew; i++){
      fTimeDomainGraph.GetY()[i] = temp.at(i);
      fTimeDomainGraph.GetX()[i] = t0 + i*dtNew;
    }

    fNeedToUpdateTimeDomain = false;
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
void SWAP(double &a, double &b){
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

