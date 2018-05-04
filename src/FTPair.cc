#include "FTPair.h"





icemc::FTPair::FTPair(int n, const double* timeDomainAmplitude,  double dt, double t0)
  : fTimeDomainGraph(n, timeDomainAmplitude, timeDomainAmplitude),
    fFreqDomain(),
    fNeedToUpdateTimeDomain(false),
    fNeedToUpdateFreqDomain(true),
    fDebug(false)
{
  for(int i=0; i < n; i++){
    fTimeDomainGraph.GetX()[i] = t0 + dt*i;
  }

  zeroPadTimeDomainToNearestPowerOf2();
  fNeedToUpdateTimeDomain = false;
  fNeedToUpdateFreqDomain = true;
}



icemc::FTPair::FTPair(const std::vector<double>& timeDomainAmplitude, double dt, double t0)
  : fTimeDomainGraph(timeDomainAmplitude.size(), &timeDomainAmplitude[0], &timeDomainAmplitude[0]),
    fFreqDomain(),
    fNeedToUpdateTimeDomain(false),
    fNeedToUpdateFreqDomain(true),
    fDebug(false)
{
  for(int i=0; i < timeDomainAmplitude.size(); i++){
    fTimeDomainGraph.GetX()[i] = t0 + dt*i;
  }
  zeroPadTimeDomainToNearestPowerOf2();
  fNeedToUpdateTimeDomain = false;
  fNeedToUpdateFreqDomain = true;
}



icemc::FTPair::FTPair(const TGraph& grTimeDomain)
  : fTimeDomainGraph(grTimeDomain),
    fFreqDomain(),
    fNeedToUpdateTimeDomain(false),
    fNeedToUpdateFreqDomain(true),
    fDebug(false)
{
  zeroPadTimeDomainToNearestPowerOf2();
}








const TGraph& icemc::FTPair::getTimeDomain() const{
  maybeUpdateTimeDomain();
  return fTimeDomainGraph;
}

TGraph& icemc::FTPair::changeTimeDomain() {
  maybeUpdateTimeDomain();
  fNeedToUpdateFreqDomain = true;
  return fTimeDomainGraph;
}

const std::vector<std::complex<double> >& icemc::FTPair::getFreqDomain() const {
  maybeUpdateFreqDomain();
  return fFreqDomain;
}

std::vector<std::complex<double> >& icemc::FTPair::changeFreqDomain() {
  maybeUpdateFreqDomain();
  fNeedToUpdateTimeDomain = true;
  return fFreqDomain;
}


void icemc::FTPair::forceUpdateTimeDomain() const {
  fNeedToUpdateTimeDomain = true;
  maybeUpdateTimeDomain();
};

void icemc::FTPair::forceUpdateFreqDomain() const {
  fNeedToUpdateFreqDomain = true;
  maybeUpdateFreqDomain();
};

TGraph icemc::FTPair::makePowerSpectrumGraph() const {
  maybeUpdateFreqDomain();
  double dt = fTimeDomainGraph.GetX()[1] - fTimeDomainGraph.GetX()[0];
  int n = fTimeDomainGraph.GetN();
  double df;
  getFreqInfo(n, dt, df);
  TGraph gr(fFreqDomain.size());
  for(int j=0; j < gr.GetN(); j++){
    gr.GetY()[j] = std::norm(fFreqDomain[j]);
    gr.GetX()[j] = df*j;
  }
  return gr;
}



















int icemc::FTPair::zeroPadTimeDomainToNearestPowerOf2() const {
  int n  = fTimeDomainGraph.GetN();
  int n2 = nextPowerOf2(n);
  if(n!=n2){
    std::cerr << "Info in " << __PRETTY_FUNCTION__ << ": " 
	      << "Padding time domain waveform to next power of 2: "
	      << n << "->" << n2 << std::endl;
    double t0 = fTimeDomainGraph.GetX()[0];
    double dt = fTimeDomainGraph.GetX()[1] - t0;
    for(int i=n; i < n2; i++){
      fTimeDomainGraph.SetPoint(i, t0 + dt*i, 0);
    }
  }
  return n2;
}





void icemc::FTPair::maybeUpdateFreqDomain() const {
  if(fNeedToUpdateFreqDomain){
    if(fDebug){
      std::cerr << __PRETTY_FUNCTION__ << ", updating frequency domain!"  << std::endl;
    }

    // time domain has been mucked around with
    // ... so it could be any length, let's sort that out
    zeroPadTimeDomainToNearestPowerOf2();
    
    double df;
    int nf = getFreqInfo(fTimeDomainGraph.GetN(), fTimeDomainGraph.GetX()[1] - fTimeDomainGraph.GetX()[0], df);
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
  if(fNeedToUpdateTimeDomain){
    if(fDebug){
      std::cerr << __PRETTY_FUNCTION__ << ", updating time domain!"  << std::endl;
    }

    int n = fTimeDomainGraph.GetN();
    double dt = fTimeDomainGraph.GetX()[1] - fTimeDomainGraph.GetX()[0];
    double df;
    int nf = getFreqInfo(n, dt,df);


    // frequency domain is dirty! it could be of any length.
    // let's make sure it's the right length!
    if(fFreqDomain.size()!=nf){
      if(fFreqDomain.size()!=0){
	std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ": unexpected number of frequency bins!";
      }
      
      if(fFreqDomain.size() < nf ){
	if(fDebug) {
	  std::cerr << " Zero padding frequency bins " << fFreqDomain.size() << "->" << nf << std::endl;
	}
	fFreqDomain.reserve(nf);
	while(fFreqDomain.size() < nf){
	  fFreqDomain.push_back(0);
	}
      }
      else {
	if(fDebug) {	
	  std::cerr << " Truncating frequency bins " << fFreqDomain.size() << "->" << nf << std::endl;
	}
	fFreqDomain.resize(nf);
      }
    }

    std::vector<double> temp;
    temp.reserve(n);
    temp.push_back(fFreqDomain.at(0).real()); // DC offset
    temp.push_back(fFreqDomain.at(nf-1).real()); // nqyuist
    for(int j=1; j < nf-1; j++){
      temp.push_back(fFreqDomain.at(j).real());
      temp.push_back(fFreqDomain.at(j).imag());
    }
    realft(&temp[0], -1, temp.size());

    for(int i=0;  i < n; i++){
      fTimeDomainGraph.GetY()[i] = temp.at(i);
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


