#include "Tools.h"
#include <iostream>
#include <cmath>
#include "TSpline.h"
#include "TH2F.h"
#include <fstream>
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Constants.h"


//using std::cout;
using namespace std;


void  Tools::MakeGraph(int index, const int n,double *time,double *volts,TGraph *&mygraph,TH2F *&h2, double scalex,double scaley,string xaxistitle,string yaxistitle) {
    
    double maxtime=-1.E20;
    double maxv=-1.E20;
    double mintime=1E20;
    double minv=1.E20;
    
    double timecopy[n];
    double voltscopy[n];
    
    for (int i=0;i<n;i++) {
        timecopy[i]=time[i];
        voltscopy[i]=volts[i];
        timecopy[i]*=scalex;
        voltscopy[i]*=scaley;
    }
    for (int i=0;i<n;i++) {
        if (timecopy[i]>maxtime)
            maxtime=timecopy[i];
        if (timecopy[i]<mintime)
            mintime=timecopy[i];
        if (voltscopy[i]>maxv)
            maxv=voltscopy[i];
        if (voltscopy[i]<minv)
            minv=voltscopy[i];
    }
    
    mygraph=new TGraph(n,timecopy,voltscopy);
    
    h2=new TH2F(Form("h2%d",index),"",10*n,mintime*1.1,maxtime*1.1,100,minv*1.1,maxv*1.1);
    h2->SetLineWidth(3);
    h2->SetXTitle(xaxistitle.c_str());
    h2->SetYTitle(yaxistitle.c_str());
}


int Tools::iSum(int* thisarray,int n) {
    
    int sum=0;
    for (int i=0;i<n;i++) {
        sum+=thisarray[i];
    } //for
    return sum;
} //iSum


double Tools::getMaxMagnitude(vector<double> v) {
    double mag=0.;
    for (int i=0;i<(int)v.size();i++) {
        if (v[i]>mag)
            mag=v[i];
    }
    return mag;
}


void Tools::ShiftLeft(double *x,const int n,int ishift) {
    
    double x_temp[n];
    // shift the x array to the left by ishift bins and fill the gap with zeroes
    for (int i=0;i<n;i++) {
      x_temp[i]=x[i];
    }
    for (int i=0;i<n-ishift;i++) {
      x[i]=x_temp[i+ishift];
    }
    for (int i=n-ishift;i<n;i++) {
      x[i]=0.;
    }
}


void Tools::ShiftRight(double *x,const int n,int ishift) {
    
    double x_temp[n];
    // shift the x array to the right by ishift bins and fill the gap with zeroes
    for (int i=0;i<n;i++) {
      x_temp[i]=x[i];
    }
    for (int i=ishift;i<n;i++) {
      x[i]=x_temp[i-ishift];
    }
    for (int i=0;i<ishift;i++) {
      x[i]=0.;
    }
}


void Tools::realft(double *data, const int isign, int nsize){
    int i, i1, i2, i3, i4;
    double c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;
    theta=3.141592653589793238/(nsize>>1);
    if (isign == 1) {
        c2 = -0.5;
        four1(data,1,nsize);
    } else {
        c2=0.5;
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
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
    } else {
        data[0]=c1*((h1r=data[0])+data[1]);
        data[1]=c1*(h1r-data[1]);
        four1(data,-1,nsize);
    }
}


void Tools::SWAP(double &a, double &b){
    double dum = a;
    a = b;
    b = dum;
}


void Tools::four1(double *data, const int isign,int nsize) {
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


double Tools::dSquare(double *p) {
    return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
} //dSquare


int Tools::WhichIsMin(double *x,int n) {
    double min=1.E22;
    int imin=0;
    for (int i=0;i<n;i++) {
        if (x[i]<min) {
            min=x[i];
            imin=i;
        }
    } //for
    return imin;
} //WhichIsMin


void Tools::Print(double *p,int i) {
    for (int j=0;j<i;j++) {
        cout << p[j] << " ";
    } //for
    cout << "\n";
} //Print (double*,int)


void Tools::Print(int *p,int i) {
    for (int j=0;j<i;j++) {
        cout << p[j] << " ";
    } //for
    cout << "\n";
} //Print (double*,int)


void Tools::GetNextNumberAsString(ifstream& fin, ofstream& fout, string& number) {
    string temp;
    getline(fin,temp); // get next line of the input file
    
    fout << temp << "\n"; // output this line to the summary file
    
    int place=0;
    place=temp.find_first_of(" \t"); // find where the first space or tab is
    
    number=temp.substr(0,place); // everything up until the first space is what we're looking for
} //GetNextNumberAsString


void Tools::GetNumbersAsStringArray(ifstream& fin, ofstream& fout,vector<string>& vnumbers, int nelements) {
    string temp;
    //getline(fin,temp);
    
    //fout << temp << "\n";
    
    //  int place_previous=0;
    //int place_next;
    vnumbers.clear();
    string s;
    for (int n=0;n<nelements;n++) {
        fin >> s;
        fout << s << "\t"; // output this line to the summary file
        vnumbers.push_back(s);
        //    place_next=temp.find_first_of("\t",place_previous+1); // find where first tab is
        //vnumbers.push_back(temp.substr(place_previous,place_next-place_previous));
        //place_previous=place_next;
    }
    getline(fin,temp);
    fout << temp << "\n";
    //  cout << "temp is " << temp << "\n";
}


void Tools::GetNext2NumbersAsString(ifstream& fin,ofstream& fout,string& number1,string& number2, string& stherest) {
    
    string temp;
    getline(fin,temp); // get next line of the input file
    
    fout << temp << "\n"; // output this line to the summary file
    
    int place=0;
    place=temp.find_first_of(" \t"); // find where the first space is
    
    number1=temp.substr(0,place); // everything up until the first space is what we're looking for
    
    temp=temp.substr(place+1,temp.size());
    
    number2=temp.substr(0,temp.find_first_of(" "));
    
    stherest=temp.substr(2,temp.size());
} //GetNext2NumbersAsString


double Tools::GetFWHM(TH1 *h1) {
    int imax=h1->GetMaximumBin();
    double max=h1->GetMaximum();
    
    //  cout << "imax, max are " << imax << " " << max << "\n";
    int ibin_plus=0;
    int ibin_minus=0;
    // now step to the right until it's half
    for (int ibin=imax;ibin<=h1->GetNbinsX();ibin++) {
        if (h1->GetBinContent(ibin)<max/2.) {
            ibin_plus=ibin;
            ibin=h1->GetNbinsX()+1;
            //  cout << "ibin_plus is " << ibin_plus << "\n";
        }
    }
    // now step to the left
    for (int ibin=imax;ibin>=1;ibin--) {
        if (h1->GetBinContent(ibin)<max/2.) {
            ibin_minus=ibin;
            ibin=0;
            //  cout << "ibin_minus is " << ibin_minus << "\n";
        }
    }

    if (ibin_plus>0 && ibin_minus==0) {
        ibin_minus=1;
        //cout << "bin_minus is " << ibin_minus << "\n";
    }
    
    if (ibin_plus==0 && ibin_minus==0) {
        cout << "Found 0 FWHM.\n";
        return 0.;
    }
    
    return (h1->GetBinCenter(ibin_plus)-h1->GetBinCenter(ibin_minus))/2.;
}


int Tools::NonZero(double *anarray,int n) { // how many bins are nonzero
  int count=0;
  for (int i=0;i<n;i++) {
    cout << "i, anarray is " << i << "\t" << anarray[i] << "\n";
    if (anarray[i]!=0)
      count++;
  }
  return count;
}


void Tools::Zero(int *anarray,int n) {
    for (int i=0;i<n;i++) {
        anarray[i]=0;
    } //for
} //Zero (int*,int)


void Tools::Zero(double *anarray,int n) {
    for (int i=0;i<n;i++) {
        anarray[i]=0.;
    } //for
} //Zero (int*,int)


double Tools::dMinNotZero(const double *x,int n) {
    double min=dMax(x,n);
    if (min==0)
        cout << "max is 0.\n";
    for (int k=1;k<n;k++) {
        if (x[k]<min && x[k]!=0)
            min=x[k];
    }
    return min;
} //dMinNotZero(double*, int)


double Tools::dMin(const double *x,int n) {
    double min=x[0];
    for (int k=1;k<n;k++) {
        if (x[k]<min)
            min=x[k];
    }
    return min;
} //dMin(double*, int)


double Tools::dMin(double x,double y) {
    double min=1.E22;
    if (x<y)
        min=x;
    else
        min=y;
    
    return min;
} //dMin(double,double)


double Tools::dMax(const double *x,int n) {
    
    double max=x[0];
    for (int k=1;k<n;k++) {
        if (x[k]>max)
            max=x[k];
    }
    return max;
} //dMax(double*, int)


double Tools::dvMax(const vector<double> x) {
    
    double max=x[0];
    for (int k=1;k<(int)x.size();k++) {
        if (x[k]>max)
            max=x[k];
    }
    return max;
} //dMax(double*, int)


double Tools::dsMax(TSpline5 *sp) {
    vector<double> y;
    double maxn;
    double blah1,blah2;
    for (int i=0;i<sp->GetNp();i++) {
        sp->GetKnot(i,blah1,blah2);
        y.push_back(blah2);
    }
    maxn=Tools::dvMax(y);
    return maxn;
}


double Tools::dMax(double a,double b) {
    if (a>b)
        return a;
    else if (a<b)
        return b;
    else if (a==b)
        return a;
    return 0;
} //dMax(double,double


int Tools::Getifreq(double freq,double freq_low,double freq_high,int n) {
    
    if (freq>=freq_high)
        return -1;
    if (freq<freq_low)
        return -1;
    
    return (int)((freq-freq_low)/(freq_high-freq_low)*(double)n);
} //Getifreq


void Tools::InterpolateReal(double* array, const unsigned n){
    // to get rid of the zero bins
    double previous_nonzero = 0.;
    double next_nonzero = 0.;
    double check;
    unsigned ifirstnonzero = 0;
    unsigned ilastnonzero = 0;
    unsigned k;
    unsigned m = 0;
    unsigned count_nonzero = 0;
    
    // find the first nonzero even element
    while (!array[m]){
        m++;
    }
    
    ifirstnonzero = m;
    
    
    // count the nonzero elements
    for (unsigned i = 0; i < n; i++) {
        if (array[i] != 0)
            count_nonzero++;
    }
    
    if (count_nonzero) {
        // loop through the elements of the array and replace the zeros with interpolated values
        for (unsigned i = ifirstnonzero; i < n; i++) {
            if (array[i]) {
                // set the lower nonzero value that we are interpolating from
                previous_nonzero = array[i];
            }
            else{
                check = 0.;
                k = i;
                while (check == 0. && k < n) {
                    check = array[k];
                    k++;
                }
                if (k < n){
                    next_nonzero = check;
                    for (unsigned j = i; j < k; j++) {
                        array[j] = previous_nonzero + (next_nonzero - previous_nonzero) * double(j - (i - 1)) / double(k - (i - 1));
                    }
                    i = k - 1;
                    previous_nonzero = next_nonzero;
                }
                else {
                    ilastnonzero = i - 1;
                    i = n;
                }
            } // end if array=0
        } // end loop over i
        for (unsigned j = 0; j < n; j++){
            array[2 * j] *= sqrt(double(count_nonzero) / double(ilastnonzero - ifirstnonzero));
        }
    }
}


void Tools::InterpolateComplex(double *array, const unsigned n) {
    // to get rid of the zero bins
    double previous_nonzero = 0.;
    double next_nonzero = 0.;
    double check;
    unsigned ifirstnonzero = 0;
    unsigned ilastnonzero = 0;
    unsigned k;
    unsigned m = 0;
    unsigned count_nonzero = 0;
    
    // find the first nonzero even element
    while (array[2 * m] == 0){
        m++;
    }
    
    ifirstnonzero = m;
    
    // count the nonzero elements
    for (unsigned i = 0; i < n; i++) {
        if (array[2 * i] != 0)
            count_nonzero++;
    }
    
    if (count_nonzero) {
        // loop through the elements of the array and replace the zeros with interpolated values
        for (unsigned i = ifirstnonzero; i < n; i++) {
            if (array[2 * i] != 0.) {
                // set the lower nonzero value that we are interpolating from
                previous_nonzero = array[2 * i];
            }
            else{
                check = 0.;
                k = i;
                while (check == 0. && k < n) {
                    check = array[2 * k];
                    k++;
                }
                if (k < n){
                    next_nonzero = check;
                    for (unsigned j = i; j < k; j++) {
                        array[2 * j] = previous_nonzero + (next_nonzero - previous_nonzero) * (double)(j - (i - 1)) / (double)(k - (i - 1));
                        array[2 * j + 1] = array[2 * j];
                    }
                    i = k - 1;
                    previous_nonzero = next_nonzero;
                }
                else {
                    ilastnonzero = i - 1;
                    i = n;
                }
            } // end if array=0
        } // end loop over i
        for (unsigned j = 0; j < n; j++){
            array[2 * j] *= sqrt((double)count_nonzero / (double)(ilastnonzero - ifirstnonzero));
            array[2 * j + 1] *= sqrt((double)count_nonzero / (double)(ilastnonzero - ifirstnonzero));
        }
    }
}


void Tools::NormalTimeOrdering(const int n, double* volts) {
    double volts_temp[n];
    for (int i=0;i<n/2;i++) {
        volts_temp[i]=volts[i+n/2];
        volts_temp[i+n/2]=volts[i];
    }
    for (int i=0;i<n;i++) {
        volts[i]=volts_temp[i];
    }
}


int Tools::findIndex(double *freqlab,double freq,int npoints,double min,double max) {
    
    //  cout << "inside findIndex, freq, min are " << freq << " " << min << "\n";
    
    if (freq<min)
        return -1;
    if (freq>=min && freq<freqlab[1])
        return 0;
    if (freq>=freqlab[npoints-1] && freq<max)
        return npoints-1;
    return ((int)((freq-freqlab[0])/(freqlab[1]-freqlab[0])))+1;
    
    return -1;
}


void Tools::get_random_rician(double signal_amplitude, double signal_phase, double sigma, double &amplitude, double &phase){
    double rand_gauss_a, rand_gauss_b;
    get_circular_bivariate_normal_random_variable(rand_gauss_a, rand_gauss_b);
    
    // Gives the gaussian-distributed random variables a standard deviation of sigma
    rand_gauss_a *= sigma;
    rand_gauss_b *= sigma;
    
    // Gives the gaussian-distributed random variables a mean of (v*cos(theta), v*sin(theta)) when v is the mean of the desired rician distribution
    rand_gauss_a += signal_amplitude * cos(signal_phase);
    rand_gauss_b += signal_amplitude * sin(signal_phase);
    
    // The Rician Distribution produces the probability of the the absolute value (radius) of a circular bivariate normal random variable:
    amplitude = sqrt(rand_gauss_a * rand_gauss_a + rand_gauss_b * rand_gauss_b);
    // Thus, the descriptor other than amplitude for the circular bivariate is given by a phase:
    phase = atan2(rand_gauss_b, rand_gauss_a);
    return;
}


void Tools::get_circular_bivariate_normal_random_variable(double& rand_gauss_a, double& rand_gauss_b){
    double rand_uni_a = gRandom->Rndm(); //gRandom->Rndm() produces uniformly-distributed floating points in ]0,1]
    double rand_uni_b = gRandom->Rndm();
    double first_term = sqrt(-2. * log(rand_uni_a));
    // Box-Muller transform from a bivariate uniform distribution from 0 to 1 to a gaussian with mean = 0 and sigma = 1
    rand_gauss_a = first_term * cos(2. * M_PI * rand_uni_b);
    rand_gauss_b = first_term * sin(2. * M_PI * rand_uni_b);
    return;
}


int Tools::round(double number){
    return (number > 0.0) ? floor(number + 0.5) : ceil(number - 0.5);
}


double Tools::AbbyPhiCalc(double x_abby, double y_abby){
    double abbyanglephi = 0;
    
    if(x_abby>=0 && y_abby>=0) //first quadrant
        abbyanglephi= atan(y_abby/x_abby)*DEGRAD;
    if(x_abby<0 && y_abby>=0) //second quadrant
        abbyanglephi= atan(y_abby/x_abby)*DEGRAD+180.;
    if(x_abby<0 && y_abby<0) //third quadrant
        abbyanglephi= atan(y_abby/x_abby)*DEGRAD+180.;
    if(x_abby>=0 && y_abby<0) //fourth quadrant
        abbyanglephi= atan(y_abby/x_abby)*DEGRAD+360.;
    //else abbyanglephi=0.;
    if(x_abby==0 && y_abby>=0)
        abbyanglephi=90.;
    if(x_abby==0 && y_abby<=0)
        abbyanglephi=270.;
    //cout<<"abby's phi is "<< abbyanglephi<<endl;
    return abbyanglephi;
}
