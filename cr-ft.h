#pragma once

# include <stdlib.h>
# include "TGraph.h"
# include "TCanvas.h"
# include "buffer.hh"
# include "TLine.h"
# include "TLegend.h"

#define ANITA_TIME_SAMPLES 512
#define ANITA_FT_SAMPLES (int(ANITA_TIME_SAMPLES) / 2 + 1) //!< Number of FFTW freq bins needed to recover ANITA_TIME_SAMPLES in time domain.
#define ANITA_FT_BINS (ANITA_FT_SAMPLES - 1)
#define ANITA_FREQ_HIGH 1300e+6 //!< Hz. Nyquist frequency corresponding to ANITA's 2.6 GSa / s.

struct cr_ft_state {
  TCanvas *cZhsEAndAlpha;
  TPad *px1;
  TGraph *grZhsTimeE;
  TLine *line_min;
  TLine *line_max;
  TPad *px2;
  TGraph *grZhsAlpha;
  TLegend *legend;
  std::unique_ptr<double[]> ZhsFftInp;
  std::unique_ptr<TGraph> grFftRe;
  std::unique_ptr<TGraph> grFftIm; 
  std::unique_ptr<TGraph> grFftRho;
  std::unique_ptr<TGraph> grFftPhi;
  double *FftRho; //!< Amplitude of ZHS FT.
  double *FftPhi; //!< Phase of ZHS FT.
  TPad *panel_ft_rho;
  //unique_ptr<TPad> panel_ft_rho;
  TPad *panel_ft_phi;
  //unique_ptr<TPad> panel_ft_phi;
  std::unique_ptr<FFTWComplex[]> ZhsFft;
  FFTWComplex AnitaFT[ANITA_FT_SAMPLES]; //!< Frequency distribution sampled on the ANITA frequency grid. 
  TCanvas *cZhsFft;
  int ind_maxval;
  double fwhm_xmin;
  double fwhm_xmax;
  double fwhm_xmaxval;
  double fwhm_ymaxval;
  double vis_xmin;
  double vis_xmax;
  int vis_xmin_bin;
  int vis_xmax_bin;
  double dfreq;
  bvv::TBuffer <int> vis_nbins;
  TCanvas *cZhsIFft;
  std::unique_ptr<TGraph> grIFft;
  std::unique_ptr<TGraph> grZhsTimeERec;
  std::unique_ptr<double[]> ZhsIFft;
  bvv::TBuffer <int> BinShift;
};
