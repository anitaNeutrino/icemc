#pragma once

# include <stdlib.h>
# include "TGraph.h"
# include "TCanvas.h"
# include "buffer.hh"
# include "TLine.h"
# include "TLegend.h"

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
  bvv::TBuffer <int> vis_nbins;
  TCanvas *cZhsIFft;
  std::unique_ptr<TGraph> grIFft;
  std::unique_ptr<TGraph> grZhsTimeERec;
  std::unique_ptr<double[]> ZhsIFft;
  bvv::TBuffer <int> BinShift;
};
