# include <iostream>
# include <stdlib.h>
# include "TCanvas.h"
# include "TSystem.h"
# include "game.h"
# include "TGraph.h"
# include "TGaxis.h"
# include "TLine.h"
# include "TH1F.h"
# include "TLegend.h"
# include "FFTtools.h"

// extern std::string some_imortant_str;
using namespace std;

struct game_state {
  TCanvas *cZhsEAndAlpha;
  TPad *px1;
  TGraph *grZhsTimeE;
  TLine *line_min;
  TLine *line_max;
  TPad *px2;
  TGraph *grZhsAlpha;
  TLegend *legend;
  double *ZhsFftInp;
  TGraph *grFft;
  FFTWComplex *ZhsFft;
  TCanvas *cZhsFft;
};

extern const double pi;
extern int ZhsTimeN;
extern double ZhsTimeStart; 
extern double ZhsTimeDelta;
extern vector<double> ZhsTimeArr;
extern vector<double> ZhsTimeE;
extern vector<double> ZhsAlpha;
bool FWHM(long int n, double *x, double *y, double &xmin, double &xmax, int &ind_maxval, double &xmaxval, double &ymaxval, double threshold_rel = 0.5);


static struct game_state *game_init()
{
  struct game_state *state = (game_state *) malloc(sizeof(*state));
  state->cZhsEAndAlpha = new TCanvas();
  // state->c = new TCanvas("c", "c", 800, 600);
  state->px1 = new TPad("px1","",0,0,1,1);
  state->px1->Draw();
  state->px1->cd();
  state->grZhsTimeE = new TGraph(ZhsTimeN, ZhsTimeArr.data(), ZhsTimeE.data());
  state->grZhsTimeE->SetLineColor(kRed);
  state->grZhsTimeE->SetLineWidth(2);
  state->grZhsTimeE->Draw("AL");
  double xmin = 183.95e+3;
  double xmax = 183.98e+3;
  double xmaxval;
  double ymaxval;
  int ind_maxval;
  bool fwhm_res = FWHM(ZhsTimeN, ZhsTimeArr.data(), ZhsTimeE.data(), xmin, xmax, ind_maxval, xmaxval, ymaxval);
  double vis_xmin = xmaxval - 20*(xmaxval - xmin);
  double vis_xmax = xmaxval + 20*(xmax - xmaxval); 
  int vis_xmin_bin = (vis_xmin - ZhsTimeStart) / ZhsTimeDelta;
  int vis_xmax_bin = (vis_xmax - ZhsTimeStart) / ZhsTimeDelta;
  int vis_nbins = vis_xmax_bin - vis_xmin_bin + 1;
  printf("Waveform xmax - xmin: %5.3f, (%5.3f%%)\n", vis_xmax - vis_xmin, (vis_xmax - vis_xmin) / ZhsTimeDelta);
  state->grZhsTimeE->GetXaxis()->SetLimits(vis_xmin, vis_xmax);
  TGaxis::SetExponentOffset(0.02, -0.04, "x");
  state->grZhsTimeE->GetXaxis()->SetTitle("Time [ns]");
  state->grZhsTimeE->GetYaxis()->SetTitle("E [V/m]");
  state->grZhsTimeE->Draw("AL");
  gPad->Update();

  state->line_min = new TLine(xmin, gPad->GetUymin(), xmin, gPad->GetUymax());
  state->line_max = new TLine(xmax, gPad->GetUymin(), xmax, gPad->GetUymax());
  state->line_min->Draw();
  state->line_max->Draw();
  gPad->Update();


  state->px2 = new TPad("px2","",0,0,1,1);
  state->px2->SetFillStyle(4000);
  state->px2->SetFrameFillStyle(0);
  state->px2->Draw();
  state->px2->cd();


  state->grZhsAlpha = new TGraph(ZhsTimeN, ZhsTimeArr.data(), ZhsAlpha.data());
  state->grZhsAlpha->SetLineColor(kBlue);
  state->grZhsAlpha->SetLineWidth(2);
  
  state->grZhsAlpha->GetXaxis()->SetLimits(vis_xmin, vis_xmax);
  state->grZhsAlpha->GetHistogram()->SetMaximum(180);
  state->grZhsAlpha->GetHistogram()->SetMinimum(0);

  state->px2->Update();
  state->grZhsAlpha->Draw("ALY+");
  state->legend = new TLegend(0.1,0.8,0.3,0.9);
  state->legend->AddEntry(state->grZhsTimeE, "E_{proj}","l");
  state->legend->AddEntry(state->grZhsAlpha, "pol. angle","l");
  state->legend->Draw();

  fflush(stdout);

  state->ZhsFftInp = new double[vis_nbins];

  printf("vis_nbins: %d, vis_xmin_bin: %d, ZhsTimeE.size(): %zu\n", vis_nbins, vis_xmin_bin, ZhsTimeE.size());
  for (int i = 0; i < vis_nbins; i++) {
    double val = *(ZhsTimeE.data() + vis_xmin_bin + i);
    state->ZhsFftInp[i] = val;
    printf("ZhsFftInp[%d]: %11.8e\n", i, val);
  }


  //  fftw_complex *in, *out;
  //  fftw_plan p;
  //  in = (fftw_complex*) malloc(sizeof(fftw_complex) * vis_nbins);
  //  out = (fftw_complex*) malloc(sizeof(fftw_complex) * vis_nbins);
  //  p = fftw_plan_dft_1d(vis_nbins, in, out, FFTW_FORWARD, FFTW_ESTIMATE); 
   double xmax_normal = 5;
   double xmin_normal = -5;
   int N = vis_nbins;
   double dx = (xmax_normal - xmin_normal) / N;
    for (int i = 0; i < vis_nbins; i++){
      double x = xmin_normal + dx * i;
      double y = exp(-0.5 * x * x) / (sqrt(2.0) * sqrt(pi));
      state->ZhsFftInp[i] = y;
  //    in[i][0] = y; // ZhsFftInp[i]; 
  //    in[i][1] = 0; 
    }
  //  fftw_execute(p); /* repeat as needed */
  state->grFft = new TGraph(vis_nbins / 2);
  state->ZhsFft = FFTtools::doFFT(vis_nbins, state->ZhsFftInp);
  for (int i = 0; i < vis_nbins / 2; i++){
  // grFft->SetPoint(i, i, ZhsFft[i].getAbsSq() /* * ZhsFft[i].re */);
    state->grFft->SetPoint(i, i / (xmax_normal - xmin_normal), state->ZhsFft[i].re * dx); // To get continuous ft values.
  //   grFft->SetPoint(i, i / (xmax_normal - xmin_normal), out[i][0] * dx); // To get continuous ft values.
  //   // printf("ZhsFft: %d, %11.4e\n", i, ZhsFft[i].re);
  //   // printf("out: %d, %11.4e\n", i, ZhsFft[i].re);
  //   printf("out: %d, %10.4e\n", i, out[i][0] * dx);
  }

  // fftw_destroy_plan(p);
  // fftw_free(in); fftw_free(out);
  // free(in); free(out);
  


  state->cZhsFft = new TCanvas();
  state->grFft->Draw("AL");
  state->cZhsFft->Update();
  printf("Done");
  

  return state;
}

static void game_finalize(struct game_state *state)
{
  delete state;
}

static void game_reload(struct game_state *state)
{
  // state->c->Clear();

  // xmin, xmax: fwhm interval endpoints.
  double xmin = 183.95e+3;
  double xmax = 183.98e+3;
  double xmaxval;
  double ymaxval;
  int ind_maxval;
  bool fwhm_res = FWHM(ZhsTimeN, ZhsTimeArr.data(), ZhsTimeE.data(), xmin, xmax, ind_maxval, xmaxval, ymaxval);
  // vis_xmin, vis_xmax: which range to visualize, should encompass (xmin, xmax).
  double vis_xmin = xmaxval - 200*(xmaxval - xmin);
  double vis_xmax = xmaxval + 200*(xmax - xmaxval); 
  int vis_xmin_bin = (vis_xmin - ZhsTimeStart) / ZhsTimeDelta;
  int vis_xmax_bin = (vis_xmax - ZhsTimeStart) / ZhsTimeDelta;
  int vis_nbins = vis_xmax_bin - vis_xmin_bin + 1;
  printf("Waveform xmax - xmin: %5.3f, (%5.3f%%)\n", vis_xmax - vis_xmin, (vis_xmax - vis_xmin) / ZhsTimeDelta);
  state->grZhsTimeE->GetXaxis()->SetLimits(vis_xmin, vis_xmax);
  state->grZhsTimeE->GetXaxis()->SetTitle("Time [ns]");
  state->grZhsTimeE->GetYaxis()->SetTitle("E [V/m]");
  // state->grZhsTimeE->Draw("AL");
  // state->px1->Update();


  state->grZhsAlpha->GetXaxis()->SetLimits(vis_xmin, vis_xmax);
  state->grZhsAlpha->GetHistogram()->SetMaximum(180);
  state->grZhsAlpha->GetHistogram()->SetMinimum(0);
  state->legend->SetX1NDC(0.12);
  state->legend->SetX2NDC(0.3);
  state->legend->SetY1NDC(0.75);
  state->legend->SetY2NDC(0.88);
  // state->px2->Update();

  state->px1->Modified();
  state->px2->Modified();
  state->cZhsEAndAlpha->Modified();
  state->cZhsEAndAlpha->Update();


  state->cZhsFft->Clear();
 
  delete[] state->ZhsFftInp;
  state->ZhsFftInp = new double[vis_nbins];

  // printf("vis_nbins: %d\n", vis_nbins);
  // state->ZhsFftInp = new double[vis_nbins];
  // for (int i = 0; i < vis_nbins; i++) {
  //   double val = *(ZhsTimeE.data() + vis_xmin_bin + i);
  //   state->ZhsFftInp[i] = val;
  //   printf("ZhsFftInp[%d]: %11.8e\n", i, val);
  // }

   double xmin_normal = -30;
   double xmax_normal = +30;
   int N = vis_nbins;
   double dx = (xmax_normal - xmin_normal) / N;
    for (int i = 0; i < vis_nbins; i++){
      double x = xmin_normal + dx * i;
      double y = exp(-0.5 * x * x) / (sqrt(2.0) * sqrt(pi));
      state->ZhsFftInp[i] = y;
  //    in[i][0] = y; // ZhsFftInp[i]; 
  //    in[i][1] = 0; 
    }

  if (state->grFft) delete state->grFft;
  state->grFft = new TGraph(vis_nbins / 2);
  printf("Before deleting ZhsFft\n");
  // delete[]: Is it a right thing to do?
  // http://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs:
  // If you allocate an array with fftw_malloc() you must deallocate it with fftw_free(). Do not use free() or, heaven forbid, _delete_. 
  // doFFT is not using fftw_malloc() though, unless "new" is overloaded.
  if (state->ZhsFft) delete[] state->ZhsFft;
  printf("After deleting ZhsFft\n");
  state->ZhsFft = FFTtools::doFFT(vis_nbins, state->ZhsFftInp);
  for (int i = 0; i < vis_nbins / 2; i++){
    state->grFft->SetPoint(i, i / (xmax_normal - xmin_normal), state->ZhsFft[i].re * dx); // To get continuous ft values.
  }
  state->cZhsFft->cd();
  state->grFft->Draw("AL");

  state->grFft->SetLineColor(kBlue);
  state->grFft->SetLineWidth(1);
  state->cZhsFft->Modified();
  state->cZhsFft->Update();
  cout << "reloaded dl" << endl;
  
}

static void game_unload(struct game_state *state)
{
  cout << "game unloaded" << endl;
}


static bool game_step(struct game_state *state)
{
  // state->c->Update();

  gSystem->ProcessEvents();
  return true;
}

__attribute__((visibility("default"))) extern const struct game_api GAME_API = {
  game_init, // game_init
  game_finalize, // game_finalize
  game_reload, // game_reload
  game_unload, // game_unload
  game_step // game_step
};
