#pragma once

# include <stdlib.h>
# include "TGraph.h"
# include "TCanvas.h"
# include "buffer.hh"
# include "TLine.h"
# include "TLegend.h"

struct hot_test_state {
  TCanvas *cHotTest;
  TGraph *grHotTest;
  std::unique_ptr<TF1> fHotTest;
};
