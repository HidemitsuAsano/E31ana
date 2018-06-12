#ifndef SLEWINGTOOLS_HH
#define SLEWINGTOOLS_HH 1

#include "ConfMan.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TF1.h"
#include "CalibTools.h"
#include "TProfile.h"

namespace SlewingTools
{
  void MakeHist();
  void Fill(int cid, int seg, ConfMan *conf, TNtuple *tup);
  void Fill(int cid, int seg, ConfMan *conf, TNtuple *tup, double mean, double thre);
  void Draw();
  bool Reset();
  TF1* GetTF1(int type);
  TF1* MakeTF1(int type);
  TF1* Fit(TH2F *h2, int type);

  void SetRange(TH2F *h2);
  TH1F* projectionX(TH2F *h2);
  TH1F* projectionY(TH2F *h2);
  double TimeReso();
  
  std::vector<double> SetParam(int cid, int seg, int ud, int ith, int type, TH2F *h2, ConfMan *conf, bool add_old=true);
  std::vector<double> SetParam(int cid, int seg, int ud, int ith, int type, TH2F *h2, ConfMan *conf, double weight, bool add_old=true);
};
#endif
