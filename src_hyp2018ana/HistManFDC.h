#ifndef HistManFDC_h
#define HistManFDC_h 1

#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "DCEffMan.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

class HistManFDC
{
 public:
  HistManFDC(TFile *f);
  virtual ~HistManFDC(){};

 private:
  TFile *fFile;

 public:
  void fill(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, DCEffMan *effMan, TKOHitCollection *tko);
  TH1* getTH(const char *name);
  TH2* getTH2(const char *name);
};
#endif
