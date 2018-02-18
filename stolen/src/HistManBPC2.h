#ifndef HistManBPC2_h
#define HistManBPC2_h 1

#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "DCEffMan.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

class HistManBPC2
{
 public:
  HistManBPC2(TFile *f);
  virtual ~HistManBPC2(){};

 private:
  TFile *fFile;

 public:
  void fill(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, DCEffMan *effMan);
  TH1* getTH(const char *name);
  TH2* getTH2(const char *name);
};
#endif
