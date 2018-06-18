#ifndef DispFDC_h
#define DispFDC_h 1

#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "ConfMan.h"

#include "TCanvas.h"
#include "TMarker.h"
#include "TLine.h"
#include "TH2.h"
#include "FDCTrackMan.h"
#include "TArc.h"

class DispFDC
{
 public:
  DispFDC(ConfMan *conf);
  virtual ~DispFDC(){};

 private:
  ConfMan *fConfMan;
  TCanvas *fCanvas;
  TH2F *fH2_U;
  TH2F *fH2_X;
  TH2F *fH2_V;

 public:
  bool draw(BeamLineHitMan *blMan, BeamLineTrackMan *trackMan);

};
#endif
