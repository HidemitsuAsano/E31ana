#ifndef DispBPC_h
#define DispBPC_h 1

#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "ConfMan.h"

#include "TCanvas.h"
#include "TMarker.h"
#include "TLine.h"
#include "TH2.h"
#include "BPCTrackMan.h"
#include "TArc.h"

class DispBPC
{
 public:
  DispBPC(ConfMan *conf);
  virtual ~DispBPC(){};

 private:
  ConfMan *fConfMan;
  TCanvas *fCanvas;
  TH2F *fH2_xz;
  TH2F *fH2_yz;

 public:
  bool draw(BeamLineHitMan *blMan, BeamLineTrackMan *trackMan, BPCTrackMan &BPCTrack);

};
#endif
