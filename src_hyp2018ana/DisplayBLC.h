#ifndef DisplayBLC_h
#define DisplayBLC_h 1

#include "BeamLineHitMan.h"
#include "ConfMan.h"
#include "BeamLineTrackMan.h"

#include "TCanvas.h"
#include "TMarker.h"
#include "TLine.h"

class DisplayBLC
{
 public:
  DisplayBLC();

 private:
  TCanvas *fCanvas;
  TH2F *fBLC2a_xz;
  TH2F *fBLC2a_yz;
  TH2F *fBLC2b_xz;
  TH2F *fBLC2b_yz;

 public:
  void Draw(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan);

  bool Wait();
};
#endif 
