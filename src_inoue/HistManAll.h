#ifndef HistManAll_h
#define HistManAll_h 1

#include "ConfMan.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "ELossTools.h"

class HistManAll
{
 public:
  HistManAll(TFile *f, ConfMan *conf);
  virtual ~HistManAll();

 private:
  TFile *fFile;
  int fBeamPID;
  BeamSpectrometer *fD5;
  std::vector<HodoscopeLikeHit*> fBHD_hit;
  std::vector<HodoscopeLikeHit*> fT0_hit;
  std::vector<HodoscopeLikeHit*> fBPD_hit;
  std::vector<HodoscopeLikeHit*> fDEF_hit;
  std::vector<HodoscopeLikeHit*> fCVC_hit;
  std::vector<HodoscopeLikeHit*> fPC_hit;
  std::vector<HodoscopeLikeHit*> fNC_hit[8];

 public:
  void fill(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan);
  TH1* getTH(const char *name);
  TH2* getTH2(const char *name);

  void clear();

 protected:
  void initBLDC();
  void initFC();
  void fillBLDC(EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bktrackMan);
  void fillFC(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bktrackMan);
};
#endif
