#ifndef HistMan_h
#define HistMan_h 1

#include "ConfMan.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "CDSTrackingMan.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "ELossTools.h"
#include "TrackTools.h"

#include "MyTools.h"
#include "SlewingData.h"

class HistMan
{
 public: 
  HistMan(TFile *f);
  virtual ~HistMan(){};

 private:
  TFile *fFile;
  int fStatus;

  std::vector<HodoscopeLikeHit*> fBHDhit;
  std::vector<HodoscopeLikeHit*> fT0hit;
  std::vector<HodoscopeLikeHit*> fBPDhit;
  std::vector<HodoscopeLikeHit*> fDEFhit;
  std::vector<HodoscopeLikeHit*> fBVChit;
  std::vector<HodoscopeLikeHit*> fCVChit;
  std::vector<HodoscopeLikeHit*> fPChit;
  std::vector<HodoscopeLikeHit*> fNChit[8];
  std::vector<HodoscopeLikeHit*> fBDhit;

 public:
  void initTrig();
  void fill(EventHeader *header, ConfMan *conf);
  void fill(EventHeader *header, ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, BeamSpectrometer *D5);

  void set(BeamLineHitMan *blMan);
  void clear();
};
#endif
