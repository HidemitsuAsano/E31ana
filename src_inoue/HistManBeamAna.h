#ifndef HistManBeamAna_h
#define HistManBeamAna_h 1

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

//static const double BEAM_MOM = 1.1;
static const double L_MIN = 1.1097;
static const double L_MAX = 1.1217;
static const double K0_MIN   = 0.4876;
static const double K0_MAX   = 0.5076;

class HistManBeamAna
{
 public:
  HistManBeamAna(TFile *f, ConfMan *conf);
  virtual ~HistManBeamAna();

 private:
  TFile *fFile;
  int fBeamPID;
  double fBeamMom;
  BeamSpectrometer *fD5;
  std::vector<HodoscopeLikeHit*> fBHD_hit;
  std::vector<HodoscopeLikeHit*> fT0_hit;
  std::vector<HodoscopeLikeHit*> fBPD_hit;
  std::vector<HodoscopeLikeHit*> fDEF_hit;
  std::vector<HodoscopeLikeHit*> fCVC_hit;
  std::vector<HodoscopeLikeHit*> fBVC_hit;
  std::vector<HodoscopeLikeHit*> fPC_hit;
  std::vector<HodoscopeLikeHit*> fNC_hit[8];
  std::vector<HodoscopeLikeHit*> fLB_hit;

  std::vector<HodoscopeLikeHit*> fCDH_hit;
  std::vector<CDSTrack*> fCDSpim;
  std::vector<CDSTrack*> fCDSkm;
  std::vector<CDSTrack*> fCDSpip;
  std::vector<CDSTrack*> fCDSp;
  std::vector<CDSTrack*> fCDSd;

  double fVtxDis;
  TVector3 fVtxBeam;
  TVector3 fVtxCDS;
  TVector3 fVtxBeam_w;
  TVector3 fVtxCDS_w;
  SlewingData *fSlewingBHDT0;
  SlewingData *fSlewingT0NC;
  SlewingData *fSlewingT0CVC;
  SlewingData *fSlewingT0PC;

  int fNumNEv;
  int fGoodNEv;

  bool fLflag;
  bool fK0flag;

 public:
  void fill(ConfMan *conf, EventHeader *header);
  void fill(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan);
  void fillFC(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bktrackMan);
  void fillFC(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bktrackMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan);
  TH1* getTH(const char *name);
  TH2* getTH2(const char *name);
  double beam_mom() const { return fBeamMom; };

  void clear();

 protected:
  void initBLDC();
  void initFC();
  void fillBLDC(EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bktrackMan);

 public:
  void initCDS();
  void fillCDS(ConfMan *conf, EventHeader *header, CDSHitMan *cdsMan, BeamLineHitMan *blMan, CDSTrackingMan *cdstrackMan, BeamLineTrackMan *bltrackMan); 
};
#endif
