#ifndef HistManwMC_h
#define HistManwMC_h 1

#include <iostream>
#include <string>
#include "KnuclRootData.h"
#include "BeamLineHitMan.h"
#include "CDSHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"
#include "TrackTools.h"
#include "ELossTools.h"

class HistManwMC
{
 public:
  HistManwMC(TFile *f, ConfMan *conf);
  virtual ~HistManwMC(){};

 private:
  // Output TFile pointer;
  TFile *rtFile;
  ConfMan *confMan;
  // Data Manager for knucl
  RunHeaderMC *runHeaderMC;
  EventHeaderMC *evHeaderMC;
  MCData *mcData;
  DetectorData *detData;
  ReactionData *reacData;
  // Data Manager for k18ana
  BeamLineHitMan *blMan;
  CDSHitMan *cdsMan;
  BeamLineTrackMan *bltrackMan;
  CDSTrackingMan *cdstrackMan;
  // for Analysis Data
  double fT0time;
  double fD5mom;
  LocalTrack *fTrackBPC;
  TLorentzVector fBeamLmom;
  TVector3 fT0pos;
  double fNCtime;
  double fNCdE;
  TVector3 fNCpos;
  TVector3 fVtxCDS; // Best Vertex Helix
  TVector3 fVtxBeam; // Best Vertex Beam
  TLorentzVector fFLmom;

  std::vector<HodoscopeLikeHit*> fT0_hit;
  std::vector<HodoscopeLikeHit*> fBPD_hit;
  std::vector<HodoscopeLikeHit*> fCDH_hit;
  std::vector<HodoscopeLikeHit*> fBVC_hit;
  std::vector<HodoscopeLikeHit*> fCVC_hit;
  std::vector<HodoscopeLikeHit*> fPC_hit;
  std::vector<HodoscopeLikeHit*> fNC_hit[8];
  std::vector<HodoscopeLikeHit*> fBD_hit;

  std::vector<CDSTrack*> fTrackPim;
  std::vector<CDSTrack*> fTrackKm;
  std::vector<CDSTrack*> fTrackPip;
  std::vector<CDSTrack*> fTrackP;
  std::vector<CDSTrack*> fTrackD;

  int fFPID;
  double fNCbeta;

 public:
  void setMC(RunHeaderMC *runhead, EventHeaderMC *evhead, MCData *mcdata, DetectorData *detdata, ReactionData *reacdata)
  {
    runHeaderMC=runhead, evHeaderMC=evhead, mcData=mcdata, detData=detdata, reacData=reacdata;
  }
  void setK18ana(BeamLineHitMan *blman, CDSHitMan *cdsman, BeamLineTrackMan *bltrackman, CDSTrackingMan *cdstrackman)
  {
    blMan=blman, cdsMan=cdsman, bltrackMan=bltrackman, cdstrackMan=cdstrackman;
  }

  void setT0time(const double &time){ fT0time=time; };
  void setD5mom(const double &mom){ fD5mom=mom; };
  void setTrackBPC(LocalTrack *track){ fTrackBPC=track; };
		 
  double T0time() const { return fT0time; };
  double D5mom() const { return fD5mom; };
  LocalTrack* trakcBPC() const { return fTrackBPC; };

  int nT0() const { return fT0_hit.size(); };
  int nBPD() const { return fT0_hit.size(); };
  int nCDH() const { return fT0_hit.size(); };
  int nBVC() const { return fT0_hit.size(); };
  int nCVC() const { return fT0_hit.size(); };
  int nPC() const { return fT0_hit.size(); };
  int nNC() const {
    int nNC = 0;
    for( int i=0; i<8; i++ ) nNC += fNC_hit[i].size();
    return nNC;
  }
  int nNC(const int &lay){ return fNC_hit[lay-1].size(); };
  int nBD() const { return fBD_hit.size(); };

  HodoscopeLikeHit *T0(const int &i){ return fT0_hit[i]; };
  HodoscopeLikeHit *BPD(const int &i){ return fBPD_hit[i]; };
  HodoscopeLikeHit *CDH(const int &i){ return fCDH_hit[i]; };
  HodoscopeLikeHit *BVC(const int &i){ return fBVC_hit[i]; };
  HodoscopeLikeHit *CVC(const int &i){ return fCVC_hit[i]; };
  HodoscopeLikeHit *PC(const int &i){ return fPC_hit[i]; };
  HodoscopeLikeHit *NC(const int &lay, const int &i){ return fNC_hit[lay-1][i]; };
  HodoscopeLikeHit *BD(const int &i){ return fBD_hit[i]; };

  void initHist();
  bool ana(); // return false if not find CDS GoodTrack w/ CDH hit
  void fill();
  void finit();

  bool Clear(const bool flag=true);

 public:
  void dump();
};
#endif
