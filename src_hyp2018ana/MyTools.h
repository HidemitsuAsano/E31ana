#ifndef MYTOOLS_H
#define MYTOOLS_H 1

#include <string>
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TMath.h"
#include "GlobalVariables.h"
#include "GeomTools.h"
#include "MathTools.h"
#include "CDSTrack.h"
#include "CDSTrackingMan.h"
#include "Particle.h"
#include "LocalTrack.h"
#include "TrackTools.h"
#include "KnuclRootData.h"
#include "MyParam.h"

#include "AnaInfo.h"
#include "EventHeader.h"

namespace MyTools
{
  double cosOP(CDS2Info *info);
  double calc_offset(const double &mm, const TLorentzVector tot_lmom, const TLorentzVector lmom, AnaInfo *anaInfo);

  bool searchCDHHit(CDSTrack *track, CDSHitMan *cdsMan, bool z_match=true);
  bool searchIHHit(CDSTrack *track, CDSHitMan *cdsMan);
  double dphiCDH(CDSTrack *track, HodoscopeLikeHit *hit);
  double dphiIH(CDSTrack *track, HodoscopeLikeHit *hit);
  bool sameCDHCluster(HodoscopeLikeHit *hit1, HodoscopeLikeHit *hit2);

  bool reanaFC(ConfMan *conf, BeamLineHitMan *blMan, AnaInfo *anaInfo, bool refit=false);
  ForwardChargeInfo makeFC(ConfMan *conf, BeamLineHitMan *blMan, AnaInfo *anaInfo);
  ForwardNeutralInfo makeFN(BeamLineHitMan *blMan, AnaInfo *anaInfo, const double &thre=8.0, bool sim=false);
  BeamInfo makeBeamInfo(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, BeamSpectrometer *spec);
  CDSInfo makeCDSInfo(int id, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, BeamInfo *beam, BeamLineTrackMan *bltrackMan, ConfMan *conf=0, bool sim=false);

  int CDSanaLevel(CDSTrack *track, CDSHitMan *cdsMan, BeamInfo *beam, BeamLineTrackMan *bltrackMan, ConfMan *conf=0, bool sim=false);
  CDSInfo makeCDSInfo_z40(int id, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, BeamInfo *beam, BeamLineTrackMan *bltrackMan, ConfMan *conf=0, bool sim=false);
  CDS2Info makeCDS2Info(CDSTrackingMan *cdstrackMan, int id1, int id2, AnaInfo *info);
  ForwardNeutralInfo makeFN(BeamLineHitMan *blMan, AnaInfo *anaInfo);

  std::vector<std::vector<HodoscopeLikeHit*> > getNChits(BeamLineHitMan *blMan);
  HodoscopeLikeHit *getCDH(CDSTrack *track, CDSHitMan *cdsMan);
  HodoscopeLikeHit *getCDH2(CDSTrack *track, CDSHitMan *cdsMan); //search true hit in cluster by Z position
  HodoscopeLikeHit *getIH(CDSTrack *track, CDSHitMan *cdsMan);

  bool isBeamFiducial(BeamLineTrackMan *bltrackMan);
  bool isBeamMom(BeamSpectrometer *beam);
  LocalTrack *trackBLC1(BeamLineTrackMan *bltrackMan);
  LocalTrack *trackBLC2(BeamLineTrackMan *bltrackMan);
  LocalTrack *trackBPC(BeamLineTrackMan *bltrackMan);
  bool BLC2BPC(BeamLineTrackMan *bltrackMan);
  bool TOF_K(BeamLineHitMan *blMan);
  bool TOF_pi(BeamLineHitMan *blMan);

  std::vector<HodoscopeLikeHit*> getCDHCluster(HodoscopeLikeHit *seed, CDSHitMan *cdsMan);
  std::vector<HodoscopeLikeHit*> getCDH(CDSHitMan *cdsMan);
  std::vector<HodoscopeLikeHit*> getIH(CDSHitMan *cdsMan);
  std::vector<HodoscopeLikeHit*> getHodo(BeamLineHitMan *blMan, const int &cid);
  bool CDCall1hit(CDSHitMan *cdsMan);

  template <class T> 
  T* get(const TString &name, TFile *f=gFile);
			 
  int checkArea(double mom, double mass2); // if mass2 between pi- and K- return -1; pi+ p return 1; else return 0;
  int PID(double mom,double mass2);
  int PIDcorr(double mom, double mass2);
  int PIDcorr2(double mom, double mass2);
  int PIDcorr_wide(double mom, double mass2);
  int PIDcorr_tigt(double mom, double mass2);

  bool IsTarget0(const TVector3 &pos, ConfMan *conf);
  bool IsTarget1(const TVector3 &pos, ConfMan *conf);
  bool IsTarget2(const TVector3 &pos, ConfMan *conf);
  bool IsTarget3(const TVector3 &pos, ConfMan *conf);
  bool IsTarget4(const TVector3 &pos, ConfMan *conf);
  bool IsTarget5(const TVector3 &pos, ConfMan *conf);

  bool IsTarget0(const std::vector<TVector3> &pos, ConfMan *conf);
  bool IsTarget1(const std::vector<TVector3> &pos, ConfMan *conf);
  bool IsTarget2(const std::vector<TVector3> &pos, ConfMan *conf);
  bool IsTarget3(const std::vector<TVector3> &pos, ConfMan *conf);
  bool IsTarget4(const std::vector<TVector3> &pos, ConfMan *conf);
  bool IsTarget5(const std::vector<TVector3> &pos, ConfMan *conf);

  bool IsElectron(const double &beta, const double &mom);

  int generation(const int &parentID, MCData *mcData);
  int trackStatus(const int &trackID, MCData *mcData);

  int nStatus(const int &trackID, MCData *mcData);
  int nStatus(const int &seg, DetectorData *detData, MCData *mcData);
  int gammaStatus(const int &trackID, MCData *mcData);

  void printCDSParams(CDSTrack *track);

  bool isShareCDH(CDSHitMan *cdsMan, CDSTrackingMan *trackMan, const int &id);

  std::string getDetName(const TVector3 &pos);

  bool isFiducial(AnaInfo *info);
  bool isFiducial(CDSInfo *info);
};
#endif

